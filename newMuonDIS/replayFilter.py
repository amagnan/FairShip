#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Skim replay cbmsim events, preserving all branches of each accepted event."""

import argparse
import json
import os
import tempfile
from itertools import combinations
from pathlib import Path

from analyseEvents import Geometry, make_vertex, selected_tracks


def selection_stage(event, geometry, pdg):
    """Return 0 (<2 selected tracks), 1 (no accepted pair), or 2 (keep event)."""
    tracks, _ = selected_tracks(event, pdg)
    if len(tracks) < 2:
        return 0
    for first, second in combinations(tracks, 2):
        if first.charge * second.charge >= 0:
            continue
        vertex = make_vertex(first, second)
        if vertex is not None and vertex.doca < 5.0 and geometry.clearance(vertex.position) >= 0:
            return 2
    return 1


def filter_file(root, input_file, output_file, geometry_file=None, start_event=0, n_events=-1):
    """Filter one file; the caller must have loaded libShipMuDIS.so once.

    Keep the original tree schema, track indices, weights and normalization.
    Other ROOT keys are copied unchanged; replayFilter records this skim's
    provenance and raw event counts. Publish only a successfully closed file.
    """
    from MuDISGenerator_postProcessing import _copy_metadata

    if start_event < 0 or n_events < -1:
        raise ValueError("Counts must be nonnegative, except n_events=-1 for all remaining entries")
    input_path = Path(input_file).expanduser().resolve()
    output_path = Path(output_file).expanduser().absolute()
    if output_path.exists() or output_path.is_symlink():
        raise FileExistsError(f"Output must not exist: {output_path}")
    if not input_path.is_file():
        raise ValueError(f"Input file does not exist: {input_path}")
    geometry_path = (
        Path(geometry_file).expanduser().resolve()
        if geometry_file
        else input_path.with_name("geo_" + input_path.name.removeprefix("sim_"))
    )
    geometry = Geometry(root, str(geometry_path))
    geometry.load(input_path)
    pdg = root.TDatabasePDG.Instance()
    with root.TFile.Open(str(input_path), "READ") as source:
        tree = source.Get("cbmsim")
        if not tree or not tree.InheritsFrom("TTree"):
            raise ValueError(f"Missing cbmsim tree: {input_path}")
        required = ("MCTrack", "strawtubesPoint", "TimeDetPoint")
        missing = [name for name in required if not tree.GetBranch(name)]
        if missing:
            raise ValueError(f"Missing branches in {input_path}: {missing}")
        total = int(tree.GetEntries())
        if start_event > total:
            raise ValueError(f"start_event {start_event} exceeds the {total} available entries")
        stop = total if n_events == -1 else min(total, start_event + n_events)
        summary = {
            "source_file": str(input_path),
            "source_uuid": str(source.GetUUID().AsString()),
            "geometry_file": str(geometry_path),
            "source_entry_start": start_event,
            "source_entry_stop": stop,
            "entries_available": total,
            "events_processed": stop - start_event,
            "events_two_tracks": 0,
            "events_selected": 0,
            "selection": "Tr1-4 + TD, p > 1 GeV/c; at least one opposite-sign pair with "
            "backward vertex, DOCA < 5 cm, inside decay helium; no clearance or IP cut",
        }
        print(f"Processing entries [{start_event}, {stop}) of {total}: {input_path}", flush=True)
        # Keep the source alive: CloneTree shares its branch addresses. All
        # branches stay enabled, including those not used for the selection.
        tree.SetBranchStatus("*", True)
        with tempfile.TemporaryDirectory(prefix=".replayFilter-", dir=output_path.parent) as temporary:
            temporary_path = Path(temporary) / "selected.root"
            with root.TFile.Open(str(temporary_path), "RECREATE") as destination:
                destination.SetCompressionSettings(source.GetCompressionSettings())
                destination.cd()
                selected = tree.CloneTree(0)
                if not selected:
                    raise OSError("Cannot clone cbmsim tree")
                for entry in range(start_event, stop):
                    if tree.GetEntry(entry) <= 0:
                        raise OSError(f"Cannot read cbmsim entry {entry}")
                    stage = selection_stage(tree, geometry, pdg)
                    summary["events_two_tracks"] += int(stage >= 1)
                    if stage == 2:
                        if selected.Fill() < 0:
                            raise OSError(f"Cannot write cbmsim entry {entry}")
                        summary["events_selected"] += 1
                    if (entry - start_event + 1) % 1000 == 0:
                        print(
                            f"Processed {entry - start_event + 1}/{stop - start_event}; "
                            f"selected {summary['events_selected']}",
                            flush=True,
                        )
                destination.cd()
                if selected.Write("cbmsim", root.TObject.kOverwrite) <= 0:
                    raise OSError("Cannot write filtered cbmsim tree")
                _copy_metadata(source, destination, skip=("cbmsim", "replayFilter"))
                previous = source.Get("replayFilter")
                if previous:
                    summary["previous_filter"] = json.loads(str(previous.GetTitle()))
                destination.cd()
                if root.TNamed("replayFilter", json.dumps(summary)).Write() <= 0:
                    raise OSError("Cannot write replayFilter metadata")
                destination.Flush()
                if destination.TestBit(root.TFile.kWriteError):
                    raise OSError("Error writing filtered replay output")
            # Reserve the destination exclusively, then replace only our own
            # placeholder. This also works on filesystems without hard links.
            with output_path.open("xb"):
                pass
            try:
                os.replace(temporary_path, output_path)
            except BaseException:
                output_path.unlink()
                raise
    print(
        f"Processed {summary['events_processed']}; >=2 tracks: {summary['events_two_tracks']}; "
        f"selected: {summary['events_selected']}. Output: {output_path}",
        flush=True,
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--inputfile", required=True, help="Replay ROOT file, after muon-hit restoration")
    parser.add_argument("-o", "--outputfile", required=True, help="New ROOT file; parent directory must exist")
    parser.add_argument("-g", "--geoFile", help="Geometry ROOT file (default: matching geo_ file beside input)")
    parser.add_argument("-s", "--start_event", type=int, default=0, help="First input entry (default: 0)")
    parser.add_argument("-n", "--n_events", type=int, default=-1, help="Input entries to process (-1: all remaining)")
    args = parser.parse_args()
    if args.start_event < 0 or args.n_events < -1:
        parser.error("Counts must be nonnegative, except -n -1 for all remaining entries")

    import ROOT

    ROOT.gROOT.SetBatch(True)
    print("Loading libShipMuDIS.so", flush=True)
    result = ROOT.gSystem.Load("libShipMuDIS.so")
    print(f"libShipMuDIS.so load result: {result}", flush=True)
    if result < 0:
        raise RuntimeError("Cannot load libShipMuDIS.so")
    filter_file(ROOT, args.inputfile, args.outputfile, args.geoFile, args.start_event, args.n_events)


if __name__ == "__main__":
    main()

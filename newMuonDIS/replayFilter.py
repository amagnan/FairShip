#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Skim replay cbmsim events, preserving all branches of each accepted event."""

import argparse
import json
import os
import tempfile
from array import array
from itertools import combinations
from pathlib import Path

from analyseEvents import Geometry, analytical_cross_section, make_vertex, selected_tracks


def clone_filtered_tree(tree, cross_section="stored"):
    """Clone the event schema, optionally replacing the scalar cross section."""
    if cross_section not in ("stored", "analytical"):
        raise ValueError("replayfilter cross_section must be stored or analytical")
    tree.SetBranchStatus("*", True)
    value = None
    if cross_section == "analytical":
        for name in ("muDIS_xsec", "muDIS_pPythia"):
            branch = tree.GetBranch(name)
            leaf = branch.GetLeaf(name) if branch else None
            if (not leaf or leaf.GetTypeName() != "Double_t"
                    or leaf.GetLenStatic() != 1 or leaf.GetLeafCount()):
                raise ValueError(f"Analytical replayfilter requires scalar Double_t {name}")
        tree.SetBranchStatus("muDIS_xsec", False)
        value = array("d", [0.0])
    selected = tree.CloneTree(0)
    if not selected:
        raise OSError("Cannot clone replay tree schema")
    if value is not None:
        selected.Branch("muDIS_xsec", value, "muDIS_xsec/D")
    return selected, value


def replace_cross_section(selected, event, value):
    if value is not None:
        value[0] = analytical_cross_section(event.muDIS_pPythia)
        # TChain may copy input addresses to clones when switching files.
        selected.GetBranch("muDIS_xsec").SetAddress(value)


def selection_stage(event, geometry, pdg):
    """Return 0 (<2 selected tracks), 1 (no accepted pair), or 2 (keep event)."""
    detail = selection_detail(event, geometry, pdg)
    return 0 if detail == 0 else (2 if detail == 5 else 1)


def selection_detail(event, geometry, pdg):
    """Furthest sequential cut passed by any one pair (0 through 5)."""
    tracks, _ = selected_tracks(event, pdg)
    if len(tracks) < 2:
        return 0
    stage = 1
    for first, second in combinations(tracks, 2):
        if first.charge * second.charge >= 0:
            continue
        stage = max(stage, 2)
        vertex = make_vertex(first, second)
        if vertex is None:
            continue
        stage = max(stage, 3)
        if not vertex.doca < 5.0:
            continue
        stage = max(stage, 4)
        if geometry.clearance(vertex.position) >= 0:
            return 5
    return stage


def filter_file(root, input_file, output_file, geometry_file=None, start_event=0, n_events=-1,
                diagnostics=False, cross_section="stored"):
    """Filter one file; the caller must have loaded libShipMuDIS.so once.

    Keep track indices, weights and normalization; optionally replace muDIS_xsec.
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
            "cross_section": cross_section,
            "selection": "(Tr1 or Tr2) + (Tr3 or Tr4) + TD, p > 1 GeV/c; at least one opposite-sign pair with "
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
                selected, xsec_value = clone_filtered_tree(tree, cross_section)
                if not selected:
                    raise OSError("Cannot clone cbmsim tree")
                diagnostic = None
                if diagnostics:
                    from replayFilterDiagnostics import ReplayFilterDiagnostics

                    diagnostic = ReplayFilterDiagnostics(root, summary)
                for entry in range(start_event, stop):
                    if tree.GetEntry(entry) <= 0:
                        raise OSError(f"Cannot read cbmsim entry {entry}")
                    detail = selection_detail(tree, geometry, pdg)
                    stage = 0 if detail == 0 else (2 if detail == 5 else 1)
                    if diagnostic is not None:
                        diagnostic.fill(tree, entry, summary["events_selected"] if stage == 2 else -1, detail)
                    summary["events_two_tracks"] += int(stage >= 1)
                    if stage == 2:
                        replace_cross_section(selected, tree, xsec_value)
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
                _copy_metadata(
                    source,
                    destination,
                    skip=(
                        "cbmsim",
                        "replayFilter",
                        "replayFilterDiagnostics",
                        "replayFilterMuons",
                        "replayFilterPlots",
                    ),
                )
                if diagnostic is not None:
                    destination.cd()
                    diagnostic.write(destination)
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
    parser.add_argument("--cross-section", choices=("stored", "analytical"), default="stored",
                        help="Preserve muDIS_xsec or replace it with sigma(muDIS_pPythia) in mb")
    parser.add_argument(
        "--diagnostics",
        action="store_true",
        help="Save all-event muon/vertex/hit diagnostics and acceptance histograms",
    )
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
    filter_file(ROOT, args.inputfile, args.outputfile, args.geoFile, args.start_event, args.n_events,
                args.diagnostics, args.cross_section)


if __name__ == "__main__":
    main()

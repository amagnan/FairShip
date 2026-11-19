#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Extract histograms and tables from events.root, without simulation files or geometry."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from analyseEvents import (
    ANALYSIS_SHA256,
    MATERIALS,
    OUTPUT_VERSION,
    STEPS,
    PairCutflow,
    ReplayDiagnostics,
    book_plots,
    replay,
    replay_tables,
    save_normalization,
    scaled_rows,
    write_table,
)
from analysisTree import CachedGeometry, event_view, read_records
from renderAnalysis import render_archive


class SavedGeometry(CachedGeometry):
    """Read saved clearances despite roundoff when recomputing pair vertices."""

    # Ten nanometres: an absolute bound, independent of the detector coordinates.
    POSITION_TOLERANCE_CM = 1.0e-6

    def _saved_clearance(self, position, values, name):
        key = tuple(map(float, position))
        if key in values:
            return values[key]
        matches = [
            value
            for saved, value in values.items()
            if all(abs(a - b) <= self.POSITION_TOLERANCE_CM for a, b in zip(key, saved))
        ]
        # A truth point and a reconstructed vertex can coincide. Multiple nearby
        # entries are safe only when they agree on the saved clearance exactly.
        if matches and all(value == matches[0] for value in matches):
            return matches[0]
        if matches:
            raise ValueError(f"Ambiguous saved {name} for {key}")
        raise ValueError(f"Missing saved {name} for {key}")

    def clearance(self, position):
        return self._saved_clearance(position, self.clearances, "clearance")

    def entrance_clearance(self, position):
        return self._saved_clearance(position, self.entrance_clearances, "entrance clearance")


def summarize(source, output, *, save_plots=True):
    source, output = Path(source).expanduser().resolve(), Path(output).expanduser().resolve()
    if output == source or (output.exists() and any(output.iterdir())):
        raise ValueError(f"Output directory must be empty and distinct from input: {output}")
    with (source / "run.json").open() as stream:
        manifest = json.load(stream)
    if manifest.get("output_version") != OUTPUT_VERSION or not manifest["analysis_signature"]["type"].startswith(
        "disreplay"
    ):
        raise ValueError(f"Expected version-{OUTPUT_VERSION} replay event tree")
    if manifest["analysis_signature"]["sha256"] != ANALYSIS_SHA256:
        raise ValueError("Analysis code differs from extraction; re-extract with the current selection")
    if "events.root" in manifest.get("omitted_files", []):
        raise ValueError("Event tree was omitted from this output")
    output.mkdir(parents=True, exist_ok=True)
    config = dict(manifest["config"], tree_only=False, output_dir=str(output), save_plots=save_plots)
    plots = book_plots(config)
    counts = [[0, 0.0, 0.0] for _ in STEPS]
    material_counts = {mat: [[0, 0.0, 0.0] for _ in STEPS] for mat in MATERIALS}
    pairs = PairCutflow()
    diagnostics = ReplayDiagnostics(plots)
    processed = 0
    factor = manifest["analysis_signature"]["weight_factor"]
    for record in read_records(source / "events.root"):
        event = event_view(record)
        charges = {track["pdg_id"]: track["charge_e"] * 3 for track in record["tracks"]}
        pdg = SimpleNamespace(GetParticle=lambda code, charges=charges: SimpleNamespace(Charge=lambda: charges[code]))
        geometry = SavedGeometry(cache=record["geometry_cache"])
        pairs.source = diagnostics.source = (record["source_file"], record["file_uuid"], record["source_entry"])
        cutflow = material_counts.setdefault(record["material"], [[0, 0.0, 0.0] for _ in STEPS])
        replay(
            event,
            plots,
            geometry,
            pdg,
            factor,
            counts,
            cutflow,
            pairs,
            diagnostics,
            cross_section=manifest["analysis_signature"].get("cross_section", "constant"),
        )
        processed += 1
        if processed % 1000 == 0:
            print(f"Summarized {processed}/{manifest['entries_processed']} events", flush=True)
    if processed != manifest["entries_processed"]:
        raise ValueError("Event tree does not match entries_processed")
    scale = manifest["normalization_scale"]
    rows = replay_tables(output, counts, material_counts, pairs, diagnostics, scale)
    write_table(output / "cutflow.csv", scaled_rows(rows, scale))
    save_normalization(output, manifest["normalization"], scale)
    plots.save_archive(output / "histograms.npz")
    if save_plots:
        print("Rendering plots", flush=True)
        render_archive(output / "histograms.npz", output, scale)
    else:
        plots.save_flows(output, scale)
    # Keep the large tree in the source directory; these are derived outputs.
    manifest.update(config=config, tree_only=False, summarized_from=str(source), omitted_files=["events.root"])
    with (output / "run.json").open("w") as stream:
        json.dump(manifest, stream, indent=2)
    print(f"Summarized {processed} events: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Analysis or merged directory containing events.root and run.json")
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument("--no-plots", action="store_true", help="Save tables and histogram arrays without rendering")
    args = parser.parse_args()
    try:
        summarize(args.input, args.output_dir, save_plots=not args.no_plots)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"Summary failed: {error}\n")


if __name__ == "__main__":
    main()

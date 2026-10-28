#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Merge analysis outputs, copying compressed event-tree baskets without decoding the payload."""

import argparse
import csv
import json
import math
from pathlib import Path

from analyseEvents import (
    OUTPUT_VERSION,
    STEPS,
    TYPES,
    PairCutflow,
    Plots,
    ReplayDiagnostics,
    combine_normalization,
    normalization_scale,
    pair_composition_rows,
    save_normalization,
    scaled_rows,
    write_table,
)
from analysisTree import TREE_NAME, validate_tree


def check_sources(manifests):
    """Reject overlapping half-open entry ranges by both path and ROOT UUID."""
    intervals, sources, file_entries, path_uuids = {}, [], {}, {}
    for manifest in manifests:
        processed = 0
        for source in manifest["source_ranges"]:
            start, stop, entries = (source[key] for key in ("start", "stop", "entries"))
            if any(type(value) is not int for value in (start, stop, entries)) or not 0 <= start < stop <= entries:
                raise ValueError(f"Invalid source entry range: {source}")
            path, file_uuid, tree = (source[key] for key in ("path", "file_uuid", "tree"))
            expected_tree = "cbmsim" if manifest["analysis_signature"]["type"].startswith("disreplay") else "MuonDIS"
            if not Path(path).is_absolute() or not file_uuid or tree != expected_tree:
                raise ValueError(f"Invalid source identity: {source}")
            if path_uuids.setdefault(path, file_uuid) != file_uuid:
                raise ValueError(f"Source file UUID changed: {path}")
            for identity in (("path", path, tree), ("uuid", file_uuid, tree)):
                intervals.setdefault(identity, []).append((start, stop))
                previous = file_entries.setdefault(identity, entries)
                if previous != entries:
                    raise ValueError(f"Source file entry count changed: {path}")
            sources.append(source)
            processed += stop - start
        if processed != manifest["entries_processed"]:
            raise ValueError("Source ranges do not match entries_processed")
    for identity, ranges in intervals.items():
        previous_stop = -1
        for start, stop in sorted(ranges):
            if start < previous_stop:
                raise ValueError(f"Overlapping source entries for {identity[1]} ({identity[2]}): [{start}, {stop})")
            previous_stop = stop
    available = sum(entries for identity, entries in file_entries.items() if identity[0] == "uuid")
    return sorted(sources, key=lambda source: (source["path"], source["tree"], source["start"])), available


def merge_table(inputs, filename, keys, integer_columns, float_columns, scaled=False):
    """Sum numeric columns; preserve blanks for undefined step-0 pair counts."""
    expected = set(keys) | set(integer_columns) | set(float_columns)
    derived = {"scaled_sum_weights", "scaled_sum_weights_squared"} if scaled else set()
    combined, fields = {}, None
    for directory in inputs:
        path = directory / filename
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            if (
                reader.fieldnames is None
                or set(reader.fieldnames) != expected | derived
                or len(reader.fieldnames) != len(expected | derived)
            ):
                raise ValueError(f"Unexpected columns: {path}")
            if fields is None:
                fields = [field for field in reader.fieldnames if field not in derived]
            seen = set()
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Malformed table row: {path}")
                for column in derived:
                    del row[column]  # Recompute using the combined preparation denominator.
                key = tuple(row[column] for column in keys)
                if key in seen:
                    raise ValueError(f"Duplicate row {key}: {path}")
                seen.add(key)
                for column in (*integer_columns, *float_columns):
                    if (
                        row[column] == ""
                        and filename == "cutflow.csv"
                        and row["step"] == "0"
                        and column.startswith("pairs")
                    ):
                        continue
                    value = int(row[column]) if column in integer_columns else float(row[column])
                    if not math.isfinite(value) or (column in integer_columns and value < 0):
                        raise ValueError(f"Invalid {column}: {path}")
                    row[column] = value
                if key not in combined:
                    combined[key] = row
                    continue
                for column in (*integer_columns, *float_columns):
                    previous, value = combined[key][column], row[column]
                    if previous == "" or value == "":
                        if previous != value:
                            raise ValueError(f"Inconsistent blank {column}: {path}")
                    else:
                        combined[key][column] += value
    return [combined[key] for key in sorted(combined)], fields


def merge_pairs(inputs, output, scale, expected):
    """Stream individual pairs and recompute their target-scaled weights."""
    fields = [*PairCutflow.detail_fields, *(["scaled_weight"] if scale is not None else [])]
    counts = {}
    with (output / "pairs.csv").open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        for directory in inputs:
            path = directory / "pairs.csv"
            with path.open(newline="") as stream:
                reader = csv.DictReader(stream)
                if reader.fieldnames != fields:
                    raise ValueError(f"Unexpected individual-pair columns: {path}")
                for row in reader:
                    if None in row or any(value is None for value in row.values()):
                        raise ValueError(f"Malformed individual pair: {path}")
                    weight = float(row["weight"])
                    if not math.isfinite(weight):
                        raise ValueError(f"Nonfinite individual-pair weight: {path}")
                    for column in ("source_entry", "track_id_1", "track_id_2"):
                        if int(row[column]) < 0:
                            raise ValueError(f"Invalid pair {column}: {path}")
                    key = (row["material"], int(row["step"]), row["composition"])
                    counts[key] = counts.get(key, 0) + 1
                    if scale is not None:
                        row["scaled_weight"] = weight * scale
                    writer.writerow(row)
    expected_counts = {
        (row["material"], row["step"], row["composition"]): row["pairs"]
        for row in expected
        if row["material"] != "ALL" and row["pairs"]
    }
    if counts != expected_counts:
        raise ValueError("Individual pair counts do not match the pair cutflow")


def merge_tagger_tracks(inputs, output, scale):
    """Keep each event/track's hits and stage flags; recompute scaled weights."""
    fields = [
        *ReplayDiagnostics.detail_fields,
        *(ReplayDiagnostics.scaled_detail_fields if scale is not None else []),
    ]
    with (output / "tagger_tracks.csv").open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        for directory in inputs:
            path = directory / "tagger_tracks.csv"
            with path.open(newline="") as stream:
                reader = csv.DictReader(stream)
                if reader.fieldnames != fields:
                    raise ValueError(f"Unexpected tagger track columns: {path}")
                for row in reader:
                    if None in row or any(value is None for value in row.values()):
                        raise ValueError(f"Malformed tagger track: {path}")
                    for column in ("source_entry", "track_id"):
                        if int(row[column]) < 0:
                            raise ValueError(f"Invalid tagger track {column}: {path}")
                    for column in ReplayDiagnostics.weight_fields:
                        if column != "weight":
                            step = column.removeprefix("weight_")
                            passed = row[f"passed_{step}"]
                            if passed not in ("True", "False") or (passed == "True") != (row[column] != ""):
                                raise ValueError(f"Invalid tagger track selection/weight for {step}: {path}")
                        if row[column] == "" and column != "weight":
                            if scale is not None:
                                row[f"scaled_{column}"] = ""
                            continue
                        weight = float(row[column])
                        if not math.isfinite(weight):
                            raise ValueError(f"Nonfinite tagger track {column}: {path}")
                        if scale is not None:
                            row[f"scaled_{column}"] = weight * scale
                    writer.writerow(row)


def merge_events(inputs, output, scale, manifests, expected=None):
    """Validate small provenance branches, then copy all compressed baskets in C++."""
    import ROOT

    counts = {}
    for directory, manifest in zip(inputs, manifests):
        path = directory / "events.root"
        if not path.is_file():
            raise FileNotFoundError(path)
        file = ROOT.TFile.Open(str(path))
        if not file or file.IsZombie():
            raise ValueError(f"Cannot open {path}")
        try:
            tree = file.Get(TREE_NAME)
            if not tree or tree.GetEntries() != manifest["entries_processed"]:
                raise ValueError(f"Event tree does not match entries_processed: {directory}")
            validate_tree(tree, path)
            ranges = {}
            for source in manifest["source_ranges"]:
                key = (source["path"], source["file_uuid"])
                ranges.setdefault(key, []).append((source["start"], source["stop"]))
            previous = {}
            tree.SetBranchStatus("*", 0)
            for name in ("source_file", "file_uuid", "source_entry", "material", "passed_steps", "stage_weights"):
                tree.SetBranchStatus(name, 1)
            for entry in range(tree.GetEntries()):
                if tree.GetEntry(entry) <= 0:
                    raise ValueError(f"Cannot read provenance: {path}, {entry}")
                key = (str(tree.source_file), str(tree.file_uuid))
                source_entry = int(tree.source_entry)
                interval = next(
                    ((start, stop) for start, stop in ranges.get(key, []) if start <= source_entry < stop), None
                )
                if interval is None:
                    raise ValueError(f"Event outside declared source ranges: {directory}, {source_entry}")
                range_key = (key, interval)
                if source_entry <= previous.get(range_key, -1):
                    raise ValueError(f"Duplicate or unordered event: {directory}, {source_entry}")
                previous[range_key] = source_entry
                flags, weights = list(tree.passed_steps), list(tree.stage_weights)
                if (
                    len(flags) != len(STEPS)
                    or len(weights) != len(STEPS)
                    or not flags[0]
                    or any(flag not in (0, 1) for flag in flags)
                    or any(flags[i] and not flags[i - 1] for i in range(1, len(flags)))
                ):
                    raise ValueError(f"Invalid event stages: {directory}, {source_entry}")
                for step, (flag, weight) in enumerate(zip(flags, weights)):
                    if (flag and not math.isfinite(weight)) or (not flag and not math.isnan(weight)):
                        raise ValueError(f"Invalid event stage weight: {directory}, {source_entry}")
                    if flag:
                        for material in ("ALL", str(tree.material)):
                            index = (material, step)
                            counts[index] = counts.get(index, 0) + 1
        finally:
            file.Close()
    if expected is not None:
        expected_counts = {(row["material"], int(row["step"])): row["events"] for row in expected if row["events"]}
        if counts != expected_counts:
            raise ValueError("Event records do not match the event cutflow")
    merger = ROOT.TFileMerger(False)
    merger.SetFastMethod(True)
    if not merger.OutputFile(str(output / "events.root"), "RECREATE", 101):
        raise OSError("Cannot create merged event tree")
    for directory in inputs:
        if not merger.AddFile(str(directory / "events.root")):
            raise OSError(f"Cannot add event tree: {directory}")
    print("Copying compressed DIS tree baskets", flush=True)
    if not merger.Merge():
        raise OSError("ROOT event-tree merge failed")
    # Weights deliberately remain unscaled; the combined scale belongs to run.json.


def merge(inputs, output, *, save_plots=True, skip_event_records=False, skip_tagger_tracks=False):
    if skip_tagger_tracks:
        raise ValueError("Output version 7 uses events.root; use --skip-event-records instead of --skip-tagger-tracks")
    inputs = [Path(directory).expanduser().resolve() for directory in inputs]
    output = Path(output).expanduser().resolve()
    if not inputs:
        raise ValueError("At least one input directory is required")
    if len(set(inputs)) != len(inputs):
        raise ValueError("An input directory was listed more than once")
    if output in inputs or (output.exists() and any(output.iterdir())):
        raise ValueError(f"Output directory must be empty and distinct from the inputs: {output}")
    manifests = []
    for directory in inputs:
        with (directory / "run.json").open() as stream:
            manifest = json.load(stream)
        if manifest.get("output_version") != OUTPUT_VERSION:
            raise ValueError(
                f"Missing or unsupported output version in {directory}; rerun the analysis with archive output"
            )
        if "events.root" in manifest.get("omitted_files", []) and not skip_event_records:
            raise ValueError(f"Event records were omitted from {directory}; use --skip-event-records")
        manifests.append(manifest)
    reference = manifests[0]
    signature = reference["analysis_signature"]
    tree_only = reference.get("tree_only", False)
    if any(manifest.get("tree_only", False) != tree_only for manifest in manifests):
        raise ValueError("Cannot mix tree-only and histogram outputs; summarize them first")
    if tree_only and skip_event_records:
        raise ValueError("Cannot omit event records from a tree-only merge")
    if signature["type"] not in TYPES:
        raise ValueError(f"Unsupported analysis type: {signature['type']}")
    for directory, manifest in zip(inputs, manifests):
        if (
            manifest["analysis_signature"] != signature
            or manifest["weight_definition"] != reference["weight_definition"]
        ):
            raise ValueError(f"Incompatible analysis code, type, weight factor, or geometry settings: {directory}")
    sources, available = check_sources(manifests)
    normalization = combine_normalization([manifest["normalization"] for manifest in manifests])
    normalization["weight_factor"] = signature["weight_factor"]
    scale = normalization_scale(signature, normalization)
    merged_tables = []
    if not tree_only:
        plots = Plots.load_archive(inputs[0] / "histograms.npz")
        for directory in inputs[1:]:
            plots.merge(Plots.load_archive(directory / "histograms.npz"))

        weights = ("sum_weights", "sum_weights_squared")
        if signature["type"].startswith("disreplay"):
            tables = [
                (
                    "cutflow.csv",
                    ("material", "step", "selection"),
                    ("events", "pairs", "pairs_DIS_DIS", "pairs_DIS_Soft", "pairs_Soft_Soft"),
                    weights,
                ),
                ("pair_cutflow.csv", PairCutflow.fields[:8], ("pairs",), weights),
                ("tagger_cutflow.csv", ReplayDiagnostics.fields[:6], ("events",), weights),
            ]
        else:
            tables = [("vertex_counts.csv", ("material",), ("vertices",), weights)]
        merged_tables = [
            (filename, *merge_table(inputs, filename, keys, integers, floats, scaled=scale is not None))
            for filename, keys, integers, floats in tables
        ]

    # Metadata and summaries are validated first; event records are checked
    # while streaming, before the completion manifest is written.
    output.mkdir(parents=True, exist_ok=True)
    for filename, rows, fields in merged_tables:
        if scale is not None:
            fields += ["scaled_sum_weights", "scaled_sum_weights_squared"]
        write_table(output / filename, scaled_rows(rows, scale), fields)
    if tree_only:
        merge_events(inputs, output, scale, manifests)
    elif signature["type"].startswith("disreplay"):
        pair_rows = next(rows for filename, rows, _ in merged_tables if filename == "pair_cutflow.csv")
        compositions = pair_composition_rows(pair_rows)
        write_table(output / "pair_composition_cutflow.csv", scaled_rows(compositions, scale))
        if not skip_event_records:
            event_rows = next(rows for filename, rows, _ in merged_tables if filename == "cutflow.csv")
            merge_events(inputs, output, scale, manifests, event_rows)
    save_normalization(output, normalization, scale)
    if not tree_only:
        plots.save_archive(output / "histograms.npz")
        if save_plots:
            print("Rendering merged plots", flush=True)
            plots.save(output, scale)
        else:
            plots.save_flows(output, scale)
    config = {
        key: reference["config"][key]
        for key in ("type", "weight_factor", "plot_axes", "geometry_file", "target_sum_weights", "failed_job_factor")
        if key in reference["config"]
    }
    config.update(output_dir=str(output), save_plots=save_plots, tree_only=tree_only)
    manifest = {
        "output_version": OUTPUT_VERSION,
        "tree_only": tree_only,
        "analysis_signature": signature,
        "config": config,
        "files": sorted({source["path"] for source in sources}),
        "source_ranges": sources,
        "normalization": normalization,
        "normalization_scale": scale,
        "entries_available": available,
        "entries_processed": sum(item["entries_processed"] for item in manifests),
        "weight_definition": reference["weight_definition"],
        "merged_from": [str(directory) for directory in inputs],
        "event_record_schema": reference.get("event_record_schema"),
    }
    if skip_event_records and signature["type"].startswith("disreplay"):
        manifest["omitted_files"] = ["events.root"]
    with (output / "run.json").open("w") as stream:
        json.dump(manifest, stream, indent=2)
    print(f"Merged {len(inputs)} outputs, {manifest['entries_processed']} entries: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Analysis or previously merged output directories")
    parser.add_argument("-o", "--output-dir", required=True, type=Path, help="New or empty output directory")
    parser.add_argument("--no-plots", action="store_true", help="Save merged data without rendering plots")
    parser.add_argument(
        "--skip-event-records",
        action="store_true",
        help="Omit events.root; merge histograms and aggregate cutflows only",
    )
    args = parser.parse_args()
    try:
        merge(args.inputs, args.output_dir, save_plots=not args.no_plots, skip_event_records=args.skip_event_records)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"Merge failed: {error}\n")


if __name__ == "__main__":
    main()

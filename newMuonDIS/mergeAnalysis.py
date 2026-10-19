#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Merge analyseEvents output directories without reading ROOT files."""

import argparse
import csv
import json
import math
from pathlib import Path

from analyseEvents import (
    OUTPUT_VERSION,
    TYPES,
    PairCutflow,
    Plots,
    combine_normalization,
    normalization_scale,
    pair_composition_rows,
    save_normalization,
    scaled_rows,
    write_table,
)


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


def merge(inputs, output, *, save_plots=True):
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
        manifests.append(manifest)
    reference = manifests[0]
    signature = reference["analysis_signature"]
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
        ]
    else:
        tables = [("vertex_counts.csv", ("material",), ("vertices",), weights)]
    merged_tables = [
        (filename, *merge_table(inputs, filename, keys, integers, floats, scaled=scale is not None))
        for filename, keys, integers, floats in tables
    ]

    # Metadata and summaries are validated first; individual pairs are checked
    # while streaming, before the completion manifest is written.
    output.mkdir(parents=True, exist_ok=True)
    for filename, rows, fields in merged_tables:
        if scale is not None:
            fields += ["scaled_sum_weights", "scaled_sum_weights_squared"]
        write_table(output / filename, scaled_rows(rows, scale), fields)
    if signature["type"].startswith("disreplay"):
        pair_rows = next(rows for filename, rows, _ in merged_tables if filename == "pair_cutflow.csv")
        compositions = pair_composition_rows(pair_rows)
        write_table(output / "pair_composition_cutflow.csv", scaled_rows(compositions, scale))
        merge_pairs(inputs, output, scale, compositions)
    save_normalization(output, normalization, scale)
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
    config.update(output_dir=str(output), save_plots=save_plots)
    manifest = {
        "output_version": OUTPUT_VERSION,
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
    }
    with (output / "run.json").open("w") as stream:
        json.dump(manifest, stream, indent=2)
    print(f"Merged {len(inputs)} outputs, {manifest['entries_processed']} entries: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Analysis or previously merged output directories")
    parser.add_argument("-o", "--output-dir", required=True, type=Path, help="New or empty output directory")
    parser.add_argument("--no-plots", action="store_true", help="Save merged data without rendering plots")
    args = parser.parse_args()
    try:
        merge(args.inputs, args.output_dir, save_plots=not args.no_plots)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"Merge failed: {error}\n")


if __name__ == "__main__":
    main()

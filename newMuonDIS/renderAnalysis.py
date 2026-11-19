#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Render saved histogram arrays without reading event trees or rerunning analysis."""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def read_histograms(path):
    """Read the common version-7/8 histogram schema, one histogram at a time."""
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata"]))
        if metadata["version"] not in (7, 8):
            raise ValueError(f"Unsupported histogram archive version: {path}")
        seen = set()
        for item in metadata["histograms"]:
            name, prefix, dimensions = item["name"], item["prefix"], item["dimensions"]
            if not name or name in seen or dimensions not in (1, 2):
                raise ValueError(f"Invalid histogram definition: {name}")
            if Path(name).is_absolute() or ".." in Path(name).parts or name == ".":
                raise ValueError(f"Invalid histogram name: {name}")
            seen.add(name)
            edges = tuple(archive[f"{prefix}_edge{axis}"] for axis in range(dimensions))
            if any(
                edge.ndim != 1 or len(edge) < 2 or not np.all(np.isfinite(edge)) or not np.all(np.diff(edge) > 0)
                for edge in edges
            ):
                raise ValueError(f"Invalid histogram edges: {name}")
            raw, weighted = archive[f"{prefix}_raw"], archive[f"{prefix}_weighted"]
            if raw.shape != tuple(len(edge) - 1 for edge in edges) or weighted.shape != raw.shape:
                raise ValueError(f"Invalid histogram shape: {name}")
            flow, wflow = float(archive[f"{prefix}_flow"]), float(archive[f"{prefix}_wflow"])
            if (
                not np.all(np.isfinite(raw))
                or not np.all(np.isfinite(weighted))
                or np.any(raw < 0)
                or np.any(raw != np.floor(raw))
                or not math.isfinite(flow)
                or flow < 0
                or flow != math.floor(flow)
                or not math.isfinite(wflow)
            ):
                raise ValueError(f"Invalid histogram contents: {name}")
            yield name, (edges, raw, weighted, int(flow), wflow)


def render_archive(source, output, scale=None):
    """Write PDF/PNG plots and flow totals; callers manage the output directory."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    if scale is not None and (isinstance(scale, bool) or not math.isfinite(scale)):
        raise ValueError("Normalization scale must be finite")
    output = Path(output)
    rows = []
    count = 0
    for name, (edges, raw, weighted, flow, wflow) in read_histograms(source):
        samples = [("unweighted", raw, flow), ("weighted", weighted, wflow)]
        if scale is not None:
            samples.append(("scaled", weighted * scale, wflow * scale))
        for label, data, outside in samples:
            destination = output / "plots" / name / label
            destination.parent.mkdir(parents=True, exist_ok=True)
            fig, ax = plt.subplots(figsize=(8, 6))
            try:
                mesh = None
                if len(edges) == 2:
                    positive = data[data > 0]
                    vmin = float(positive.min()) if positive.size else 1.0
                    vmax = float(positive.max()) if positive.size else 10.0
                    if vmax <= vmin:
                        vmax = 10.0 * vmin
                    mesh = ax.pcolormesh(
                        edges[0],
                        edges[1],
                        np.ma.masked_less_equal(data.T, 0),
                        shading="auto",
                        norm=LogNorm(vmin=vmin, vmax=vmax),
                    )
                    fig.colorbar(mesh, ax=ax, label="Entries" if label == "unweighted" else "Sum of weights")
                    vertical, horizontal = name.split("/")[-1].split("_vs_")
                    ax.set_xlabel(f"{horizontal} [cm]")
                    ax.set_ylabel("Number of hits" if vertical == "n_hits" else f"{vertical} [cm]")
                else:
                    ax.stairs(data, edges[0])
                    ax.set_xlabel(name.split("/")[-1])
                    ax.set_ylabel("Entries" if label == "unweighted" else "Sum of weights")
                ax.set_title(f"{name} ({label})\nOutside plotted range: {outside:g}")
                fig.tight_layout()
                if mesh is not None:
                    mesh.set_rasterized(True)
                fig.savefig(destination.with_suffix(".pdf"), dpi=300)
                if mesh is not None:
                    mesh.set_rasterized(False)
                fig.savefig(destination.with_suffix(".png"))
            finally:
                plt.close(fig)
        row = {"histogram": name, "outside_range": flow, "weighted_outside_range": wflow}
        if scale is not None:
            row["scaled_weighted_outside_range"] = wflow * scale
        rows.append(row)
        count += 1
        if count % 25 == 0:
            print(f"Rendered {count} histograms", flush=True)
    if rows:
        with (output / "histogram_flows.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(f"Rendered {count} histograms: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="histograms.npz or its containing summary directory")
    parser.add_argument("-o", "--output-dir", required=True, type=Path, help="New or empty plot output directory")
    scaling = parser.add_mutually_exclusive_group()
    scaling.add_argument("--scale", type=float, help="Normalization scale; defaults to the adjacent run.json")
    scaling.add_argument("--unscaled", action="store_true", help="Render only unweighted and weighted plots")
    args = parser.parse_args()
    try:
        source = args.input.expanduser().resolve()
        if source.is_dir():
            source /= "histograms.npz"
        if not source.is_file():
            raise ValueError(f"Histogram archive not found: {source}")
        output = args.output_dir.expanduser().resolve()
        if output == source.parent or (output.exists() and (not output.is_dir() or any(output.iterdir()))):
            raise ValueError(f"Output directory must be empty and distinct from input: {output}")
        scale = args.scale
        if not args.unscaled and scale is None:
            manifest = source.with_name("run.json")
            if not manifest.is_file():
                raise ValueError("No adjacent run.json; specify --scale or --unscaled")
            with manifest.open() as stream:
                scale = json.load(stream)["normalization_scale"]
        if scale is not None and (isinstance(scale, bool) or not math.isfinite(scale)):
            raise ValueError("Normalization scale must be finite")
        output.mkdir(parents=True, exist_ok=True)
        print(f"Rendering {source} (normalization scale: {scale})", flush=True)
        render_archive(source, output, scale)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f"Rendering failed: {error}\n")


if __name__ == "__main__":
    main()

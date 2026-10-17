#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Plot prepared DIS vertices or analyse replay truth tracks; see analyseEvents.yaml."""

import argparse
import csv
import hashlib
import json
import math
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import yaml

MATERIALS = ("MS", "UBT", "SBTsens", "SBTfr", "SSTsens", "SSTfr", "HE", "AIR", "HeBalloon", "REST")
TYPES = ("disgen", "disfilter", "disreplay_gen", "disreplay_filter")
OUTPUT_VERSION = 3
ANALYSIS_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
STEPS = (
    "all",
    "at least two tracks: Tr1-4 + TD, p > 1 GeV/c",
    "opposite sign, backward vertex, DOCA < 1 cm",
    "inside decay helium, wall clearance >= 5 cm",
    "target-plane IP < 250 cm",
    "target-plane IP < 10 cm",
)


def combine_normalization(inputs):
    records = {}
    complete = bool(inputs)
    for metadata in inputs:
        complete = complete and metadata["preparation_complete"]
        for record in metadata["records"]:
            identity = record["id"]
            if identity in records and records[identity] != record:
                raise ValueError(f"Conflicting normalization record: {identity}")
            records[identity] = record
    return {"preparation_complete": complete, "records": [records[key] for key in sorted(records)]}


def read_normalization(root, files, tree_name):
    inputs = []
    for path in files:
        source = root.TFile.Open(str(path), "READ")
        if not source or source.IsZombie():
            raise ValueError(f"Cannot read normalization metadata: {path}")
        try:
            tree = source.Get(tree_name)
            if not tree:
                raise ValueError(f"Missing {tree_name} tree: {path}")
            metadata = tree.GetUserInfo().FindObject("muDIS_normalization")
            records = []
            complete = False
            if metadata:
                version = metadata.FindObject("version")
                if not version or int(version.GetVal()) != 1:
                    raise ValueError(f"Unsupported normalization metadata: {path}")
                flag = metadata.FindObject("preparation_complete")
                complete = bool(flag and flag.GetVal())
                for item in metadata:
                    if not item.InheritsFrom("TList"):
                        continue
                    stage = item.FindObject("stage")
                    record = {"id": str(item.GetName()), "stage": str(stage.GetTitle()) if stage else ""}
                    for name in ("input_muons", "sum_weights", "sum_weights_squared"):
                        value = item.FindObject(name)
                        if not value:
                            raise ValueError(f"Missing normalization field {name}: {path}")
                        record[name] = int(value.GetVal()) if name == "input_muons" else float(value.GetVal())
                    if (
                        record["stage"] not in ("prepareEvents", "filterEvents", "NewMuDISGenerator")
                        or record["input_muons"] < 0
                        or record["sum_weights_squared"] < 0
                        or not all(math.isfinite(record[key]) for key in ("sum_weights", "sum_weights_squared"))
                    ):
                        raise ValueError(f"Invalid normalization record: {path}")
                    records.append(record)
            inputs.append({"preparation_complete": complete, "records": records})
        finally:
            source.Close()
    return combine_normalization(inputs)


def normalization_scale(config, metadata):
    target = config.get("target_sum_weights")
    failed = float(config.get("failed_job_factor", 1.0))
    if not math.isfinite(failed) or failed <= 0:
        raise ValueError("failed_job_factor must be finite and positive")
    if target is None:
        return None
    target = float(target)
    if not math.isfinite(target) or target <= 0:
        raise ValueError("target_sum_weights must be finite and positive")
    if not metadata["preparation_complete"]:
        raise ValueError("Target normalization requires complete prepareEvents metadata in every input file")
    denominator = math.fsum(
        record["sum_weights"] for record in metadata["records"] if record["stage"] == "prepareEvents"
    )
    if denominator <= 0 or not math.isfinite(denominator):
        raise ValueError("Preparation sum_weights must be finite and positive for target normalization")
    return target / denominator * failed


def scaled_rows(rows, scale):
    if scale is None:
        return rows
    return [
        dict(
            row,
            scaled_sum_weights=row["sum_weights"] * scale,
            scaled_sum_weights_squared=row["sum_weights_squared"] * scale * scale,
        )
        for row in rows
    ]


def save_normalization(output, metadata, scale):
    rows = []
    for stage in ("prepareEvents", "filterEvents", "NewMuDISGenerator"):
        records = [record for record in metadata["records"] if record["stage"] == stage]
        if records:
            rows.append(
                {
                    "stage": stage,
                    "input_muons": sum(record["input_muons"] for record in records),
                    "sum_weights": math.fsum(record["sum_weights"] for record in records),
                    "sum_weights_squared": math.fsum(record["sum_weights_squared"] for record in records),
                }
            )
    fields = ["stage", "input_muons", "sum_weights", "sum_weights_squared"]
    if scale is not None:
        fields += ["scaled_sum_weights", "scaled_sum_weights_squared"]
    write_table(output / "normalization.csv", scaled_rows(rows, scale), fields)


def discover(path, kind):
    """Select event files only, including the legacy unsuffixed generation name."""
    if kind not in TYPES:
        raise ValueError(f"Unknown type {kind!r}; choose from {TYPES}")
    path = Path(path).expanduser()
    if not path.is_dir():
        raise ValueError(f"Production path is not a directory: {path}")
    result = []
    for file in path.rglob("*.root"):
        name = file.stem
        if kind.startswith("disreplay"):
            match = name.startswith("sim_") and name.endswith("_" + kind)
        elif kind == "disfilter":
            match = name.startswith("muonDis_") and name.endswith("_disfilter")
        else:
            match = name.startswith("muonDis_") and (
                name.endswith("_disgen") or not any(token in name for token in ("disfilter", "disreplay"))
            )
        if match:
            result.append(file)
    if not result:
        raise ValueError(f"No {kind} event files found under {path}")
    return sorted(result)


@dataclass
class Track:
    position: np.ndarray
    momentum: np.ndarray
    charge: float
    weight: float = 1.0
    pdg_id: int = 0
    track_id: int = -1


@dataclass
class Vertex:
    position: np.ndarray
    momentum: np.ndarray
    doca: float
    opposite: bool
    weight: float = 1.0
    tracks: tuple = ()


def shared_weight(weights, muon_weight, *, event=False):
    """Count shared factors once, for a pair (AND) or surviving pairs (OR).

    Soft tracks/pairs have the original muon weight. All other tracks/pairs
    must share one DIS weight. For an OR, a soft pair already suffices to pass.
    A relative tolerance accommodates fW's Double32_t storage.
    """
    weights = list(weights)
    if not math.isfinite(muon_weight) or not all(math.isfinite(w) for w in weights):
        raise ValueError("Nonfinite replay weight")
    soft = [math.isclose(w, muon_weight, rel_tol=1.0e-6, abs_tol=0.0) for w in weights]
    dis = [w for w, is_soft in zip(weights, soft) if not is_soft]
    if dis and not all(math.isclose(w, dis[0], rel_tol=1.0e-6, abs_tol=0.0) for w in dis):
        raise ValueError(f"Inconsistent DIS track/pair weights: {dis}; muon weight: {muon_weight}")
    if not dis or (event and any(soft)):
        return muon_weight
    return dis[0]


def make_vertex(a, b):
    """Closest points of two lines; accept only backward extrapolations.

    Parallel lines have no unique vertex and are rejected. The vertex is the
    midpoint of the closest points, with the sum of momenta measured at Tr1.
    """
    u = a.momentum / np.linalg.norm(a.momentum)
    v = b.momentum / np.linalg.norm(b.momentum)
    w = a.position - b.position
    uv = float(u @ v)
    denominator = 1.0 - uv * uv
    if denominator <= 1.0e-12:
        return None
    s = (uv * float(v @ w) - float(u @ w)) / denominator
    t = (float(v @ w) - uv * float(u @ w)) / denominator
    if s > 0 or t > 0:
        return None
    pa, pb = a.position + s * u, b.position + t * v
    return Vertex((pa + pb) / 2.0, a.momentum + b.momentum, float(np.linalg.norm(pa - pb)), a.charge * b.charge < 0)


def target_ip(vertex):
    p = vertex.momentum
    if abs(p[2]) < 1.0e-12:
        return math.inf
    at_target = vertex.position - vertex.position[2] * p / p[2]
    return float(np.linalg.norm(at_target[:2]))


def selected_tracks(event, pdg):
    """Associate hits by MCTrack index; use the smallest-z Tr1 hit per track."""
    stations, first = {}, {}
    count = len(event.MCTrack)
    for hit in event.strawtubesPoint:
        index, station = int(hit.GetTrackID()), int(hit.GetDetectorID()) // 1000000
        if not 0 <= index < count or station not in (1, 2, 3, 4):
            continue
        stations.setdefault(index, set()).add(station)
        if station == 1 and (index not in first or hit.GetZ() < first[index][0][2]):
            first[index] = (
                np.array([hit.GetX(), hit.GetY(), hit.GetZ()]),
                np.array([hit.GetPx(), hit.GetPy(), hit.GetPz()]),
            )
    timing = {int(hit.GetTrackID()) for hit in event.TimeDetPoint}
    tracks, momenta = [], []
    for index, seen in stations.items():
        if seen != {1, 2, 3, 4} or index not in timing:
            continue
        position, momentum = first[index]
        if not np.all(np.isfinite(position)) or not np.all(np.isfinite(momentum)):
            raise ValueError("Nonfinite Tr1 measurement")
        p = float(np.linalg.norm(momentum))
        weight = float(event.MCTrack[index].GetWeight())
        if not math.isfinite(weight):
            raise ValueError(f"Nonfinite MCTrack weight for track {index}")
        momenta.append((p, weight))
        if p <= 1.0:
            continue
        code = int(event.MCTrack[index].GetPdgCode())
        particle = pdg.GetParticle(code)
        if particle:
            charge = float(particle.Charge())
        elif 1000000000 <= abs(code) < 1100000000:
            # Nuclear codes are +/-10LZZZAAAI; ROOT charges are in units of e/3.
            z = (abs(code) // 10000) % 1000
            charge = 3.0 * z * (1 if code > 0 else -1)
        else:
            raise ValueError(f"Unknown PDG code: {code}")
        tracks.append(Track(position, momentum, charge, weight, code, index))
    return tracks, momenta


class Geometry:
    def __init__(self, root, override=None):
        self.root, self.override, self.path, self.manager = root, override, None, None

    def load(self, event_file):
        path = (
            Path(self.override).expanduser()
            if self.override
            else event_file.with_name("geo_" + event_file.name.removeprefix("sim_"))
        )
        if path == self.path:
            return
        if not path.is_file():
            raise ValueError(f"Missing replay geometry: {path}; set geometry_file in YAML")
        self.manager = self.root.TGeoManager.Import(str(path))
        if not self.manager:
            raise ValueError(f"Cannot import geometry: {path}")
        self.path = path

    def clearance(self, position):
        node = self.manager.FindNode(*map(float, position))
        if not node:
            return -1.0
        volume = node.GetVolume()
        material = volume.GetMaterial()
        if (
            not material
            or str(material.GetName()).lower() != "helium"
            or "DecayVolume" not in str(self.manager.GetPath())
        ):
            return -1.0
        # ROOT navigator safety is the distance to the nearest boundary in
        # any direction, including endcaps and any daughter-volume surfaces.
        return float(self.manager.GetCurrentNavigator().Safety())


class Plots:
    """Accumulate bounded-memory histograms, retaining under/overflow counts."""

    def __init__(self, config):
        self.histograms = {}
        self.axes = {
            axis: np.linspace(*config.get(axis, default))
            for axis, default in {"x": [-600, 600, 121], "y": [-600, 600, 121], "z": [0, 14000, 141]}.items()
        }

    def book(self, name, edges):
        edges = tuple(np.asarray(edge) for edge in edges)
        if name not in self.histograms:
            shape = tuple(len(edge) - 1 for edge in edges)
            self.histograms[name] = [edges, np.zeros(shape), np.zeros(shape), 0, 0.0]

    def fill(self, name, values, edges, weight):
        self.book(name, edges)
        stored, raw, weighted, _, _ = self.histograms[name]
        indices = tuple(int(np.searchsorted(edge, value, side="right") - 1) for edge, value in zip(stored, values))
        if any(i < 0 or i >= len(edge) - 1 for i, edge in zip(indices, stored)):
            self.histograms[name][3] += 1
            self.histograms[name][4] += weight
        else:
            raw[indices] += 1
            weighted[indices] += weight

    def one(self, name, value, low, high, weight, bins=100):
        self.fill(name, (value,), (np.linspace(low, high, bins + 1),), weight)

    def vertex(self, prefix, position, weight):
        coordinates = dict(zip("xyz", position))
        for vertical, horizontal in [("x", "z"), ("y", "z"), ("y", "x")]:
            self.fill(
                f"{prefix}/{vertical}_vs_{horizontal}",
                (coordinates[horizontal], coordinates[vertical]),
                (self.axes[horizontal], self.axes[vertical]),
                weight,
            )

    def save_archive(self, path):
        arrays, histograms = {}, []
        for i, (name, (edges, raw, weighted, flow, wflow)) in enumerate(self.histograms.items()):
            prefix = f"h{i}"
            histograms.append({"name": name, "prefix": prefix, "dimensions": len(edges)})
            for axis, edge in enumerate(edges):
                arrays[f"{prefix}_edge{axis}"] = edge
            arrays[f"{prefix}_raw"] = raw
            arrays[f"{prefix}_weighted"] = weighted
            arrays[f"{prefix}_flow"] = np.array(flow, dtype=np.int64)
            arrays[f"{prefix}_wflow"] = np.array(wflow)
        arrays["metadata"] = np.array(json.dumps({"version": OUTPUT_VERSION, "histograms": histograms}))
        np.savez_compressed(path, **arrays)

    @classmethod
    def load_archive(cls, path):
        plots = cls({})
        with np.load(path, allow_pickle=False) as archive:
            metadata = json.loads(str(archive["metadata"]))
            if metadata["version"] != OUTPUT_VERSION:
                raise ValueError(f"Unsupported histogram archive version: {path}")
            for item in metadata["histograms"]:
                name, prefix, dimensions = item["name"], item["prefix"], item["dimensions"]
                if name in plots.histograms or dimensions not in (1, 2):
                    raise ValueError(f"Invalid histogram definition: {name}")
                if Path(name).is_absolute() or ".." in Path(name).parts:
                    raise ValueError(f"Invalid histogram name: {name}")
                edges = tuple(archive[f"{prefix}_edge{axis}"] for axis in range(dimensions))
                if any(
                    edge.ndim != 1 or len(edge) < 2 or not np.all(np.isfinite(edge)) or not np.all(np.diff(edge) > 0)
                    for edge in edges
                ):
                    raise ValueError(f"Invalid histogram edges: {name}")
                plots.book(name, edges)
                raw, weighted = archive[f"{prefix}_raw"], archive[f"{prefix}_weighted"]
                if raw.shape != plots.histograms[name][1].shape or weighted.shape != raw.shape:
                    raise ValueError(f"Invalid histogram shape: {name}")
                flow, wflow = int(archive[f"{prefix}_flow"]), float(archive[f"{prefix}_wflow"])
                if (
                    not np.all(np.isfinite(raw))
                    or not np.all(np.isfinite(weighted))
                    or np.any(raw < 0)
                    or np.any(raw != np.floor(raw))
                    or flow < 0
                    or not math.isfinite(wflow)
                ):
                    raise ValueError(f"Invalid histogram contents: {name}")
                plots.histograms[name] = [edges, raw.astype(float), weighted.astype(float), flow, wflow]
        return plots

    def merge(self, other):
        if self.histograms.keys() != other.histograms.keys():
            raise ValueError("Histogram names differ between jobs")
        for name, (edges, raw, weighted, flow, wflow) in other.histograms.items():
            stored = self.histograms[name]
            if len(edges) != len(stored[0]) or any(not np.array_equal(a, b) for a, b in zip(edges, stored[0])):
                raise ValueError(f"Histogram bin edges differ: {name}")
            stored[1] += raw
            stored[2] += weighted
            stored[3] += flow
            stored[4] += wflow

    def save_flows(self, output, scale=None):
        rows = [
            {"histogram": name, "outside_range": values[3], "weighted_outside_range": values[4]}
            for name, values in self.histograms.items()
        ]
        if scale is not None:
            for row in rows:
                row["scaled_weighted_outside_range"] = row["weighted_outside_range"] * scale
        write_table(output / "histogram_flows.csv", rows)

    def save(self, output, scale=None):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        for name, (edges, raw, weighted, flow, wflow) in self.histograms.items():
            samples = [("unweighted", raw, flow), ("weighted", weighted, wflow)]
            if scale is not None:
                samples.append(("scaled", weighted * scale, wflow * scale))
            for label, data, outside in samples:
                destination = output / "plots" / name / label
                destination.parent.mkdir(parents=True, exist_ok=True)
                fig, ax = plt.subplots(figsize=(8, 6))
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
                    ax.set_ylabel(f"{vertical} [cm]")
                else:
                    ax.stairs(data, edges[0])
                    ax.set_xlabel(name.split("/")[-1])
                    ax.set_ylabel("Entries" if label == "unweighted" else "Sum of weights")
                ax.set_title(f"{name} ({label})\nOutside plotted range: {outside:g}")
                fig.tight_layout()
                for extension in ("pdf", "png"):
                    fig.savefig(destination.with_suffix("." + extension))
                plt.close(fig)
        self.save_flows(output, scale)


def write_table(path, rows, fieldnames=None):
    if not rows and fieldnames is None:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class PairCutflow:
    """Count unordered pairs, retaining the PDG ID associated with each origin."""

    fields = (
        "material",
        "step",
        "selection",
        "composition",
        "pdg_id_1",
        "pdg_id_2",
        "origin_1",
        "origin_2",
        "pairs",
        "sum_weights",
        "sum_weights_squared",
    )
    compositions = ("DIS-DIS", "DIS-Soft", "Soft-Soft")
    detail_fields = (
        "source_file",
        "file_uuid",
        "source_entry",
        "material",
        "step",
        "selection",
        "composition",
        "track_id_1",
        "track_id_2",
        "pdg_id_1",
        "pdg_id_2",
        "origin_1",
        "origin_2",
        "weight",
    )

    def __init__(self, writer=None, scale=None):
        self.counts = {}
        self.totals = {}
        self.writer = writer
        self.scale = scale
        self.source = None

    def add(self, material, step, tracks, muon_weight, weight):
        particles = sorted(
            (
                "Soft" if math.isclose(track.weight, muon_weight, rel_tol=1.0e-6, abs_tol=0.0) else "DIS",
                track.pdg_id,
                track.track_id,
            )
            for track in tracks
        )
        (origin1, pdg1, track1), (origin2, pdg2, track2) = particles
        composition = f"{origin1}-{origin2}"
        if self.writer is not None:
            if self.source is None:
                raise ValueError("Individual pair output requires source event identity")
            row = dict(
                zip(
                    self.detail_fields,
                    (
                        *self.source,
                        material,
                        step,
                        STEPS[step],
                        composition,
                        track1,
                        track2,
                        pdg1,
                        pdg2,
                        origin1,
                        origin2,
                        weight,
                    ),
                )
            )
            if self.scale is not None:
                row["scaled_weight"] = weight * self.scale
            self.writer.writerow(row)
        for mat in ("ALL", material):
            key = (mat, step, composition, pdg1, pdg2, origin1, origin2)
            count = self.counts.setdefault(key, [0, 0.0, 0.0])
            count[0] += 1
            count[1] += weight
            count[2] += weight * weight
            total_key = (mat, step, composition)
            self.totals[total_key] = self.totals.get(total_key, 0) + 1

    def summary(self, material, step):
        counts = [self.totals.get((material, step, composition), 0) for composition in self.compositions]
        # Step 0 counts input events, before candidate pairs are defined.
        return dict(
            zip(
                ("pairs", "pairs_DIS_DIS", "pairs_DIS_Soft", "pairs_Soft_Soft"),
                [sum(counts), *counts] if step else [""] * 4,
            )
        )

    def rows(self):
        return [
            dict(zip(self.fields, (mat, step, STEPS[step], composition, pdg1, pdg2, origin1, origin2, *count)))
            for (mat, step, composition, pdg1, pdg2, origin1, origin2), count in sorted(self.counts.items())
        ]


def pair_composition_rows(rows):
    """Include all three compositions, even when no pairs survive a cut."""
    materials = list(dict.fromkeys(["ALL", *MATERIALS, *(row["material"] for row in rows)]))
    totals = {
        (mat, step, composition): [0, 0.0, 0.0]
        for mat in materials
        for step in range(1, len(STEPS))
        for composition in PairCutflow.compositions
    }
    for row in rows:
        count = totals[row["material"], int(row["step"]), row["composition"]]
        for i, column in enumerate(("pairs", "sum_weights", "sum_weights_squared")):
            count[i] += row[column]
    return [
        {
            "material": mat,
            "step": step,
            "selection": STEPS[step],
            "composition": composition,
            "pairs": count[0],
            "sum_weights": count[1],
            "sum_weights_squared": count[2],
        }
        for (mat, step, composition), count in totals.items()
    ]


def replay(event, plots, geometry, pdg, factor, counts, material_counts=None, pair_counts=None):
    if not len(event.MCTrack):
        raise ValueError("Replay event has no initial MCTrack[0]")
    muon_weight = float(event.MCTrack[0].GetWeight())
    base_weight = shared_weight([], muon_weight) * factor

    def event_weight(weights):
        return shared_weight(weights, muon_weight, event=True) * factor

    def passed(step, weight):
        for cutflow in (counts,) if material_counts is None else (counts, material_counts):
            cutflow[step][0] += 1
            cutflow[step][1] += weight
            cutflow[step][2] += weight * weight

    def pair_passed(step, tracks, weight):
        if pair_counts is not None:
            pair_counts.add(str(event.muDIS_material), step, tracks, muon_weight, weight)

    passed(0, base_weight)
    tracks, momenta = selected_tracks(event, pdg)
    plots.one("step0/n_tracks_all_detectors", len(momenta), -0.5, 20.5, base_weight, 21)
    for p, track_weight in momenta:
        plots.one("step0/track_p_GeV", p, 0, 100, track_weight * factor)
    plots.one("step0/n_tracks_p_gt_1", len(tracks), -0.5, 20.5, base_weight, 21)
    if len(tracks) < 2:
        return
    # Validate the common DIS factor, even when no vertex can be fitted.
    shared_weight((track.weight for track in tracks), muon_weight)
    pair_weights = [shared_weight((a.weight, b.weight), muon_weight) for a, b in combinations(tracks, 2)]
    stage1_weight = event_weight(pair_weights)
    passed(1, stage1_weight)
    candidates = []
    for (a, b), pair_weight in zip(combinations(tracks, 2), pair_weights):
        weight = pair_weight * factor
        pair_passed(1, (a, b), weight)
        plots.one("step1/opposite_sign", int(a.charge * b.charge < 0), -0.5, 1.5, weight, 2)
        vertex = make_vertex(a, b)
        if vertex is None:
            continue
        vertex.weight = pair_weight
        vertex.tracks = (a, b)
        plots.one("step1/doca_cm", vertex.doca, 0, 10, weight)
        if vertex.opposite:
            plots.one("step1/opposite_sign_doca_cm", vertex.doca, 0, 10, weight)
        if vertex.opposite and vertex.doca < 1.0:
            pair_passed(2, vertex.tracks, weight)
            candidates.append(vertex)
    plots.one("step1/n_good_pairs", len(candidates), -0.5, 20.5, stage1_weight, 21)
    if not candidates:
        return
    passed(2, event_weight(v.weight for v in candidates))
    contained = []
    for vertex in candidates:
        weight = vertex.weight * factor
        plots.vertex("step2/vertices", vertex.position, weight)
        clearance = geometry.clearance(vertex.position)
        plots.one("step2/helium_clearance_cm", clearance, -2, 300, weight)
        if clearance >= 5.0:
            pair_passed(3, vertex.tracks, weight)
            contained.append(vertex)
    if not contained:
        return
    passed(3, event_weight(v.weight for v in contained))
    loose = []
    for vertex in contained:
        weight = vertex.weight * factor
        ip = target_ip(vertex)
        plots.one("step3/target_ip_cm", ip, 0, 500, weight)
        if ip < 250.0:
            pair_passed(4, vertex.tracks, weight)
            loose.append(vertex)
    if not loose:
        return
    passed(4, event_weight(v.weight for v in loose))
    tight = []
    for vertex in loose:
        weight = vertex.weight * factor
        ip = target_ip(vertex)
        plots.one("step4/target_ip_cm", ip, 0, 250, weight)
        plots.one("step4/target_ip_zoom_cm", ip, 0, 20, weight)
        if ip < 10.0:
            pair_passed(5, vertex.tracks, weight)
            tight.append(vertex)
    if tight:
        passed(5, event_weight(v.weight for v in tight))
        for vertex in tight:
            plots.vertex("step5/vertices", vertex.position, vertex.weight * factor)


def run(config):
    with ExitStack() as output_streams:
        _run(config, output_streams)


def _run(config, output_streams):
    kind = config["type"]
    files = discover(config["production_path"], kind)
    limit = config.get("n_events", -1)
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < -1:
        raise ValueError("n_events must be -1 or a nonnegative integer")
    first = config.get("first_entry", 0)
    if isinstance(first, bool) or not isinstance(first, int) or first < 0:
        raise ValueError("first_entry must be a nonnegative integer")
    save_plots = config.get("save_plots", True)
    if not isinstance(save_plots, bool):
        raise ValueError("save_plots must be true or false")
    factor = float(config.get("weight_factor", 1.0))
    if not math.isfinite(factor):
        raise ValueError("weight_factor must be finite")
    output = Path(config.get("output_dir", f"analyseEvents_{kind}")).expanduser()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output directory must be empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    import ROOT

    ROOT.gROOT.SetBatch(True)
    print("Loading libShipMuDIS.so", flush=True)
    result = ROOT.gSystem.Load("libShipMuDIS.so")
    print(f"libShipMuDIS.so load result: {result}", flush=True)
    if result < 0:
        raise RuntimeError("Cannot load libShipMuDIS.so")
    is_replay = kind.startswith("disreplay")
    tree_name = "cbmsim" if is_replay else "MuonDIS"
    normalization = read_normalization(ROOT, files, tree_name)
    scale = normalization_scale(config, normalization)
    chain = ROOT.TChain(tree_name)
    for file in files:
        if chain.Add(str(file)) != 1:
            raise ValueError(f"Cannot add {file}")
    total = int(chain.GetEntries())
    if first > total:
        raise ValueError(f"first_entry {first} exceeds the {total} available entries")
    n = total - first if limit == -1 else min(limit, total - first)
    plots = Plots(config.get("plot_axes", {}))
    # Book even empty stages/materials so every requested output is present.
    prefixes = ["step2/vertices", "step5/vertices"] if is_replay else MATERIALS
    for prefix in prefixes:
        for vertical, horizontal in [("x", "z"), ("y", "z"), ("y", "x")]:
            plots.book(f"{prefix}/{vertical}_vs_{horizontal}", (plots.axes[horizontal], plots.axes[vertical]))
    if is_replay:
        for name, low, high, bins in [
            ("step0/n_tracks_all_detectors", -0.5, 20.5, 21),
            ("step0/track_p_GeV", 0, 100, 100),
            ("step0/n_tracks_p_gt_1", -0.5, 20.5, 21),
            ("step1/opposite_sign", -0.5, 1.5, 2),
            ("step1/doca_cm", 0, 10, 100),
            ("step1/opposite_sign_doca_cm", 0, 10, 100),
            ("step1/n_good_pairs", -0.5, 20.5, 21),
            ("step2/helium_clearance_cm", -2, 300, 100),
            ("step3/target_ip_cm", 0, 500, 100),
            ("step4/target_ip_cm", 0, 250, 100),
            ("step4/target_ip_zoom_cm", 0, 20, 100),
        ]:
            plots.book(name, (np.linspace(low, high, bins + 1),))
    geometry = Geometry(ROOT, config.get("geometry_file"))
    counts = [[0, 0.0, 0.0] for _ in STEPS]
    replay_counts = {mat: [[0, 0.0, 0.0] for _ in STEPS] for mat in MATERIALS}
    pair_writer = None
    if is_replay:
        pair_stream = output_streams.enter_context((output / "pairs.csv").open("w", newline=""))
        detail_fields = [*PairCutflow.detail_fields, *(["scaled_weight"] if scale is not None else [])]
        pair_writer = csv.DictWriter(pair_stream, fieldnames=detail_fields)
        pair_writer.writeheader()
    pair_counts = PairCutflow(pair_writer, scale)
    material_counts = {mat: [0, 0.0, 0.0] for mat in MATERIALS}
    previous_tree = -1
    source_ranges = []
    print(f"{len(files)} files; processing entries [{first}, {first + n}) of {total}", flush=True)
    for entry in range(first, first + n):
        local_entry = int(chain.LoadTree(entry))
        if local_entry < 0:
            raise RuntimeError(f"Cannot load entry {entry}")
        if chain.GetTreeNumber() != previous_tree:
            previous_tree = chain.GetTreeNumber()
            required = (
                ["MCTrack", "strawtubesPoint", "TimeDetPoint", "muDIS_material"]
                if is_replay
                else [
                    f"{prefix}_{mat}"
                    for mat in MATERIALS
                    for prefix in ("mudis_DISvx", "mudis_DISvy", "mudis_DISvz", "muon_wDIS")
                ]
            )
            missing = [name for name in required if not chain.GetTree().GetBranch(name)]
            if missing:
                raise ValueError(f"Missing branches in {chain.GetFile().GetName()}: {missing}")
            if is_replay:
                geometry.load(Path(chain.GetFile().GetName()))
            source_ranges.append(
                {
                    "path": str(Path(chain.GetFile().GetName()).resolve()),
                    "file_uuid": str(chain.GetFile().GetUUID().AsString()),
                    "tree": tree_name,
                    "entries": int(chain.GetTree().GetEntries()),
                    "start": local_entry,
                    "stop": local_entry,
                }
            )
        if chain.GetEntry(entry) <= 0:
            raise RuntimeError(f"Cannot read entry {entry}")
        source_ranges[-1]["stop"] = local_entry + 1
        if is_replay:
            material = str(chain.muDIS_material)
            pair_counts.source = (source_ranges[-1]["path"], source_ranges[-1]["file_uuid"], local_entry)
            cutflow = replay_counts.setdefault(material, [[0, 0.0, 0.0] for _ in STEPS])
            replay(chain, plots, geometry, ROOT.TDatabasePDG.Instance(), factor, counts, cutflow, pair_counts)
        else:
            for mat in MATERIALS:
                xyz = [getattr(chain, f"mudis_DISv{axis}_{mat}") for axis in "xyz"]
                if len({len(values) for values in xyz}) != 1:
                    raise ValueError(f"Inconsistent vertex arrays at entry {entry}, {mat}")
                weight = float(getattr(chain, f"muon_wDIS_{mat}")) * factor
                if not math.isfinite(weight):
                    raise ValueError(f"Nonfinite wDIS at entry {entry}, {mat}")
                for position in zip(*xyz):
                    if not all(math.isfinite(value) for value in position):
                        raise ValueError(f"Nonfinite DIS vertex at entry {entry}, {mat}")
                    plots.vertex(mat, position, weight)
                    material_counts[mat][0] += 1
                    material_counts[mat][1] += weight
                    material_counts[mat][2] += weight * weight
        if (entry - first + 1) % 1000 == 0:
            print(f"Processed {entry - first + 1}/{n}", flush=True)
    output_streams.close()
    if is_replay:
        rows = [
            {
                "material": material,
                "step": i,
                "selection": label,
                "events": count[0],
                "sum_weights": count[1],
                "sum_weights_squared": count[2],
                **pair_counts.summary(material, i),
            }
            for material, cutflow in [("ALL", counts), *replay_counts.items()]
            for i, (label, count) in enumerate(zip(STEPS, cutflow))
        ]
        pair_fields = list(PairCutflow.fields)
        if scale is not None:
            pair_fields += ["scaled_sum_weights", "scaled_sum_weights_squared"]
        write_table(output / "pair_cutflow.csv", scaled_rows(pair_counts.rows(), scale), pair_fields)
        write_table(
            output / "pair_composition_cutflow.csv", scaled_rows(pair_composition_rows(pair_counts.rows()), scale)
        )
    else:
        rows = [
            {"material": mat, "vertices": count[0], "sum_weights": count[1], "sum_weights_squared": count[2]}
            for mat, count in material_counts.items()
        ]
    rows = scaled_rows(rows, scale)
    write_table(output / ("cutflow.csv" if is_replay else "vertex_counts.csv"), rows)
    save_normalization(output, normalization, scale)
    plots.save_archive(output / "histograms.npz")
    if save_plots:
        print("Rendering plots", flush=True)
        plots.save(output, scale)
    else:
        plots.save_flows(output, scale)
    # Write the completion manifest last, after every requested output succeeds.
    with (output / "run.json").open("w") as stream:
        json.dump(
            {
                "output_version": OUTPUT_VERSION,
                "analysis_signature": {
                    "sha256": ANALYSIS_SHA256,
                    "type": kind,
                    "weight_factor": factor,
                    "target_sum_weights": float(config["target_sum_weights"])
                    if config.get("target_sum_weights") is not None
                    else None,
                    "failed_job_factor": float(config.get("failed_job_factor", 1.0)),
                    "geometry_file": str(Path(config["geometry_file"]).expanduser().resolve())
                    if config.get("geometry_file")
                    else None,
                },
                "config": config,
                "files": [str(file.resolve()) for file in files],
                "source_ranges": source_ranges,
                "normalization": normalization,
                "normalization_scale": scale,
                "chain_entry_start": first,
                "chain_entry_stop": first + n,
                "entries_available": total,
                "entries_processed": n,
                "weight_definition": (
                    "track fW; pairs count the common muon/DIS factor once; "
                    "event OR uses muon fW if any soft-soft pair survives, otherwise DIS fW; "
                    "step0 uses muon fW; all multiplied by weight_factor"
                    if is_replay
                    else "wDIS * weight_factor"
                ),
            },
            stream,
            indent=2,
        )
    for row in rows:
        print(row, flush=True)
    print(f"Output: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="Analysis YAML file")
    parser.add_argument("--first-entry", type=int, help="First TChain entry (zero-based)")
    parser.add_argument("--n-events", type=int, help="Number of entries, or -1 for all remaining")
    parser.add_argument("--output-dir", help="Unique, empty output directory for this job")
    parser.add_argument("--no-plots", action="store_true", help="Save tables and histograms without rendering")
    args = parser.parse_args()
    with args.config.open() as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        parser.error("YAML must contain a mapping")
    for key in ("first_entry", "n_events", "output_dir"):
        if getattr(args, key) is not None:
            config[key] = getattr(args, key)
    if args.no_plots:
        config["save_plots"] = False
    run(config)


if __name__ == "__main__":
    main()

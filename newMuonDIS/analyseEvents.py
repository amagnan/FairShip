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
from analysisTree import CachedGeometry, NullPlots, TreeWriter

MATERIALS = ("MS", "UBT", "SBTsens", "SBTfr", "SSTsens", "SSTfr", "HE", "AIR", "HeBalloon", "REST")
TYPES = ("disgen", "disfilter", "disreplay_gen", "disreplay_filter")
OUTPUT_VERSION = 7
ANALYSIS_SHA256 = hashlib.sha256(
    Path(__file__).read_bytes() + Path(__file__).with_name("analysisTree.py").read_bytes()
).hexdigest()
STEPS = (
    "filtered input: Tr1-4 + TD, p > 1 GeV/c; opposite sign, backward vertex, DOCA < 5 cm, inside decay helium",
    "DOCA < 1 cm and helium wall clearance >= 5 cm",
    "target-plane IP < 250 cm",
    "target-plane IP < 10 cm",
)
UBT_WINDOWS = {"full": None, "400x600": (400, 600), "360x300": (360, 300), "140x310": (140, 310)}
PAIR_MOMENTUM_HISTOGRAMS = (
    ("p1_GeV", 0, 100),
    ("p2_GeV", 0, 100),
)
HIT_CATEGORIES = ("neither", "UBT_only", "SBT_only", "both")
SBT_WALLS = {"top": (1, 1), "bottom": (1, -1), "left": (0, -1), "right": (0, 1)}


def particle_charge(code, pdg):
    particle = pdg.GetParticle(code)
    if particle:
        return float(particle.Charge())
    if 1000000000 <= abs(code) < 1100000000:
        return 3.0 * ((abs(code) // 10000) % 1000) * (1 if code > 0 else -1)
    raise ValueError(f"Unknown PDG code: {code}")


def is_soft(weight, muon_weight):
    return math.isclose(weight, muon_weight, rel_tol=1.0e-6, abs_tol=0.0)


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
    """Return normalization and usable inputs, excluding keyless ROOT files."""
    inputs = []
    usable_files = []
    for path in files:
        source = root.TFile.Open(str(path), "READ")
        if not source or source.IsZombie():
            raise ValueError(f"Cannot read normalization metadata: {path}")
        try:
            if source.GetNkeys() == 0:
                print(f"WARNING: Skipping empty ROOT file (no keys): {path}", flush=True)
                continue
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
            usable_files.append(path)
        finally:
            source.Close()
    if not usable_files:
        raise ValueError("No usable analysis inputs remain: all ROOT files are empty")
    return combine_normalization(inputs), usable_files


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


def detector_track_states(event):
    """Associate detector coverage and the smallest-z Tr1 state by MCTrack ID."""
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
    return stations, first, timing


def selected_tracks(event, pdg):
    """Associate hits by MCTrack index; use the smallest-z Tr1 hit per track."""
    stations, first, timing = detector_track_states(event)
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
        charge = particle_charge(code, pdg)
        tracks.append(Track(position, momentum, charge, weight, code, index))
    return sorted(tracks, key=lambda track: track.track_id), momenta


class Geometry:
    def __init__(self, root, override=None):
        self.root, self.override, self.path, self.manager = root, override, None, None
        self.sbt_cells = {}

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
        self.sbt_cells = self._sbt_cells()
        self.path = path

    @staticmethod
    def _cell_wall(name, position):
        if name.startswith(("LiScX_", "LiSc_S3_", "LiSc_S5_")):
            return "top" if position[1] > 0 else "bottom"
        if name.startswith(("LiScY_", "LiSc_S4_", "LiSc_S6_")):
            return "right" if position[0] > 0 else "left"
        return None

    def _sbt_cells(self):
        """Map sensitive copy numbers to walls using their placed cell centres."""
        cells = {}

        def visit(node, parent):
            matrix = self.root.TGeoHMatrix(parent)
            matrix.Multiply(node.GetMatrix())
            volume = node.GetVolume()
            name = str(volume.GetName())
            if self._cell_wall(name, (1, 1, 0)) is not None:
                center = np.empty(3, dtype=np.float64)
                matrix.LocalToMaster(volume.GetShape().GetOrigin(), center)
                wall = self._cell_wall(name, center)
                detector_id = int(node.GetNumber())
                if detector_id in cells and cells[detector_id] != wall:
                    raise ValueError(f"Ambiguous SBT detector ID {detector_id} in replay geometry")
                cells[detector_id] = wall
            for index in range(node.GetNdaughters()):
                visit(node.GetDaughter(index), matrix)

        visit(self.manager.GetTopNode(), self.root.TGeoHMatrix())
        return cells

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

    def sbt_wall(self, position, detector_id=None):
        """Identify the actual scintillator cell, including tapered corner cells.

        veto::ProcessHits records the sensitive cell's copy number, but its
        stored entry/exit midpoint need not locate that cell with FindNode.
        Prefer the cell placement over the measured position when available.
        For other geometry versions, use transverse helium boundaries at hit z.
        """
        if detector_id in self.sbt_cells:
            return self.sbt_cells[detector_id]
        x, y, z = map(float, position)
        node = self.manager.FindNode(x, y, z)
        name = str(node.GetVolume().GetName()) if node else ""
        wall = self._cell_wall(name, position)
        if wall is not None:
            return wall
        # Find the decay-medium cross section at this z, in world coordinates.
        # DistFromInside on the actual shape accounts for its taper, rather
        # than using the bounding box's maximum width/height at every z.
        if self.clearance((0, 0, z)) < 0:
            raise ValueError(f"Cannot locate decay helium at SBT hit z={z}: {position}")
        navigator = self.manager.GetCurrentNavigator()
        shape = navigator.GetCurrentNode().GetVolume().GetShape()
        center = np.array([0.0, 0.0, z], dtype=np.float64)
        local = np.empty(3, dtype=np.float64)
        navigator.MasterToLocal(center, local)
        bounds = {}
        for wall, (axis, sign) in SBT_WALLS.items():
            direction = np.zeros(3, dtype=np.float64)
            direction[axis] = sign
            local_direction = np.empty(3, dtype=np.float64)
            navigator.MasterToLocalVect(direction, local_direction)
            distance = float(shape.DistFromInside(local, local_direction, 3))
            if not math.isfinite(distance) or not 0 < distance < 1.0e10:
                raise ValueError(f"Invalid decay-volume boundary for SBT hit: {position}")
            bounds[wall] = sign * distance
        # Distance to the finite transverse side segment resolves corner hits
        # uniquely. Exact ties go to top/bottom (dictionary order).
        distances = {}
        for wall, (axis, _) in SBT_WALLS.items():
            along = x if axis == 1 else y
            low, high = (bounds["left"], bounds["right"]) if axis == 1 else (bounds["bottom"], bounds["top"])
            beyond_end = max(low - along, 0.0, along - high)
            across = (y if axis == 1 else x) - bounds[wall]
            distances[wall] = math.hypot(across, beyond_end)
        return min(distances, key=distances.get)


class Plots:
    """Accumulate bounded-memory histograms, retaining under/overflow counts."""

    def __init__(self, config):
        self.histograms = {}
        self.material = None
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
        self._fill(name, values, edges, weight)
        if self.material is not None:
            self._fill(f"materials/{self.material}/{name}", values, edges, weight)

    def _fill(self, name, values, edges, weight):
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
        # Material histograms are allocated only when filled. Split jobs may
        # see different materials/stages; their common aggregate schema must match.
        def base(plots):
            return {name for name in plots.histograms if not name.startswith("materials/")}

        if base(self) != base(other):
            raise ValueError("Histogram names differ between jobs")
        for plots in (self, other):
            for name, (edges, *_) in plots.histograms.items():
                if name.startswith("materials/"):
                    parts = name.split("/", 2)
                    template = self.histograms.get(parts[2]) if len(parts) == 3 else None
                    if (
                        template is None
                        or len(edges) != len(template[0])
                        or any(not np.array_equal(a, b) for a, b in zip(edges, template[0]))
                    ):
                        raise ValueError(f"Invalid material histogram: {name}")
        for name, (edges, raw, weighted, flow, wflow) in other.histograms.items():
            self.book(name, edges)
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
                    ax.set_ylabel("Number of hits" if vertical == "n_hits" else f"{vertical} [cm]")
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
        "truth_dis_vertices",
        "truth_dis_in_decay_volume",
        "truth_dis_all_in_decay_volume",
    )

    def __init__(self, writer=None, scale=None):
        self.counts = {}
        self.totals = {}
        self.writer = writer
        self.scale = scale
        self.source = None

    def add(self, material, step, tracks, muon_weight, weight, truth=None):
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
            row.update(truth or dict.fromkeys(ReplayDiagnostics.truth_fields, ""))
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
        return dict(
            zip(
                ("pairs", "pairs_DIS_DIS", "pairs_DIS_Soft", "pairs_Soft_Soft"),
                [sum(counts), *counts],
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
        for step in range(len(STEPS))
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


def particle_counts(tracks, excluded=()):
    """Count stored MCTracks, not unique physical particles across a shower."""
    charges = [track["charge_e"] for track in tracks if track["track_id"] not in excluded]
    charged = sum(charge != 0 for charge in charges)
    return {"charged": charged, "neutral": len(charges) - charged}


def initial_dis_products(event):
    """Generator-injected DIS products, excluding soft primaries and showers."""
    muon_weight = float(event.MCTrack[0].GetWeight())
    return [
        i
        for i, track in enumerate(event.MCTrack)
        if i and track.GetMotherId() == 0 and not is_soft(float(track.GetWeight()), muon_weight)
    ]


def selected_pair_info(event, tracks, vertex):
    """Order by MCTrack ID; find the nearest shared ancestor beyond input muon 0."""
    tracks = sorted(tracks, key=lambda track: track.track_id)

    def ancestry(index):
        chain = []
        while index > 0:
            if index >= len(event.MCTrack) or index in chain:
                raise ValueError("Invalid MCTrack ancestry")
            chain.append(index)
            index = int(event.MCTrack[index].GetMotherId())
        return chain

    first, second = (ancestry(track.track_id) for track in tracks)
    common = set(first).intersection(second)
    ancestor = min(common, key=lambda i: (first.index(i) + second.index(i), i)) if common else -1
    p1, p2 = (float(np.linalg.norm(track.momentum)) for track in tracks)
    return {
        "track_ids": [track.track_id for track in tracks],
        "p1_GeV": p1,
        "p2_GeV": p2,
        "doca_cm": vertex.doca,
        "vertex_cm": vertex.position.tolist(),
        "target_ip_cm": target_ip(vertex),
        "common_ancestor_id": ancestor,
        "common_ancestor_pdg": int(event.MCTrack[ancestor].GetPdgCode()) if ancestor >= 0 else 0,
    }


class ReplayRecords:
    """Stream one self-contained record per event; never expand pairs by stage."""

    def __init__(self, stream, scale=None):
        self.stream, self.scale = stream, scale
        self.source = None
        self.record = None

    @staticmethod
    def optional_branch(event, name):
        # Inspect the current tree, not TChain's first-file branch list: older
        # replay files may lack an optional branch present in another file.
        tree = event.GetTree() if hasattr(event, "GetTree") else event
        if hasattr(tree, "GetBranch") and not tree.GetBranch(name):
            return None
        return getattr(event, name, None)

    def start_event(self, event, geometry, pdg, factor):
        if self.source is None:
            raise ValueError("Event records require source identity")
        self.event, self.geometry = event, geometry
        stations, first, timing = detector_track_states(event)
        muon_weight = float(event.MCTrack[0].GetWeight())
        tracks = []
        for index, mc in enumerate(event.MCTrack):
            weight = float(mc.GetWeight())
            code = int(mc.GetPdgCode())
            state = first.get(index)
            tracks.append(
                {
                    "track_id": index,
                    "pdg_id": code,
                    "mother_id": int(mc.GetMotherId()),
                    "charge_e": particle_charge(code, pdg) / 3.0,
                    "origin": "Soft" if is_soft(weight, muon_weight) else "DIS",
                    "mc_weight": weight,
                    "weight": weight * factor,
                    "start_position_cm": [float(getattr(mc, f"GetStart{axis}")()) for axis in "XYZ"],
                    "start_momentum_GeV": [float(getattr(mc, f"GetP{axis}")()) for axis in "xyz"],
                    "straw_stations": sorted(stations.get(index, ())),
                    "has_timing_hit": index in timing,
                    "tr1_position_cm": state[0].tolist() if state else None,
                    "tr1_momentum_GeV": state[1].tolist() if state else None,
                    "selected": bool(
                        state is not None
                        and stations[index] == {1, 2, 3, 4}
                        and index in timing
                        and np.linalg.norm(state[1]) > 1.0
                    ),
                    "stage_mask": 0,
                }
            )
        hits = {}
        for detector, branch in (("UBT", "UpstreamTaggerPoint"), ("SBT", "vetoPoint")):
            hits[detector] = [
                [
                    int(hit.GetTrackID()),
                    float(hit.GetX()),
                    float(hit.GetY()),
                    float(hit.GetZ()),
                    float(hit.GetTime()),
                    int(hit.GetDetectorID()),
                    float(hit.GetPx()),
                    float(hit.GetPy()),
                    float(hit.GetPz()),
                    float(hit.GetEnergyLoss()),
                ]
                for hit in getattr(event, branch)
            ]
        lengths = {}
        for material in ("total", *MATERIALS):
            branch = "muon_path_length" + ("" if material == "total" else f"_{material}")
            value = self.optional_branch(event, branch)
            if value is not None:
                value = float(value)
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f"Invalid {branch}: {value}")
            lengths[material] = value
        identity = {}
        for name in ("muDIS_muEntry", "muDIS_disIndex"):
            value = self.optional_branch(event, name)
            identity[name] = int(value) if value is not None else None
        self.record = {
            "source_file": self.source[0],
            "file_uuid": self.source[1],
            "source_entry": self.source[2],
            "material": str(event.muDIS_material),
            "geometry_file": str(geometry.path) if getattr(geometry, "path", None) else None,
            "replay_identity": identity,
            "passed_steps": [False] * len(STEPS),
            "stage_weights": [None] * len(STEPS),
            "particle_counts": {
                "initial_DIS": particle_counts([tracks[i] for i in initial_dis_products(event)]),
                "excluding_input_muon": particle_counts(tracks, excluded={0}),
                "all": particle_counts(tracks),
                "direct_daughters": particle_counts([track for track in tracks if track["mother_id"] == 0]),
            },
            "leading_pair": None,
            "selected_pairs": [None] * len(STEPS),
            "tracks": tracks,
            "hits": hits,
            "muon_path_length_cm": lengths,
        }
        self.leading_rank = None
        self.stage_ranks = {}

    def stage(self, step, weight, tracks=()):
        self.record["passed_steps"][step] = True
        self.record["stage_weights"][step] = weight
        for track in tracks:
            self.record["tracks"][track.track_id]["stage_mask"] |= 1 << step

    def pair(self, step, tracks, vertex):
        if vertex is None:
            return
        indices = sorted(track.track_id for track in tracks)
        rank = (-step, vertex.doca, *indices)
        if step not in self.stage_ranks or rank < self.stage_ranks[step]:
            self.stage_ranks[step] = rank
            self.record["selected_pairs"][step] = selected_pair_info(self.event, tracks, vertex)
        if self.leading_rank is None or rank < self.leading_rank:
            self.leading_rank = rank
            self.record["leading_pair"] = {"track_ids": indices, "step": step}

    @staticmethod
    def rescale(record, scale):
        record.pop("scaled_stage_weights", None)
        for track in record["tracks"]:
            track.pop("scaled_weight", None)
        if scale is not None:
            record["scaled_stage_weights"] = [None if w is None else w * scale for w in record["stage_weights"]]
            for track in record["tracks"]:
                track["scaled_weight"] = track["weight"] * scale

    def finish_event(self):
        record = self.record
        pair = record["leading_pair"]
        record["particle_counts"]["excluding_input_muon_and_leading_pair"] = particle_counts(
            record["tracks"], excluded={0, *(pair["track_ids"] if pair else [])}
        )
        record["geometry_cache"] = self.geometry.export()
        if hasattr(self.stream, "write_record"):
            self.stream.write_record(record)
        else:
            # In-memory diagnostic clients can still inspect a record as JSON.
            self.rescale(record, self.scale)
            self.stream.write(json.dumps(record, allow_nan=False, separators=(",", ":")) + "\n")
        self.record = None


class ReplayDiagnostics:
    """Stage-wise event categories and one lossless hit record per event/track."""

    fields = (
        "material",
        "step",
        "selection",
        "sample",
        "window",
        "category",
        "events",
        "sum_weights",
        "sum_weights_squared",
    )
    detail_fields = (
        "source_file",
        "file_uuid",
        "source_entry",
        "material",
        "sample",
        "track_id",
        "pdg_id",
        "weight",
        "categories",
        "ubt_hits",
        "sbt_hits",
        "sbt_walls",
        *(f"passed_step{step}" for step in range(len(STEPS))),
        *(f"weight_step{step}" for step in range(len(STEPS))),
    )
    weight_fields = ("weight", *(f"weight_step{step}" for step in range(len(STEPS))))
    scaled_detail_fields = tuple(f"scaled_{column}" for column in weight_fields)
    truth_fields = ("truth_dis_vertices", "truth_dis_in_decay_volume", "truth_dis_all_in_decay_volume")

    def __init__(self, plots, writer=None, scale=None):
        self.plots, self.writer, self.scale = plots, writer, scale
        self.source = None
        self.counts = {}
        self.materials = set()
        for material in ("ALL", *MATERIALS):
            self.book_counts(material)
        for step in range(len(STEPS)):
            for charge in ("charged", "neutral"):
                plots.book(f"step{step}/particles/n_{charge}", (np.arange(-0.5, 1001.5),))
            for axis in ("x", "y"):
                plots.book(f"step{step}/tagger/all/SBT/{axis}_vs_z", (plots.axes["z"], plots.axes[axis]))
            for sample in ("muon", "daughter"):
                prefix = f"step{step}/tagger/{sample}"
                for window in UBT_WINDOWS:
                    plots.book(f"{prefix}/UBT/{window}/y_vs_x", (plots.axes["x"], plots.axes["y"]))
                for wall, (axis, _) in SBT_WALLS.items():
                    coordinate = "x" if axis == 1 else "y"
                    plots.book(f"{prefix}/SBT/{wall}/{coordinate}_vs_z", (plots.axes["z"], plots.axes[coordinate]))
                    plots.book(f"{prefix}/SBT/{wall}/n_hits_vs_entry_z", (plots.axes["z"], np.arange(-0.5, 101.5)))
            for name in ("distance_xy_cm", "distance_z_cm", "distance_3d_cm"):
                plots.book(f"step{step}/truth_DIS/{name}", (np.linspace(0, 15000, 301),))
                plots.book(f"step{step}/truth_DIS/{name.removesuffix('_cm')}_zoom_cm", (np.linspace(0, 100, 201),))
            plots.book(f"step{step}/truth_DIS/in_decay_volume", (np.arange(-0.5, 2.5),))

    def book_counts(self, material):
        if material in self.materials:
            return
        self.materials.add(material)
        for step in range(len(STEPS)):
            for sample in ("muon", "daughter_any", "daughter_all"):
                for window in UBT_WINDOWS:
                    for category in (*HIT_CATEGORIES, "no_daughters"):
                        if sample == "muon" and category == "no_daughters":
                            continue
                        self.counts.setdefault((material, step, sample, window, category), [0, 0.0, 0.0])

    def start_event(self, event, pdg, geometry, factor):
        self.event, self.geometry, self.factor = event, geometry, factor
        self.stage_weights = {}
        self.muon_weight = float(event.MCTrack[0].GetWeight())
        self.material = str(event.muDIS_material)
        charges = [particle_charge(int(event.MCTrack[i].GetPdgCode()), pdg) for i in initial_dis_products(event)]
        self.particle_counts = {
            "charged": sum(charge != 0 for charge in charges),
            "neutral": sum(charge == 0 for charge in charges),
        }
        self.book_counts(self.material)
        self.daughters = {
            i
            for i, track in enumerate(event.MCTrack)
            if i
            and track.GetMotherId() == 0
            and not is_soft(float(track.GetWeight()), self.muon_weight)
            and particle_charge(int(track.GetPdgCode()), pdg) != 0
        }
        self.hits = {i: {"UBT": [], "SBT": []} for i in (0, *sorted(self.daughters))}
        for branch, detector in (("UpstreamTaggerPoint", "UBT"), ("vetoPoint", "SBT")):
            for hit in getattr(event, branch):
                index = int(hit.GetTrackID())
                if index not in self.hits:
                    continue
                point = [
                    float(hit.GetX()),
                    float(hit.GetY()),
                    float(hit.GetZ()),
                    float(hit.GetTime()),
                    int(hit.GetDetectorID()),
                ]
                if not all(math.isfinite(v) for v in point):
                    raise ValueError(f"Nonfinite {branch} hit for track {index}")
                self.hits[index][detector].append(point)
        self.truth_cache = {}
        self.wall_hits = {
            index: [geometry.sbt_wall(point[:3], point[4]) for point in hits["SBT"]]
            for index, hits in self.hits.items()
        }

    def window_hits(self, index, window):
        hits = self.hits[index]["UBT"]
        widths = UBT_WINDOWS[window]
        return (
            hits if widths is None else [p for p in hits if abs(p[0]) <= widths[0] / 2 and abs(p[1]) <= widths[1] / 2]
        )

    def categories(self, index):
        return {
            window: HIT_CATEGORIES[int(bool(self.window_hits(index, window))) + 2 * int(bool(self.hits[index]["SBT"]))]
            for window in UBT_WINDOWS
        }

    def stage(self, step, weight, tracks=()):
        for charge, count in self.particle_counts.items():
            self.plots.one(f"step{step}/particles/n_{charge}", count, -0.5, 1000.5, weight, 1001)
        # Include all veto points, even those with an absent/invalid MCTrack ID.
        for hit in self.event.vetoPoint:
            for axis in ("x", "y"):
                self.plots.fill(
                    f"step{step}/tagger/all/SBT/{axis}_vs_z",
                    (float(hit.GetZ()), float(getattr(hit, f"Get{axis.upper()}")())),
                    (self.plots.axes["z"], self.plots.axes[axis]),
                    weight,
                )
        daughters = self.daughters if step == 0 else self.daughters.intersection(t.track_id for t in tracks)
        categories = {i: self.categories(i) for i in (0, *sorted(daughters))}
        for window in UBT_WINDOWS:
            selected = [categories[i][window] for i in daughters]
            samples = {
                "muon": [categories[0][window]],
                "daughter_any": sorted(set(selected)) if selected else ["no_daughters"],
                "daughter_all": [selected[0]]
                if selected and len(set(selected)) == 1
                else ([] if selected else ["no_daughters"]),
            }
            for sample, labels in samples.items():
                for label in labels:
                    for mat in ("ALL", self.material):
                        count = self.counts[mat, step, sample, window, label]
                        count[0] += 1
                        count[1] += weight
                        count[2] += weight * weight
        for index in categories:
            sample = "muon" if index == 0 else "daughter"
            track = self.event.MCTrack[index]
            track_weight = weight if index == 0 else float(track.GetWeight()) * self.factor
            if not math.isfinite(track_weight):
                raise ValueError(f"Nonfinite diagnostic weight for track {index}")
            if self.writer is not None:
                self.stage_weights.setdefault(index, [None] * len(STEPS))[step] = track_weight
            prefix = f"step{step}/tagger/{sample}"
            for window in UBT_WINDOWS:
                for x, y, *_ in self.window_hits(index, window):
                    self.plots.fill(
                        f"{prefix}/UBT/{window}/y_vs_x",
                        (x, y),
                        (self.plots.axes["x"], self.plots.axes["y"]),
                        track_weight,
                    )
            sbt = self.hits[index]["SBT"]
            for wall, (axis, _) in SBT_WALLS.items():
                hits = [p for p, assigned in zip(sbt, self.wall_hits[index]) if assigned == wall]
                coordinate = "x" if axis == 1 else "y"
                for point in hits:
                    self.plots.fill(
                        f"{prefix}/SBT/{wall}/{coordinate}_vs_z",
                        (point[2], point[1 - axis]),
                        (self.plots.axes["z"], self.plots.axes[coordinate]),
                        track_weight,
                    )
                if sbt:
                    # Time ordering also handles backward-going daughters. Use the
                    # first SBT hit overall, with zero counts for unvisited walls.
                    entry_z = min(sbt, key=lambda p: p[3])[2]
                    self.plots.fill(
                        f"{prefix}/SBT/{wall}/n_hits_vs_entry_z",
                        (entry_z, len(hits)),
                        (self.plots.axes["z"], np.arange(-0.5, 101.5)),
                        track_weight,
                    )

    def finish_event(self):
        """Flush only this event's tracks after all selection stages have run."""
        if self.writer is None:
            return
        if self.source is None:
            raise ValueError("Track diagnostics require source event identity")
        for index, weights in sorted(self.stage_weights.items()):
            track = self.event.MCTrack[index]
            row = dict(
                zip(
                    self.detail_fields,
                    (
                        *self.source,
                        self.material,
                        "muon" if index == 0 else "daughter",
                        index,
                        int(track.GetPdgCode()),
                        float(track.GetWeight()) * self.factor,
                        json.dumps(self.categories(index)),
                        json.dumps(self.hits[index]["UBT"]),
                        json.dumps(self.hits[index]["SBT"]),
                        json.dumps(self.wall_hits[index]),
                        *(weight is not None for weight in weights),
                        *(weight if weight is not None else "" for weight in weights),
                    ),
                )
            )
            if self.scale is not None:
                for column in self.weight_fields:
                    row[f"scaled_{column}"] = row[column] * self.scale if row[column] != "" else ""
            self.writer.writerow(row)
        self.stage_weights.clear()

    def truth_vertices(self, tracks):
        vertices = set()
        for track in tracks:
            if is_soft(track.weight, self.muon_weight):
                continue
            index, visited = track.track_id, set()
            while index > 0:
                if index in visited or index >= len(self.event.MCTrack):
                    raise ValueError("Invalid MCTrack ancestry")
                visited.add(index)
                mc = self.event.MCTrack[index]
                mother = int(mc.GetMotherId())
                if mother == 0:
                    if not is_soft(float(mc.GetWeight()), self.muon_weight):
                        position = tuple(float(getattr(mc, f"GetStart{axis}")()) for axis in "XYZ")
                        if not all(math.isfinite(v) for v in position):
                            raise ValueError("Nonfinite truth DIS vertex")
                        vertices.add(position)
                    break
                index = mother
        return sorted(vertices)

    def pair_truth(self, step, tracks, weight, vertex=None):
        positions = self.truth_vertices(tracks)
        contained = []
        for position in positions:
            if position not in self.truth_cache:
                self.truth_cache[position] = self.geometry.clearance(position) >= 0
            contained.append(self.truth_cache[position])
            if vertex is not None:
                delta = vertex.position - np.asarray(position)
                for name, distance in (
                    ("xy", np.linalg.norm(delta[:2])),
                    ("z", abs(delta[2])),
                    ("3d", np.linalg.norm(delta)),
                ):
                    self.plots.one(f"step{step}/truth_DIS/distance_{name}_cm", distance, 0, 15000, weight, 300)
                    self.plots.one(f"step{step}/truth_DIS/distance_{name}_zoom_cm", distance, 0, 100, weight, 200)
        if positions:
            self.plots.one(f"step{step}/truth_DIS/in_decay_volume", int(any(contained)), -0.5, 1.5, weight, 2)
        return dict(
            zip(
                self.truth_fields,
                (
                    json.dumps(positions),
                    int(any(contained)) if positions else "",
                    int(all(contained)) if positions else "",
                ),
            )
        )

    def rows(self):
        return [
            dict(zip(self.fields, (mat, step, STEPS[step], sample, window, category, *count)))
            for (mat, step, sample, window, category), count in sorted(self.counts.items())
        ]


def replay(
    event, plots, geometry, pdg, factor, counts, material_counts=None, pair_counts=None, diagnostics=None, records=None
):
    if records is not None:
        geometry = CachedGeometry(geometry)
    if diagnostics is None:
        diagnostics = ReplayDiagnostics(plots)
    material = str(event.muDIS_material)
    if not material or Path(material).name != material or material in (".", ".."):
        raise ValueError(f"Invalid material label: {material!r}")
    plots.material = material
    try:
        _replay(event, plots, geometry, pdg, factor, counts, material_counts, pair_counts, diagnostics, records)
    finally:
        plots.material = None
    # _replay returns early when an event fails a cut; flush those events too.
    diagnostics.finish_event()
    if records is not None:
        records.finish_event()


def _replay(event, plots, geometry, pdg, factor, counts, material_counts, pair_counts, diagnostics, records=None):
    if not len(event.MCTrack):
        raise ValueError("Replay event has no initial MCTrack[0]")
    muon_weight = float(event.MCTrack[0].GetWeight())
    shared_weight([], muon_weight)  # Validate the incoming-muon weight.
    diagnostics.start_event(event, pdg, geometry, factor)
    if records is not None:
        records.start_event(event, geometry, pdg, factor)

    def event_weight(weights):
        return shared_weight(weights, muon_weight, event=True) * factor

    def passed(step, weight, tracks=()):
        for cutflow in (counts,) if material_counts is None else (counts, material_counts):
            cutflow[step][0] += 1
            cutflow[step][1] += weight
            cutflow[step][2] += weight * weight
        diagnostics.stage(step, weight, tracks)
        if records is not None:
            records.stage(step, weight, tracks)

    def pair_passed(step, tracks, weight, vertex=None):
        truth = diagnostics.pair_truth(step, tracks, weight, vertex)
        if pair_counts is not None:
            pair_counts.add(str(event.muDIS_material), step, tracks, muon_weight, weight, truth)
        if records is not None:
            records.pair(step, tracks, vertex)

    def surviving_tracks(vertices):
        return [track for vertex in vertices for track in vertex.tracks]

    tracks, momenta = selected_tracks(event, pdg)
    shared_weight((track.weight for track in tracks), muon_weight)
    candidates = []
    for a, b in combinations(tracks, 2):
        vertex = make_vertex(a, b)
        if vertex is None or not vertex.opposite or vertex.doca >= 5.0:
            continue
        if geometry.clearance(vertex.position) < 0:
            continue
        vertex.weight = shared_weight((a.weight, b.weight), muon_weight)
        vertex.tracks = (a, b)
        candidates.append(vertex)
    if not candidates:
        raise ValueError("Replay input does not pass replayFilter selection; run replayFilter.py first")

    for step in range(len(STEPS)):
        weight = event_weight(v.weight for v in candidates)
        passed(step, weight, surviving_tracks(candidates))
        leading = min(candidates, key=lambda v: (v.doca, *sorted(t.track_id for t in v.tracks)))
        info = selected_pair_info(event, leading.tracks, leading)
        for name, low, high in PAIR_MOMENTUM_HISTOGRAMS:
            plots.one(f"step{step}/selected_pair/{name}", info[name], low, high, weight)
        daughters = {t.track_id: t for t in surviving_tracks(candidates) if t.track_id != 0}
        for track in daughters.values():
            plots.one(
                f"step{step}/daughters/p_GeV", float(np.linalg.norm(track.momentum)), 0, 100, track.weight * factor
            )
        plots.one(f"step{step}/n_pairs", len(candidates), -0.5, 20.5, weight, 21)
        for vertex in candidates:
            pair_weight = vertex.weight * factor
            pair_passed(step, vertex.tracks, pair_weight, vertex)
            plots.vertex(f"step{step}/vertices", vertex.position, pair_weight)
            plots.one(f"step{step}/doca_cm", vertex.doca, 0, 5, pair_weight)
            plots.one(f"step{step}/helium_clearance_cm", geometry.clearance(vertex.position), 0, 300, pair_weight)
            ip = target_ip(vertex)
            plots.one(f"step{step}/target_ip_cm", ip, 0, 500, pair_weight)
            plots.one(f"step{step}/target_ip_zoom_cm", ip, 0, 20, pair_weight)
        if step == 0:
            plots.one("step0/n_tracks_all_detectors", len(momenta), -0.5, 20.5, weight, 21)
            plots.one("step0/n_tracks_p_gt_1", len(tracks), -0.5, 20.5, weight, 21)
            for momentum, track_weight in momenta:
                plots.one("step0/track_p_GeV", momentum, 0, 100, track_weight * factor)
            candidates = [v for v in candidates if v.doca < 1.0 and geometry.clearance(v.position) >= 5.0]
        elif step == 1:
            candidates = [v for v in candidates if target_ip(v) < 250.0]
        elif step == 2:
            candidates = [v for v in candidates if target_ip(v) < 10.0]
        if not candidates:
            break


def book_plots(config, is_replay=True, tree_only=False):
    plots = (NullPlots if tree_only else Plots)(config.get("plot_axes", {}))
    # Book aggregate templates, including empty stages. Material copies are lazy.
    prefixes = [f"step{step}/vertices" for step in range(len(STEPS))] if is_replay else MATERIALS
    for prefix in prefixes:
        for vertical, horizontal in [("x", "z"), ("y", "z"), ("y", "x")]:
            plots.book(f"{prefix}/{vertical}_vs_{horizontal}", (plots.axes[horizontal], plots.axes[vertical]))
    if is_replay:
        for name, low, high, bins in [
            ("step0/n_tracks_all_detectors", -0.5, 20.5, 21),
            ("step0/track_p_GeV", 0, 100, 100),
            ("step0/n_tracks_p_gt_1", -0.5, 20.5, 21),
        ]:
            plots.book(name, (np.linspace(low, high, bins + 1),))
        for step in range(len(STEPS)):
            for name, low, high, bins in [
                ("n_pairs", -0.5, 20.5, 21),
                ("doca_cm", 0, 5, 100),
                ("helium_clearance_cm", 0, 300, 100),
                ("target_ip_cm", 0, 500, 100),
                ("target_ip_zoom_cm", 0, 20, 100),
            ]:
                plots.book(f"step{step}/{name}", (np.linspace(low, high, bins + 1),))
    if is_replay:
        for step in range(len(STEPS)):
            for name, low, high in PAIR_MOMENTUM_HISTOGRAMS:
                plots.book(f"step{step}/selected_pair/{name}", (np.linspace(low, high, 101),))
            plots.book(f"step{step}/daughters/p_GeV", (np.linspace(0, 100, 101),))
    return plots


def replay_tables(output, counts, replay_counts, pair_counts, diagnostics, scale):
    tagger_fields = [
        *ReplayDiagnostics.fields,
        *(["scaled_sum_weights", "scaled_sum_weights_squared"] if scale is not None else []),
    ]
    write_table(output / "tagger_cutflow.csv", scaled_rows(diagnostics.rows(), scale), tagger_fields)
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
    write_table(output / "pair_composition_cutflow.csv", scaled_rows(pair_composition_rows(pair_counts.rows()), scale))
    return rows


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
    tree_only = config.get("tree_only", False)
    if not isinstance(tree_only, bool) or (tree_only and not is_replay):
        raise ValueError("tree_only must be a boolean and is supported only for replay")
    tree_name = "cbmsim" if is_replay else "MuonDIS"
    normalization, files = read_normalization(ROOT, files, tree_name)
    normalization["weight_factor"] = factor
    scale = normalization_scale(config, normalization)
    chain = ROOT.TChain(tree_name)
    for file in files:
        if chain.Add(str(file)) != 1:
            raise ValueError(f"Cannot add {file}")
    total = int(chain.GetEntries())
    if first > total:
        raise ValueError(f"first_entry {first} exceeds the {total} available entries")
    n = total - first if limit == -1 else min(limit, total - first)
    plots = book_plots(config, is_replay, tree_only)
    geometry = Geometry(ROOT, config.get("geometry_file"))
    counts = [[0, 0.0, 0.0] for _ in STEPS]
    replay_counts = {mat: [[0, 0.0, 0.0] for _ in STEPS] for mat in MATERIALS}
    pair_counts = PairCutflow()
    diagnostics = records = None
    if is_replay:
        event_stream = TreeWriter(ROOT, output / "events.root")
        output_streams.callback(event_stream.close)
        records = ReplayRecords(event_stream)
        diagnostics = ReplayDiagnostics(plots)
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
                ["MCTrack", "strawtubesPoint", "TimeDetPoint", "UpstreamTaggerPoint", "vetoPoint", "muDIS_material"]
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
            diagnostics.source = pair_counts.source
            records.source = pair_counts.source
            cutflow = replay_counts.setdefault(material, [[0, 0.0, 0.0] for _ in STEPS])
            replay(
                chain,
                plots,
                geometry,
                ROOT.TDatabasePDG.Instance(),
                factor,
                counts,
                cutflow,
                pair_counts,
                diagnostics,
                records,
            )
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
        rows = replay_tables(output, counts, replay_counts, pair_counts, diagnostics, scale) if not tree_only else []
    else:
        rows = [
            {"material": mat, "vertices": count[0], "sum_weights": count[1], "sum_weights_squared": count[2]}
            for mat, count in material_counts.items()
        ]
    rows = scaled_rows(rows, scale)
    if not tree_only:
        write_table(output / ("cutflow.csv" if is_replay else "vertex_counts.csv"), rows)
    save_normalization(output, normalization, scale)
    if not tree_only:
        plots.save_archive(output / "histograms.npz")
    if save_plots and not tree_only:
        print("Rendering plots", flush=True)
        plots.save(output, scale)
    elif not tree_only:
        plots.save_flows(output, scale)
    # Write the completion manifest last, after every requested output succeeds.
    with (output / "run.json").open("w") as stream:
        json.dump(
            {
                "output_version": OUTPUT_VERSION,
                "tree_only": tree_only,
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
                    "all stages including filtered baseline use the surviving-pair OR; all multiplied by weight_factor"
                    if is_replay
                    else "wDIS * weight_factor"
                ),
                "event_record_schema": {
                    "stages": list(STEPS),
                    "leading_pair": "smallest DOCA among built vertices at the furthest passed stage; "
                    "ties ordered by MCTrack IDs",
                    "particle_count_default": "initial DIS products: mother ID 0 and weight different from input muon",
                    "tree": "events.root:DIS; typed branches, one entry per event",
                    "selected_pairs": "smallest DOCA at each passed stage, ties by MCTrack IDs",
                    "common_ancestor": "nearest shared ancestor, including a selected parent; excludes input muon 0; -1/0 ID/PDG when absent",
                    "scaled_weights": "multiply saved weights by normalization_scale from run.json",
                    "hit_columns": [
                        "track_id",
                        "x_cm",
                        "y_cm",
                        "z_cm",
                        "time_ns",
                        "detector_id",
                        "px_GeV",
                        "py_GeV",
                        "pz_GeV",
                        "energy_loss_GeV",
                    ],
                    "missing_path_length": "NaN; zero means an available branch with zero path length",
                }
                if is_replay
                else None,
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
    parser.add_argument(
        "--tree-only", action="store_true", help="Extract replay event tree; defer histograms and tables"
    )
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
    if args.tree_only:
        config["tree_only"] = True
    run(config)


if __name__ == "__main__":
    main()

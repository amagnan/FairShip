#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Plot prepared DIS vertices or analyse replay truth tracks; see analyseEvents.yaml."""

import argparse
import csv
import json
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import yaml

MATERIALS = ("MS", "UBT", "SBTsens", "SBTfr", "SSTsens", "SSTfr", "HE", "AIR", "HeBalloon", "REST")
TYPES = ("disgen", "disfilter", "disreplay_gen", "disreplay_filter")
STEPS = (
    "all",
    "at least two tracks: Tr1-4 + TD, p > 1 GeV/c",
    "opposite sign, backward vertex, DOCA < 1 cm",
    "inside decay helium, wall clearance >= 5 cm",
    "target-plane IP < 250 cm",
    "target-plane IP < 10 cm",
)


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
            match = name.startswith("muonDis_") and not any(token in name for token in ("disfilter", "disreplay"))
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


@dataclass
class Vertex:
    position: np.ndarray
    momentum: np.ndarray
    doca: float
    opposite: bool
    weight: float = 1.0


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
        tracks.append(Track(position, momentum, charge, weight))
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

    def save(self, output):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        flows = []
        for name, (edges, raw, weighted, flow, wflow) in self.histograms.items():
            for label, data in [("unweighted", raw), ("weighted", weighted)]:
                destination = output / "plots" / name / label
                destination.parent.mkdir(parents=True, exist_ok=True)
                fig, ax = plt.subplots(figsize=(8, 6))
                if len(edges) == 2:
                    mesh = ax.pcolormesh(edges[0], edges[1], data.T, shading="auto")
                    fig.colorbar(mesh, ax=ax, label="Sum of weights" if label == "weighted" else "Entries")
                    vertical, horizontal = name.split("/")[-1].split("_vs_")
                    ax.set_xlabel(f"{horizontal} [cm]")
                    ax.set_ylabel(f"{vertical} [cm]")
                else:
                    ax.stairs(data, edges[0])
                    ax.set_xlabel(name.split("/")[-1])
                    ax.set_ylabel("Sum of weights" if label == "weighted" else "Entries")
                ax.set_title(f"{name} ({label})\nOutside plotted range: {wflow if label == 'weighted' else flow:g}")
                fig.tight_layout()
                for extension in ("pdf", "png"):
                    fig.savefig(destination.with_suffix("." + extension))
                plt.close(fig)
            flows.append({"histogram": name, "outside_range": flow, "weighted_outside_range": wflow})
        write_table(output / "histogram_flows.csv", flows)


def write_table(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def replay(event, plots, geometry, pdg, factor, counts):
    if not len(event.MCTrack):
        raise ValueError("Replay event has no initial MCTrack[0]")
    muon_weight = float(event.MCTrack[0].GetWeight())
    base_weight = shared_weight([], muon_weight) * factor

    def event_weight(weights):
        return shared_weight(weights, muon_weight, event=True) * factor

    def passed(step, weight):
        counts[step][0] += 1
        counts[step][1] += weight
        counts[step][2] += weight * weight

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
        plots.one("step1/opposite_sign", int(a.charge * b.charge < 0), -0.5, 1.5, weight, 2)
        vertex = make_vertex(a, b)
        if vertex is None:
            continue
        vertex.weight = pair_weight
        plots.one("step1/doca_cm", vertex.doca, 0, 10, weight)
        if vertex.opposite:
            plots.one("step1/opposite_sign_doca_cm", vertex.doca, 0, 10, weight)
        if vertex.opposite and vertex.doca < 1.0:
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
            tight.append(vertex)
    if tight:
        passed(5, event_weight(v.weight for v in tight))
        for vertex in tight:
            plots.vertex("step5/vertices", vertex.position, vertex.weight * factor)


def run(config):
    kind = config["type"]
    files = discover(config["production_path"], kind)
    limit = config.get("n_events", -1)
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < -1:
        raise ValueError("n_events must be -1 or a nonnegative integer")
    factor = float(config.get("weight_factor", 1.0))
    if not math.isfinite(factor):
        raise ValueError("weight_factor must be finite")
    output = Path(config.get("output_dir", f"analyseEvents_{kind}")).expanduser()
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
    chain = ROOT.TChain(tree_name)
    for file in files:
        if chain.Add(str(file)) != 1:
            raise ValueError(f"Cannot add {file}")
    total = int(chain.GetEntries())
    n = total if limit == -1 else min(limit, total)
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
    material_counts = {mat: [0, 0.0, 0.0] for mat in MATERIALS}
    previous_tree = -1
    print(f"{len(files)} files; processing {n}/{total} tree entries", flush=True)
    for entry in range(n):
        if chain.LoadTree(entry) < 0:
            raise RuntimeError(f"Cannot load entry {entry}")
        if chain.GetTreeNumber() != previous_tree:
            previous_tree = chain.GetTreeNumber()
            required = (
                ["MCTrack", "strawtubesPoint", "TimeDetPoint"]
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
        if chain.GetEntry(entry) <= 0:
            raise RuntimeError(f"Cannot read entry {entry}")
        if is_replay:
            replay(chain, plots, geometry, ROOT.TDatabasePDG.Instance(), factor, counts)
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
        if (entry + 1) % 1000 == 0:
            print(f"Processed {entry + 1}/{n}", flush=True)
    if is_replay:
        rows = [
            {
                "step": i,
                "selection": label,
                "events": count[0],
                "sum_weights": count[1],
                "sum_weights_squared": count[2],
            }
            for i, (label, count) in enumerate(zip(STEPS, counts))
        ]
    else:
        rows = [
            {"material": mat, "vertices": count[0], "sum_weights": count[1], "sum_weights_squared": count[2]}
            for mat, count in material_counts.items()
        ]
    write_table(output / ("cutflow.csv" if is_replay else "vertex_counts.csv"), rows)
    with (output / "run.json").open("w") as stream:
        json.dump(
            {
                "config": config,
                "files": [str(file) for file in files],
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
    print("Rendering plots", flush=True)
    plots.save(output)
    for row in rows:
        print(row, flush=True)
    print(f"Output: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="Analysis YAML file")
    args = parser.parse_args()
    with args.config.open() as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        parser.error("YAML must contain a mapping")
    run(config)


if __name__ == "__main__":
    main()

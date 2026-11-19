# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Typed, self-contained replay summaries. No FairShip dictionaries are needed to read them."""

import math
from array import array
from types import SimpleNamespace

import numpy as np

TREE_NAME = "DIS"
N_STEPS = 4
MATERIALS = ("MS", "UBT", "SBTsens", "SBTfr", "SSTsens", "SSTfr", "HE", "AIR", "HeBalloon", "REST")
PAIR_FLOATS = ("p1_GeV", "p2_GeV", "doca_cm", "target_ip_cm")
TRACK_INTS = ("pdg_id", "mother_id", "has_timing_hit", "selected", "stage_mask")
TRACK_FLOATS = ("charge_e", "mc_weight", "weight")
TRACK_XYZ = ("start_position_cm", "start_momentum_GeV", "tr1_position_cm", "tr1_momentum_GeV")
HIT_COLUMNS = (
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
)
WALLS = (None, "top", "bottom", "left", "right")

# S: std::string, L: Long64_t, D: double, vi/vd: vector<int/double>.
SCHEMA = dict.fromkeys(("source_file", "file_uuid", "material", "geometry_file"), "S")
SCHEMA.update(dict.fromkeys(("source_entry", "muDIS_muEntry", "muDIS_disIndex", "n_charged", "n_neutral"), "L"))
SCHEMA.update(muDIS_pPythia="D", muDIS_xsec="D")
SCHEMA.update(
    passed_steps="vi",
    stage_weights="vd",
    selected_track_1="vi",
    selected_track_2="vi",
    common_ancestor_id="vi",
    common_ancestor_pdg="vi",
)
SCHEMA.update(dict.fromkeys(PAIR_FLOATS, "vd"))
SCHEMA.update({f"vertex_{axis}_cm": "vd" for axis in "xyz"})
SCHEMA.update({f"track_{name}": "vi" for name in (*TRACK_INTS, "soft", "station_mask")})
SCHEMA.update({f"track_{name}": "vd" for name in TRACK_FLOATS})
SCHEMA.update({f"track_{name}_{axis}": "vd" for name in TRACK_XYZ for axis in "xyz"})
for detector in ("UBT", "SBT"):
    SCHEMA.update({f"{detector}_{name}": "vi" if i in (0, 5) else "vd" for i, name in enumerate(HIT_COLUMNS)})
SCHEMA.update({f"muon_path_length_{name}": "D" for name in ("total", *MATERIALS)})
SCHEMA.update({f"clearance_{name}": "vd" for name in ("x", "y", "z", "value")})
SCHEMA.update({f"entrance_clearance_{name}": "vd" for name in ("x", "y", "z", "value")})
SCHEMA.update({f"wall_{name}": "vd" for name in ("x", "y", "z")})
SCHEMA.update(wall_detector_id="vi", wall_value="vi")


def nullable(value):
    return float("nan") if value is None else value


def flatten(record):
    values = {name: record[name] or "" for name in ("source_file", "file_uuid", "material", "geometry_file")}
    for name in ("muDIS_pPythia", "muDIS_xsec"):
        values[name] = nullable(record[name])
    values["source_entry"] = record["source_entry"]
    for name in ("muDIS_muEntry", "muDIS_disIndex"):
        value = record["replay_identity"][name]
        values[name] = -1 if value is None else value
    for charge in ("charged", "neutral"):
        values[f"n_{charge}"] = record["particle_counts"]["initial_DIS"][charge]
    values["passed_steps"] = list(map(int, record["passed_steps"]))
    values["stage_weights"] = list(map(nullable, record["stage_weights"]))
    pairs = record["selected_pairs"]
    for i in range(2):
        values[f"selected_track_{i + 1}"] = [pair["track_ids"][i] if pair else -1 for pair in pairs]
    for name in ("common_ancestor_id", "common_ancestor_pdg", *PAIR_FLOATS):
        missing = (-1 if name.endswith("_id") else 0) if name.startswith("common_") else float("nan")
        values[name] = [pair[name] if pair else missing for pair in pairs]
    for i, axis in enumerate("xyz"):
        values[f"vertex_{axis}_cm"] = [pair["vertex_cm"][i] if pair else float("nan") for pair in pairs]
    tracks = record["tracks"]
    for name in (*TRACK_INTS, *TRACK_FLOATS):
        values[f"track_{name}"] = [int(t[name]) if name in TRACK_INTS else t[name] for t in tracks]
    values["track_soft"] = [int(t["origin"] == "Soft") for t in tracks]
    values["track_station_mask"] = [sum(1 << s for s in t["straw_stations"]) for t in tracks]
    for name in TRACK_XYZ:
        for i, axis in enumerate("xyz"):
            values[f"track_{name}_{axis}"] = [t[name][i] if t[name] is not None else float("nan") for t in tracks]
    for detector, hits in record["hits"].items():
        for i, name in enumerate(HIT_COLUMNS):
            values[f"{detector}_{name}"] = [h[i] for h in hits]
    for name, value in record["muon_path_length_cm"].items():
        values[f"muon_path_length_{name}"] = nullable(value)
    for prefix, names in (
        ("clearance", ("x", "y", "z", "value")),
        ("entrance_clearance", ("x", "y", "z", "value")),
        ("wall", ("x", "y", "z", "detector_id", "value")),
    ):
        rows = record["geometry_cache"][prefix]
        for i, name in enumerate(names):
            values[f"{prefix}_{name}"] = [r[i] for r in rows]
    return values


def unflatten(values):
    """Reconstruct the small Python view used by the existing selection/diagnostics."""
    record = {name: values[name] for name in ("source_file", "file_uuid", "source_entry", "material", "geometry_file")}
    for name in ("muDIS_pPythia", "muDIS_xsec"):
        record[name] = None if math.isnan(values[name]) else values[name]
    record["geometry_file"] = record["geometry_file"] or None
    record["replay_identity"] = {
        name: None if values[name] < 0 else values[name] for name in ("muDIS_muEntry", "muDIS_disIndex")
    }
    record["passed_steps"] = list(map(bool, values["passed_steps"]))
    record["stage_weights"] = [None if math.isnan(w) else w for w in values["stage_weights"]]
    tracks = []
    for i in range(len(values["track_pdg_id"])):
        track = {name: values[f"track_{name}"][i] for name in (*TRACK_INTS, *TRACK_FLOATS)}
        track.update(
            track_id=i,
            origin="Soft" if values["track_soft"][i] else "DIS",
            straw_stations=[s for s in range(1, 5) if values["track_station_mask"][i] & (1 << s)],
        )
        for name in ("selected", "has_timing_hit"):
            track[name] = bool(track[name])
        for name in TRACK_XYZ:
            xyz = [values[f"track_{name}_{axis}"][i] for axis in "xyz"]
            track[name] = None if name.startswith("tr1_") and math.isnan(xyz[0]) else xyz
        tracks.append(track)
    record["tracks"] = tracks
    record["selected_pairs"] = []
    record["leading_pair"] = None
    for step, passed in enumerate(record["passed_steps"]):
        pair = None
        if passed:
            pair = {name: values[name][step] for name in ("common_ancestor_id", "common_ancestor_pdg", *PAIR_FLOATS)}
            pair["track_ids"] = [values[f"selected_track_{i}"][step] for i in (1, 2)]
            pair["vertex_cm"] = [values[f"vertex_{axis}_cm"][step] for axis in "xyz"]
            record["leading_pair"] = {"track_ids": pair["track_ids"], "step": step}
        record["selected_pairs"].append(pair)
    record["hits"] = {
        detector: list(map(list, zip(*(values[f"{detector}_{name}"] for name in HIT_COLUMNS))))
        for detector in ("UBT", "SBT")
    }
    record["muon_path_length_cm"] = {
        name: None if math.isnan(values[f"muon_path_length_{name}"]) else values[f"muon_path_length_{name}"]
        for name in ("total", *MATERIALS)
    }
    record["geometry_cache"] = {
        prefix: list(map(list, zip(*(values[f"{prefix}_{name}"] for name in names))))
        for prefix, names in (
            ("clearance", ("x", "y", "z", "value")),
            ("entrance_clearance", ("x", "y", "z", "value")),
            ("wall", ("x", "y", "z", "detector_id", "value")),
        )
    }
    from analyseEvents import particle_counts

    leading = record["leading_pair"]
    record["particle_counts"] = {
        "initial_DIS": {charge: values[f"n_{charge}"] for charge in ("charged", "neutral")},
        "all": particle_counts(tracks),
        "excluding_input_muon": particle_counts(tracks, {0}),
        "direct_daughters": particle_counts([t for t in tracks if t["mother_id"] == 0]),
        "excluding_input_muon_and_leading_pair": particle_counts(
            tracks, {0, *(leading["track_ids"] if leading else [])}
        ),
    }
    return record


class TreeWriter:
    def __init__(self, root, path):
        self.file = root.TFile.Open(str(path), "RECREATE")
        if not self.file or self.file.IsZombie():
            raise OSError(f"Cannot create {path}")
        # ZLIB level 1: modest compression CPU cost and broad ROOT compatibility.
        self.file.SetCompressionSettings(101)
        self.tree = root.TTree(TREE_NAME, "One entry per filtered DIS event; cm and GeV/c")
        self.tree.SetAutoFlush(-32 * 1024 * 1024)
        # Keep a single output file, including productions larger than ROOT's default 100 GB.
        self.tree.SetMaxTreeSize(1_000_000_000_000_000)
        self.buffers = {}
        for name, kind in SCHEMA.items():
            if kind == "S":
                buffer = root.std.string()
            elif kind.startswith("v"):
                buffer = root.std.vector("int" if kind == "vi" else "double")()
            else:
                buffer = array("q" if kind == "L" else "d", [0])
            self.buffers[name] = buffer
            if kind in ("L", "D"):
                self.tree.Branch(name, buffer, f"{name}/{kind}")
            else:
                self.tree.Branch(name, buffer)

    def write_record(self, record):
        values = flatten(record)
        for name, kind in SCHEMA.items():
            buffer, value = self.buffers[name], values[name]
            if kind == "S":
                buffer.assign(value)
            elif kind.startswith("v"):
                buffer.clear()
                # PyROOT converts the contiguous numpy array in C++, avoiding a
                # Python/C++ call per hit or track component.
                buffer.assign(np.asarray(value, dtype=np.int32 if kind == "vi" else np.float64))
            else:
                buffer[0] = value
        if self.tree.Fill() < 0:
            raise OSError("Failed to fill DIS tree")

    def close(self):
        if self.file.IsOpen():
            self.file.cd()
            try:
                if self.tree.Write() <= 0:
                    raise OSError("Failed to write DIS tree")
            finally:
                self.file.Close()


def validate_tree(tree, path):
    if not tree or {b.GetName() for b in tree.GetListOfBranches()} != set(SCHEMA):
        raise ValueError(f"Invalid DIS tree schema: {path}")
    for name, kind in SCHEMA.items():
        branch = tree.GetBranch(name)
        if kind in ("L", "D"):
            leaf = branch.GetLeaf(name)
            valid = bool(leaf) and leaf.GetTypeName() == ("Long64_t" if kind == "L" else "Double_t")
        else:
            expected = {"S": "string", "vi": "vector<int>", "vd": "vector<double>"}[kind]
            valid = branch.GetClassName().replace("std::", "") == expected
        if not valid:
            raise ValueError(f"Invalid DIS branch type for {name}: {path}")


def read_records(path):
    import ROOT

    file = ROOT.TFile.Open(str(path))
    if not file or file.IsZombie():
        raise OSError(f"Cannot open {path}")
    try:
        tree = file.Get(TREE_NAME)
        validate_tree(tree, path)
        for entry in range(tree.GetEntries()):
            if tree.GetEntry(entry) <= 0:
                raise OSError(f"Cannot read {path}, entry {entry}")
            values = {}
            for name, kind in SCHEMA.items():
                value = getattr(tree, name)
                values[name] = str(value) if kind == "S" else list(value) if kind.startswith("v") else value
            yield unflatten(values)
    finally:
        file.Close()


class CachedGeometry:
    """Persist geometry queries so histogramming needs neither TGeo nor the source files."""

    def __init__(self, geometry=None, cache=None):
        self.geometry = geometry
        self.path = getattr(geometry, "path", None)
        self.clearances = {tuple(row[:3]): row[3] for row in cache["clearance"]} if cache else {}
        self.entrance_clearances = {tuple(row[:3]): row[3] for row in cache["entrance_clearance"]} if cache else {}
        self.walls = {tuple(row[:4]): WALLS[row[4]] for row in cache["wall"]} if cache else {}

    def clearance(self, position):
        key = tuple(map(float, position))
        if key not in self.clearances:
            if self.geometry is None:
                raise ValueError(f"Missing saved clearance for {key}")
            self.clearances[key] = self.geometry.clearance(position)
        return self.clearances[key]

    def entrance_clearance(self, position):
        key = tuple(map(float, position))
        if key not in self.entrance_clearances:
            if self.geometry is None:
                raise ValueError(f"Missing saved entrance clearance for {key}")
            self.entrance_clearances[key] = self.geometry.entrance_clearance(position)
        return self.entrance_clearances[key]

    def sbt_wall(self, position, detector_id=None):
        key = (*map(float, position), -1 if detector_id is None else int(detector_id))
        if key not in self.walls:
            if self.geometry is None:
                raise ValueError(f"Missing saved SBT wall for {key}")
            self.walls[key] = self.geometry.sbt_wall(position, detector_id)
        return self.walls[key]

    def export(self):
        return {
            "clearance": [[*key, value] for key, value in self.clearances.items()],
            "entrance_clearance": [[*key, value] for key, value in self.entrance_clearances.items()],
            "wall": [[*key, WALLS.index(value)] for key, value in self.walls.items()],
        }


class NullPlots:
    """Skip histogram allocation and filling during tree-only extraction."""

    def __init__(self, config):
        self.axes = {
            axis: np.linspace(*config.get(axis, default))
            for axis, default in {"x": [-600, 600, 121], "y": [-600, 600, 121], "z": [0, 14000, 141]}.items()
        }
        self.material = None

    def book(self, *args):
        pass

    fill = one = vertex = book


def event_view(record):
    """Expose only the MC/hit getters consumed by replay(), using saved double precision states."""

    def getters(values):
        return SimpleNamespace(**{name: lambda value=value: value for name, value in values.items()})

    tracks, straws, timing = [], [], []
    for t in record["tracks"]:
        values = {"GetPdgCode": t["pdg_id"], "GetMotherId": t["mother_id"], "GetWeight": t["mc_weight"]}
        values.update({f"GetStart{axis.upper()}": value for axis, value in zip("xyz", t["start_position_cm"])})
        values.update({f"GetP{axis}": value for axis, value in zip("xyz", t["start_momentum_GeV"])})
        tracks.append(getters(values))
        for station in t["straw_stations"]:
            values = {"GetTrackID": t["track_id"], "GetDetectorID": station * 1000000}
            values.update(
                {f"Get{axis.upper()}": value for axis, value in zip("xyz", t["tr1_position_cm"] or [0, 0, 0])}
            )
            values.update({f"GetP{axis}": value for axis, value in zip("xyz", t["tr1_momentum_GeV"] or [0, 0, 0])})
            straws.append(getters(values))
        if t["has_timing_hit"]:
            timing.append(getters({"GetTrackID": t["track_id"]}))
    branches = {}
    for detector, branch in (("UBT", "UpstreamTaggerPoint"), ("SBT", "vetoPoint")):
        names = (
            "GetTrackID",
            "GetX",
            "GetY",
            "GetZ",
            "GetTime",
            "GetDetectorID",
            "GetPx",
            "GetPy",
            "GetPz",
            "GetEnergyLoss",
        )
        branches[branch] = [getters(dict(zip(names, hit))) for hit in record["hits"][detector]]
    return SimpleNamespace(
        MCTrack=tracks,
        strawtubesPoint=straws,
        TimeDetPoint=timing,
        muDIS_material=record["material"],
        muDIS_pPythia=record["muDIS_pPythia"],
        **branches,
    )

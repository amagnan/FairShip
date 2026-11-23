# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Optional, raw-count diagnostics for replayFilter; no custom ROOT dictionary."""

import math
from array import array

from analyseEvents import initial_dis_products


def optional(event, name):
    if hasattr(event, "GetBranch") and not event.GetBranch(name):
        return None
    return getattr(event, name, None)


def dis_track_ids(event, primaries):
    """Follow ancestry, including neutral products and all shower descendants."""
    origins = {0: False, -1: False, **dict.fromkeys(primaries, True)}
    for index in range(1, len(event.MCTrack)):
        path, seen = [], set()
        current = index
        while current not in origins:
            if current in seen or not 0 <= current < len(event.MCTrack):
                raise ValueError("Invalid MCTrack ancestry in replay diagnostics")
            seen.add(current)
            path.append(current)
            current = int(event.MCTrack[current].GetMotherId())
        for ancestor in path:
            origins[ancestor] = origins[current]
    return {index for index, is_dis in origins.items() if is_dis}


def event_values(event):
    """Extract truth and hit counts independently of the filter outcome."""
    if not len(event.MCTrack):
        raise ValueError("Replay diagnostics require MCTrack[0]")
    muon = event.MCTrack[0]
    values = {f"muon_{axis}_cm": float(getattr(muon, f"GetStart{axis.upper()}")()) for axis in "xyz"}
    values.update({f"muon_p{axis}_GeV": float(getattr(muon, f"GetP{axis}")()) for axis in "xyz"})
    px, py, pz = (values[f"muon_p{axis}_GeV"] for axis in "xyz")
    values.update(
        muon_p_GeV=math.sqrt(px * px + py * py + pz * pz),
        muon_pt_GeV=math.hypot(px, py),
        muon_theta_rad=math.atan2(math.hypot(px, py), pz),
        muon_r_cm=math.hypot(values["muon_x_cm"], values["muon_y_cm"]),
        muon_energy_GeV=float(muon.GetEnergy()),
        muon_weight=float(muon.GetWeight()),
        muon_pdg=int(muon.GetPdgCode()),
    )
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("Nonfinite incoming-muon state in replay diagnostics")
    for name in ("muDIS_muEntry", "muDIS_disIndex"):
        value = optional(event, name)
        values[name] = int(value) if value is not None else -1
    material = optional(event, "muDIS_material")
    values["muDIS_material"] = str(material) if material is not None else "UNKNOWN"
    primaries = set(initial_dis_products(event))
    descendants = dis_track_ids(event, primaries)
    positions = [
        tuple(float(getattr(event.MCTrack[index], f"GetStart{axis}")()) for axis in "XYZ")
        for index in sorted(primaries)
    ]
    # 0: missing; 1: valid common position; 2: disagreeing positions; 3: nonfinite.
    status = 0
    if positions:
        status = 1
        if not all(math.isfinite(v) for position in positions for v in position):
            status = 3
        elif not all(
            all(math.isclose(a, b, rel_tol=0, abs_tol=1e-4) for a, b in zip(positions[0], position))
            for position in positions[1:]
        ):
            status = 2
    values["dis_vertex_status"] = status
    values["n_dis_primary"] = len(primaries)
    values.update(zip((f"dis_{axis}_cm" for axis in "xyz"), positions[0] if status == 1 else (math.nan,) * 3))
    for detector, branch in (("ubt", "UpstreamTaggerPoint"), ("sbt", "vetoPoint")):
        hits = optional(event, branch)
        counts = dict.fromkeys(("muon", "dis_direct", "dis", "muon_or_dis", "other", "unassigned"), 0)
        if hits is None:
            counts = dict.fromkeys(counts, -1)
        else:
            for hit in hits:
                index = int(hit.GetTrackID())
                if index == 0:
                    counts["muon"] += 1
                elif index in descendants:
                    counts["dis"] += 1
                    counts["dis_direct"] += int(index in primaries)
                elif 0 <= index < len(event.MCTrack):
                    counts["other"] += 1
                else:
                    counts["unassigned"] += 1
            counts["muon_or_dis"] = counts["muon"] + counts["dis"]
        values.update({f"n_{detector}_{origin}": count for origin, count in counts.items()})
    return values


class ScalarTree:
    """Own branch buffers for a compact tree containing only built-in types."""

    def __init__(self, root, name, integers, doubles, strings):
        self.tree = root.TTree(name, name)
        self.buffers = {}
        for names, code, leaf in ((integers, "q", "L"), (doubles, "d", "D")):
            for name in names:
                self.buffers[name] = array(code, [0])
                self.tree.Branch(name, self.buffers[name], f"{name}/{leaf}")
        for name in strings:
            self.buffers[name] = root.std.string()
            self.tree.Branch(name, self.buffers[name])
        self.strings = set(strings)

    def fill(self, values):
        for name, buffer in self.buffers.items():
            if name in self.strings:
                buffer.assign(values[name])
            else:
                buffer[0] = values[name]
        if self.tree.Fill() < 0:
            raise OSError(f"Cannot fill {self.tree.GetName()}")

    def write(self):
        if self.tree.Write() <= 0:
            raise OSError(f"Cannot write {self.tree.GetName()}")


class ReplayFilterDiagnostics:
    def __init__(self, root, summary):
        self.root, self.summary = root, summary
        integers = [
            "source_entry",
            "selected_entry",
            "accepted",
            "selection_stage",
            "selection_detail",
            "muDIS_muEntry",
            "muDIS_disIndex",
            "muon_pdg",
            "dis_vertex_status",
            "n_dis_primary",
        ]
        integers.extend(
            f"n_{detector}_{origin}"
            for detector in ("ubt", "sbt")
            for origin in ("muon", "dis_direct", "dis", "muon_or_dis", "other", "unassigned")
        )
        doubles = [f"muon_{axis}_cm" for axis in "xyz"] + [f"muon_p{axis}_GeV" for axis in "xyz"]
        doubles += [
            "muon_p_GeV",
            "muon_pt_GeV",
            "muon_theta_rad",
            "muon_r_cm",
            "muon_energy_GeV",
            "muon_weight",
            "dis_x_cm",
            "dis_y_cm",
            "dis_z_cm",
        ]
        self.events = ScalarTree(root, "replayFilterDiagnostics", integers, doubles, ("source_uuid", "muDIS_material"))
        self.muon_fields = tuple(name for name in doubles if name.startswith("muon_"))
        self.axes, self.specs = self.plot_specs()
        self.muons = {}
        self.histograms = {}
        self.book("ALL")
        summary["diagnostics"] = {
            "version": 1,
            "selection_detail": [
                "fewer_than_two_tracks",
                "no_opposite_sign_pair",
                "no_backward_vertex",
                "no_pair_with_doca_below_5_cm",
                "no_pair_inside_helium",
                "accepted",
            ],
            "dis_vertex_status": {"0": "missing", "1": "valid", "2": "ambiguous", "3": "nonfinite"},
            "dis_origin": "initial_dis_products ancestry/weight convention; DIS includes all descendants",
            "vertex_tolerance_cm": 1e-4,
            "hit_counts": "stored points, no energy/time cut; UBT=UpstreamTaggerPoint, SBT=vetoPoint; -1 if branch missing",
            "muon_grouping": "source UUID + muDIS_muEntry, within this file and inspected range only; missing IDs use source_entry",
            "efficiency": "unweighted replay-event accepted/all; bins with zero denominator have efficiency 0; retain all histogram",
            "histogram_ranges": "fixed ranges; underflow/overflow retained, including in efficiency division",
        }

    @staticmethod
    def plot_specs():
        axes = {
            "muon_p_GeV": (100, 0, 1000),
            "muon_pt_GeV": (100, 0, 20),
            "muon_theta_rad": (100, 0, math.pi),
            "muon_r_cm": (100, 0, 1000),
            "muon_x_cm": (100, -1000, 1000),
            "muon_y_cm": (100, -1000, 1000),
            "muon_z_cm": (150, -10000, 20000),
            "dis_x_cm": (100, -1000, 1000),
            "dis_y_cm": (100, -1000, 1000),
            "dis_z_cm": (150, -10000, 20000),
        }
        specs = {name: (name,) for name in axes}
        for prefix in ("muon", "dis"):
            for first, second in (("x", "y"), ("z", "x"), ("z", "y")):
                specs[f"{prefix}_{second}_vs_{first}"] = (f"{prefix}_{first}_cm", f"{prefix}_{second}_cm")
        specs["muon_p_vs_r"] = ("muon_r_cm", "muon_p_GeV")
        specs["muon_p_vs_dis_z"] = ("dis_z_cm", "muon_p_GeV")
        return axes, specs

    def book(self, material):
        if material in self.histograms:
            return
        histograms = {}
        for name, fields in self.specs.items():
            for sample in ("all", "accepted", "rejected"):
                hist_name = f"{name}_{sample}"
                title = ";".join([f"{material}: {sample}", *fields, "Replay events"])
                dimensions = [value for field in fields for value in self.axes[field]]
                histogram = (self.root.TH1D if len(fields) == 1 else self.root.TH2D)(hist_name, title, *dimensions)
                histogram.SetDirectory(0)
                histogram.Sumw2()
                histograms[name, sample] = histogram
        self.histograms[material] = histograms

    def fill(self, event, entry, selected_entry, detail):
        values = event_values(event)
        accepted = int(detail == 5)
        values.update(
            source_entry=entry,
            selected_entry=selected_entry,
            accepted=accepted,
            selection_stage=0 if detail == 0 else (2 if accepted else 1),
            selection_detail=detail,
            source_uuid=self.summary["source_uuid"],
        )
        self.events.fill(values)
        muon_entry = values["muDIS_muEntry"]
        key = (muon_entry, -1 if muon_entry >= 0 else entry)
        if key not in self.muons:
            self.muons[key] = dict(
                {name: values[name] for name in (*self.muon_fields, "muon_pdg", "muDIS_muEntry", "source_uuid")},
                first_source_entry=entry,
                events_processed=0,
                events_selected=0,
                identity_valid=int(muon_entry >= 0),
            )
        muon = self.muons[key]
        muon["events_processed"] += 1
        muon["events_selected"] += accepted
        material = values["muDIS_material"]
        # Material is a directory label, never a path supplied by the input.
        if not material or not all(c.isalnum() or c == "_" for c in material) or material == "ALL":
            material = "UNKNOWN"
        self.book(material)
        for group in ("ALL", material):
            for name, fields in self.specs.items():
                coordinates = [values[field] for field in fields]
                if not all(math.isfinite(v) for v in coordinates):
                    continue
                for sample in ("all", "accepted" if accepted else "rejected"):
                    self.histograms[group][name, sample].Fill(*coordinates)

    def write(self, destination):
        self.events.write()
        muons = ScalarTree(
            self.root,
            "replayFilterMuons",
            (
                "muDIS_muEntry",
                "identity_valid",
                "first_source_entry",
                "events_processed",
                "events_selected",
                "any_accepted",
                "muon_pdg",
            ),
            self.muon_fields + ("accepted_fraction",),
            ("source_uuid",),
        )
        for values in self.muons.values():
            muons.fill(
                dict(
                    values,
                    any_accepted=int(values["events_selected"] > 0),
                    accepted_fraction=values["events_selected"] / values["events_processed"],
                )
            )
        muons.write()
        directory = destination.mkdir("replayFilterPlots")
        for material, histograms in self.histograms.items():
            directory.mkdir(material).cd()
            for histogram in histograms.values():
                if histogram.Write() <= 0:
                    raise OSError("Cannot write replay diagnostic histogram")
            for name, sample in histograms:
                if sample != "all":
                    continue
                efficiency = histograms[name, "accepted"].Clone(f"{name}_efficiency")
                efficiency.SetDirectory(0)
                efficiency.SetTitle(
                    efficiency.GetTitle().replace("accepted", "efficiency").replace("Replay events", "Accepted / all")
                )
                efficiency.Divide(histograms[name, "accepted"], histograms[name, "all"], 1, 1, "B")
                if efficiency.Write() <= 0:
                    raise OSError("Cannot write replay diagnostic efficiency")
        destination.cd()

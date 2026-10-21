#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration
"""Plot simulation coordinates using only Python and PyROOT (no FairShip setup).

Example:
    python3 validateLayout.py -f sim.root -g geo.root -o layout
    python3 validateLayout.py -g geo.root --geometry-only -o layout

Also writes geometry_dimensions.md and geometry_dimensions.csv: main volume
envelopes from the target through the ECAL scoring plane, and cavern dimensions.

Writes one PNG per detector point collection and one for primary muon starts,
plus a multipage PDF and a ROOT file containing the histograms and canvases.
Each canvas shows a 1D Z histogram and a 2D Y-versus-X histogram.
Veto points also have X-versus-Z and Y-versus-Z panels. SST Y-versus-X
panels are separated into the four tracking stations using detector IDs.
Primary muon ranges are X=[-300,300], Y=[-400,400], Z=[3199,3201] cm,
covering a centred 600 x 800 cm gun plane at Z=32 m.
Reads split data members, so FairShip dictionaries are not needed. Unsplit
collections cannot be read this way and produce an explicit error.

Geometry annotations and dashed boundaries are global axis-aligned bounding boxes,
not sensitive-area boundaries. Dimensions are full widths, in cm. Override a
geometry association with --volume 'branch=volume1,volume2' (repeatable; exact
TGeoVolume names). The associated geometry paths are printed to the terminal.
"""

import argparse
import csv
import itertools
import re
import sys
from array import array
from contextlib import ExitStack
from functools import cache
from pathlib import Path

import ROOT

# Detector envelopes, rather than individual bars, straws or silicon sensors.
DETECTOR_VOLUMES = {
    "strawtubesPoint": ("Tr1", "Tr2", "Tr3", "Tr4"),
    "TimeDetPoint": ("Timing Detector",),
    "UpstreamTaggerPoint": ("Upstream_Tagger",),
    "vetoPoint": ("DecayVolume",),
    "CaloScoringPlanePoint": ("CaloScoringPlane",),
    "TargetPoint": ("volTarget",),
    "TTPoint": ("TT",),
    "SiliconTargetPoint": ("SiliconTarget",),
    "MTCDetPoint": ("MTC",),
}

MC_START_RANGES = ((-300.0, 300.0), (-400.0, 400.0), (3199.0, 3201.0))

CAVERN_PARTS = {
    "muon_shield_cavern": "TCC8 tunnel interior",
    "experiment_cavern": "ECN3 hall interior",
    "stair_step": "Intermediate step interior",
    "yoke_pit": "Spectrometer pit interior",
}


def shape_bounds(shape, matrix):
    """Transform all eight local bounding-box corners to global coordinates."""
    shape.ComputeBBox()
    origin = shape.GetOrigin()
    half = (shape.GetDX(), shape.GetDY(), shape.GetDZ())
    corners = []
    for signs in itertools.product((-1, 1), repeat=3):
        local = array("d", (origin[i] + signs[i] * half[i] for i in range(3)))
        world = array("d", [0, 0, 0])
        matrix.LocalToMaster(local, world)
        corners.append(tuple(world))
    return tuple((min(p[i] for p in corners), max(p[i] for p in corners)) for i in range(3))


def cavern_parts(shape, matrix, path):
    """Read the named cavern cut-outs from the persisted Boolean shape tree."""
    if shape.GetName() in CAVERN_PARTS:
        yield CAVERN_PARTS[shape.GetName()], path, shape_bounds(shape, matrix)
    if shape.InheritsFrom("TGeoCompositeShape"):
        boolean = shape.GetBoolNode()
        for child, placement in (
            (boolean.GetLeftShape(), boolean.GetLeftMatrix()),
            (boolean.GetRightShape(), boolean.GetRightMatrix()),
        ):
            transform = ROOT.TGeoHMatrix(matrix)
            transform.Multiply(placement)
            yield from cavern_parts(child, transform, path + "/" + child.GetName())


def split_fields(branch):
    """Collect member names, including members of split base classes."""
    fields = {leaf.GetName().rsplit(".", 1)[-1] for leaf in branch.GetListOfLeaves()}
    for child in branch.GetListOfBranches():
        fields.update(split_fields(child))
    return fields


def geometry_boxes(geo, wanted, include_cavern=False):
    """Find every placement of the requested volumes, in global coordinates."""
    volumes = {v.GetNumber(): v for v in geo.GetListOfVolumes()}

    @cache
    def contains_match(number):
        volume = volumes[number]
        return volume.GetName() in wanted or any(
            contains_match(volume.GetNode(i).GetVolume().GetNumber()) for i in range(volume.GetNdaughters())
        )

    boxes = {name: [] for name in wanted}

    def visit(node, parent_matrix, parent_path):
        volume = node.GetVolume()
        if not contains_match(volume.GetNumber()):
            return
        matrix = ROOT.TGeoHMatrix(parent_matrix)
        matrix.Multiply(node.GetMatrix())
        path = parent_path + "/" + node.GetName()
        name = volume.GetName()
        if name in wanted:
            shape = volume.GetShape()
            boxes[name].append((path, shape_bounds(shape, matrix)))
            if include_cavern and name == "Cavern":
                for label, part_path, bounds in cavern_parts(shape, matrix, path):
                    boxes.setdefault(label, []).append((part_path, bounds))
        for i in range(volume.GetNdaughters()):
            child = volume.GetNode(i)
            if contains_match(child.GetVolume().GetNumber()):
                visit(child, matrix, path)

    visit(geo.GetTopNode(), ROOT.TGeoHMatrix(), "")
    return boxes


def write_geometry_table(geo, output_dir):
    """Summarise main envelopes, excluding repetitive detector internals."""
    top = geo.GetTopVolume()
    wanted = {top.GetNode(i).GetVolume().GetName() for i in range(top.GetNdaughters())}
    wanted.update(volume for names in DETECTOR_VOLUMES.values() for volume in names)
    wanted.update(
        (
            "TargetArea",
            "target_vacuum_box",
            "proximity_shielding",
            "AbsorberVol",
            "MuonShieldArea",
            "SHiPMagnet",
            "Cavern",
        )
    )
    boxes = geometry_boxes(geo, wanted, include_cavern=True)
    target = [bounds for name in ("TargetArea", "target_vacuum_box") for _, bounds in boxes[name]]
    calo = [bounds for _, bounds in boxes["CaloScoringPlane"]]
    start = min(b[2][0] for b in target) if target else float("-inf")
    end = max(b[2][1] for b in calo) if calo else float("inf")
    notes = [
        "All values are in cm. Dimensions are full global axis-aligned bounding-box widths, not half-lengths.",
        "Main envelopes are sorted by Z centre; nested envelopes can overlap. Detector internals are omitted.",
        "Cavern interior rows describe the modelled cut-outs; Cavern itself is the outer concrete envelope.",
    ]
    if not target or not calo:
        missing = [label for label, found in (("target", target), ("ECAL scoring plane", calo)) if not found]
        notes.append(f"Missing {' and '.join(missing)}: the corresponding Z cutoff is not applied.")
    if not boxes["Cavern"]:
        notes.append("No Cavern volume found in this geometry.")
    elif not any(boxes.get(label) for label in CAVERN_PARTS.values()):
        notes.append("No recognised cavern cut-out shapes found; only the outer concrete bounds are available.")
    records = []
    for name, placements in boxes.items():
        section = "Cavern" if name == "Cavern" or name in CAVERN_PARTS.values() else "Beamline"
        for path, bounds in placements:
            if section == "Beamline" and (bounds[2][1] < start or bounds[2][0] > end):
                continue
            label = path.rsplit("/", 1)[-1] if section == "Beamline" else name
            records.append((section, label, path, bounds))
    records.sort(key=lambda row: (row[0] == "Cavern", sum(row[3][2]) / 2, row[1], row[2]))
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "geometry_dimensions.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "section",
                "volume",
                "geometry_path",
                "width_x_cm",
                "height_y_cm",
                "length_z_cm",
                "x_min_cm",
                "x_max_cm",
                "y_min_cm",
                "y_max_cm",
                "z_min_cm",
                "z_max_cm",
                "z_centre_cm",
            )
        )
        for section, label, path, bounds in records:
            writer.writerow(
                (
                    section,
                    label,
                    path,
                    *(hi - lo for lo, hi in bounds),
                    *(edge for pair in bounds for edge in pair),
                    sum(bounds[2]) / 2,
                )
            )
    lines = ["# Geometry dimensions", "", *notes, ""]
    for section in ("Beamline", "Cavern"):
        lines.extend(
            (
                f"## {section}",
                "",
                "| Volume | Width X | Height Y | Length Z | X bounds | Y bounds | Z bounds | Z centre |",
                "| :--- | ---: | ---: | ---: | :--- | :--- | :--- | ---: |",
            )
        )
        for group, label, _, bounds in records:
            if group != section:
                continue
            cells = [
                label,
                *(f"{hi - lo:.6g}" for lo, hi in bounds),
                *(f"[{lo:.6g}, {hi:.6g}]" for lo, hi in bounds),
                f"{sum(bounds[2]) / 2:.6g}",
            ]
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    lines.append("Exact geometry paths and unrounded values are included in geometry_dimensions.csv.\n")
    markdown = "\n".join(lines)
    (output_dir / "geometry_dimensions.md").write_text(markdown)
    print(markdown)


def draw_checked(tree, expression, selection, entries):
    count = tree.Draw(expression, selection, "goff", entries)
    if count < 0:
        raise RuntimeError(f"Cannot read split fields: {expression}")
    return count


def axis_range(tree, expression, selection, entries, boxes, axis):
    """Use ROOT's automatic data range, expanded to include the geometry."""
    count = draw_checked(tree, f"{expression}>>layout_range(100)", selection, entries)
    hist = ROOT.gDirectory.Get("layout_range")
    limits = [edge for _, bounds in boxes for edge in bounds[axis]]
    if count:
        limits.extend((hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax()))
    hist.Delete()
    low, high = (min(limits), max(limits)) if limits else (-1.0, 1.0)
    margin = max((high - low) * 0.05, 0.01)
    return low - margin, high + margin


def box_description(path, bounds):
    x, y, z = bounds
    return (
        f"{path.rsplit('/', 1)[-1]}: "
        f"size XYZ=({x[1] - x[0]:.4g}, {y[1] - y[0]:.4g}, {z[1] - z[0]:.4g}) cm; "
        f"X=[{x[0]:.4g}, {x[1]:.4g}], Y=[{y[0]:.4g}, {y[1]:.4g}], "
        f"Z=[{z[0]:.4g}, {z[1]:.4g}] cm; z centre={(z[0] + z[1]) / 2:.4g} cm"
    )


def plot_specs(branch, selection, boxes):
    """Histogram suffix, horizontal/vertical axes, selection and envelopes."""
    specs = [("Z", 2, None, selection, boxes)]
    if branch == "strawtubesPoint":
        for station in range(1, 5):
            # StrawDecode: station = detector ID / 1e6 (integer division).
            cut = f"{branch}.fDetectorID>={station * 1000000} && {branch}.fDetectorID<{(station + 1) * 1000000}"
            station_boxes = [box for box in boxes if box[0].rsplit("/", 1)[-1].startswith(f"Tr{station}_")]
            specs.append((f"YX_station{station}", 0, 1, cut, station_boxes))
    else:
        specs.append(("YX", 0, 1, selection, boxes))
    if branch == "vetoPoint":
        specs.extend((("XZ", 2, 0, selection, boxes), ("YZ", 2, 1, selection, boxes)))
    return specs


def plot_collection(tree, branch, fields, selection, boxes, notes, args, output, pdf):
    ROOT.gROOT.cd()
    coordinates = [f"{branch}.{field}" for field in fields]
    ranges = (
        MC_START_RANGES
        if branch == "MCTrack"
        else [axis_range(tree, expr, selection, args.events, boxes, i) for i, expr in enumerate(coordinates)]
    )
    label = "MCTrack_primary_muons" if branch == "MCTrack" else branch
    name = re.sub(r"[^A-Za-z0-9_]", "_", label)
    specs = plot_specs(branch, selection, boxes)
    rows = (len(specs) + 1) // 2
    plot_height = 650 if rows == 1 else 500 * rows
    descriptions = [box_description(path, bounds) for path, bounds in boxes]
    header = [label, *notes, *descriptions]
    header_height = 35 + 24 * len(header)
    height = plot_height + header_height
    canvas = ROOT.TCanvas(name, label, 1600, height)
    panels = ROOT.TPad(f"{name}_panels", "", 0, 0, 1, plot_height / height)
    panels.Draw()
    panels.Divide(2, rows)
    # Keep objects alive until the canvas has been written and exported.
    objects = [panels]
    histograms = []
    counts = []
    for panel, (suffix, horizontal, vertical, cut, panel_boxes) in enumerate(specs, start=1):
        pad = panels.cd(panel)
        pad.SetRightMargin(0.15)
        pad.SetLeftMargin(0.13)
        pad.SetBottomMargin(0.12)
        title = f"{label} station {suffix[-1]}" if suffix.startswith("YX_station") else label
        if vertical is None:
            hist = ROOT.TH1D(f"{name}_Z", f"{label};{fields[2]} [cm];Points", args.bins, *ranges[2])
            expression = coordinates[2]
        else:
            hist = ROOT.TH2D(
                f"{name}_{suffix}",
                f"{title};{fields[horizontal]} [cm];{fields[vertical]} [cm];Points",
                args.bins,
                *ranges[horizontal],
                args.bins,
                *ranges[vertical],
            )
            expression = f"{coordinates[vertical]}:{coordinates[horizontal]}"
        count = draw_checked(tree, f"{expression}>>{hist.GetName()}", cut, args.events)
        counts.append(count)
        hist.SetStats(False)
        if vertical is None:
            hist.SetMinimum(0)
            hist.SetMaximum(max(1.0, 1.1 * hist.GetMaximum()))
            hist.Draw("HIST")
        else:
            hist.Draw("COLZ")
        histograms.append(hist)
        for _, bounds in panel_boxes:
            if vertical is None:
                outlines = [ROOT.TLine(z, 0, z, hist.GetMaximum()) for z in bounds[2]]
            else:
                outline = ROOT.TBox(
                    bounds[horizontal][0], bounds[vertical][0], bounds[horizontal][1], bounds[vertical][1]
                )
                outline.SetFillStyle(0)
                outlines = [outline]
            for outline in outlines:
                outline.SetLineColor(ROOT.kMagenta + 2)
                outline.SetLineStyle(2)
                outline.SetLineWidth(2)
                outline.Draw("same")
                objects.append(outline)
        if count == 0:
            empty = ROOT.TLatex()
            empty.SetNDC()
            empty.DrawLatex(0.3, 0.5, "No selected points")
            objects.append(empty)
    canvas.cd()
    latex = ROOT.TLatex()
    latex.SetNDC()
    # Font precision 1 displays literal underscores in branch and node names.
    latex.SetTextFont(41)
    latex.SetTextSize(18 / height)
    header[0] += f"  ({counts[0]} points, {args.events} events)"
    for i, line in enumerate(header):
        latex.DrawLatex(0.02, 1 - (26 + 24 * i) / height, line)
    canvas.Update()
    canvas.SaveAs(str(args.output / f"{name}.png"))
    canvas.Print(str(pdf))
    output.cd()
    for hist in histograms:
        hist.Write()
        hist.SetDirectory(0)
    canvas.Write()
    canvas.Close()
    print(f"{label}: {counts[0]} points in {args.events} events")
    if branch == "strawtubesPoint":
        print("  Station counts: " + ", ".join(f"Tr{i}={n}" for i, n in enumerate(counts[1:], start=1)))
        if sum(counts[1:]) != counts[0]:
            print("  WARNING: Some detector IDs are outside stations 1-4; those points appear only in Z.")
    for path, bounds in boxes:
        print(f"  {path}\n    {box_description(path, bounds)}")
    for note in notes:
        print(f"  {note}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "-f", "--input", help="Simulation ROOT file containing cbmsim (required unless --geometry-only)"
    )
    parser.add_argument("-g", "--geometry", required=True, help="Geometry ROOT file containing FAIRGeom")
    parser.add_argument(
        "--geometry-only", action="store_true", help="Write geometry tables without reading simulation data"
    )
    parser.add_argument("-o", "--output", type=Path, default=Path("layout_plots"), help="Output directory")
    parser.add_argument("-n", "--events", type=int, default=-1, help="First N events; -1 means all (default)")
    parser.add_argument("--bins", type=int, default=150, help="Bins per axis (default: 150)")
    parser.add_argument(
        "--volume",
        action="append",
        default=[],
        metavar="BRANCH=VOLUME[,VOLUME]",
        help="Override the detector envelope volume names for a branch",
    )
    args = parser.parse_args()
    if not args.geometry_only and not args.input:
        parser.error("-f/--input is required unless --geometry-only is used")
    if args.bins <= 0 or args.events < -1:
        parser.error("--bins must be positive and --events must be -1 or nonnegative")
    overrides = {}
    for value in args.volume:
        branch, separator, names = value.partition("=")
        if not separator or not branch or not all(names.split(",")):
            parser.error("--volume requires BRANCH=VOLUME[,VOLUME]")
        overrides[branch] = tuple(names.split(","))

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    with ExitStack() as stack:
        geometry = stack.enter_context(ROOT.TFile.Open(args.geometry))
        geo = geometry.Get("FAIRGeom")
        if not geo or not geo.InheritsFrom("TGeoManager"):
            raise RuntimeError("Geometry file does not contain a FAIRGeom TGeoManager")
        write_geometry_table(geo, args.output)
        if args.geometry_only:
            print(f"Geometry tables saved in {args.output.resolve()}")
            return 0
        simulation = stack.enter_context(ROOT.TFile.Open(args.input))
        tree = simulation.Get("cbmsim")
        if not tree or not tree.InheritsFrom("TTree"):
            raise RuntimeError("Simulation file does not contain a cbmsim TTree")
        branches = {b.GetName().rstrip("."): split_fields(b) for b in tree.GetListOfBranches()}
        points = []
        for name, members in branches.items():
            # The event header has vertex coordinates but is not a detector collection.
            if name == "MCEventHeader":
                continue
            if {"fX", "fY", "fZ"} <= members:
                points.append(name)
            elif "point" in name.lower():
                raise RuntimeError(f"{name} has no split fX/fY/fZ fields; regenerate the file with branch splitting")
        if not {"fStartX", "fStartY", "fStartZ", "fPdgCode", "fMotherId"} <= branches.get("MCTrack", set()):
            raise RuntimeError("MCTrack is missing the split start-position, PDG or mother-ID fields")
        if "strawtubesPoint" in points and "fDetectorID" not in branches["strawtubesPoint"]:
            raise RuntimeError("strawtubesPoint needs a split fDetectorID field to separate the four stations")
        if unknown := overrides.keys() - set(points):
            parser.error(f"Unknown point branches in --volume: {', '.join(sorted(unknown))}")
        associations = {
            name: overrides.get(name, DETECTOR_VOLUMES.get(name, (re.sub(r"Points?$", "", name),))) for name in points
        }
        boxes = geometry_boxes(geo, {v for names in associations.values() for v in names})
        args.events = tree.GetEntries() if args.events == -1 else min(args.events, tree.GetEntries())
        args.output.mkdir(parents=True, exist_ok=True)
        pdf = args.output / "validateLayout.pdf"
        print(f"Reading {args.events} events; found point collections: {', '.join(points) or '(none)'}")
        with ROOT.TFile.Open(str(args.output / "validateLayout.root"), "RECREATE") as output:
            book = ROOT.TCanvas("layout_book", "")
            book.Print(str(pdf) + "[")
            try:
                for name in points:
                    detector_boxes = [box for volume in associations[name] for box in boxes[volume]]
                    missing = [v for v in associations[name] if not boxes[v]]
                    notes = ["Dashed lines: envelope limits in the plotted coordinates; all dimensions in cm."]
                    if missing:
                        notes.append(f"Geometry not found: {', '.join(missing)}. Use --volume to set the association.")
                    plot_collection(tree, name, ("fX", "fY", "fZ"), "", detector_boxes, notes, args, output, pdf)
                plot_collection(
                    tree,
                    "MCTrack",
                    ("fStartX", "fStartY", "fStartZ"),
                    "abs(MCTrack.fPdgCode)==13 && MCTrack.fMotherId<0",
                    [],
                    [
                        "Input muons: |PDG| = 13, mother ID < 0 (all primary muons).",
                        "Gun plane: X=[-300, 300] cm, Y=[-400, 400] cm; Z=3200 cm (32 m).",
                    ],
                    args,
                    output,
                    pdf,
                )
            finally:
                book.Print(str(pdf) + "]")
                book.Close()
    print(f"Plots saved in {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)

#!/usr/bin/env python
"""Plot filterEvents.py histograms in a multipage PDF using YAML settings."""

import argparse
import math
from pathlib import Path

import ROOT
import yaml

ROOT.gROOT.SetBatch(True)

DEFAULTS = {
    "x_range": None, "y_range": None, "z_range": None,
    "log_x": False, "log_y": False, "log_z": False,
    "rebin": 1, "enabled": True,
}


def histogram_paths(directory, prefix=""):
    """Discover latest histogram cycles without keeping all histograms in memory."""
    paths = []
    for name in sorted({key.GetName() for key in directory.GetListOfKeys()}):
        key = directory.GetKey(name)
        cls = ROOT.TClass.GetClass(key.GetClassName())
        path = f"{prefix}{name}"
        if cls and cls.InheritsFrom("TDirectory"):
            paths.extend(histogram_paths(directory.GetDirectory(name), path + "/"))
        elif cls and cls.InheritsFrom("TH1") and not cls.InheritsFrom("TH3"):
            if name != "filter_counts":
                paths.append(path)
    return paths


def read_histogram(source, path):
    histogram = source.Get(path)
    if not histogram or not histogram.InheritsFrom("TH1"):
        raise ValueError(f"Histogram not found: {path}")
    histogram.SetDirectory(0)
    ROOT.SetOwnership(histogram, True)
    return histogram


def settings_for(config, path):
    settings = DEFAULTS | config.get("defaults", {}) | config.get("histograms", {}).get(path, {})
    unknown = settings.keys() - DEFAULTS.keys()
    if unknown:
        raise ValueError(f"{path}: unknown settings: {', '.join(sorted(unknown))}")
    for key in ("enabled", "log_x", "log_y", "log_z"):
        if not isinstance(settings[key], bool):
            raise ValueError(f"{path}: {key} must be true or false")
    for axis in "xyz":
        value = settings[f"{axis}_range"]
        if value is not None:
            if (not isinstance(value, list) or len(value) != 2
                    or any(isinstance(v, bool) or not isinstance(v, (int, float))
                           or not math.isfinite(v) for v in value)
                    or value[0] >= value[1]):
                raise ValueError(f"{path}: {axis}_range must be null or [minimum, maximum]")
            if settings[f"log_{axis}"] and value[0] <= 0:
                raise ValueError(f"{path}: a logarithmic {axis} range must start above zero")
    return settings


def rebin_factors(histogram, value):
    factors = value if isinstance(value, list) else [value] * histogram.GetDimension()
    if (len(factors) != histogram.GetDimension()
            or any(type(factor) is not int or factor < 1 for factor in factors)):
        raise ValueError("rebin must be a positive integer, or [x_factor, y_factor] for 2D")
    axes = [histogram.GetXaxis(), histogram.GetYaxis()]
    for axis, factor in zip(axes, factors):
        if axis.GetNbins() % factor:
            raise ValueError("Rebin factors must divide the bin counts to preserve the full range")
        if factor != 1 and any(axis.GetBinLabel(i) for i in range(1, axis.GetNbins() + 1)):
            raise ValueError("Rebinning labelled category axes would merge different species")
    return factors


def validate_plot(histogram, settings):
    rebin_factors(histogram, settings["rebin"])
    for axis, getter in (("x", histogram.GetXaxis), ("y", histogram.GetYaxis)):
        if axis == "y" and histogram.GetDimension() == 1:
            continue
        if settings[f"log_{axis}"]:
            limits = settings[f"{axis}_range"]
            low = limits[0] if limits else getter().GetXmin()
            if low <= 0:
                raise ValueError(f"log_{axis} requires an explicit positive {axis}_range for this histogram")
    if histogram.GetDimension() == 1 and (settings["log_z"] or settings["z_range"] is not None):
        raise ValueError("z_range/log_z apply only to 2D histograms")


def draw_histogram(canvas, histogram, settings, path):
    dimension = histogram.GetDimension()
    factors = rebin_factors(histogram, settings["rebin"])
    if dimension == 2:
        histogram.Rebin2D(*factors)
    elif factors[0] != 1:
        histogram.Rebin(factors[0])
    canvas.Clear()
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.22 if dimension == 2 else 0.06)
    labelled = any(histogram.GetXaxis().GetBinLabel(i) for i in range(1, histogram.GetNbinsX() + 1))
    canvas.SetBottomMargin(0.20 if labelled else 0.14)
    canvas.SetTopMargin(0.09)
    canvas.SetLogx(settings["log_x"])
    canvas.SetLogy(settings["log_y"])
    canvas.SetLogz(settings["log_z"])
    histogram.SetStats(False)
    histogram.SetTitle(f"{path}: {histogram.GetTitle()}")
    histogram.GetXaxis().SetTitleOffset(1.2)
    histogram.GetYaxis().SetTitleOffset(1.45)
    if settings["x_range"] is not None:
        histogram.GetXaxis().SetRangeUser(*settings["x_range"])
    if dimension == 2:
        if settings["y_range"] is not None:
            histogram.GetYaxis().SetRangeUser(*settings["y_range"])
        if not histogram.GetZaxis().GetTitle():
            histogram.GetZaxis().SetTitle("Entries")
        histogram.GetZaxis().SetTitleOffset(1.35)
        content_axis = "z"
    else:
        content_axis = "y"
    limits = settings[f"{content_axis}_range"]
    if limits is not None:
        histogram.SetMinimum(limits[0])
        histogram.SetMaximum(limits[1])
    # Give empty histograms a drawable positive scale on log-count axes.
    if settings[f"log_{content_axis}"] and limits is None and histogram.GetMaximum() <= 0:
        histogram.SetMinimum(0.1)
        histogram.SetMaximum(1.)
    histogram.Draw("COLZ" if dimension == 2 else "HIST")
    canvas.Update()
    if dimension == 2:
        palette = histogram.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.SetX1NDC(0.80)
            palette.SetX2NDC(0.84)
    canvas.Modified()
    canvas.Update()


def write_default_config(source, paths, destination):
    histograms = {}
    for path in paths:
        histogram = read_histogram(source, path)
        limits = {"x_range": [float(histogram.GetXaxis().GetXmin()), float(histogram.GetXaxis().GetXmax())]}
        if histogram.GetDimension() == 2:
            limits["y_range"] = [float(histogram.GetYaxis().GetXmin()), float(histogram.GetYaxis().GetXmax())]
        histograms[path] = limits
    config = {"defaults": DEFAULTS.copy(), "histograms": histograms}
    with destination.open("x") as output:
        output.write("# null uses the full stored axis range, or automatic limits for bin contents.\n")
        yaml.safe_dump(config, output, sort_keys=False)


def latex_escape(value):
    replacements = {"\\": r"\textbackslash{}", "_": r"\_", "%": r"\%", "&": r"\&",
                    "#": r"\#", "$": r"\$", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
                    "^": r"\textasciicircum{}"}
    return "".join(replacements.get(char, char) for char in value)


def write_table(source, paths, destination):
    materials = sorted({path.rsplit("/", 1)[0] for path in paths if "/" in path})
    rows = []
    for material in materials:
        summary = source.Get(f"{material}/filter_counts")
        counts = {}
        if summary:
            counts = {summary.GetXaxis().GetBinLabel(i): summary.GetBinContent(i)
                      for i in range(1, summary.GetNbinsX() + 1)}
        else:
            # Old files have exact raw counts, but no joint multiplicity/weight sums.
            for sample, suffix in (("processed", ""), ("selected", "filtered")):
                for kind, name in (("muons", "muon_p_"), ("dis", "n_daughters_")):
                    histogram = source.Get(f"{material}/{name}{suffix}")
                    if histogram:
                        counts[f"{kind}_{sample}_raw"] = histogram.GetEntries()
        if not counts:
            continue
        for kind, label in (("muons", "Muons"), ("dis", "DIS")):
            values = []
            for sample in ("processed", "selected"):
                for weighting in ("raw", "weighted"):
                    value = counts.get(f"{kind}_{sample}_{weighting}")
                    values.append("---" if value is None else
                                  (f"{value:.0f}" if weighting == "raw" else f"{value:.8g}"))
            rows.append(" & ".join([latex_escape(material), label, *values]) + r" \\")
    if not rows:
        raise ValueError("No per-material filter counts or validation histograms found")
    lines = [r"\begin{tabular}{llrrrr}", r"\hline",
             r"Material & Sample & Processed raw & Processed weighted & Selected raw & Selected weighted \\",
             r"\hline", *rows, r"\hline", r"\end{tabular}",
             "% Processed counts exclude unreadable/malformed entries. Selected muons are counted per material.",
             "% Muon weights sum wDIS; DIS weights sum nDIS * wDIS. --- means unavailable in older files."]
    destination.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--inputfile", required=True, type=Path)
    parser.add_argument("-c", "--config", type=Path, default=Path(__file__).with_name("plotFilterEvents.yaml"))
    parser.add_argument("-o", "--output", type=Path, help="Multipage PDF (default: input filename with .pdf)")
    parser.add_argument("--latex", type=Path, help="Also write a LaTeX count table")
    parser.add_argument("--write-config", type=Path, help="Write YAML with actual histogram ranges and exit")
    args = parser.parse_args()
    output = args.output or args.inputfile.with_suffix(".pdf")
    destinations = [args.write_config] if args.write_config else [output, args.latex]
    resolved = [path.resolve() for path in destinations if path is not None]
    if (len(set(resolved)) != len(resolved)
            or any(path in (args.inputfile.resolve(), args.config.resolve()) for path in resolved)):
        parser.error("Output paths must differ from each other and from the input/config files")
    if not args.write_config and output.suffix.lower() != ".pdf":
        parser.error("The plot output must have a .pdf extension")
    source = ROOT.TFile.Open(str(args.inputfile), "READ")
    if not source or source.IsZombie():
        parser.error(f"Cannot open {args.inputfile}")
    try:
        paths = histogram_paths(source)
        if not paths:
            raise ValueError("No 1D or 2D histograms found")
        if args.write_config:
            write_default_config(source, paths, args.write_config)
            print(f"Wrote settings for {len(paths)} histograms to {args.write_config}")
            return
        config = yaml.safe_load(args.config.read_text()) or {}
        if not isinstance(config, dict) or config.keys() - {"defaults", "histograms"}:
            raise ValueError("YAML must contain only defaults and histograms mappings")
        if not isinstance(config.get("defaults", {}), dict) or not isinstance(config.get("histograms", {}), dict):
            raise ValueError("defaults and histograms must be mappings")
        if any(not isinstance(value, dict) for value in config.get("histograms", {}).values()):
            raise ValueError("Each histogram override must be a mapping")
        unknown = config.get("histograms", {}).keys() - set(paths)
        if unknown:
            raise ValueError(f"Unknown histogram paths: {', '.join(sorted(unknown))}")
        plots = []
        for path in paths:
            settings = settings_for(config, path)
            if settings["enabled"]:
                histogram = read_histogram(source, path)
                try:
                    validate_plot(histogram, settings)
                except ValueError as error:
                    raise ValueError(f"{path}: {error}") from error
                plots.append((path, settings))
        if not plots:
            raise ValueError("No histograms enabled")
        if args.latex:
            write_table(source, paths, args.latex)
        canvas = ROOT.TCanvas("filter_plots", "Filter validation", 1000, 800)
        canvas.Print(str(output) + "[")
        try:
            for path, settings in plots:
                histogram = read_histogram(source, path)
                draw_histogram(canvas, histogram, settings, path)
                canvas.Print(str(output), f"Title:{path}")
                canvas.Clear()
        finally:
            canvas.Print(str(output) + "]")
            canvas.Close()
        print(f"Saved {len(plots)} pages to {output}")
        if args.latex:
            print(f"Saved count table to {args.latex}")
    except (OSError, ValueError, yaml.YAMLError) as error:
        parser.error(str(error))
    finally:
        source.Close()


if __name__ == "__main__":
    main()

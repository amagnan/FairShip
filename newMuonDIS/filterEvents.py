#!/usr/bin/env python
"""Filter prepareEvents.py output, retaining selected DIS interactions."""

import argparse
import time
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-f", "--inputfile", nargs="+", required=True,
        help="Input file(s) or directories to search recursively for .root files",
    )
    parser.add_argument("-o", "--outputfile", required=True, help="New ROOT file (must not exist)")
    parser.add_argument("-n", "--n_events", type=int, default=-1)
    parser.add_argument("-s", "--start_event", type=int, default=0)
    parser.add_argument("--min-charged", type=int, default=2)
    parser.add_argument("--exclude-muons", action="store_true")
    args = parser.parse_args()
    if args.min_charged < 0 or args.start_event < 0 or args.n_events < -1:
        parser.error("Counts must be nonnegative, except -n -1 for all entries")
    file_names = []
    seen = set()
    output_path = Path(args.outputfile).resolve()
    for name in args.inputfile:
        path = Path(name).expanduser()
        if path.is_dir():
            candidates = sorted(p for p in path.rglob("*.root") if p.is_file())
        elif path.is_file():
            candidates = [path]
        else:
            parser.error(f"Input does not exist or is not a file/directory: {name}")
        for candidate in candidates:
            resolved = candidate.resolve()
            if resolved == output_path:
                parser.error(f"Output file is also an input: {candidate}")
            if resolved not in seen:
                seen.add(resolved)
                file_names.append(str(candidate))
    if not file_names:
        parser.error("No ROOT input files found")
    print(f"Found {len(file_names)} input file(s).", flush=True)

    # Measure library initialization and filtering, excluding input discovery.
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    if ROOT.gSystem.Load("libShipMuDIS.so") < 0:
        raise RuntimeError("Cannot load libShipMuDIS.so")
    selection = ROOT.MuDISFilter()
    selection.init(args.n_events, args.start_event)
    selection.SetMinChargedDaughters(args.min_charged)
    selection.SetIncludeMuons(not args.exclude_muons)
    inputs = ROOT.std.vector("string")()
    for name in file_names:
        inputs.push_back(name)
    selection.process_file(inputs, args.outputfile)
    wall_seconds = time.perf_counter() - wall_start
    cpu_seconds = time.process_time() - cpu_start
    cpu_percent = 100.0 * cpu_seconds / wall_seconds if wall_seconds else 0.0
    print(
        f"Performance (initialization + filtering): wall time {wall_seconds:.3f} s, "
        f"CPU time {cpu_seconds:.3f} s, average CPU utilization {cpu_percent:.1f}%"
    )


if __name__ == "__main__":
    main()

# Replay event filter

`replayFilter.py` writes a smaller replay ROOT file for subsequent analysis.
Run it **after** the incoming-muon hit restoration performed by
`MuDISGenerator_postProcessing` (automatically called by `run_simScript.py`).

```sh
pixi cvmfs python -u newMuonDIS/replayFilter.py \
    -f production/sim_example_disreplay_gen.root \
    -o selected/sim_example_disreplay_gen.root
```

Create the output directory first. Each invocation processes one input file;
the output file must not already exist. Both replay types are supported.
Geometry defaults to the matching `geo_` file beside the input; use
`-g /path/to/geometry.root` to override it. External geometry and parameter files
are not copied. When running `analyseEvents.py` on the selected output, supply
`geometry_file` or place the matching geometry beside it. Preserve the usual
`sim_*_disreplay_gen.root` / `sim_*_disreplay_filter.root` naming for automatic
analysis discovery, and keep original and selected productions separate.

The selection reuses the helpers in `analyseEvents.py`:

1. Require at least two distinct MC tracks with hits in all four straw stations
   and the timing detector, each with momentum **strictly greater than 1 GeV/c**.
   Position and momentum come from the smallest-z Tr1 hit. This is the analysis
   original step1 requirement, now part of the filtered analysis baseline.
2. Form unordered pairs of these tracks. Keep the event if at least one pair
   has opposite charges, a valid backward straight-line vertex, **DOCA < 5 cm**,
   and its reconstructed vertex inside helium under `DecayVolume`.
   There is no 5 cm wall-clearance requirement and no target-IP cut.

The pair is used only to decide whether to keep the event. **Every branch and
every track/hit of an accepted event is copied**, including those not involved
in an accepted pair. Track indices, original event IDs, replay identity branches
and weights are unchanged. Events retain their input order; output tree entry
numbers are compacted. No pair branches or derived analysis variables are added.

The output retains the `cbmsim` schema and tree user information, including
`muDIS_normalization`. Other file objects and nested directories are copied
(latest key cycles), including the hit-restoration marker. Existing summaries
and auxiliary trees describe the original production and are copied unchanged;
normalization is not recomputed from the selected sample. A `replayFilter`
`TNamed` stores JSON with source path/UUID, geometry path, processed entry range,
selection description, and raw counts before selection, after the two-track
requirement, and after pair selection. Refiltering retains the earlier filter
record inside this metadata. An empty selection still writes an empty `cbmsim`
tree with its schema and metadata.

Optional `-s/--start_event` and `-n/--n_events` select an input entry range;
defaults are entry 0 and all remaining entries. They limit entries inspected,
not the number of accepted events. The input is opened read-only, and the output
is published only after successful filtering and writing to a temporary file.

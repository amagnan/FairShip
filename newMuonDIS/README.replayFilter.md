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

1. Require at least two distinct MC tracks with hits in `(Tr1 or Tr2)` and `(Tr3 or Tr4)`
   and the timing detector, each with momentum **strictly greater than 1 GeV/c**.
   Position and momentum come from the smallest-z hit in Tr1 or Tr2.
   Tr2 supplies the state when no Tr1 hit is present.
2. Form unordered pairs of these tracks. Keep the event if at least one pair
   has opposite charges, a valid backward straight-line vertex, **DOCA < 5 cm**,
   and its reconstructed vertex inside helium under `DecayVolume`.
   There is no 5 cm wall-clearance requirement and no target-IP cut.

Use `--cross-section analytical` to replace `muDIS_xsec` with
`0.0000599 + 0.0002839 * log(muDIS_pPythia / 2)**2` in mb (natural log,
momentum in GeV/c). Both inputs must be scalar double branches, and selected
events require finite positive Pythia momentum. The default `stored` mode
preserves the input cross section, including any legacy placeholder. The
chosen mode is recorded in `replayFilter` metadata, also for empty outputs.
Analyse analytical filter outputs with `cross_section: stored` and
`weight_factor: 1.0`; the interaction conversion is applied once during analysis.

The pair is used only to decide whether to keep the event. **Every branch and
every track/hit of an accepted event is copied**, including those not involved
in an accepted pair. Track indices, original event IDs, replay identity branches
and weights are unchanged; analytical mode replaces only `muDIS_xsec`.
Events retain their input order; output tree entry
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

## Optional acceptance snapshot

Add `--diagnostics` to capture accepted **and rejected** replay events in the
same output file. This is Python-only and uses the existing library:

```sh
pixi cvmfs python -u newMuonDIS/replayFilter.py \
    -f production/sim_example_disreplay_gen.root \
    -o selected/sim_example_disreplay_gen.root --diagnostics
```

`replayFilterDiagnostics` has one row per inspected input entry, even when no
events survive. It does not have the same entry numbering as the skimmed
`cbmsim`. `source_entry` and `source_uuid` identify the input; `selected_entry`
is the corresponding output `cbmsim` entry or `-1` for a rejection. Source path,
geometry, range, units/conventions and schema version are recorded in the
`replayFilter` JSON metadata (units also appear in branch names).

The diagnostic branches include:

- `muDIS_muEntry`, `muDIS_disIndex`, `muDIS_material`: replay identity; absent
  integer IDs are `-1`, absent material is `UNKNOWN`.
- `muon_x_cm`, `muon_y_cm`, `muon_z_cm`, `muon_px_GeV`, `muon_py_GeV`,
  `muon_pz_GeV`, `muon_energy_GeV`, `muon_pdg`, `muon_weight`: original
  `MCTrack[0]` state, plus `muon_p_GeV`, `muon_pt_GeV`, `muon_theta_rad`,
  and `muon_r_cm` (radius about the global z axis).
- `dis_x_cm`, `dis_y_cm`, `dis_z_cm`: common starting position of the direct
  DIS products identified by `analyseEvents.initial_dis_products`. This uses
  mother ID 0 and weight distinct from the input-muon weight, under the existing
  analysis tolerance. `n_dis_primary` counts these products.
  `dis_vertex_status` is 0 for missing, 1 for valid, 2 for disagreeing positions
  (absolute tolerance `1e-4 cm`), or 3 for nonfinite positions. Invalid vertices
  are NaN and excluded from vertex histograms. This is the truth vertex, not
  the reconstructed pair vertex. The classification inherits the existing
  weight convention; equal soft/DIS weights cannot be distinguished by it.
- `accepted`, `selection_stage` (the original 0/1/2), and `selection_detail`:
  0 = fewer than two selected tracks; 1 = no opposite-sign pair;
  2 = no valid backward vertex; 3 = no pair with DOCA < 5 cm;
  4 = no surviving pair inside helium; 5 = accepted. Successive cuts must
  be passed by the same pair; the furthest surviving pair sets the detail.
- For each detector `ubt` and `sbt`, `n_<detector>_muon` counts points with
  track ID 0; `n_<detector>_dis_direct` counts points from direct DIS products;
  `n_<detector>_dis` includes direct products **and all their descendants**,
  including neutral products. `n_<detector>_muon_or_dis` is the sum of muon and
  DIS counts, without double counting direct products. Soft daughters and
  companion muons/descendants are excluded from this sum and counted in
  `n_<detector>_other`. Invalid track IDs are counted in
  `n_<detector>_unassigned`. Points whose tracks were not retained (for example,
  UBT secondaries with track ID `-2`) cannot be attributed from the replay
  file and remain unassigned. These are numbers of stored points, not unique
  tracks or fired detector cells; there is no additional energy or timing cut.
  UBT uses **`UpstreamTaggerPoint` only**, never `UpstreamTaggerScoringPoint`.
  SBT uses `vetoPoint`. All counts for a missing branch are `-1`, whereas an
  available empty branch gives zero. Incoming-muon counts use the restored
  hits upstream of the DIS vertex, as supplied by the replay input.

`replayFilterMuons` groups inspected events by `muDIS_muEntry` within the source
file, across materials, saving the muon state, `events_processed`,
`events_selected`, `accepted_fraction`, `any_accepted`, and `first_source_entry`.
`source_uuid` scopes the identity. If the muon ID is missing, each event gets
its own row with `identity_valid=0`. This is a summary of the inspected range,
not necessarily all DIS events for that muon. Do not merge input-muon IDs
across production files without checking their preparation-input provenance;
the utility does not perform cross-file deduplication.

`replayFilterPlots/ALL` and one directory per material contain ROOT histograms
for all, accepted, rejected, and accepted/all replay events. Plots cover muon
position, radius, momentum, transverse momentum, angle, position maps, DIS
position/maps, momentum versus starting radius, and momentum versus DIS z.
They use raw counts; `muon_weight` is preserved for later weighted studies.
Efficiencies have binomial errors; bins with no denominator are zero, so
always consult the matching `_all` histogram. Events from the same muon are
correlated: these errors are nominal event-level errors, not independent-muon
uncertainties. Fixed ranges are 0–1000 GeV/c for momentum, 0–20 GeV/c for pT,
0–pi for angle, 0–1000 cm for radius, ±1000 cm for x/y, and -10000–20000 cm for
z. Underflow/overflow bins are retained. The diagnostic tree permits rebinning.
These efficiencies describe the replay sample entering this filter, including
any earlier selections, not the original full muon flux.

On refiltering, these three diagnostic objects are replaced when diagnostics
are enabled and omitted otherwise, so an earlier snapshot is never mistaken
for the current inspected sample. Earlier filter JSON remains in
`previous_filter`. All diagnostics are published atomically with the skim;
`cbmsim` branches, event selection and normalization are unchanged.

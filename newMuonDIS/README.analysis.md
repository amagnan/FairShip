# DIS event analysis

Run in the FairShip environment (on lxplus, `pixi cvmfs`):

```sh
pixi cvmfs python -u newMuonDIS/analyseEvents.py newMuonDIS/analyseEvents.yaml
```

Set `production_path`, `type`, `n_events`, `weight_factor`, and `output_dir` in
YAML. Input discovery recurses below the production directory. `disgen` and
`disfilter` read `MuonDIS`; replay types read `cbmsim`. Keep the normal
`sim_*_disreplay_gen.root` / `sim_*_disreplay_filter.root` names for replay
inputs, in a separate production directory containing only selected files.
`n_events: -1` reads all entries. `first_entry` and `n_events` address the sorted
input-file TChain; use disjoint ranges for split jobs.

Keyless ROOT files are skipped; nonempty files missing the required tree or
branches fail. Zero-entry trees retain their normalization. Outputs must be
new or empty directories. `run.json` is written last to mark successful completion.

## Prepared vertices

Generation/filter analysis is unchanged: read `mudis_DISvx/y/z_<material>` and
`muon_wDIS_<material>`. Each vertex fills x-versus-z, y-versus-z and y-versus-x
histograms once, with weight `wDIS * weight_factor`. `vertex_counts.csv` records
counts, sums of weights and sums of squared weights per material.

## Replay input and four-stage selection

**Replay analysis now starts from `replayFilter.py` output**, after incoming-muon
hit restoration. Unfiltered events that fail the baseline cause an error.
The original six-stage numbering has been replaced:

| Stage | Cumulative requirement |
| --- | --- |
| 0 | Filtered input: at least two tracks with Tr1–4 + TD hits and p > 1 GeV/c; at least one opposite-sign backward pair with DOCA < 5 cm and reconstructed vertex inside decay helium |
| 1 | At least one of those pairs also has DOCA < 1 cm and helium clearance >= 5 cm |
| 2 | At least one surviving pair has target-plane IP < 250 cm |
| 3 | At least one surviving pair has target-plane IP < 10 cm |

Hits are associated by MCTrack index. Straw stations are decoded as
`GetDetectorID() // 1000000`. Position and momentum are taken from the
smallest-z Tr1 hit. Pairs use the same backward straight-line closest-approach
calculation as the filter. Containment requires material `helium` under
`DecayVolume`; ROOT navigator safety supplies clearance, including endcaps.
IP is the transverse distance from the target axis at z=0, using the summed
Tr1 momentum. All distances are in cm and momenta in GeV/c.

An event passes a stage when any pair survives; the leading pair does not
control selection. The **leading pair** is the pair with the smallest DOCA
at the furthest stage reached by the event. Ties are resolved by ascending
MCTrack IDs. The best pair is stored independently at each passed stage, so changes in pair choice are retained.
Only this best pair is saved per stage, not an individual record for every
passing pair. All passing pairs contribute to pair cutflows and general pair
histograms. The saved tracks retain their Tr1 states and stage participation
masks, but those masks do not encode which tracks were paired together.

Geometry defaults to the matching `geo_` file beside each `sim_` file and is
updated when the TChain changes files. Set `geometry_file` for an explicit
shared geometry, e.g. when selected files are stored separately from geometry.

## One ROOT tree entry per event

Output version **7** writes `events.root`, containing the `DIS` TTree, with one
entry per filtered replay event (including events failing stages 1–3).
Branches are typed scalars, strings and vectors; there is no JSON event payload.
Integers retain track/PDG IDs; double precision retains momenta, positions and
weights. ZLIB level 1 compression and buffered baskets bound writer memory.
The small completion/normalization manifest remains `run.json`.

| Branches | Content |
| --- | --- |
| `source_file`, `file_uuid`, `source_entry` | Original file identity and local tree entry |
| `material`, `geometry_file`, `muDIS_muEntry`, `muDIS_disIndex` | Material, geometry reference and replay identity (missing IDs are -1) |
| `passed_steps`, `stage_weights` | Four stage flags and event weights; unreached weights are NaN |
| `n_charged`, `n_neutral` | Initial DIS product multiplicities |
| `selected_track_1`, `selected_track_2` | Four best-pair MCTrack IDs, ordered ascending; -1 for unreached stages |
| `p1_GeV`, `p2_GeV` | Four best-pair momentum magnitudes, measured at Tr1 |
| `doca_cm`, `target_ip_cm`, `vertex_{x,y,z}_cm` | Best-pair reconstruction at each stage |
| `common_ancestor_id`, `common_ancestor_pdg` | Nearest shared ancestor for each best pair; absent is -1/0 |
| `track_*` | All stored MCTracks, indexed by original MCTrack ID |
| `UBT_*`, `SBT_*` | All tagger hits, once per event |
| `muon_path_length_*` | Available input-muon path lengths, otherwise NaN |
| `clearance_*`, `wall_*` | Geometry-query results needed to reproduce the diagnostics |

Initial DIS products have mother ID 0 and raw weight different from the incoming
muon's weight, following `NewMuDISGenerator` and the existing Soft/DIS classification
(relative weight tolerance 1e-6). This excludes injected soft tracks, companion
muons and Geant4 shower/decay descendants. Counts include charged and neutral
initial products regardless of detector acceptance. They do not count every
generation of a shower. If the generator assigns identical muon and DIS weights,
the existing weight-based origin convention cannot distinguish them.

Shared ancestry follows mother IDs, not weights. It includes a selected track
when that track is the other's ancestor, but excludes incoming muon 0 (which is
otherwise common to unrelated DIS products). The nearest shared ancestor
minimizes the sum of parent-link distances; ties use the smaller MCTrack ID.
The signed PDG ID identifies the particle species, while the ancestor ID
identifies the particular MCTrack. Invalid cyclic/out-of-range ancestry fails.

`track_pdg_id`, `track_mother_id`, `track_charge_e`, `track_soft`, `track_mc_weight`
and `track_weight` describe each track. Position/momentum components are
`track_{start_position_cm,start_momentum_GeV,tr1_position_cm,tr1_momentum_GeV}_{x,y,z}`.
Missing Tr1 states use NaN. `track_station_mask` has bit N for straw station N,
`track_has_timing_hit` and `track_selected` retain detector selection, and
`track_stage_mask` has bit N when the track participates in any surviving pair
at stage N. Thus selected daughter momenta can be read directly by masking the
track vectors, with no stage duplication. Low-momentum and unselected tracks
remain available too.

Hit columns are `track_id`, `x_cm`, `y_cm`, `z_cm`, `time_ns`, `detector_id`,
`px_GeV`, `py_GeV`, `pz_GeV`, `energy_loss_GeV`, each prefixed by `UBT_` or `SBT_`.
Neutral-particle hits and invalid/unassociated track IDs remain in the tree.
There is no histogram-range truncation of the saved data.

The geometry cache retains exactly the helium-clearance and SBT-wall queries
used during extraction. The separate summary command reuses those results;
it needs neither the simulation nor geometry files and loads no FairShip
library. Changing the geometric selection to require additional queries needs
re-extraction. The summary is not a replacement for the complete simulation.

## Extract, merge, then summarize

For replay productions, defer histogram filling, large histogram archives and
plot rendering until after merging:

```sh
pixi cvmfs python -u newMuonDIS/analyseEvents.py replay.yaml \
    --tree-only --output-dir analysis/job_0
# Run other disjoint file/entry ranges into analysis/job_1, etc.
pixi cvmfs python -u newMuonDIS/mergeAnalysis.py analysis/job_* \
    --output-dir analysis/combined
pixi cvmfs python -u newMuonDIS/summarizeAnalysis.py analysis/combined \
    --output-dir analysis/summary --no-plots
```

Use a YAML with `type: disreplay_gen` or `disreplay_filter` and filtered inputs.
`tree_only: true` is the equivalent YAML option; it is invalid for gen/filter
vertex analysis. Without that option the existing combined extraction and
histogram workflow is retained, with ROOT replacing the JSON event file.

The summary command writes the existing tables and histogram archive. Omit
`--no-plots` to render PDF/PNG plots as well. It leaves the large tree in the
input directory and marks `events.root` omitted in the summary manifest. To
merge summary directories, use `--skip-event-records`. Tree-only and combined
histogram outputs cannot be mixed in one merge. Use the same analysis code
for extraction and summary; selection-code mismatches are rejected.

Read individual branches with ROOT/RDataFrame, or use the Python iterator
`analysisTree.read_records(path)` to reconstruct the dictionary view of an
event. The latter is for inspection; merging does not decode this view.

## Input-muon path lengths

Read `muon_path_length` and `muon_path_length_<material>` when those branches
exist in the current input tree. Zero is a measured zero path length;
`NaN` means unavailable. Missing values are never inferred from weights or
from the straight-line reconstructed track. Per-material keys are MS, UBT,
SBTsens, SBTfr, SSTsens, SSTfr, HE, AIR, HeBalloon and REST; `total` is the
upstream total excluding REST.

These branches currently exist in prepared `MuonDIS` files but are **not
propagated by NewMuDISGenerator into its standard replay output**. Such replay
files therefore produce NaN values. The analysis does not look up external
prepared files, whose identity cannot be safely inferred from a replay entry
number alone. `replayFilter.py` preserves the branches if an input has them.
The MS branch represents the preparation's restricted sampling path (last
20 cm in accumulated z, converted to trajectory length), not the full shield.

## Histograms and summary tables

Existing aggregate replay paths are `plots/step0/` through `plots/step3/`.
Every filled replay histogram also has a material-separated copy under
`plots/materials/<material>/stepN/`. The archive uses the same names without
`plots/`. Material copies are allocated only when filled, to avoid allocating
all the large empty 2D histograms for every material. Missing material/stage
histograms mean zero entries; merging handles disjoint material populations.

At every stage, pair plots contain vertex projections, DOCA, clearance, IP and
zoomed IP, and truth-DIS transverse/longitudinal/3D distances plus truth
containment. Pair multiplicity counts one entry per event. Particle histograms
`particles/n_charged` and `particles/n_neutral` count only initial DIS products,
with the stage event weight.

`selected_pair/{p1_GeV,p2_GeV}`
contains one entry per passed event/stage, weighted by the stage event weight.
The pair is chosen by minimum DOCA at that stage, with MCTrack-ID tie breaking.
`daughters/p_GeV` contains each non-input-muon track participating in any
surviving pair once per stage, weighted by its processed track weight. These
plots also have the material-separated copies described above. The daughter
momentum sample includes selected descendants and soft tracks; initial-DIS
multiplicity and the direct-daughter tagger sample have narrower definitions.

`tagger/all/SBT/{x,y}_vs_z` includes **all SBT hits**, weighted by the stage's
event weight, independent of pair membership or valid track association.
The existing `tagger/{muon,daughter}/` diagnostics remain: UBT windows are
full plane and centered rectangles of full width 400x600, 360x300 and 140x310 cm;
SBT maps split top/bottom/left/right and retain hit multiplicity versus earliest
SBT-hit z. Geometry cell IDs assign each associated hit to exactly one wall.

For the daughter diagnostic sample, require mother ID 0, nonzero charge and
weight different from the incoming muon's. Stage 0 includes all such daughters;
later stages include only daughters in surviving pairs, deduplicated within
the stage. `daughter_any` event categories can overlap; `daughter_all` requires
a nonempty sample all in one category. Empty samples use `no_daughters`.

Summary tables retain counts, sum of weights and sum of squared weights:

- `cutflow.csv`: events and pair counts per material/stage, including `ALL`.
- `pair_cutflow.csv`: pairs per material/stage, composition and signed PDG IDs.
- `pair_composition_cutflow.csv`: DIS-DIS, DIS-Soft and Soft-Soft totals.
- `tagger_cutflow.csv`: material/stage/sample/window/hit category.
- `normalization.csv`: original input-muon normalization per production stage.

`histograms.npz` stores bin edges, raw/weighted contents and outside-range totals.
`histogram_flows.csv` tabulates those flow totals. Histograms do not retain
per-bin sum of squared weights. `plot_axes` changes position ranges, not cuts.
PDF/PNG rendering is optional; `--no-plots` still writes all data and tables.

## Weights, normalization and merging

Track weights are raw fW times `weight_factor`. Pair weights count the shared
DIS/muon factor once. At **every stage, including stage 0**, the event weight
uses the incoming muon's fW if any surviving pair is Soft-Soft; otherwise it
uses the common DIS weight. The factor is applied once. Weighted cutflows can
therefore change when the last Soft-Soft pair fails, while raw event counts
remain nested. Original unselected preparation normalization remains intact.

The event tree stores processed weights before target-normalization scaling.
Multiply them by `normalization_scale` in `run.json` when that value is not null.
The merger recomputes this factor from the combined preparation denominator;
no per-event or per-track scaled branches need to be rewritten. Tables and
rendered plots retain their existing scaled outputs.

```sh
pixi cvmfs python -u newMuonDIS/mergeAnalysis.py analysis/job_* \
    --output-dir analysis/combined --no-plots
```

Merging validates source ranges and scans only the small provenance/stage
branches to reject duplicate, unordered or out-of-range entries and malformed
stage weights. For combined outputs it checks stage counts against cutflow
rows. It checks tree schemas and entry counts even for empty trees, then uses
ROOT's fast `TFileMerger` to copy compressed baskets without decoding/re-encoding
the track/hit payload. See [ROOT's fast tree-copy documentation](https://root.cern.ch/doc/master/classTTree.html).
Material histograms and aggregate tables are merged when present.

Use `--skip-event-records` to retain only histograms and aggregate tables;
`run.json` records the omission. A further merge of a reduced result requires
the same explicit option. This option is invalid for tree-only inputs.
Version-5/6 outputs must be regenerated: version 7 changes the multiplicity
meaning and adds per-stage pair information, so old JSON outputs cannot be
mixed with the new trees. No speedup is assumed from JSON compression alone;
the key change is avoiding Python decoding/re-encoding during tree merging.

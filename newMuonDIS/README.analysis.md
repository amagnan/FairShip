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
MCTrack IDs. It is stored once as two IDs plus the stage number.

Geometry defaults to the matching `geo_` file beside each `sim_` file and is
updated when the TChain changes files. Set `geometry_file` for an explicit
shared geometry, e.g. when selected files are stored separately from geometry.

## One detailed record per event

Output version **6** writes `events.jsonl`: one JSON object per input replay
event, streamed as it completes, including events failing stages 1–3. It replaces
`pairs.csv` and `tagger_tracks.csv`; neither detailed CSV is produced.

| Field | Content |
| --- | --- |
| `source_file`, `file_uuid`, `source_entry` | Source identity and local ROOT tree entry |
| `material`, `geometry_file`, `replay_identity` | DIS material, geometry reference and available `muDIS_muEntry` / `muDIS_disIndex` |
| `passed_steps` | Four event-level Boolean flags in stage order |
| `stage_weights` | Four processed event weights; `null` for failed/unreached stages |
| `particle_counts` | Charged/neutral counts for the populations below |
| `leading_pair` | `track_ids` and `step`, selected by smallest DOCA at the furthest stage |
| `tracks` | One record for **every stored MCTrack**, in original index order |
| `hits.UBT`, `hits.SBT` | All upstream-tagger and veto points, once per event |
| `muon_path_length_cm` | Available total and per-material input-muon path lengths, otherwise `null` |

The primary particle population is **all stored MCTracks except MCTrack[0]**,
including Geant4 descendants. Counts are not detector-acceptance counts, and
they do not deduplicate different generations of a shower. The named populations
are `excluding_input_muon`, `all`, `direct_daughters` (mother ID 0, including
injected soft tracks), and `excluding_input_muon_and_leading_pair`.

Each track saves its original index, signed PDG ID, mother index, charge in
units of e, Soft/DIS label, raw `mc_weight`, processed `weight`, start position
and momentum, straw-station coverage, timing-hit presence, smallest-z Tr1
position and momentum (or `null` if absent), and the detector/momentum
`selected` flag. Both low-momentum tracks and tracks missing detectors remain
in the record. The Soft/DIS label uses the same comparison with the incoming
muon's raw weight as the existing analysis; it is not a replacement for ancestry.

Each hit is
`[track_id, x_cm, y_cm, z_cm, time_ns, detector_id, px_GeV, py_GeV, pz_GeV, energy_loss_GeV]`.
All hits are kept, including neutral-particle hits, tracks outside the leading
pair, and invalid/unassociated track IDs. There is no stage duplication or
histogram-range truncation in these lists. The column order is also recorded
in `run.json` under `event_record_schema`.

Any pair's vertex, DOCA, IP and weight can be recomputed from the saved Tr1
states and weights. Ancestry and start positions retain DIS-truth associations.
Helium containment and clearance still require the referenced geometry file.
SBT wall assignments and alternative tagger windows can likewise be recovered
from the saved hit coordinates/IDs and geometry. This record is an analysis
summary, not a replacement for the complete ROOT simulation file.

## Input-muon path lengths

Read `muon_path_length` and `muon_path_length_<material>` when those branches
exist in the current input tree. Zero is a measured zero path length;
`null` means unavailable. Missing values are never inferred from weights or
from the straight-line reconstructed track. Per-material keys are MS, UBT,
SBTsens, SBTfr, SSTsens, SSTfr, HE, AIR, HeBalloon and REST; `total` is the
upstream total excluding REST.

These branches currently exist in prepared `MuonDIS` files but are **not
propagated by NewMuDISGenerator into its standard replay output**. Such replay
files therefore produce null values. The analysis does not look up external
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
`particles/n_charged` and `particles/n_neutral` count all stored tracks except
the incoming muon, with the stage event weight.

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

If target normalization is enabled, records also contain `scaled_stage_weights`
and each track has `scaled_weight`. The merger recalculates these from the
combined preparation denominator, rather than summing job-specific scale factors.

```sh
pixi cvmfs python -u newMuonDIS/mergeAnalysis.py analysis/job_* \
    --output-dir analysis/combined --no-plots
```

Merging streams event records, verifies ordered unique source entries against
the manifests, checks event stage counts against cutflow totals, merges material
histograms and preserves the event schema. Use `--skip-event-records` to keep
only histograms and aggregate tables; `run.json` records that omission. Merging
such a reduced result again requires the same explicit option.
`--skip-tagger-tracks` is replaced by `--skip-event-records` in version 6.
Version-5 outputs must be regenerated and cannot be mixed with version 6.

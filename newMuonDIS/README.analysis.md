# DIS event analysis

Run from an initialized FairShip environment:

```sh
pixi run --as-is python -u newMuonDIS/analyseEvents.py newMuonDIS/analyseEvents.yaml
```

Copy the example YAML and set `production_path`, `type`, `n_events`,
`weight_factor`, and `output_dir`. Discovery recurses below the production path.
It selects `muonDis_*.root` for `disgen` (excluding filter/replay names),
`muonDis_*_disfilter.root` for `disfilter`, and
`sim_*_disreplay_gen.root` or `sim_*_disreplay_filter.root` for replay.
Geometry and parameter files are excluded. `n_events: -1` processes all entries;
for prepared files each tree entry is an input muon and can contain multiple DIS
vertices per material. The limit is global across the sorted input files.

Empty ROOT files with no keys are skipped with a warning before normalization
and chain construction. They are excluded from the completion manifest's input
list. If all inputs are empty, analysis fails with a clear error. Unreadable
files and nonempty files missing the expected tree still fail. A valid tree with
zero entries is retained, including its normalization metadata.

The `MuonDIS` tree supplies `mudis_DISvx/y/z_<material>` and
`muon_wDIS_<material>`. Every DIS vertex contributes once to the raw plots and
with `wDIS * weight_factor` to the weighted plots. There is no additional
normalization by the number of generated or retained vertices.

Replay uses `cbmsim`, `MCTrack`, `strawtubesPoint`, `TimeDetPoint`,
`UpstreamTaggerPoint`, and `vetoPoint`. Both tagger branches are required;
missing branches are not interpreted as zero hits. Run replay post-processing
first to restore the input muon's upstream hits into `MCTrack[0]`.
Hits are associated by MCTrack index, with straw stations decoded as
`GetDetectorID() // 1000000`. Each distinct track must have hits in all four
stations and TD. Momentum and position come from the smallest-z Tr1 hit;
the momentum threshold is strictly greater than 1 GeV/c. Charge comes from the
MCTrack PDG identity. No restriction is imposed on the track's ancestry.

For each pair, the closest points of the straight lines through the Tr1 states
are found analytically. Both extrapolations must run backward from those states.
The vertex is their midpoint and DOCA is their separation. Parallel lines do not
define a unique vertex and are rejected. Successive cuts retain the same
candidates: opposite charges and DOCA < 1 cm, helium clearance >= 5 cm,
IP < 250 cm, then IP < 10 cm. The event passes if at least one candidate survives.
The IP is the transverse distance to (0,0) at z=0 after extrapolating with the
sum of the two Tr1 momenta. All geometry distances are in cm.

Replay weights use `MCTrack.GetWeight()` (`fW`). Each track plot receives its
own `fW * weight_factor`. For a pair, weights are compared with the original
muon's `MCTrack[0].fW`: two soft tracks use that muon weight; a mixed soft/DIS
pair uses the DIS track's full weight; two DIS tracks use their common full
weight once. The YAML factor is applied once, after selecting this weight.
Neither the muon factor nor the DIS factor is multiplied twice. Matching uses
relative tolerance 1e-6 (no absolute tolerance) for `Double32_t` storage;
inconsistent non-muon weights raise an error rather than silently choosing one.
A DIS weight numerically equal to the muon weight gives the same result either way.

The event cutflow counts each TChain event once per passed stage, regardless of
the number of surviving pairs. Step0 and its multiplicity plots use the initial
muon weight. Step1 uses all pairs of tracks passing the detector/momentum cuts;
subsequent steps use the pairs surviving that step. If any such pair is
soft–soft, it already satisfies the selection without requiring DIS, so the
event receives the muon weight. Otherwise it receives the common DIS weight.
The pair-count plot before step2 uses the step1 event weight, including events
with zero good pairs. Weighted sums can change when the last soft–soft pair
fails a cut; the raw event counts always form a nested cutflow. Pair plots still
contain one entry per pair with that pair's own weight.

The geometry defaults to the matching `geo_` file beside each replay `sim_`
file, and is updated as the TChain changes files. `geometry_file` can override
this when all inputs share one geometry. Containment requires material `helium`
under `DecayVolume`; ROOT navigator safety supplies clearance to the closest
boundary in any direction, including endcaps. Outside vertices are assigned
clearance -1 cm in the diagnostic histogram.

`plots/<material>/` holds prepared vertex plots. Replay uses `plots/step0/`
through `plots/step5/`: step0 shows track counts and momenta; step1 shows pair
charges, DOCA and accepted pair count; step2 shows vertices and wall clearance;
steps 3 and 4 show IP before their respective cuts; step5 shows final vertices.
Track and pair plots have one entry per track or candidate, respectively, while
multiplicity plots have one entry per event. Multiple candidates do not multiply
the event count. All plots are saved as both PDF and PNG, weighted and unweighted.
Empty stages and materials still get plots.

`cutflow.csv` (replay) or `vertex_counts.csv` (prepared files) contains counts,
sums of weights and sums of squared weights. `run.json` records the configuration,
input file list, entry counts and weighting convention. `histogram_flows.csv`
records entries outside the plotted ranges; those entries still contribute to
the cutflow. `plot_axes` adjusts the vertex histogram ranges without changing
selection. One-dimensional histogram ranges are fixed in the script.

## Replay tagger and DIS truth diagnostics

Every passed stage also records tagger diagnostics. The input muon is always
`MCTrack[0]`. The daughter sample requires `GetMotherId() == 0`, nonzero charge,
and a weight different from the input muon's weight, using the same tolerance
as the pair composition. The weight requirement excludes injected soft tracks,
which also have mother ID 0 in the replay generator. At step0 all such daughters
are considered; at step1 only daughters passing the detector/momentum selection
are considered; at steps 2–5 only daughters belonging to a surviving pair are
considered. A daughter shared by several pairs contributes once per stage.

`tagger_cutflow.csv` counts events by material (including `ALL`), stage, sample,
UBT window, and hit category: `UBT_only`, `SBT_only`, `both`, or `neither`.
The UBT windows are the full plane and centered rectangles of **full width**
400×600, 360×300, and 140×310 cm. Boundaries are included; for example the
400×600 window requires |x| ≤ 200 and |y| ≤ 300 cm. SBT presence always uses
all associated `vetoPoint` hits, independent of the UBT window.

The `muon` sample partitions the events into four categories. `daughter_any`
counts an event in a category when at least one eligible daughter has that
category; these counts can overlap. `daughter_all` requires a nonempty daughter
sample with every daughter in that same category. Mixed categories therefore
do not contribute to `daughter_all`. Events with no eligible daughters are
recorded as `no_daughters` for both daughter summaries, never as `neither` or
as satisfying an empty "all" condition. Summary weights are the existing
event weights for that stage.

`tagger_tracks.csv` retains one row per input muon or eligible daughter per
passed stage, including source file/UUID/entry, MCTrack index, PDG ID, weight,
and JSON columns for the four window categories and all associated UBT/SBT
points. Each point is `[x_cm, y_cm, z_cm, time_ns, detector_id]`; `sbt_walls`
contains the corresponding wall assignment for each SBT point. This preserves
hit coordinates and multiplicities, including hits outside histogram ranges,
and permits alternative per-daughter combinations without rerunning ROOT.
Muon hit plots use the stage's event weight; daughter hit plots use the
daughter's weight times `weight_factor`. The row's `weight` is the plot weight.

Under `stepN/tagger/{muon,daughter}/`, UBT plots show y versus x for every
window. SBT plots show x versus z for top (y>0) and bottom (y<0), and y versus z
for left (x<0) and right (x>0). **Each hit belongs to exactly one wall.**
The hit's detector ID identifies the sensitive cell copy in the loaded geometry;
its placed centre determines the side of the vessel. This remains valid when
the stored hit midpoint lies outside the cell or outside the helium z range.
`LiScX`, `LiSc_S3`, and `LiSc_S5` cells belong to top/bottom;
`LiScY`, `LiSc_S4`, and `LiSc_S6` belong to left/right. These volume families,
including the corner cells, are defined in `veto::add_block`. The copy numbers
are read from geometry nodes, without assuming a detector-ID encoding. If the
ID does not match a known cell, the hit position is used to locate one. For
other volume names, the fallback finds the helium
cross section at the hit's z using the actual ROOT shape, then assigns the
nearest transverse wall segment (exact corner ties go to top/bottom). This
accounts for the changing vessel dimensions; it does not use fixed x/y limits.
Failure to locate the helium cross section is reported explicitly.
For each wall, `n_hits_vs_entry_z` counts one entry per track with
at least one SBT hit, including zero hits in an unvisited wall.
The entry z is the z position of the earliest-time hit **anywhere in SBT**;
tracks without an SBT hit have no defined entry z and do not fill this plot.
Multiplicity bins cover 0–100 hits; overflow is retained in the usual flow table.

Under `stepN/truth_DIS/` (steps 1–5), `distance_xy_cm`, `distance_z_cm`, and
`distance_3d_cm` compare each built pair vertex with its associated truth DIS
vertex. The z distance is absolute. The broad plots cover 0–15000 cm and
`*_zoom_cm` plots cover 0–100 cm. A DIS secondary is traced to its mother-0
DIS ancestor, whose `GetStartX/Y/Z()` supplies the interaction position; soft
tracks are excluded. A shared truth position contributes once per pair, not
once per daughter. If distinct truth positions occur, each is retained.
Step1 includes all successfully built vertices before the opposite-sign/DOCA
cut; later steps use only pairs surviving that stage. Pure soft pairs do not
fill these plots.

`pairs.csv` includes JSON `truth_dis_vertices` and the flags
`truth_dis_in_decay_volume` (any associated truth vertex inside) and
`truth_dis_all_in_decay_volume` (all inside). Containment means helium under
`DecayVolume`, without the reconstructed vertex's 5 cm clearance cut. Pairs
without an associated DIS vertex have an empty vertex list and blank flags.
`truth_DIS/in_decay_volume` counts one entry per pair with associated DIS truth,
including step1 pairs for which a reconstructed vertex could not be built.

The histogram archive and both tagger tables are supported by `mergeAnalysis.py`.
Output version 4 requires rerunning older analysis outputs before merging with
these diagnostics. Use `--no-plots` for split jobs to defer rendering the
additional plots until the merge.

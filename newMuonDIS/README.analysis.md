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

The `MuonDIS` tree supplies `mudis_DISvx/y/z_<material>` and
`muon_wDIS_<material>`. Every DIS vertex contributes once to the raw plots and
with `wDIS * weight_factor` to the weighted plots. There is no additional
normalization by the number of generated or retained vertices.

Replay uses `cbmsim`, `MCTrack`, `strawtubesPoint`, and `TimeDetPoint`.
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

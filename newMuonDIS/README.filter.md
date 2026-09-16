Filtering prepared DIS events
============================

After building ShipMuDIS, run in the FairShip environment:

```sh
python newMuonDIS/filterEvents.py -f prepared.root -o selected.root \
    -g geometry.root --min-charged 2
```

Multiple input files are accepted after `-f`. The output must not already exist.
`-s` sets the first muon entry; `-n` limits the number of input entries (zero
processes none; -1 processes all remaining entries).

The filter cut applies independently to each DIS interaction in each material. By
default it requires at least two charged daughters whose extrapolated positions
are within a centred 4 m by 6 m plane: `|x| <= 200 cm`, `|y| <= 300 cm`.
All daughters start at their DIS vertex. Only crossings in the particle's
direction of motion count. Outgoing muons are included by default. Use
`--exclude-muons` to count only charged non-muons. ROOT's PDG database supplies
charges, with a nuclear-PDG-code fallback; unknown non-nuclear codes do not
count as charged.

The plane defaults to the global z coordinate of the `Tr1` station centre in
the imported ROOT geometry, including parent placements. `--detector-z` overrides
it in cm. Like `macro/ShipReco.py`, the filter loads `ShipGeo` from the geometry
file and calls `geomGeant4.addVMCFields(..., withVirtualMC=False)`. Only the SST
`MainSpecMap` is configured; the field maker remains alive throughout filtering.
The filename and translation default to `ShipGeo.Bfield.fieldMap` (relative to
`VMCWORKDIR`) and `ShipGeo.Bfield.z`. `--field-map` overrides the filename;
explicit paths are resolved relative to the current directory.
`--field-z` supplies an explicit map offset (cm), for example `0` for a map
already expressed in global coordinates. Map coordinates are in cm and map
file field values in Tesla, as required by `ShipBFieldMap`.

`--min-charged` sets the minimum number of accepted charged daughters.
`--no-detector-acceptance` explicitly restores multiplicity-only filtering;
geometry and field inputs are then optional. Geometry/map files are excluded
from recursive event-file discovery.

Only passing interactions are copied, together with their complete daughter
lists, cross sections, targets, vertices and times. Muons with no passing DIS
interaction are omitted. MC tracks and detector hits are preserved, as are
per-material `wDIS` weights (already normalized during generation). The
`MuonDIS` tree remains compatible with the prepared input schema.

Each material directory contains unweighted validation histograms before the
cut and with suffix `filtered` after the cut. Muon histograms before the cut
include zero-DIS entries; filtered muon histograms include only muons with a
passing interaction in that material. `muon_wDIS_` and `muon_wDIS_filtered`
record the per-material `wDIS` once per input muon in their respective samples,
with unit entry weight. Their initial range is [0, 1], extending automatically
to accommodate weights outside that range. Momentum fractions use sums of momentum
magnitudes; the muon fraction includes all muon daughters. Charged multiplicity
histograms always count all charged species, independent of `--exclude-muons`.
Malformed entries are logged and skipped before histogram filling.

`daughter_pdg_grouped_` and `daughter_pdg_grouped_filtered` group daughters
into 11 labelled bins: e+, e-, mu+, mu-, gamma, all neutrinos/antineutrinos,
pi+/-, other charged hadrons, pi0, other neutral hadrons, and other particles.
Hadron categories use ROOT's meson/baryon classification and charge; nuclei
and unrecognized PDG codes are included in the last bin.

For other selections, supply a C++ predicate (it replaces the default cut):

```cpp
MuDISFilter selection;
selection.SetFilter([](const std::vector<DISparticle>& daughters) {
  unsigned energetic = 0;
  for (const auto& p : daughters)
    if (std::abs(p.pid) == 211 && p.E > 1.) ++energetic;
  return energetic >= 2;
});
selection.process_file("prepared.root", "selected.root");
```

`SetFilter({})` restores the default multiplicity and detector-acceptance cut. The predicate
receives all daughters of a single interaction. It can also be declared via
`ROOT.gInterpreter.Declare` when steering from Python.

Reusable magnetic propagation
-----------------------------

`MagneticTrackPropagator` is independent of DIS event data. It takes a
`ShipBFieldMap*` and an optional `TGeoManager*`. Positions are in cm, momentum
vectors in GeV/c (magnitude and direction together), and charge in units of e.
Both objects are borrowed: keep them alive and unchanged during propagation.
A null field explicitly means zero field for generic transport; the filter
requires a supplied map when detector acceptance is enabled.

```cpp
ShipBFieldMap field("SST", "spectrometer.root", 0., 0., fieldZ);
MagneticTrackPropagator transport(&field, gGeoManager);
TVector3 start(0., 0., 4000.), momentum(0.1, 0., 10.);
TVector3 hit, momentumAtPlane;
double planeZ = transport.GetPlaneZ();  // Tr1, or pass another volume name
bool ok = transport.Extrapolate(-1., start, momentum, planeZ, hit, momentumAtPlane);

// Cache once, then query many z positions without evaluating the field again.
if (transport.BuildTrajectory(-1., start, momentum, planeZ)) {
  transport.PositionAt((start.Z() + planeZ) / 2., hit);
}

MuDISFilter selection;
selection.SetDetectorAcceptance(&field, gGeoManager);  // optional third argument: plane z
selection.process_file("prepared.root", "selected.root");
```

The map is scanned once for nonzero z layers. Neighbouring interpolation cells,
map placement, rotations and quadrant symmetry determine conservative global
field regions. Transport across field-free gaps is exactly linear. Inside the
field regions, adaptive RK4 with step doubling follows the Lorentz force; steps
are bounded by half the smallest grid spacing and a 5 cm maximum path step.
No threshold removes weak fringe fields. Grids with nonzero noise throughout
their volume consequently require integration throughout that volume.

`SetAccuracy(positionTolerance, relativeMomentumTolerance, maxStep)` controls
the local step comparison tolerances (defaults: 0.001 cm and 1e-6) and maximum
step in cm. These are local controls, not a bound on accumulated global error.
Cached positions use cubic Hermite interpolation between accepted steps and
linear segments; `GetTrajectorySize()` reports the number of cached points.
Queries outside the cached interval return false. Changing accuracy clears the
cache; failed trajectory construction also clears it.

Transport includes magnetic bending only, with constant momentum magnitude:
no energy loss, scattering, decays or geometry/material navigation. Generic
extrapolation can run in either z direction. Tracks with zero momentum, nearly
zero pz, a reversal in pz, or failed numerical convergence return false and do
not count towards filter acceptance. The cached path represents one monotonic
z crossing, not multiple crossings of a curling track. This implementation
does not change `MuDISProcessor` or `MuonPath`.

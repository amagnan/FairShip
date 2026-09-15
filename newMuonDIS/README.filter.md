Filtering prepared DIS events
============================

After building ShipMuDIS, run in the FairShip environment:

```sh
python newMuonDIS/filterEvents.py -f prepared.root -o selected.root --min-charged 2
```

Multiple input files are accepted after `-f`. The output must not already exist.
`-s` sets the first muon entry; `-n` limits the number of input entries (zero
processes none; -1 processes all remaining entries).

The filter cut applies independently to each DIS interaction in each material. By
default it requires two charged daughters, including outgoing muons. Use
`--exclude-muons` to count only charged non-muons. ROOT's PDG database supplies
charges, with a nuclear-PDG-code fallback; unknown non-nuclear codes do not
count as charged.

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

`SetFilter({})` restores the configurable multiplicity cut. The predicate
receives all daughters of a single interaction. It can also be declared via
`ROOT.gInterpreter.Declare` when steering from Python.

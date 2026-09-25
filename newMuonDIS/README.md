# Muon DIS simulation
## Introduction

This folder contains the necessary processors to create DIS events from input muons.

- Input: root files with cbmsim tree with muons from FairShip Particle Gun or MuonBack generators.
  - assumes muons are available at about the z position of the last interaction length of the muon shield.
  - To do 1: provide this position from MuonBack generator.
  - To do 2: process several muons from the same "event" = pot interaction.
- Output: tree "MuonDIS" with the initial muon information, hits in veto detectors and SST, soft particles emitted by the muon along its initial path, and daughters from DIS events in several volumes with random vertex positions within each volume, and associated probability weight.
  - Soft tracks are for all processes except destructive "Muon nuclear interaction". The first one is the initial input muon.
  - UBT, SBT, SST and TD hits are all hits with a GetTrackID() equal to the input muon track ID. TD hits are persisted in `muon_TDPoints`, including in filtered output. Older inputs without TD branches are supported.
  - The separate volumes are:
    - muon shield "MS"
    - UBT detector "UBT" (now just a dummy plane)
    - SBT detector "SBTsens" for liquid scintillator and "SBTfr" for frame material
    - SST detector "SSTsens" for liquid scintillator and "SSTfr" for frame material
    - Helium "HE"
    - Air "AIR" from before and after UBT, after balloon, cavern, SST.
    - The rest "REST" anything not taken into account in the previous categories.
    - to add material or change categories, edit top part of MuDISDefs.h file, and method GetLabel() in MuonPath.cxx

## Quick recipe

```bash
cd FairShip
pixi run build
pixi shell
cd <your_workdir>
python3 <relative_path_to_FairShip>/FairShip/newMuonDIS/prepareEvents.py -f <your_sim_input_root_file.root> -o <your_output_root_file_name>.root -n <number_of_initial_p.o.t./muon_events_to_process> -d <number_of_DIS_per_muon_per_volume> -g <your_geometry_file.root>
python3 <path_to_FairShip>/newMuonDIS/run_simScript_newDIS.py -f <muonDis_output.root> --tag <your_tag> --MuDIS --debug 1 -n 10
```

## Overview of classes:


- prepareEvents.py: python macro with argument parameters to pass to the main processor.
- class MuDISProcessor: main processor, reading input and creating output, and calling the others.
- class MuGeoProcessor: class defining the interface with the geometry. A map of objects of type "MuonPath" is filled, for each input muons, with the specific volumes traversed by the muon. This info will be used to generate vertices for DIS in each specific material, separately, and associate a random vertex position within each of these volumes.
- class MuonPath: to go along the input muon trajectory, with magnetic field extrapolation (to do), and distribute the path along the list of separate materials defined in MuDISDefs.h
- class DISparticle: DIS particle pid and 4-vector momentum.
- header file MuDISDefs.h: all helper classes and struct being used.
- class NewMuDISGenerator: FairGenerator to read again the MuonDIS tree and process particles again through Geant4 to give cbmsim tree.


## Path to volumes

The path uses the MC start and the first recorded muon hit in UBT, SBT,
each of Tr1--Tr4, and TD, omitting missing measurements. The available hits
are ordered in z. A non-forward (`pz <= 0`) selected measurement rejects the
muon. POCAs between adjacent measurements set the z planes at which the
momentum direction changes; each segment remains anchored to its own measured
position, momentum and time. Geometry navigation and DIS vertex positions use
these same straight segments. Times before a reference measurement are
extrapolated backwards. DIS vertices are sampled uniformly in segment length,
using a cached cumulative length and binary search. The existing shield
restriction to the last 20 cm in accumulated z is retained.

Transverse gaps at the switching planes are excluded from the path length.
At the end of processing, `LOG(info)` reports rejected muons, the number of
switches and muons with a transverse gap above the configured threshold, and
the largest gap. Set the threshold with `--poca-jump-threshold <cm>` in
`prepareEvents.py`, or `MuDISProcessor::SetPocaJumpThreshold(cm)` in C++.
The default is 1 cm, and a gap equal to the threshold is not counted.

## DIS events settings

- Pythia6 initialised with
  ```bash
  //set process 1=QCD, 2=DY/others
  fPythia->SetMSEL(2);
  //set min hard scale: 2 GeV --->try 1.5 for soft muons ?
  fPythia->SetPARP(2, 2);
  ```

- first nDIS/2 events are with proton target, second half with neutron target.
  - To do: understand why no difference between p and n target at the moment.


## Structure of output tree:

- branches muon_* : input muon information. Size: number of entries processed.
- branches muon_nDISevt_* : number of DIS interactions generated for each volume. Should be input parameter nDIS, but real number generated (in case some Pythia6 evt fail). Size: number of entries processed.
- branches mudis_*: DIS events information for each material. The vector size is variable and is determined by the actual per-entry `muon_nDISevt_<VOL>` count.
- branches mudis_DISproducts_* : all DIS daughters in each material <VOL>, with all DIS events stored together. Use `mudis_nDISdaughters_<VOL>` to determine the actual daughter ranges per DIS event.

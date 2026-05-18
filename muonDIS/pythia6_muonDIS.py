#!/usr/bin/env python
"""Script to generate DIS events for muons in Pythia6, and save them to a ROOT file (along with the original muon's soft interactions)."""

import argparse
import logging
import os
import time
from array import array

import ROOT as r
from tabulate import tabulate

PDG = r.TDatabasePDG.Instance()
PDG.AddParticle("C12", "Carbon-12", 12.0, True, 0, 6.0, "nucleus", 1000060120)
PDG.AddParticle("C13", "Carbon-13", 13.003355, True, 0, 6.0, "nucleus", 1000060130)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-o", "--outputFile", help="Output file to use", required=False, default="muonDis.root")
parser.add_argument(
    "-n",
    "--nEvents",
    dest="n_events",
    help="Number of muons to generate DIS for",
    required=False,
    default=10,
    type=int,
)
parser.add_argument(
    "-nDISPerMuon",
    "--nDIS",
    help="Number of DIS per muon to generate",
    required=False,
    default=10,
    type=int,
)
parser.add_argument(
    "-pmin",
    "--pMin",
    help="min p of muon to generate",
    required=False,
    default=2,
    type=float,
)
parser.add_argument(
    "-pmax",
    "--pMax",
    help="max p of muon to generate",
    required=False,
    default=400,
    type=float,
)
parser.add_argument(
    "-qmin",
    "--qMin",
    help="min q of pythia6 hard scale",
    required=False,
    default=2,
    type=float,
)


args = parser.parse_args()
n_events = args.n_events


def rotate(px, py, pz, theta, phi):
    """Rotate the daughter particle momentum to align with respect to the muon's momentum."""
    momentum = r.TVector3(px, py, pz)

    rotation = r.TRotation()
    rotation.RotateY(theta)  # Rotate around the Y-axis
    rotation.RotateZ(phi)  # Rotate around the Z-axis

    # Apply the rotation to the momentum vector
    rotated_momentum = rotation * momentum

    return rotated_momentum.X(), rotated_momentum.Y(), rotated_momentum.Z()


headers = [
    "DIS_index",
    "target type",
    "muon p (GeV)",
    "nParticles in event",
    "xsec (mb)"
]
Fixtarget = {1: "p+", 0: "n0"}


def inspect_file(filename):
    """Inspect the contents of muonDis file."""
    file = r.TFile.Open(filename, "READ")
    tree = file.DIS

    table_rows = []

    for i, event in enumerate(tree):

        nParticles = event.nDIS

        table_rows.append(
            [i, nParticles]
        )

    file.Close()
    logging.info(f"\n" + tabulate(table_rows, headers=headers, tablefmt="grid"))


def pythia6_muonDIS():
    """Generate DIS events."""

    logging.info(f"Creating output file: {args.outputFile}")

    outputFile = r.TFile.Open(args.outputFile, "recreate")
    output_tree = r.TTree("DIS", "muon DIS")

    # --- Size variable ---
    nDIS = array('i', [0])   # int branch

    # --- Vector branch ---
    brpid = r.std.vector('int')()
    brtid = r.std.vector('int')()
    brpx = r.std.vector('double')()
    brpy = r.std.vector('double')()
    brpz = r.std.vector('double')()
    brE = r.std.vector('double')()
    brxsec = r.std.vector('double')()

    # --- Create branches ---
    output_tree.Branch("nDIS", nDIS, "n/I")
    output_tree.Branch("pid", brpid)
    output_tree.Branch("target", brtid)
    output_tree.Branch("px", brpx)
    output_tree.Branch("py", brpy)
    output_tree.Branch("pz", brpz)
    output_tree.Branch("E", brE)
    output_tree.Branch("xsec", brxsec)


    myPythia = r.TPythia6()
    #set process 1=QCD, 2=DY/others
    myPythia.SetMSEL(2)
    #set min hard scale: 2 GeV --->try 1.5 for soft muons ?
    myPythia.SetPARP(2, args.qMin)
    #disable decay for those PDGID
    for kf in [211, 321, 130, 310, 3112, 3122, 3222, 3312, 3322, 3334]:
        kc = myPythia.Pycomp(kf)
        myPythia.SetMDCY(kc, 1, 0)

    seed = int(time.time() % 900000000)
    myPythia.SetMRPY(1, seed)
    #dictionary: pythia beam definition to enable gamma radiations
    mutype = {-13: "gamma/mu+", 13: "gamma/mu-"}
    #for pythia6 output, 11=set to output file. 6=stdout.
    myPythia.SetMSTU(11, 11)
    logging.info(
        f"Processing muon events"
    )

    nMade = 0
    R = r.TRandom(seed)
    
    for k in range(0,n_events):
        DIS_table = []  # debug
        cross_sections = []

        pid = 13
        px = 0
        py = 0
        pz = R.Uniform(args.pMin,args.pMax)

        p = r.TMath.Sqrt(px**2 + py**2 + pz**2)
        mass = PDG.GetParticle(abs(int(pid))).Mass()
        E = r.TMath.Sqrt(mass**2 + p**2)

        theta = r.TMath.ACos(pz / p)
        #returns phi between -pi and pi
        phi = r.TMath.ATan2(py, px)

        isProton = 1
        xsec = 0

        myPythia.Initialize("FIXT", mutype[pid], "p+", p)  # target = "p+"
        #print summary of initialisation params
        myPythia.Pylist(1)

        for a in range(args.nDIS):
            #half-way through, we change to neutron target with 50-50 : ---> update to real material ??
            if a == args.nDIS // 2:
                myPythia.Initialize("FIXT", mutype[pid], "n0", p)  # target = "n0"
                isProton = 0
                # logging.debug("Switching to neutron interaction")



            myPythia.GenerateEvent()
            #clean all but final stable particles
            myPythia.Pyedit(1)

            xsec = myPythia.GetPARI(1) #in mb
            ndaugh =  myPythia.GetN()
            nDIS[0] = ndaugh+1

            brpid.resize(ndaugh+1,0)
            brtid.resize(ndaugh+1,0)
            brpx.resize(ndaugh+1,0)
            brpy.resize(ndaugh+1,0)
            brpz.resize(ndaugh+1,0)
            brE.resize(ndaugh+1,0)
            brxsec.resize(ndaugh+1,0)

            brpid[0] = pid
            brtid[0] = isProton
            brpx[0] = px
            brpy[0] = py
            brpz[0] = pz
            brE[0] = E
            brxsec[0] = xsec
            
            #loop over daughters and rotate in muon input direction
            for itrk in range(1, ndaugh + 1):
                brpid[itrk] = myPythia.GetK(itrk, 2)
                brpx[itrk], brpy[itrk], brpz[itrk] = rotate(
                    myPythia.GetP(itrk, 1),
                    myPythia.GetP(itrk, 2),
                    myPythia.GetP(itrk, 3),
                    theta,
                    phi,
                )
                dpsq = brpx[itrk]**2 + brpy[itrk]**2 + brpz[itrk]**2
                dmasssq = PDG.GetParticle(brpid[itrk]).Mass() ** 2
                brE[itrk] = r.TMath.Sqrt(dmasssq + dpsq)
                brtid[itrk] = isProton
                brxsec[itrk] = xsec


            cross_sections.append(xsec)

            output_tree.Fill()
            #DIS_table.append(
            #    [
            #        a,
            #        Fixtarget[isProton],
            #        p,
            #        myPythia.GetN(),
            #        xsec,
            #    ]
            #)

        nMade += 1
        #logging.debug(
        #    f"\nMuon index:{k} \n\tPID = {pid}, weight = {w}\n\tpx = {px}, py = {py}, pz = {pz}, E = {E},\n\tx = {x}, y = {y}, z = {z}\nAvg xsec = {sum(cross_sections)/len(cross_sections)}\n\n\tDIS Events Summary\n{tabulate(DIS_table, headers=headers, tablefmt='grid')}"
        #)
        if nMade % 10 == 0:
            logging.info(f"Muons processed: {nMade} a={a}") #\n Avg xsec = {sum(cross_sections)/len(cross_sections)}\n \tDIS Events Summary\n{tabulate(DIS_table, headers=headers, tablefmt='grid')}")

    outputFile.cd()
    output_tree.Write()
    myPythia.SetMSTU(11, 6)
    logging.info(
        f"DIS generated for muons, output saved in {args.outputFile}, nDISPerMuon = {args.nDIS}"
    )
    outputFile.Close()



if __name__ == "__main__":
    logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
    )
    pythia6_muonDIS()
    #inspect_file(args.outputFile)

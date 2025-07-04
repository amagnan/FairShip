#include "vetoHit.h"

#include "FairLogger.h"   // for FairLogger, etc
#include "FairRun.h"
#include "FairRunSim.h"
#include "TGeoArb8.h"
#include "TGeoManager.h"
#include "TGeoPhysicalNode.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "veto.h"
#include "vetoPoint.h"

#include <iostream>
#include <math.h>

// -----   Default constructor   -------------------------------------------
vetoHit::vetoHit()
    : ShipHit()
{
    flag = true;
}
// -----   Standard constructor   ------------------------------------------
vetoHit::vetoHit(Int_t detID, Float_t adc)
    : ShipHit(detID, adc)
{
    ft = -1;
    flag = true;
}

// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
vetoHit::~vetoHit() {}
// -------------------------------------------------------------------------

TVector3 vetoHit::GetXYZ()
{
    TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
    TGeoNode* node = GetNode();
    TGeoVolume* volume = node->GetVolume();
    TGeoBBox* shape = dynamic_cast<TGeoBBox*>(volume->GetShape());
    Double_t origin[3] = {shape->GetOrigin()[0], shape->GetOrigin()[1], shape->GetOrigin()[2]};
    Double_t master[3] = {0, 0, 0};
    nav->LocalToMaster(origin, master);
    TVector3 pos = TVector3(master[0], master[1], master[2]);
    return pos;
}

Double_t vetoHit::GetX()
{
    TVector3 pos = GetXYZ();
    return pos.X();
}
Double_t vetoHit::GetY()
{
    TVector3 pos = GetXYZ();
    return pos.Y();
}
Double_t vetoHit::GetZ()
{
    TVector3 pos = GetXYZ();
    return pos.Z();
}
TGeoNode* vetoHit::GetNode()
{
    TGeoNode* node;
    TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
    TString path = "cave/DecayVolume_1/T2_1/VetoLiSc_0/";
    // id = ShapeType*100000 + blockNr*10000 + Zlayer*100 + number*10 + position;

    Int_t ShapeType = fDetectorID / 100000;
    Int_t blockNr = (fDetectorID % 100000) / 10000;
    Int_t Zlayer = (fDetectorID % 10000) / 100;
    Int_t number = (fDetectorID % 100) / 10;
    // Int_t position    = (fDetectorID % 10);

    if (ShapeType == 1)
        path += "LiScX_";
    else if (ShapeType == 2)
        path += "LiScY_";
    else if (ShapeType == 3)
        path += "LiSc_S3_";
    else if (ShapeType == 4)
        path += "LiSc_S4_";
    else if (ShapeType == 5)
        path += "LiSc_S5_";
    else if (ShapeType == 6)
        path += "LiSc_S6_";

    path += ShapeType;
    path += blockNr;
    if (Zlayer < 10)
        path += "0";
    path += Zlayer;
    path += number;
    path += "0_";
    path += fDetectorID;

    nav->cd(path);
    node = nav->GetCurrentNode();
    return node;
}

// -----   Public method Print   -------------------------------------------
void vetoHit::Print(Int_t detID) const
{
    LOG(INFO) << "vetoHit: veto hit "
              << " in detector " << fDetectorID;
    LOG(INFO) << "  ADC " << fdigi << " ns";
}

// -------------------------------------------------------------------------

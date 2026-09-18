#include "MuDISFilter.h"
#include "MagneticTrackPropagator.h"

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <iostream>

#include "FairLogger.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TVectorD.h"

using namespace ShipMuDIS;

namespace {
// Keep the nominal bin width and reserve one visible bin at each end.
template <typename Histogram>
Histogram* WithFlowBins(const char* name, const char* title, int bins, double low, double high) {
  const double width = (high - low) / bins;
  return new Histogram(name, title, bins + 2, low - width, high + width);
}

template <typename Histogram>
Histogram* WithFlowBins(const char* name, const char* title,
                        int nx, double xmin, double xmax, int ny, double ymin, double ymax) {
  const double dx = (xmax - xmin) / nx, dy = (ymax - ymin) / ny;
  return new Histogram(name, title, nx + 2, xmin - dx, xmax + dx, ny + 2, ymin - dy, ymax + dy);
}

double VisibleFlowValue(const TAxis* axis, double value) {
  if (value < axis->GetBinLowEdge(2)) return axis->GetBinCenter(1);
  if (value >= axis->GetBinLowEdge(axis->GetNbins())) return axis->GetBinCenter(axis->GetNbins());
  return value;
}

void FillWithFlow(TH1* histogram, double value) {
  histogram->Fill(VisibleFlowValue(histogram->GetXaxis(), value));
}

void FillWithFlow(TH2* histogram, double x, double y) {
  histogram->Fill(VisibleFlowValue(histogram->GetXaxis(), x), VisibleFlowValue(histogram->GetYaxis(), y));
}
}


// -----   Default constructor   -------------------------------------------
MuDISFilter::MuDISFilter() {
  ftree = 0;
  fouttree = 0;

  fLogger = FairLogger::GetLogger();
  fnEvts = -1;
  fstartEvt = 0;
  fPDG = TDatabasePDG::Instance();

}

Histograms MuDISFilter::BookHistograms(TDirectory* dir, const TString& mat, const TString & label)
{
  dir->cd();

  Histograms h;

  h.dis_vxz = WithFlowBins<TH2D>(
		       Form("vertex_x_vs_z_%s",label.Data()),
		       Form("DIS vertex x vs z - %s;z [cm];x [cm]", mat.Data()),
		       3300, 2700., 9000.,
		       600, -300., 300.
		       );
  h.dis_vyz = WithFlowBins<TH2D>(
		       Form("vertex_y_vs_z_%s",label.Data()),
		       Form("DIS vertex y vs z - %s;z [cm];y [cm]", mat.Data()),
		       3300, 2700., 9000.,
		       600, -300., 300.
		       );
  h.dis_vxy = WithFlowBins<TH2D>(
		       Form("vertex_y_vs_x_%s",label.Data()),
		       Form("DIS vertex y vs x - %s;x [cm];y [cm]", mat.Data()),
		       600, -300., 300.,
		       600, -300., 300.
		       );
  h.dis_vr = WithFlowBins<TH1D>(
      Form("vertex_r_%s", label.Data()),
      Form("DIS vertex radius - %s;r [cm];Events", mat.Data()),
      600, 0., 600.);
  h.dis_vz = WithFlowBins<TH1D>(
      Form("vertex_z_%s", label.Data()),
      Form("DIS vertex z - %s;z [cm];Events", mat.Data()),
      3300, 2700., 9000.);
  h.dis_pdg = WithFlowBins<TH1I>(
		       Form("daughter_pdg_%s",label.Data()),
		       Form("PDG code of DIS daughters - %s;PDG code;Particles", mat.Data()),
		       12001, -6000.5, 6000.5
		       );
  h.dis_pdgGrouped = WithFlowBins<TH1I>(
                             Form("daughter_pdg_grouped_%s",label.Data()),
                             Form("DIS daughter species - %s;Species;Particles", mat.Data()),
                             11, 0.5, 11.5
                             );
  const char* labels[] = {"e^{+}", "e^{-}", "#mu^{+}", "#mu^{-}", "#gamma",
                          "#nu/#bar{#nu}", "#pi^{#pm}", "h^{#pm}", "#pi^{0}", "h^{0}", "other"};
  for (int bin = 1; bin <= 11; ++bin)
    h.dis_pdgGrouped->GetXaxis()->SetBinLabel(bin + 1, labels[bin - 1]);
  h.dis_n = WithFlowBins<TH1I>(
		     Form("n_daughters_%s",label.Data()),
		     Form("Number of DIS daughters - %s;multiplicity;Event", mat.Data()),
		     50,0,50
		     );
  h.dis_nCharged = WithFlowBins<TH1I>(
			    Form("n_daughters_charged_%s",label.Data()),
			    Form("Number of charged DIS daughters - %s;charged multiplicity;Event", mat.Data()),
			    50,0,50
			    );
  h.dis_pChargedFrac = WithFlowBins<TH1D>(
				Form("pfrac_charged_%s",label.Data()),
				Form("Charged fraction - %s;p(charged)/p(all);Event", mat.Data()),
				101,0,1.01
				);
  h.dis_pMuFrac = WithFlowBins<TH1D>(
			   Form("pfrac_mu_%s",label.Data()),
			   Form("p fraction of outgoing mu - %s;p(mu)/p(all);Event", mat.Data()),
			   101,0,1.01
			   );
  h.mu_p = WithFlowBins<TH1D>(
		    Form("muon_p_%s",label.Data()),
		    Form("muon momentum - %s;p_{#mu,in} [GeV]; Input mu events", mat.Data()),
		    400, 0., 400.
		    );
  h.mu_pt = WithFlowBins<TH1D>(
		     Form("muon_pt_%s",label.Data()),
		     Form("muon p_{T} - %s;p_{T,#mu,in} [GeV];Input mu events", mat.Data()),
		     100, 0., 10.
		     );
  h.mu_ppt = WithFlowBins<TH2D>(
		      Form("muon_pt_vs_p_%s",label.Data()),
		      Form("muon p_{T} vs p - %s;p_{#mu,in} [GeV]; p_{T,#mu,in} [GeV];Input mu events", mat.Data()),
		      400,0,400,
		      100, 0., 10.
		      );
  h.mu_ndis = WithFlowBins<TH1I>(
		       Form("muon_n_dis_%s",label.Data()),
		       Form("Number of DIS events - %s;DIS multiplicity;Input mu events", mat.Data()),
		       1001,0,1001
		       );
  h.mu_wdis = WithFlowBins<TH1D>(
		       Form("muon_vtx_weight_%s",label.Data()),
		       Form("DIS vertex weight - %s;wDIS;Input mu events", mat.Data()),
		       100, 0., 1.
		       );
  return h;
}

void MuDISFilter::init(const int& aEvts, const int& aStart){
  if (aEvts < -1 || aStart < 0)
    throw std::invalid_argument("Invalid event range");
  fnEvts = aEvts;
  fstartEvt = aStart;
}

Bool_t MuDISFilter::InitFile(const char* fileName) {
  return InitFile(fileName,0);
}

Bool_t MuDISFilter::InitFiles(const std::vector<std::string>& fileNames) {
  return InitFiles(fileNames, 0);
}

// -----   Default constructor   -------------------------------------------
Bool_t MuDISFilter::InitFile(const char* fileName, const int startEvent) {
  std::vector<std::string> fileNames = {fileName};
  return InitFiles(fileNames, startEvent);
}

Bool_t MuDISFilter::InitFiles(const std::vector<std::string>& fileNames,
			    const int startEvent) {
  
  if (startEvent < 0) return kFALSE;
  fstartEvt = startEvent;
  if (fileNames.empty()) {
    LOG(error) << "MuDISFilter: no input files provided. "
               << "Check the -f/--inputFile argument or input file glob.";
    return kFALSE;
  }
  for (const auto& fileName : fileNames) {
    if (fileName.empty()) {
      LOG(error) << "MuDISFilter: received an empty input file name. "
                 << "Check the -f/--inputFile argument.";
      return kFALSE;
    }
  }
  
  for (const auto& name : fileNames) {
    std::unique_ptr<TFile> file(TFile::Open(name.c_str(), "READ"));
    TTree* tree = nullptr;
    if (file && !file->IsZombie()) file->GetObject("MuonDIS", tree);
    if (!tree) {
      LOG(error) << "MuDISFilter: missing MuonDIS tree in " << name;
      return kFALSE;
    }
    MuonInBranches check;
    if (!check.Setup(tree)) return kFALSE;
    tree->ResetBranchAddresses();
  }

  {
    delete ftree;
    ftree = new TChain("MuonDIS");
    for (const auto& name : fileNames) {
      if (!ftree->Add(name.c_str())) return kFALSE;
    }
    Long64_t treeEvts = ftree->GetEntries();
    LOG(info) << "Reading " << treeEvts << " entries.";
    bool ok = finEv.Setup(ftree);
    
    if (!ok) {
      LOG(error)
        << "MuDISFilter: failed to bind one or more required branches";
      return kFALSE;
    }    
    LOG(info) << "MuDISFilter: Initialization successful.";
    return kTRUE;
  }
  return kFALSE;
}

void MuDISFilter::process_file(const std::string& input,
			       const std::string& output) {
  std::vector<std::string> fileNames = {input};
  return process_file(fileNames,output);
}

void MuDISFilter::process_file(const std::vector<std::string>& input,
			       const std::string& output) {

  if (!fFilter && fFilterOption != 0 && !fUseDetectorAcceptance)
    throw std::runtime_error("Filter options 1 and 2 require detector acceptance");
  if (!fFilter && fUseDetectorAcceptance && !fPropagator)
    throw std::runtime_error("Configure detector acceptance with geometry and field map first");

  Bool_t treeOK = InitFiles(input, fstartEvt);

  if (!treeOK) {
    throw std::runtime_error("MuDISFilter: failed to initialize input files");
  }

  if (!fFilter && fUseDetectorAcceptance && fPropagator->HasMuonShieldField()) {
    // A separate, vertex-only reader preserves the filtering chain's branch
    // addresses and cache learning. Restrict the scan to the requested entries.
    TChain vertices("MuonDIS");
    for (const auto& name : input) vertices.Add(name.c_str());
    TTreeReader reader(&vertices);
    std::vector<std::unique_ptr<TTreeReaderValue<std::vector<double>>>> vz;
    for (unsigned imat = 0; imat < nMats; ++imat)
      vz.emplace_back(new TTreeReaderValue<std::vector<double>>(reader, "mudis_DISvz_" + MatTypeStr[imat]));
    const Long64_t end = fnEvts >= 0
        ? std::min(Long64_t(fstartEvt) + fnEvts, ftree->GetEntries()) : ftree->GetEntries();
    double minimum = std::numeric_limits<double>::infinity(), msMinimum = minimum;
    bool valid = true;
    if (fstartEvt < end) {
      reader.SetEntriesRange(fstartEvt, end);
      Long64_t scanned = 0;
      while (reader.Next()) {
        ++scanned;
        for (unsigned imat = 0; imat < nMats; ++imat) {
          const auto* values = vz[imat]->Get();
          if (!values || vz[imat]->GetSetupStatus() < 0
              || vz[imat]->GetReadStatus() != TTreeReaderValue<std::vector<double>>::kReadSuccess) {
            valid = false; continue;
          }
          for (double z : *values) {
            if (!std::isfinite(z)) { valid = false; continue; }
            minimum = std::min(minimum, z);
            if (MatTypeStr[imat] == "MS") msMinimum = std::min(msMinimum, z);
          }
        }
      }
      valid = valid && scanned == end - fstartEvt;
    }
    if (valid && std::isfinite(minimum)) {
      const double detectorMinimum = fFilterOption == 2 ? fDetectorVolumeZ.first
          : fFilterOption == 1 ? *std::min_element(fStationZ.begin(), fStationZ.end()) : fDetectorZ;
      minimum = std::min(minimum, detectorMinimum);
      LOG(info) << "Muon shield field bounds: input MS minimum z = " << msMinimum
                << " cm; conservative scan minimum z = " << minimum << " cm";
      fPropagator->SetMuonShieldMinZ(minimum);
    } else {
      // No usable vertices or failed reads: retain the complete field map.
      fPropagator->SetMuonShieldMinZ(-std::numeric_limits<double>::infinity());
      if (!valid) LOG(warn) << "Vertex range scan failed; using full muon shield field bounds";
    }
  }

  std::unique_ptr<TFile> outfile(TFile::Open(output.c_str(), "CREATE"));
  if (!outfile || outfile->IsZombie()) {
    throw std::runtime_error("MuDISFilter: cannot create output file " + output);
  }
  outfile->cd();

  for (unsigned imat = 0; imat < nMats; ++imat) {
    TDirectory* dir =
      outfile->mkdir(MatTypeStr[imat].Data());
    
    hist_all[imat] = BookHistograms(dir, MatTypeStr[imat]);
    hist_filt[imat] = BookHistograms(dir, MatTypeStr[imat], "filtered");
  }

  
  outfile->cd();
  fouttree = new TTree("MuonDIS",
		       "Muon information, DIS products and soft interaction tracks");
  foutEv.InitTree(fouttree);
  
  Long64_t n = ftree->GetEntries();
  LOG(info) << " * input tree with " << n << " entries" << std::endl;

  ProcessEvents();
  
  outfile->cd();
  if (outfile->Write() <= 0 || outfile->TestBit(TFile::kWriteError))
    throw std::runtime_error("MuDISFilter: failed to write output file");
  outfile->Close();
  fouttree = nullptr;

}

MuDISFilter::~MuDISFilter() {
  delete ftree;
}

void MuDISFilter::initEvent() {
  foutEv.initEvent();
  for (unsigned i = 0; i < nMats; ++i) foutEv.br[i].initEvent(0);
}

bool MuDISFilter::IsCharged(const DISparticle& particle) const {
  return Charge(particle) != 0.;
}

double MuDISFilter::Charge(const DISparticle& particle) const {
  const auto* pdg = fPDG->GetParticle(particle.pid);
  // Nuclear PDG codes encode Z in digits 5--7 (10LZZZAAAI).
  if (!pdg && std::abs(particle.pid) >= 1000000000)
    return ((std::abs(particle.pid) / 10000) % 1000) * (particle.pid > 0 ? 1. : -1.);
  return pdg ? pdg->Charge() / 3. : 0.;
}

int MuDISFilter::DaughterCategory(const DISparticle& particle) const {
  switch (particle.pid) {
    case -11: return 1;
    case 11: return 2;
    case -13: return 3;
    case 13: return 4;
    case 22: return 5;
    case 12: case -12: case 14: case -14: case 16: case -16: return 6;
    case 211: case -211: return 7;
    case 111: return 9;
  }
  const auto* pdg = fPDG->GetParticle(particle.pid);
  if (pdg) {
    const TString particleClass = pdg->ParticleClass();
    if (particleClass == "Meson" || particleClass == "Baryon")
      return pdg->Charge() != 0. ? 8 : 10;
  }
  return 11;
}

void MuDISFilter::ConfigureFilterGeometry(unsigned option, const MagneticTrackPropagator& propagator) {
  if (option == 1) {
    const std::array<double, 5> planes = {propagator.GetPlaneZ("Tr1"), propagator.GetPlaneZ("Tr2"),
                                         propagator.GetPlaneZ("Tr3"), propagator.GetPlaneZ("Tr4"),
                                         propagator.GetPlaneZ("Timing Detector")};
    fStationZ = planes;
  } else if (option == 2) {
    const double front = propagator.GetVolumeZRange("Tr1").first;
    const double end = propagator.GetVolumeZRange("SplitCalDetector").second;
    if (front >= end) throw std::runtime_error("Calorimeter end must be downstream of the Tr1 front");
    fDetectorVolumeZ = {front, end};
  }
}

void MuDISFilter::SetFilterOption(unsigned option) {
  if (option > 2) throw std::invalid_argument("Filter option must be 0, 1 or 2");
  if (fPropagator) ConfigureFilterGeometry(option, *fPropagator);
  fFilterOption = option;
}

bool MuDISFilter::HitsTrackingAndTiming(double charge, const DISparticle& particle, const TVector3& vertex) const {
  TVector3 position = vertex, momentum(particle.px, particle.py, particle.pz);
  std::array<unsigned, 5> order = {0, 1, 2, 3, 4};
  std::sort(order.begin(), order.end(), [this, &particle](unsigned a, unsigned b) {
    return particle.pz >= 0. ? fStationZ[a] < fStationZ[b] : fStationZ[a] > fStationZ[b];
  });
  bool firstPair = false, secondPair = false, timing = false;
  for (unsigned station : order) {
    if ((fStationZ[station] - vertex.Z()) * particle.pz < 0.) continue;
    TVector3 hit, nextMomentum;
    if (!fPropagator->Extrapolate(charge, position, momentum, fStationZ[station], hit, nextMomentum)) return false;
    position = hit;
    momentum = nextMomentum;
    if (std::abs(hit.X()) > 200. || std::abs(hit.Y()) > 300.) continue;
    if (station < 2) firstPair = true;
    else if (station < 4) secondPair = true;
    else timing = true;
    if (firstPair && secondPair && timing) return true;
  }
  return false;
}

bool MuDISFilter::PassDetectorFilter(const std::vector<DISparticle>& daughters, const TVector3& vertex) const {
  if (fFilterOption == 2) {
    const TVector3 minimum(-200., -300., fDetectorVolumeZ.first);
    const TVector3 maximum(200., 300., fDetectorVolumeZ.second);
    for (const auto& particle : daughters) {
      if (!fIncludeMuons && std::abs(particle.pid) == 13) continue;
      if (fPropagator->IntersectsBox(Charge(particle), vertex,
                                   TVector3(particle.px, particle.py, particle.pz), minimum, maximum)) return true;
    }
    return false;
  }
  // Count candidates before expensive transport; option 1 always requires two.
  unsigned remaining = 0, accepted = 0;
  std::vector<double> charges;
  charges.reserve(daughters.size());
  for (const auto& particle : daughters) {
    const double charge = !fIncludeMuons && std::abs(particle.pid) == 13 ? 0. : Charge(particle);
    charges.push_back(charge);
    if (charge != 0.) ++remaining;
  }
  for (std::size_t i = 0; i < daughters.size(); ++i) {
    if (accepted + remaining < 2) return false;
    if (charges[i] == 0.) continue;
    --remaining;
    if (HitsTrackingAndTiming(charges[i], daughters[i], vertex) && ++accepted == 2) return true;
  }
  return false;
}

void MuDISFilter::SetDetectorAcceptance(ShipBFieldMap* field, TGeoManager* geometry, double z,
                                      ShipBFieldMap* muonShieldField) {
  if (!field) throw std::invalid_argument("Detector acceptance requires a field map");
  auto propagator = std::make_unique<MagneticTrackPropagator>(field, geometry, muonShieldField);
  if (muonShieldField) {
    const auto range = propagator->GetMuonShieldZRange();
    LOG(info) << "MuDISFilter: muon shield geometry z range [" << range.first << ", "
              << range.second << "] cm; propagating with shield map including nonzero fringe fields";
  }
  const double planeZ = std::isnan(z) ? propagator->GetPlaneZ() : z;
  if (!std::isfinite(planeZ)) throw std::invalid_argument("Detector z must be finite");
  ConfigureFilterGeometry(fFilterOption, *propagator);
  fPropagator = std::move(propagator);
  fDetectorZ = planeZ;
  fUseDetectorAcceptance = true;
}

bool MuDISFilter::PassFilter(const std::vector<DISparticle>& daughters) const {
  if (!fFilter && (fUseDetectorAcceptance || fFilterOption != 0))
    throw std::runtime_error("Detector acceptance requires the DIS vertex; use PassFilter(daughters, vertex)");
  return PassFilter(daughters, TVector3());
}

bool MuDISFilter::PassFilter(const std::vector<DISparticle>& daughters, const TVector3& vertex) const {
  if (fFilter) return fFilter(daughters);
  if (fFilterOption != 0 && !fUseDetectorAcceptance)
    throw std::runtime_error("Filter options 1 and 2 require detector acceptance");
  if (fUseDetectorAcceptance && !fPropagator)
    throw std::runtime_error("Configure detector acceptance with geometry and field map first");
  if (fFilterOption != 0) return PassDetectorFilter(daughters, vertex);
  if (fMinChargedDaughters == 0) return true;
  if (daughters.size() < fMinChargedDaughters) return false;
  std::vector<double> charges;
  unsigned remaining = 0;
  if (fUseDetectorAcceptance) {
    charges.reserve(daughters.size());
    for (const auto& p : daughters) {
      // Count only candidates that can make a forward physical crossing.
      const double charge = (!fIncludeMuons && std::abs(p.pid) == 13)
          || (fDetectorZ - vertex.Z()) * p.pz < 0. ? 0. : Charge(p);
      charges.push_back(charge);
      if (charge != 0.) ++remaining;
    }
    if (remaining < fMinChargedDaughters) return false;
  }
  unsigned charged = 0;
  for (std::size_t i = 0; i < daughters.size(); ++i) {
    const auto& p = daughters[i];
    if (!fIncludeMuons && std::abs(p.pid) == 13) continue;
    const double charge = fUseDetectorAcceptance ? charges[i] : Charge(p);
    if (charge == 0.) continue;
    if (fUseDetectorAcceptance) {
      if (charged + remaining < fMinChargedDaughters) return false;
      --remaining;
      TVector3 hit, momentum;
      if (!fPropagator->Extrapolate(charge, vertex, TVector3(p.px, p.py, p.pz),
                                    fDetectorZ, hit, momentum)
          || std::abs(hit.X()) > 200. || std::abs(hit.Y()) > 300.) continue;
    }
    ++charged;
    if (charged >= fMinChargedDaughters) return true;
  }
  return charged >= fMinChargedDaughters;
}

void MuDISFilter::FillDIS(Histograms& h, const MuonDISInBranches& br,
                         int idis, const std::vector<DISparticle>& daughters,
                         Histograms* filtered) {
  unsigned charged = 0;
  double totalP = 0., chargedP = 0., muonP = 0.;
  for (const auto& p : daughters) {
    const int category = DaughterCategory(p);
    for (auto* target : {&h, filtered}) {
      if (!target) continue;
      FillWithFlow(target->dis_pdg, p.pid);
      FillWithFlow(target->dis_pdgGrouped, category);
    }
    const double momentum = std::sqrt(p.px*p.px + p.py*p.py + p.pz*p.pz);
    totalP += momentum;
    if (IsCharged(p)) { ++charged; chargedP += momentum; }
    if (std::abs(p.pid) == 13) muonP += momentum;
  }
  for (auto* target : {&h, filtered}) {
    if (!target) continue;
    FillWithFlow(target->dis_vxz, br.DISvz->at(idis), br.DISvx->at(idis));
    FillWithFlow(target->dis_vyz, br.DISvz->at(idis), br.DISvy->at(idis));
    FillWithFlow(target->dis_vxy, br.DISvx->at(idis), br.DISvy->at(idis));
    FillWithFlow(target->dis_vr, std::hypot(br.DISvx->at(idis), br.DISvy->at(idis)));
    FillWithFlow(target->dis_vz, br.DISvz->at(idis));
    FillWithFlow(target->dis_n, daughters.size());
    FillWithFlow(target->dis_nCharged, charged);
    if (totalP > 0.) {
      FillWithFlow(target->dis_pChargedFrac, chargedP / totalP);
      FillWithFlow(target->dis_pMuFrac, muonP / totalP);
    }
  }
}

void MuDISFilter::ProcessEvents() {
  if (!ftree || !fouttree) throw std::runtime_error("Initialize input and output first");
  if (fstartEvt < 0) throw std::invalid_argument("Negative start event");
  const Long64_t end = fnEvts >= 0
      ? std::min(Long64_t(fstartEvt) + fnEvts, ftree->GetEntries())
      : ftree->GetEntries();
  Long64_t selected = 0, skipped = 0;
  Long64_t selectedDIS[nMats] = {};
  double weightedDIS[nMats] = {};
  Long64_t processedDIS[nMats] = {};
  double weightedProcessedDIS[nMats] = {};
  double weightedMuons[nMats] = {}, weightedSelectedMuons[nMats] = {};
  for (Long64_t event = fstartEvt; event < end; ++event) {
    if ((event - fstartEvt) % 100 == 0)
      LOG(info) << "MuDISFilter: processing entry " << event;
    if (ftree->GetEntry(event) <= 0 || !finEv.mcTrks || finEv.mcTrks->empty()
        || !finEv.sbtPt || !finEv.ubtPt || !finEv.sstPt) {
      ++skipped;
      continue;
    }
    // Check all vector lengths before filling histograms or copying any data.
    bool valid = true;
    for (const auto& br : finEv.br) {
      const int n = br.nDISevts;
      if (n < 0 || !br.DISxsec || !br.DIStarget || !br.DISvx || !br.DISvy
          || !br.DISvz || !br.DISvt || !br.nDISdau || !br.DISparticles) {
        valid = false;
        break;
      }
      const auto size = static_cast<std::size_t>(n);
      if (br.DISxsec->size() != size || br.DIStarget->size() != size
          || br.DISvx->size() != size || br.DISvy->size() != size
          || br.DISvz->size() != size || br.DISvt->size() != size
          || br.nDISdau->size() != size) { valid = false; break; }
      std::size_t total = 0;
      for (int count : *br.nDISdau) {
        if (count < 0 || static_cast<std::size_t>(count) > br.DISparticles->size() - total) {
          valid = false;
          break;
        }
        total += count;
      }
      if (total != br.DISparticles->size()) valid = false;
    }
    if (!valid) {
      LOG(error) << "MuDISFilter: malformed entry " << event << "; skipping";
      ++skipped;
      continue;
    }
    initEvent();
    const auto& muon = finEv.mcTrks->at(0);
    const double muonP = muon.GetP(), muonPt = muon.GetPt();
    bool keep = false;
    for (unsigned imat = 0; imat < nMats; ++imat) {
      const auto& in = finEv.br[imat];
      auto& out = foutEv.br[imat];
      out.wDIS = in.wDIS;
      std::size_t offset = 0;
      for (int idis = 0; idis < in.nDISevts; ++idis) {
        const auto endOffset = offset + in.nDISdau->at(idis);
        std::vector<DISparticle> daughters(in.DISparticles->begin() + offset,
                                           in.DISparticles->begin() + endOffset);
        offset = endOffset;
        const bool accepted = PassFilter(daughters, TVector3(in.DISvx->at(idis), in.DISvy->at(idis),
                                                            in.DISvz->at(idis)));
        FillDIS(hist_all[imat], in, idis, daughters, accepted ? &hist_filt[imat] : nullptr);
        if (!accepted) continue;
        ++out.nDISevts;
        out.DISxsec.push_back(in.DISxsec->at(idis));
        out.DIStarget.push_back(in.DIStarget->at(idis));
        out.DISvx.push_back(in.DISvx->at(idis));
        out.DISvy.push_back(in.DISvy->at(idis));
        out.DISvz.push_back(in.DISvz->at(idis));
        out.DISvt.push_back(in.DISvt->at(idis));
        out.nDISdau.push_back(daughters.size());
        out.DISparticles.insert(out.DISparticles.end(), daughters.begin(), daughters.end());
        keep = true;
      }
      auto fillMuon = [muonP, muonPt, &in](Histograms& h, int count) {
        FillWithFlow(h.mu_p, muonP);
        FillWithFlow(h.mu_pt, muonPt);
        FillWithFlow(h.mu_ppt, muonP, muonPt);
        FillWithFlow(h.mu_ndis, count);
        FillWithFlow(h.mu_wdis, in.wDIS);
      };
      fillMuon(hist_all[imat], in.nDISevts);
      if (out.nDISevts > 0) fillMuon(hist_filt[imat], out.nDISevts);
      selectedDIS[imat] += out.nDISevts;
      if (out.nDISevts > 0) weightedDIS[imat] += out.nDISevts * in.wDIS;
      processedDIS[imat] += in.nDISevts;
      if (in.nDISevts > 0) weightedProcessedDIS[imat] += in.nDISevts * in.wDIS;
      weightedMuons[imat] += in.wDIS;
      if (out.nDISevts > 0) weightedSelectedMuons[imat] += in.wDIS;
    }
    if (keep) {
      foutEv.mcTrks = *finEv.mcTrks;
      foutEv.sbtPt = *finEv.sbtPt;
      foutEv.ubtPt = *finEv.ubtPt;
      foutEv.sstPt = *finEv.sstPt;
      if (fouttree->Fill() < 0) throw std::runtime_error("Failed writing MuonDIS entry");
      ++selected;
    }
  }
  LOG(info) << "MuDISFilter: saved " << selected << " muon entries; skipped "
            << skipped << " unreadable or malformed entries.";
  for (unsigned imat = 0; imat < nMats; ++imat) {
    LOG(info) << "MuDISFilter: selected DIS events in " << MatTypeStr[imat].Data()
              << ": raw = " << selectedDIS[imat] << ", weighted = " << weightedDIS[imat];
    auto* counts = new TH1D("filter_counts", "Filter counts;Sample;Count", 8, 0., 8.);
    counts->SetDirectory(fouttree->GetDirectory()->GetDirectory(MatTypeStr[imat].Data()));
    const char* labels[] = {"muons_processed_raw", "muons_processed_weighted",
                           "muons_selected_raw", "muons_selected_weighted",
                           "dis_processed_raw", "dis_processed_weighted",
                           "dis_selected_raw", "dis_selected_weighted"};
    const double values[] = {hist_all[imat].mu_p->GetEntries(), weightedMuons[imat],
                             hist_filt[imat].mu_p->GetEntries(), weightedSelectedMuons[imat],
                             double(processedDIS[imat]), weightedProcessedDIS[imat],
                             double(selectedDIS[imat]), weightedDIS[imat]};
    for (int bin = 1; bin <= 8; ++bin) {
      counts->GetXaxis()->SetBinLabel(bin, labels[bin - 1]);
      counts->SetBinContent(bin, values[bin - 1]);
    }
  }
}

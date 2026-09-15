#include "MuDISFilter.h"

#include <TFile.h>
#include <TTree.h>

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

  h.dis_vxz = new TH2D(
		       Form("vertex_x_vs_z_%s",label.Data()),
		       Form("DIS vertex x vs z - %s;z [cm];x [cm]", mat.Data()),
		       3000, 3000., 9000.,
		       600, -300., 300.
		       );
  h.dis_vyz = new TH2D(
		       Form("vertex_y_vs_z_%s",label.Data()),
		       Form("DIS vertex y vs z - %s;z [cm];y [cm]", mat.Data()),
		       3000, 3000., 9000.,
		       600, -300., 300.
		       );
  h.dis_vxy = new TH2D(
		       Form("vertex_y_vs_x_%s",label.Data()),
		       Form("DIS vertex y vs x - %s;x [cm];y [cm]", mat.Data()),
		       600, -300., 300.,
		       600, -300., 300.
		       );
  h.dis_pdg = new TH1I(
		       Form("daughter_pdg_%s",label.Data()),
		       Form("PDG code of DIS daughters - %s;PDG code;Particles", mat.Data()),
		       12001, -6000.5, 6000.5
		       );
  h.dis_pdgGrouped = new TH1I(
                             Form("daughter_pdg_grouped_%s",label.Data()),
                             Form("DIS daughter species - %s;Species;Particles", mat.Data()),
                             11, 0.5, 11.5
                             );
  const char* labels[] = {"e^{+}", "e^{-}", "#mu^{+}", "#mu^{-}", "#gamma",
                          "#nu/#bar{#nu}", "#pi^{#pm}", "h^{#pm}", "#pi^{0}", "h^{0}", "other"};
  for (int bin = 1; bin <= 11; ++bin)
    h.dis_pdgGrouped->GetXaxis()->SetBinLabel(bin, labels[bin - 1]);
  h.dis_n = new TH1I(
		     Form("n_daughters_%s",label.Data()),
		     Form("Number of DIS daughters - %s;multiplicity;Event", mat.Data()),
		     50,0,50
		     );
  h.dis_nCharged = new TH1I(
			    Form("n_daughters_charged_%s",label.Data()),
			    Form("Number of charged DIS daughters - %s;charged multiplicity;Event", mat.Data()),
			    50,0,50
			    );
  h.dis_pChargedFrac = new TH1D(
				Form("pfrac_charged_%s",label.Data()),
				Form("Charged fraction - %s;p(charged)/p(all);Event", mat.Data()),
				101,0,1.01
				);
  h.dis_pMuFrac = new TH1D(
			   Form("pfrac_mu_%s",label.Data()),
			   Form("p fraction of outgoing mu - %s;p(mu)/p(all);Event", mat.Data()),
			   101,0,1.01
			   );
  h.mu_p = new TH1D(
		    Form("muon_p_%s",label.Data()),
		    Form("muon momentum - %s;p_{#mu,in} [GeV]; Input mu events", mat.Data()),
		    400, 0., 400.
		    );
  h.mu_pt = new TH1D(
		     Form("muon_pt_%s",label.Data()),
		     Form("muon p_{T} - %s;p_{T,#mu,in} [GeV];Input mu events", mat.Data()),
		     100, 0., 10.
		     );
  h.mu_ppt = new TH2D(
		      Form("muon_pt_vs_p_%s",label.Data()),
		      Form("muon p_{T} vs p - %s;p_{#mu,in} [GeV]; p_{T,#mu,in} [GeV];Input mu events", mat.Data()),
		      400,0,400,
		      100, 0., 10.
		      );
  h.mu_ndis = new TH1I(
		       Form("muon_n_dis_%s",label.Data()),
		       Form("Number of DIS events - %s;DIS multiplicity;Input mu events", mat.Data()),
		       1001,0,1001
		       );
  h.mu_wdis = new TH1D(
		       Form("muon_vtx_weight_%s",label.Data()),
		       Form("DIS vertex weight - %s;wDIS;Input mu events", mat.Data()),
		       100, 0., 1.
		       );
  h.mu_wdis->SetCanExtend(TH1::kXaxis);
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

  Bool_t treeOK = InitFiles(input, fstartEvt);

  if (!treeOK) {
    throw std::runtime_error("MuDISFilter: failed to initialize input files");
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
  const auto* pdg = fPDG->GetParticle(particle.pid);
  // Nuclear PDG codes encode Z in digits 5--7 (10LZZZAAAI).
  if (!pdg && std::abs(particle.pid) >= 1000000000)
    return (std::abs(particle.pid) / 10000) % 1000 != 0;
  return pdg && pdg->Charge() != 0;
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

bool MuDISFilter::PassFilter(const std::vector<DISparticle>& daughters) const {
  if (fFilter) return fFilter(daughters);
  unsigned charged = 0;
  for (const auto& p : daughters) {
    if ((!fIncludeMuons && std::abs(p.pid) == 13) || !IsCharged(p)) continue;
    ++charged;
  }
  return charged >= fMinChargedDaughters;
}

void MuDISFilter::FillDIS(Histograms& h, const MuonDISInBranches& br,
                         int idis, const std::vector<DISparticle>& daughters) {
  h.dis_vxz->Fill(br.DISvz->at(idis), br.DISvx->at(idis));
  h.dis_vyz->Fill(br.DISvz->at(idis), br.DISvy->at(idis));
  h.dis_vxy->Fill(br.DISvx->at(idis), br.DISvy->at(idis));
  unsigned charged = 0;
  double totalP = 0., chargedP = 0., muonP = 0.;
  for (const auto& p : daughters) {
    h.dis_pdg->Fill(p.pid);
    h.dis_pdgGrouped->Fill(DaughterCategory(p));
    const double momentum = std::sqrt(p.px*p.px + p.py*p.py + p.pz*p.pz);
    totalP += momentum;
    if (IsCharged(p)) { ++charged; chargedP += momentum; }
    if (std::abs(p.pid) == 13) muonP += momentum;
  }
  h.dis_n->Fill(daughters.size());
  h.dis_nCharged->Fill(charged);
  if (totalP > 0.) {
    h.dis_pChargedFrac->Fill(chargedP / totalP);
    h.dis_pMuFrac->Fill(muonP / totalP);
  }
}

void MuDISFilter::ProcessEvents() {
  if (!ftree || !fouttree) throw std::runtime_error("Initialize input and output first");
  if (fstartEvt < 0) throw std::invalid_argument("Negative start event");
  const Long64_t end = fnEvts >= 0
      ? std::min(Long64_t(fstartEvt) + fnEvts, ftree->GetEntries())
      : ftree->GetEntries();
  Long64_t selected = 0, skipped = 0;
  for (Long64_t event = fstartEvt; event < end; ++event) {
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
    foutEv.mcTrks = *finEv.mcTrks;
    foutEv.sbtPt = *finEv.sbtPt;
    foutEv.ubtPt = *finEv.ubtPt;
    foutEv.sstPt = *finEv.sstPt;
    const auto& muon = finEv.mcTrks->at(0);
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
        FillDIS(hist_all[imat], in, idis, daughters);
        if (!PassFilter(daughters)) continue;
        FillDIS(hist_filt[imat], in, idis, daughters);
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
      auto fillMuon = [&muon, &in](Histograms& h, int count) {
        h.mu_p->Fill(muon.GetP());
        h.mu_pt->Fill(muon.GetPt());
        h.mu_ppt->Fill(muon.GetP(), muon.GetPt());
        h.mu_ndis->Fill(count);
        h.mu_wdis->Fill(in.wDIS);
      };
      fillMuon(hist_all[imat], in.nDISevts);
      if (out.nDISevts > 0) fillMuon(hist_filt[imat], out.nDISevts);
    }
    if (keep) {
      if (fouttree->Fill() < 0) throw std::runtime_error("Failed writing MuonDIS entry");
      ++selected;
    }
  }
  LOG(info) << "MuDISFilter: saved " << selected << " muon entries; skipped "
            << skipped << " unreadable or malformed entries.";
}

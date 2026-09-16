#ifndef SHIPMuDIS_MUDISFILTER_H_
#define SHIPMuDIS_MUDISFILTER_H_

#include "FairLogger.h"  // for FairLogger, MESSAGE_ORIGIN
#include "MuDISDefs.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TChain.h"  // for TTree
#include "TVector3.h"
#include <vector>
#include <functional>
#include <utility>
#include <limits>
#include <memory>
#include "TH1.h"
#include "TH2.h"

using namespace ShipMuDIS;

class MagneticTrackPropagator;
class ShipBFieldMap;
class TGeoManager;

struct Histograms {
  TH2D* dis_vxz = nullptr;
  TH2D* dis_vyz = nullptr;
  TH2D* dis_vxy = nullptr;

  TH1I* dis_pdg = nullptr;
  TH1I* dis_pdgGrouped = nullptr;
  TH1I* dis_n = nullptr;
  TH1I* dis_nCharged = nullptr;
  TH1D* dis_pChargedFrac = nullptr;
  TH1D* dis_pMuFrac = nullptr;

  TH2D* mu_ppt = nullptr;
  TH1D* mu_p = nullptr;
  TH1D* mu_pt = nullptr;
  TH1I* mu_ndis = nullptr;
  TH1D* mu_wdis = nullptr;

};

class MuDISFilter {
 public:
  /** default constructor **/
  MuDISFilter();

  /** destructor **/
  ~MuDISFilter();
  MuDISFilter(const MuDISFilter&) = delete;
  MuDISFilter& operator=(const MuDISFilter&) = delete;

  // Generic filter on daughter particles
  using Filter = std::function<bool(const std::vector<DISparticle>&)>;
  void SetFilter(Filter filter) { fFilter = std::move(filter); }
  //for default filter
  void SetMinChargedDaughters(unsigned minimum) { fMinChargedDaughters = minimum; }
  void SetIncludeMuons(bool include) { fIncludeMuons = include; }
  void SetDetectorAcceptance(ShipBFieldMap* field, TGeoManager* geometry,
                            double z = std::numeric_limits<double>::quiet_NaN());
  void SetUseDetectorAcceptance(bool enabled) { fUseDetectorAcceptance = enabled; }
  double GetDetectorZ() const { return fDetectorZ; }

  bool PassFilter(const std::vector<DISparticle>& daughters) const;
  bool PassFilter(const std::vector<DISparticle>& daughters, const TVector3& vertex) const;

  Histograms BookHistograms(TDirectory* dir, const TString& mat, const TString & label="");
  void init(const int& aEvts, const int& aStart);

  Bool_t InitFile(const char*, int);
  Bool_t InitFile(const char*);
  Bool_t InitFiles(const std::vector<std::string>&, int);
  Bool_t InitFiles(const std::vector<std::string>&);
  void process_file(const std::string& input, const std::string& output);
  void process_file(const std::vector<std::string>& input, const std::string& output);
  void initEvent();
  void ProcessEvents();



 private:
  bool IsCharged(const DISparticle& particle) const;
  double Charge(const DISparticle& particle) const;
  int DaughterCategory(const DISparticle& particle) const;
  void FillDIS(Histograms& h, const ShipMuDIS::MuonDISInBranches& br, int idis,
               const std::vector<DISparticle>& daughters, Histograms* filtered);
  unsigned fMinChargedDaughters = 2;
  bool fIncludeMuons = true;
  bool fUseDetectorAcceptance = true;
  double fDetectorZ = std::numeric_limits<double>::quiet_NaN();
  std::unique_ptr<MagneticTrackPropagator> fPropagator;  //! Runtime transport
  Filter fFilter;  //! User-supplied runtime predicate
  TChain* ftree;
  ShipMuDIS::MuonInBranches finEv;

  TTree* fouttree;
  ShipMuDIS::MuonBranches foutEv;

  FairLogger* fLogger;  //!   don't make it persistent, magic ROOT command
  int fnEvts;
  int fstartEvt;

  TDatabasePDG* fPDG;
  Histograms hist_all[nMats];
  Histograms hist_filt[nMats];

};

#endif

//
//  Zinv background prediction for RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#ifndef RA2BZINVANALYSIS_H
#define RA2BZINVANALYSIS_H

#define _CRT_SECURE_NO_WARNINGS

#define ISSKIM
#define ISV12
/* #define ISMC */

#include <TString.h>
#include <TChain.h>
#include <TTreeReaderValue.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTreeFormula.h>
#include <TChainElement.h>
#include "../../Analysis/btag/BTagCorrector.h"

// members needed by nested class cutHistos
static TString HTcut_;
static TString MHTcut_;
static TString NJetscut_;
static TString objcut_;
static TString minDphicut_;
static TString commonCuts_;
static TString trigCuts_;
static TString ptCut_;
static TString massCut_;
static TString photonDeltaRcut_;

class RA2bZinvAnalysis {

public:
  RA2bZinvAnalysis();
  RA2bZinvAnalysis(const std::string& cfg_filename);
  virtual ~RA2bZinvAnalysis() {};

  TChain* getChain(const char* sample, Int_t* fCurrent = nullptr, bool setBrAddr = true);
  std::vector<TString> fileList(TString sampleKey);
  std::vector<TH1F*> makeHistograms(const char* sample);
  TH1F* makeCChist(const char* sample);
  TCut getCuts(const TString sampleKey);
  int kinBin(double& ht, double& mht);
  void checkTrigPrescales(const char* sample);
  void checkActiveTrigPrescales(const char* sample);
  void runMakeClass(const char* sample, const char* ext);

  struct hist1D {
    TH1F* hist;
    TString name;
    const char* title;
    std::pair<const char*, const char*> axisTitles;
    Int_t Nbins;
    Double_t lowEdge;
    Double_t highEdge;
    Double_t* dvalue;
    Int_t* ivalue;
    void (RA2bZinvAnalysis::*filler)(TH1F* h, double wt);
    std::vector<TString*> omitCuts;
    const char* addCuts;
    TString NminusOneCuts;
    TTreeFormula* NminusOneFormula;
    hist1D() : dvalue(nullptr), ivalue(nullptr), filler(nullptr), addCuts("") {}
  };

  class cutHistos {
  public:
    cutHistos(TChain* chain, TObjArray* forNotify);
    ~cutHistos() {};
    void fill(TH1F* hcf, Double_t wt);
  private:
    TObjArray* forNotify_;
    TTreeFormula* HTcutf_;
    TTreeFormula* MHTcutf_;
    TTreeFormula* NJetscutf_;
    TTreeFormula* minDphicutf_;
    TTreeFormula* objcutf_;
    TTreeFormula* ptcutf_;
    TTreeFormula* masscutf_;
    TTreeFormula* trigcutf_;
    TTreeFormula* commoncutf_;
  };

private:
#ifdef ISV12
  const TString ntupleVersion_ = "V12";
#else
  const TString ntupleVersion_ = "V15";
#endif
#ifdef ISMC
  const bool isMC_ = true;
#else
  const bool isMC_ = false;
#endif
#ifdef ISSKIM
  const bool isSkim_ = true;
#else
  const bool isSkim_ = false;
#endif
  std::string treeName_;
  std::string treeLoc_;
  std::string fileListsFile_;
  std::string era_;  // "2016", ...
  double intLumi_;
  std::string deltaPhi_;  // "nominal", "hdp", "ldp"
  bool applyMassCut_;
  bool applyPtCut_;
  bool applyMinDeltaRCut_;
  // bool applySF_;
  // bool njSplit_;
  bool useTreeCCbin_;
  bool applyBTagSF_;
  bool applyPuWeight_;
  bool customPuWeight_;
  TH1* puHist_;
  const char* BTagSFfile_;
  std::vector< std::vector<double> > kinThresholds_;
  std::vector<int> nJetThresholds_;
  std::vector<int> nbThresholds_;
  unsigned kinSize_;

#ifdef ISMC

#ifdef ISV12
#include "LeafDeclaration_MC_V12.h"
#endif

#else  // ISMC

#ifdef ISV12
#include "LeafDeclaration_data_V12.h"
#else
  // V15
#ifndef ISSKIM
#include "LeafDeclaration_unskimmed_data_V15.h"
#endif

  // Declare dummy tree variables missing in some versions
#endif  // !ISMC
  Double_t        puWeight;
  Double_t        Weight;
  Double_t        TrueNumInteractions;
#endif

#ifndef ISSKIM
  UInt_t          RA2bin;
#endif

#ifdef ISV12
  Int_t           NElectrons;
  Int_t           NMuons;
#endif

  typedef std::map<TString, std::vector<TString> > vstring_map;
  typedef std::map<TString, TString> string_map;
  typedef std::map<std::vector<int>, Int_t> ivector_map;
  vstring_map triggerMap_;
  string_map objCutMap_;
  string_map minDphiCutMap_;
  string_map MHTCutMap_;
  string_map sampleKeyMap_;
  ivector_map toCCbin_;
  std::vector<const char*> activeBranches_;

  void Init(const std::string& cfg_filename="");
  void fillCutMaps();
  void bookAndFillHistograms(const char* sample, std::vector<hist1D*>& histograms);
  void fillCutFlow(TH1F* hcf, Double_t wt);

  void cleanVars() {
#ifndef ISSKIM
    NJets = NJetsclean;
    BTags = BTagsclean;
    HT = HTclean;
    MHT = MHTclean;
    JetID = JetIDclean;
    Jets = Jetsclean;
    Jets_hadronFlavor = Jetsclean_hadronFlavor;
    Jets_HTMask = Jetsclean_HTMask;
    isoElectronTracks = isoElectronTracksclean;
    isoMuonTracks = isoMuonTracksclean;
    isoPionTracks = isoPionTracksclean;
    DeltaPhi1 = DeltaPhi1clean;
    DeltaPhi2 = DeltaPhi2clean;
    DeltaPhi3 = DeltaPhi3clean;
    DeltaPhi4 = DeltaPhi4clean;
#endif
  };

  // Functions to fill histograms with non-double, non-int types
  void fillZmass(TH1F* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.M(), wt);}
  void fillZpt(TH1F* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.Pt(), wt);}

  ClassDef(RA2bZinvAnalysis, 1) // 2nd arg is ClassVersionID
};

#endif

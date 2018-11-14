//
//  Zinv background prediction for RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#ifndef RA2BZINVANALYSIS_H
#define RA2BZINVANALYSIS_H

#define _CRT_SECURE_NO_WARNINGS

#define VERSION 15
#define ISSKIM
/* #define ISMC */

#include "CCbinning.h"
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
  RA2bZinvAnalysis(const std::string& cfg_filename, const std::string& runBlock = "");
  virtual ~RA2bZinvAnalysis() {};

  TChain* getChain(const char* sample, Int_t* fCurrent = nullptr, bool makeClass = false);
  std::vector<TString> fileList(TString sampleKey);
  std::vector<TH1*> makeHistograms(const char* sample);
  void dumpSelEvIDs(const char* sample, const char* idFileName);
  TCut getCuts(const TString sampleKey);
  int kinBin(double& ht, double& mht);
  void checkTrigPrescales(const char* sample);
  void checkActiveTrigPrescales(const char* sample);
  void runMakeClass(const std::string& sample);

  struct histConfig {
    // 1D or 2D histogram; select by value of NbinsY, = 0 for 1D.
    TH1* hist;
    TString name;
    const char* title;
    std::pair<const char*, const char*> axisTitles;
    Int_t NbinsX;
    std::pair<Double_t, Double_t> rangeX;
    Int_t NbinsY;
    std::pair<Double_t, Double_t> rangeY;
    Double_t* dvalue;
    Int_t* ivalue;
    void (RA2bZinvAnalysis::*filler1D)(TH1F* h, double wt);
    void (RA2bZinvAnalysis::*filler2D)(TH2F* h, double wt);
    std::vector<TString*> omitCuts;
    const char* addCuts;
    TString NminusOneCuts;
    TTreeFormula* NminusOneFormula;
  histConfig() : dvalue(nullptr), ivalue(nullptr), filler1D(nullptr), filler2D(nullptr), addCuts(""), NbinsY(0) {}
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
#if VERSION == 12
  const TString ntupleVersion_ = "V12";
#elif  VERSION == 15
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
  int verbosity_;
  std::string treeName_;
  std::string treeLoc_;
  std::string fileListsFile_;
  std::string runBlock_;
  std::string era_;  // "2016", ...
  double intLumi_;
  std::string deltaPhi_;  // "nominal", "hdp", "ldp"
  bool applyMassCut_;
  bool applyPtCut_;
  bool applyMinDeltaRCut_;
  // bool applySF_;
  // bool njSplit_;
  bool useTreeCCbin_;
  bool useDeepCSV_;
  bool applyBTagSF_;
  bool applyPuWeight_;
  bool customPuWeight_;
  TH1* puHist_;
  const char* BTagSFfile_;
  BTagCorrector* btagcorr_;
  std::vector< std::vector<double> > kinThresholds_;
  std::vector<int> nJetThresholds_;
  std::vector<int> nJet1Thresholds_;
  std::vector<int> nbThresholds_;
  unsigned kinSize_;
  TString isoSFlepTksVeto_;
  TString isoSFlepTksCut_;
  TString photonVeto_;
  TString photonCut_;
  double csvMthreshold_;

#ifdef ISMC

  #if VERSION == 12
  #include "LeafDeclaration_MC_V12.h"
  #endif

#else  // ISMC

  #if VERSION == 12
    #include "LeafDeclaration_data_V12.h"
  #elif VERSION == 15
    #ifdef ISSKIM
      #include "LeafDeclaration_data_V15.h"
    #else
      #include "LeafDeclaration_unskimmed_data_V15.h"
    #endif

  #endif  // VERSION

  // Declare dummy tree variables missing in some versions

  Double_t        puWeight;
  Double_t        Weight;
  Double_t        TrueNumInteractions;

#endif  // !ISMC

#ifndef ISSKIM
  UInt_t          RA2bin;
#endif

#if VERSION == 12
  Int_t           NElectrons;
  Int_t           NMuons;
  Int_t           BTagsDeepCSV;
  Int_t           ecalBadCalibFilter;
#endif

  typedef std::map<TString, std::vector<TString> > vstring_map;
  typedef std::map<TString, TString> string_map;
  vstring_map triggerMap_;
  string_map objCutMap_;
  string_map minDphiCutMap_;
  string_map MHTCutMap_;
  string_map sampleKeyMap_;
  CCbinning::ivector_map toCCbin_, toCCbinjb_, toCCbinSpl_, toCCbinJb_;
  std::vector<const char*> activeBranches_;

  void Init(const std::string& cfg_filename="");
  void fillCutMaps();
  void bookAndFillHistograms(const char* sample, std::vector<histConfig*>& histograms, TCut baselineCuts);
  void fillCutFlow(TH1F* hcf, Double_t wt);
  Int_t setBTags();

  void cleanVars() {
#ifndef ISSKIM
    NJets = NJetsclean;
    BTags = BTagsclean;
    BTagsDeepCSV = BTagsDeepCSVclean;
    HT = HTclean;
    HT5 = HT5clean;
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
  void fillCC(TH1F* h, double wt);
  void fillnZcand(TH1F* h, double wt) {h->Fill(ZCandidates->size(), wt);}
  void fillZmass(TH1F* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.M(), wt);}
  void fillZmassjb(TH1F* h, double wt);
  void fillZpt(TH1F* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.Pt(), wt);}
  void fillGpt(TH1F* h, double wt) {for (auto & theG : *Photons) h->Fill(theG.Pt(), wt);}
  void fillFilterCuts(TH1F* h, double wt);
  void fillZGmass(TH1F* h, double wt);
  void fillGJdR(TH1F* h, double wt);
  void fillZGdRvsM(TH2F* h, double wt);
  void fillGLdRnoPixelSeed(TH1F* h, double wt);
  void fillGLdRpixelSeed(TH1F* h, double wt);

  ClassDef(RA2bZinvAnalysis, 1) // 2nd arg is ClassVersionID
};

#endif

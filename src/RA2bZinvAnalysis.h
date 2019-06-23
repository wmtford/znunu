//
//  Make histograms for Zinv background prediction in the RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#ifndef RA2BZINVANALYSIS_H
#define RA2BZINVANALYSIS_H

#ifndef NtupleClass_cxx
#define NtupleClass_cxx
#endif

#define _CRT_SECURE_NO_WARNINGS

#include "NtupleClass.h"
#include "CCbinning.h"
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>
#include <TTreeFormula.h>
#include "../../Analysis/btag/BTagCorrector.h"

#include <TMath.h>
using namespace TMath;

void NtupleClass::Loop() {}

// members needed by nested classes
static TString HTcut_;
static TString MHTcut_;
static TString NJetscut_;
static TString objcut_;
static TString minDphicut_;
static TString commonCuts_;
static TString ptCut_;
static TString massCut_;
static TString photonDeltaRcut_;

class RA2bZinvAnalysis : public NtupleClass {

public:
  RA2bZinvAnalysis();
  RA2bZinvAnalysis(const std::string& cfg_filename, const std::string& runBlock = "");
  virtual ~RA2bZinvAnalysis() {};

  Bool_t Notify() override {newFileInChain_ = kTRUE;  return(kTRUE);};
  std::vector<TH1*> makeHistograms(const char* sample);
  void dumpSelEvIDs(const char* sample, const char* idFileName);
  TCut getCuts(const TString sampleKey);
  void setTriggerIndexList(const char* sample);
  double getPtZ() {
    if (!isMC_) return -1;
    for (int iGen = 0, nGen =  GenParticles_PdgId->size(); iGen < nGen; ++iGen) {
      if (GenParticles_PdgId->at(iGen)==23 && GenParticles_Status->at(iGen) == 62)
	return GenParticles->at(iGen).Pt();
    }
    return -1;
  };
  void checkTrigPrescales(const char* sample);

  enum yearFirstRun {Start2016 = 271036, Start2017 = 294645, Start2018 = 315252, StartHEM = 319077, Start2018C = 319313};
  // First HEM run 319077
  enum runYear{Year2016 = 0, Year2017 = 1, Year2018 = 2, Year2018HEP = 3, Year2018HEM = 4};

  struct histConfig {
    // 1D or 2D histogram; select by value of NbinsY, = 0 for 1D.
    TH1* hist;
    TString name;
    const char* title;
    std::pair<const char*, const char*> axisTitles;
    Int_t NbinsX;
    std::pair<Double_t, Double_t> rangeX;
    Double_t* binsX;
    Int_t NbinsY;
    std::pair<Double_t, Double_t> rangeY;
    Double_t* binsY;
    Double_t* dvalue;
    Int_t* ivalue;
    void (RA2bZinvAnalysis::*filler1D)(TH1D* h, double wt);
    void (RA2bZinvAnalysis::*filler2D)(TH2F* h, double wt);
    std::vector<TString*> omitCuts;
    const char* addCuts;
    TString NminusOneCuts;
    TTreeFormula* NminusOneFormula;
  histConfig() : binsX(nullptr), binsY(nullptr), dvalue(nullptr), ivalue(nullptr),
      filler1D(nullptr), filler2D(nullptr), addCuts(""), NbinsY(0) {}
  };

  class cutHistos {
  public:
    cutHistos(TChain* chain, TObjArray* forNotify);
    ~cutHistos() {};
    void setAxisLabels(TH1D* hcf);
    void fill(TH1D* hcf, Double_t wt, bool passTrg, bool passHEM);
  private:
    TObjArray* forNotify_;
    TTreeFormula* HTcutf_;
    TTreeFormula* MHTcutf_;
    TTreeFormula* NJetscutf_;
    TTreeFormula* minDphicutf_;
    TTreeFormula* objcutf_;
    TTreeFormula* ptcutf_;
    TTreeFormula* masscutf_;
    TTreeFormula* photonDeltaRcutf_;
    TTreeFormula* commoncutf_;
  };

  class efficiencyAndPurity {
  public:
    efficiencyAndPurity() : deltaPhi_("nominal") {};
    efficiencyAndPurity(string deltaPhi) : deltaPhi_(deltaPhi) {};
    ~efficiencyAndPurity() {};
    void openFiles();
    void getHistos(const char* sample, int currentYear);
    pair<double, double> weight(CCbinning* CCbins,
				Int_t NJets, Int_t BTags, Double_t MHT, Double_t HT,
				vector<TLorentzVector> ZCandidates,
				vector<TLorentzVector> Photons,
				vector<TLorentzVector> Electrons,
				vector<TLorentzVector> Muons,
				vector<double> EBphoton,
				bool applyDRfitWt,
				int currentYear);

    double prefiring_weight_photon(vector<TLorentzVector>* Photons, unsigned p){
    /* double prefiring_weight_photon(unsigned p){ */
      if (hPrefiring_photon_ == nullptr) return 1;
      return (1 - hPrefiring_photon_->GetBinContent(hPrefiring_photon_->GetXaxis()->FindBin(Photons->at(p).Eta()),
						    hPrefiring_photon_->GetYaxis()->FindBin(Photons->at(p).Pt())));
    };
    double prefiring_weight_electron(vector<TLorentzVector>* Electrons, unsigned e){
    /* double prefiring_weight_electron(unsigned e){ */
      if (hPrefiring_photon_ == nullptr) return 1;
      return (1 - hPrefiring_photon_->GetBinContent(hPrefiring_photon_->GetXaxis()->FindBin(Electrons->at(e).Eta()),
						    hPrefiring_photon_->GetYaxis()->FindBin(Electrons->at(e).Pt())));
    };
    double prefiring_weight_jet(vector<TLorentzVector>* Jets, unsigned j){
    /* double prefiring_weight_jet(unsigned j){ */
      if (hPrefiring_jet_ == nullptr) return 1;
      return (1 - hPrefiring_jet_->GetBinContent(hPrefiring_jet_->GetXaxis()->FindBin(Jets->at(j).Eta()),
						 hPrefiring_jet_->GetYaxis()->FindBin(Jets->at(j).Pt()))) ;
    };
    double quadSum(double x, double y) {return Sqrt(Power(x,2) + Power(y,2));};
    double quadSum(double x, double y, double z) {return Sqrt(Power(x,2) + Power(y,2) + Power(z,2));};
  private:
    std::vector<TFile*> purityTrigEffFile_;
    std::vector<TFile*> photonTrigEffFile_;
    std::vector<TFile*> photonSFFile_;
    //std::vector<TFile*> elecSFFile_;
    std::vector<TFile*> elecIDandIsoSFFile_;
    std::vector<TFile*> elecRecoLowSFFile_;
    std::vector<TFile*> elecRecoHighSFFile_;
    std::vector<TFile*> muonIDSFFile_;
    std::vector<TFile*> muonIsoSFFile_;
    TFile* prefiringWeightFile_;
    TH2F *hPrefiring_photon_, *hPrefiring_jet_;
    TString theSample_;
    string deltaPhi_;
    std::vector<TH1F*> hPurity_, hTrigEff_;
    std::vector<TF1*> fTrigEff_;
    std::vector<TEfficiency*> eTrigEff_;
    TF1* DRfun_ = new TF1;
    std::vector< std::vector<Double_t> > DRpars_;

    /* TH1F* hSFeff_; */
    std::vector<TH2F*> hSFeff_;
    /* TH1D* FdirHist_; */
    TGraphErrors* FdirGraph_;
  };

private:
  std::string era_;  // "2016", "Run2"
  TString ntupleVersion_;
  bool isSkim_;
  bool isMC_;
  string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  int verbosity_;
  string rootVerbosity_;
  string treeName_;
  string treeLoc_;
  string fileListsFile_;
  string runBlock_;
  double intLumi_;
  bool applyMassCut_;
  bool applyPtCut_;
  bool useDeepCSV_;
  bool applyBTagSF_;
  bool applyPuWeight_;
  bool customPuWeight_;
  bool applyHEMjetVeto_;
  bool applyZptWt_;
  bool applyDRfitWt_;
  bool applySFwtToMC_;
  Bool_t newFileInChain_;
  TH1* puHist_;
  CCbinning* CCbins_;
  BTagCorrector* btagcorr_;
  const char* BTagSFfile_;
  TString isoSFlepTksVeto_;
  TString isoSFlepTksCut_;
  TString photonVeto_;
  TString photonCut_;
  double csvMthreshold_;
  double effWt_, effSys_;

  typedef std::map<TString, std::vector<TString> > vstring_map;
  typedef std::map<TString, TString> string_map;
  std::vector<unsigned> triggerIndexList_;
  vstring_map triggerMapByName_;
  string_map objCutMap_;
  string_map minDphiCutMap_;
  string_map MHTCutMap_;
  string_map sampleKeyMap_;

  void Config(const std::string& cfg_filename="");
  void getChain(const char* dataSet);
  void setActiveBranches(const bool activateAll = false);
  std::vector<TString> fileList(TString sampleKey);
  void fillCutMaps();
  void bookAndFillHistograms(const char* sample, std::vector<histConfig*>& histograms, TCut baselineCuts);
  void fillCutFlow(TH1D* hcf, Double_t wt);
  Int_t setBTags(int runYear);
  efficiencyAndPurity* effPurCorr_;

  bool testHEM() {
    // HEM veto for data depending on RunNum; for MC weight by lumi unless
    // forced by a substring of runBlock_ (HE Present or Missing)
    if (!isMC_ && RunNum < StartHEM) return false;
    if (isMC_ && runBlock_.find("2018") == std::string::npos) return false;
    if (isMC_ && runBlock_.find("HEP") != std::string::npos) return false;
    if (isMC_ && runBlock_.find("HEM") == std::string::npos
	&& EvtNum % 1000 < 1000*21.0/59.6) return false;
    return true;
  }

  bool passHEMobjVeto(TLorentzVector& obj, double ptThresh = 0, bool extendedHEM = true) {
    if (!testHEM()) return true;
    double etalo, etahi, philo, phihi;
    if (extendedHEM) {etalo = -3.2; etahi = -1.2; philo = -1.77; phihi = -0.67;}
    else {etalo = -3.0; etahi = -1.4; philo = -1.57; phihi = -0.87;}
    if (etalo <= obj.Eta() && obj.Eta() <= etahi &&
	philo <= obj.Phi() && obj.Phi() <= phihi &&
	obj.Pt() > ptThresh)
      return false;
    else return true;
  };

  bool passHEMjetVeto(double ptThresh = 30, double dPhiThresh = 0.5, TH1D* hg = nullptr, double wt = 1) {
    if (!testHEM()) return true;
    for (auto & jet : *Jets) {
      if (!passHEMobjVeto(jet, ptThresh)) {
	Double_t dPhi = TVector2::Phi_mpi_pi(jet.Phi() - MHTPhi);
	if (hg != nullptr) hg->Fill(dPhi, wt);
	if (abs(dPhi) < dPhiThresh) return false;
      }
    }
    return true;
  };

  void cleanVars() {
    if (isSkim_) return;
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
  };

  // Functions to fill histograms with non-double, non-int types
  void fillFilterCuts(TH1D* h, double wt);
  void fillCC(TH1D* h, double wt);
  void fillnZcand(TH1D* h, double wt) {h->Fill(ZCandidates->size(), wt);}
  void fillZmass(TH1D* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.M(), wt);}
  void fillZmassjb(TH1D* h, double wt);
  void fillPUwtvsNint(TH2F* h, double wt) {h->Fill(TrueNumInteractions, puWeight, wt);}
  void fillHT_DR_xWt(TH1D* h, double wt) {h->Fill(HT, wt*HT);}
  void fillMHT_DR_xWt(TH1D* h, double wt) {h->Fill(MHT, wt*MHT);}
  void fillNJets_DR_xWt(TH1D* h, double wt) {h->Fill(Double_t(NJets), wt*NJets);}
  void fillSFwt_DR(TH1D* h, double wt) {double wtt = effWt_ > 0 ? wt/effWt_ : wt;  h->Fill(effWt_, wtt);}
  void fillSFsys_DR(TH1D* h, double wt) {double wtt = effWt_ > 0 ? wt/effWt_ : wt;  h->Fill(effSys_, wtt);}
  void fillgenZpt(TH1D* h, double wt) {h->Fill(getPtZ(), wt);}
  void fillZpt(TH1D* h, double wt) {for (auto & theZ : *ZCandidates) h->Fill(theZ.Pt(), wt);}
  void fillPhotonPt(TH1D* h, double wt) {for (auto & theG : *Photons) h->Fill(theG.Pt(), wt);}
  void fillMuonEta(TH1D* h, double wt) {for (auto & theMu : *Muons) h->Fill(theMu.Eta(), wt);}
  void fillElectronEta(TH1D* h, double wt) {for (auto & theE : *Electrons) h->Fill(theE.Eta(), wt);}
  void fillJetDphi(TH1D* h, double wt) {bool passHEM = passHEMjetVeto(30, 1, h, wt);}
  void fillPhotonEta(TH1D* h, double wt) {for (auto & theG : *Photons) h->Fill(theG.Eta(), wt);}
  void fillGpt(TH1D* h, double wt) {for (auto & theG : *Photons) h->Fill(theG.Pt(), wt);}
  void fillZGmass(TH1D* h, double wt);
  void fillGJdR(TH1D* h, double wt);
  void fillZGdRvsM(TH2F* h, double wt);
  void fillGLdRnoPixelSeed(TH1D* h, double wt);
  void fillGLdRpixelSeed(TH1D* h, double wt);

  /* ClassDef(RA2bZinvAnalysis, 1) // 2nd arg is ClassVersionID */
};

#endif

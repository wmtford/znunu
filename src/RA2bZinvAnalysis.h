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

#include "Ntuple.h"
#include "TreeConfig.h"
#include "CutManager.h"
#include "CCbinning.h"
#include "EfficWt.h"
#include "BTagSF.h"
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include "../../Analysis/btag/BTagCorrector.h"

#include <TMath.h>
using namespace TMath;

#include <time.h>
#include <sys/time.h>

class RA2bZinvAnalysis {

public:
  RA2bZinvAnalysis();
  RA2bZinvAnalysis(const std::string& cfg_filename, const std::string& runBlock = "");
  virtual ~RA2bZinvAnalysis() {
    delete treeConfig_;
    delete CCbins_;
    delete effPurCorr_;
  };

  std::vector<TH1*> makeHistograms(const char* sample);
  void dumpSelEvIDs(const char* sample, const char* idFileName);
  void checkTrigPrescales(const char* sample);

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
    std::vector<const TString*> omitCuts;
    const char* addCuts;
    TString NminusOneCuts;
    histConfig() : binsX(nullptr), binsY(nullptr), dvalue(nullptr), ivalue(nullptr),
      filler1D(nullptr), filler2D(nullptr), addCuts(""), NbinsY(0) {}
  };

  class cutHistos {
  public:
    cutHistos(Ntuple* tuple, CutManager* selector);
    ~cutHistos() {};
    void setAxisLabels(TH1D* hcf);
    void fill(TH1D* hcf, Double_t wt, bool zeroChg, bool passTrg, bool passHEM, bool passEcalNoiseJetFilter);
    void fill(TH1D* hcf, Double_t wt, bool DiLepton, bool passPhotonVeto);
  private:
    Ntuple* tuple_;
    CutManager* selector_;
    vector<TString> cutNames_, skimCutNames_;
  };

private:
  string era_;  // "2016", "Run2"
  TString ntupleVersion_;
  bool isSkim_;
  bool isMC_;
  string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  int verbosity_;
  string rootVerbosity_;
  string runBlock_;
  double intLumi_;
  bool applyMassCut_;
  bool applyPtCut_;
  int minNbCut_;
  bool useDeepCSV_;
  bool applyBTagSF_;
  bool applyPuWeight_;
  bool customPuWeight_;
  bool applyHEMjetVeto_;
  bool applyZptWt_;
  bool applyDRfitWt_;
  bool applySFwtToMC_;
  bool restrictClean_;
  Ntuple* Tupl;
  TreeConfig* treeConfig_;
  CCbinning* CCbins_;
  CutManager* evSelector_;
  EfficWt* effPurCorr_;
  BTagSF* btagsf_;
  double csvMthreshold_;
  TH1* puHist_;
  double effWt_, effSys_;

  void Config(const std::string& cfg_filename="");
  void bookAndFillHistograms(const char* sample,
			     vector<histConfig*>& histograms,
			     vector<histConfig*>& cutHistograms);
  void bookHistograms(vector<histConfig*>& histoList, TCut baselineCuts, cutHistos cutHistFiller);
  Int_t setBTags(int runYear);
  bool diLepton(Int_t NLeptons, vector<TLorentzVector>* Leptons, vector<bool>* Leptons_passIso,
		vector<int>* Leptons_charge, vector<bool>* Leptons_mediumID = nullptr);
  bool ecalNoiseJetFilter();
  void fillCutFlow(TH1D* hcf, Double_t wt);
  double get_cpu_time() {return (double)clock() / CLOCKS_PER_SEC;};

  double getGenPtZ() {
    if (!isMC_) return -1;
    for (int iGen = 0, nGen =  Tupl->GenParticles_PdgId->size(); iGen < nGen; ++iGen) {
      if (Tupl->GenParticles_PdgId->at(iGen)==23 && Tupl->GenParticles_Status->at(iGen) == 62)
	return Tupl->GenParticles->at(iGen).Pt();
    }
    return -1;
  };

  bool testHEM() {
    // HEM veto for data depending on RunNum; for MC weight by lumi unless
    // forced by a substring of runBlock_ (HE Present or Missing)
    if (!isMC_ && Tupl->RunNum < CutManager::StartHEM) return false;
    if (isMC_ && runBlock_.find("2018") == std::string::npos) return false;
    if (isMC_ && runBlock_.find("HEP") != std::string::npos) return false;
    if (isMC_ && runBlock_.find("HEM") == std::string::npos
	&& Tupl->EvtNum % 1000 < 1000*21.0/59.6) return false;
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
    for (auto & jet : *(Tupl->Jets)) {
      if (!passHEMobjVeto(jet, ptThresh)) {
	Double_t dPhi = TVector2::Phi_mpi_pi(jet.Phi() - Tupl->MHTPhi);
	if (hg != nullptr) hg->Fill(dPhi, wt);
	if (abs(dPhi) < dPhiThresh) return false;
      }
    }
    return true;
  };

  void cleanVars() {
    // For TreeMaker these replacements are made when skims are produced,
    // so this applies only to raw ntuple data sets.
    if (isSkim_) return;
    if (restrictClean_ && Tupl->Jetsclean->size() == 0) return;
    Tupl->NJets = Tupl->NJetsclean;
    Tupl->BTags = Tupl->BTagsclean;
    Tupl->BTagsDeepCSV = Tupl->BTagsDeepCSVclean;
    Tupl->HT = Tupl->HTclean;
    Tupl->HT5 = Tupl->HT5clean;
    Tupl->MHT = Tupl->MHTclean;
    Tupl->MHTPhi = Tupl->MHTPhiclean;
    Tupl->JetID = Tupl->JetIDclean;
    Tupl->Jets = Tupl->Jetsclean;
    Tupl->Jets_hadronFlavor = Tupl->Jetsclean_hadronFlavor;
    Tupl->Jets_HTMask = Tupl->Jetsclean_HTMask;
    Tupl->Jets_bDiscriminatorCSV = Tupl->Jetsclean_bDiscriminatorCSV;
    Tupl->isoElectronTracks = Tupl->isoElectronTracksclean;
    Tupl->isoMuonTracks = Tupl->isoMuonTracksclean;
    Tupl->isoPionTracks = Tupl->isoPionTracksclean;
    Tupl->DeltaPhi1 = Tupl->DeltaPhi1clean;
    Tupl->DeltaPhi2 = Tupl->DeltaPhi2clean;
    Tupl->DeltaPhi3 = Tupl->DeltaPhi3clean;
    Tupl->DeltaPhi4 = Tupl->DeltaPhi4clean;
  };

  // Functions to fill histograms with non-double, non-int types
  void fillFilterCuts(TH1D* h, double wgt);
  void fillCC(TH1D* h, double wt);
  void fillnZcand(TH1D* h, double wt) {h->Fill(Tupl->ZCandidates->size(), wt);}
  void fillZmass(TH1D* h, double wt) {for (auto & theZ : *(Tupl->ZCandidates)) h->Fill(theZ.M(), wt);}
  void fillZmassjb(TH1D* h, double wt);
  void fillPUwtvsNint(TH2F* h, double wt) {h->Fill(Tupl->TrueNumInteractions, Tupl->puWeight, wt);}
  void fillHT_DR_xWt(TH1D* h, double wt) {h->Fill(Tupl->HT, wt*Tupl->HT);}
  void fillMHT_DR_xWt(TH1D* h, double wt) {h->Fill(Tupl->MHT, wt*Tupl->MHT);}
  void fillNJets_DR_xWt(TH1D* h, double wt) {h->Fill(Double_t(Tupl->NJets), wt*Tupl->NJets);}
  void fillSFwt_DR(TH1D* h, double wt) {double wtt = effWt_ > 0 ? wt/effWt_ : wt;  h->Fill(effWt_, wtt);}
  void fillSFsys_DR(TH1D* h, double wt) {double wtt = effWt_ > 0 ? wt/effWt_ : wt;  h->Fill(effSys_, wtt);}
  void fillgenZpt(TH1D* h, double wt) {h->Fill(getGenPtZ(), wt);}
  void fillZpt(TH1D* h, double wt) {for (auto & theZ : *(Tupl->ZCandidates)) h->Fill(theZ.Pt(), wt);}
  void fillPhotonPt(TH1D* h, double wt) {for (auto & theG : *(Tupl->Photons)) h->Fill(theG.Pt(), wt);}
  void fillMuonEta(TH1D* h, double wt) {for (auto & theMu : *(Tupl->Muons)) h->Fill(theMu.Eta(), wt);}
  void fillElectronEta(TH1D* h, double wt) {for (auto & theE : *(Tupl->Electrons)) h->Fill(theE.Eta(), wt);}
  void fillJetDphi(TH1D* h, double wt) {bool passHEM = passHEMjetVeto(30, 1, h, wt);}
  void fillPhotonEta(TH1D* h, double wt) {for (auto & theG : *(Tupl->Photons)) h->Fill(theG.Eta(), wt);}
  void fillttZpt(TH1D* h, double wt) {for (auto & theZ : *(Tupl->ZCandidates))
      if (Tupl->NJets >= 4 && Tupl->BTags >= 2) h->Fill(theZ.Pt(), wt);}
  void fillGpt(TH1D* h, double wt) {for (auto & theG : *(Tupl->Photons)) h->Fill(theG.Pt(), wt);}
  void fillZGmass(TH1D* h, double wt);
  void fillGJdR(TH1D* h, double wt);
  void fillZGdRvsM(TH2F* h, double wt);
  void fillGLdRnoPixelSeed(TH1D* h, double wt);
  void fillGLdRpixelSeed(TH1D* h, double wt);

  /* ClassDef(RA2bZinvAnalysis, 1) // 2nd arg is ClassVersionID */
};

#endif

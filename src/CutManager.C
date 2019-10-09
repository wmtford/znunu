//
//  Define cuts for RA2b analysis
//

#include "CutManager.h"

// ClassImp(CutManager)

CutManager::CutManager(const TString sample, const TString ntupleVersion, bool isSkim, bool isMC, int verbosity,
		       std::string era, string deltaPhi, bool applyMassCut, bool applyPtCut, CCbinning* CCbins) :
  ntupleVersion_(ntupleVersion), isSkim_(isSkim), verbosity_(verbosity), isMC_(isMC), era_(era), deltaPhi_(deltaPhi),
  applyMassCut_(applyMassCut), applyPtCut_(applyPtCut), CCbins_(CCbins) {

  fillCutMaps();  // Depends on isMC_
  string_map::iterator iter = sampleKeyMap_.find(sample);
  if (iter == sampleKeyMap_.end()) {
    cout << "CutManager instantiated with invalid sample = " << sample << endl;
    cuts_ = "0";
    return;
  }
  TString sampleKey = iter->second;

  if (ntupleVersion_ == "V12" || ntupleVersion_ == "V15")
    commonCuts_ =    "globalTightHalo2016Filter==1";
  else
    commonCuts_ =    "globalSuperTightHalo2016Filter==1";
  commonCuts_ += " && HBHENoiseFilter==1";
  commonCuts_ += " && HBHEIsoNoiseFilter==1";
  commonCuts_ += " && EcalDeadCellTriggerPrimitiveFilter==1";
  commonCuts_ += " && BadChargedCandidateFilter";
  commonCuts_ += " && BadPFMuonFilter";
  commonCuts_ += " && NVtx > 0";
  if (ntupleVersion_ != "V12") {
    // commonCuts_ += " && ecalBadCalibFilter==1";  // Added for 94X
  }
  if (!isMC_) {
    commonCuts_ += " && eeBadScFilter==1";
  }
  // Kevin Pedro, re V16:  Please note that I have included only the JetID "event cleaning" cut in these skims.
  // The PFCaloMETRatio, HT5/HT, and noMuonJet cuts are left out,
  // so the impact of these cuts can be tested and refined if necessary.
  TString EcalNoiseJetFilterCut = "1";
  if (isSkim_ && !(ntupleVersion_ == "V12" || ntupleVersion_ == "V15")) {
    commonCuts_ += " && METRatioFilter";
    if (!isMC_) {
      char cutbuf[256];
      sprintf(cutbuf, "(RunNum < %7d || RunNum >= %7d || EcalNoiseJetFilter)", Start2017, Start2018);
      EcalNoiseJetFilterCut = cutbuf;
    }
    if (era_ == "2016")
      commonCuts_ += " && HTRatioFilter";
    else
      commonCuts_ += " && HTRatioDPhiFilter";
      // Commoncuts_ += " && HT5/HT <= (DeltaPhi1 - (-0.5875))/1.025";
      // above implements "DeltaPhi1 >= 1.025*HT5/HT - 0.5875", avoiding "+" and "*" for regexp;
      // Revised version is "(htratio < 1.2 ? true : (looper->DeltaPhi1 >= 
      //                     (tight ? 5.3*htratio - 4.78 : 1.025*htratio - 0.5875) ) )"
      // https://github.com/kpedro88/Analysis/blob/SUSY2017/KCode/KCommonSelectors.h
    // commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming (xV16), single lepton
  }
  commonCuts_ += " && " + EcalNoiseJetFilterCut;
  if (!isSkim_) {
    commonCuts_ += " && JetIDclean";
    commonCuts_ += " && PFCaloMETRatio < 5";  // METRatioFilter in skims
    if (era_ == "2016")
      commonCuts_ += " && HT5clean/HTclean <= 2";  // HTRatioFilter
    else
      // HTRatioDPhiFilter:
      commonCuts_ += " && ((HT5clean/HTclean < 1.2) || (HT5clean/HTclean <= (DeltaPhi1clean - (-0.5875))/1.025))";
    // FIXME:  For !isSkim && !isMC_, need to define EcalNoiseJetFilter
    // commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming (xV16), single lepton
    // commonCuts_ += " && noFakeJet";  // Defined in loop, applied in skimming FastSim
  }

  massCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyMassCut_)
    massCut_ = "@ZCandidates.size()==1 && ZCandidates.M()>=76.188 && ZCandidates.M()<=106.188";

  ptCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyPtCut_)
    ptCut_ = "@ZCandidates.size()==1 && ZCandidates.Pt()>=200.";

  photonDeltaRcut_ = "1";  // For cut histograms only; incorporated into objCutMap
  if (sampleKey == "photon") {
    if (applyPtCut_) ptCut_ = "@Photons.size()==1 && Photons.Pt()>=200.";
    if (isMC_) photonDeltaRcut_ = "madMinPhotonDeltaR>=0.4";
  }

  if (isSkim_) {
    HTcut_ = string("HT>=") + std::to_string(CCbins_->htThreshold(0, 0));
    NJetscut_ = string("NJets>=") + std::to_string(CCbins_->nJetThreshold(0));
  } else {
    HTcut_ = string("HTclean>=") + std::to_string(CCbins_->htThreshold(0, 0));
    NJetscut_ = string("NJetsclean>=") + std::to_string(CCbins_->nJetThreshold(0));
  }
  MHTcut_ = MHTCutMap_.at(deltaPhi_);
  objcut_ = objCutMap_.at(sampleKey);
  minDphicut_ = minDphiCutMap_.at(deltaPhi_);
  // HEMgeomcut_ = "(-3.0 < Electrons[0].Eta() && Electrons[0].Eta() < -1.4 && -1.57 < Electrons[0].Phi() && Electrons[0].Phi() < -0.87)
  //             || (-3.0 < Electrons[1].Eta() && Electrons[1].Eta() < -1.4 && -1.57 < Electrons[1].Phi() && Electrons[1].Phi() < -0.87)";
  // Extended HEM cut [-3.2, -1.2] [-1.77, -0.67] 


  // For test of same-flavor isoLeptonTracks veto
  isoSFlepTksVeto_ = "1";
  if (isSkim_) {
    if (sampleKey == "zmm") {isoSFlepTksVeto_ = "isoMuonTracks!=0";  isoSFlepTksCut_ = "isoMuonTracks==0";}
    if (sampleKey == "zee") {isoSFlepTksVeto_ = "isoElectronTracks!=0";  isoSFlepTksCut_ = "isoElectronTracks==0";}
  } else {
    if (sampleKey == "zmm") {isoSFlepTksVeto_ = "isoMuonTracksclean!=0";  isoSFlepTksCut_ = "isoMuonTracksclean==0";}
    if (sampleKey == "zee") {isoSFlepTksVeto_ = "isoElectronTracksclean!=0";  isoSFlepTksCut_ = "isoElectronTracksclean==0";}
  }
  photonVeto_ = "@Photons.size()==0";
  photonCut_ = "@Photons.size()!=0";

  cuts_ += objcut_;
  cuts_ += HTcut_;
  cuts_ += NJetscut_;
  cuts_ += MHTcut_;
  if (!isSkim_) cuts_ += minDphicut_;
  cuts_ += ptCut_;
  cuts_ += massCut_;
  cuts_ += commonCuts_;

}  // ======================================================================================

void
CutManager::fillCutMaps() {

  sampleKeyMap_["sig"] = "sig";
  sampleKeyMap_["T1bbbb"] = "sig";
  sampleKeyMap_["T1qqqq"] = "sig";
  sampleKeyMap_["T1tttt"] = "sig";
  sampleKeyMap_["zinv"] = "sig";
  sampleKeyMap_["topW"] = "sig";
  sampleKeyMap_["topWsle"] = "sle";
  sampleKeyMap_["topWslm"] = "slm";
  sampleKeyMap_["topsle"] = "sle";
  sampleKeyMap_["topslm"] = "slm";
  sampleKeyMap_["wjetssle"] = "sle";
  sampleKeyMap_["wjetsslm"] = "slm";
  sampleKeyMap_["qcd"] = "sig";
  sampleKeyMap_["qcdIDP"] = "sig";
  sampleKeyMap_["zinvIDP"] = "sig";
  sampleKeyMap_["wjets"] = "sig";
  sampleKeyMap_["other"] = "sig";
  sampleKeyMap_["ttzvv"] = "ttz";
  sampleKeyMap_["zmm"] = "zmm";
  sampleKeyMap_["zmm20"] = "zmm";
  sampleKeyMap_["dymm"] = "zmm";
  sampleKeyMap_["ttmm"] = "zmm";
  sampleKeyMap_["ttzmm"] = "zmm";
  sampleKeyMap_["VVmm"] = "zmm";
  sampleKeyMap_["tribosonmm"] = "zmm";
  sampleKeyMap_["zee"] = "zee";
  sampleKeyMap_["zee20"] = "zee";
  sampleKeyMap_["dyee"] = "zee";
  sampleKeyMap_["ttzee"] = "zee";
  sampleKeyMap_["ttee"] = "zee";
  sampleKeyMap_["VVee"] = "zee";
  sampleKeyMap_["tribosonee"] = "zee";
  sampleKeyMap_["zll"] = "zll";
  sampleKeyMap_["zll20"] = "zll";
  sampleKeyMap_["dyll"] = "zll";
  sampleKeyMap_["ttzll"] = "zll";
  sampleKeyMap_["ttll"] = "zll";
  sampleKeyMap_["VVll"] = "zll";
  sampleKeyMap_["photon"] = "photon";
  sampleKeyMap_["photon20"] = "photon";
  sampleKeyMap_["gjets"] = "photon";
  sampleKeyMap_["gjetsold"] = "photon";
  sampleKeyMap_["ttgjets"] = "photon";
  sampleKeyMap_["gjetsqcd"] = "photonqcd";

  if (isSkim_) {
    if (ntupleVersion_ == "V12") {
      objCutMap_["sig"] = "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // objCutMap_["zmm"] = "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoMuonTracks==0";
      objCutMap_["zmm"] = "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0";  // Troy mod+
      // objCutMap_["zee"] = "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoElectronTracks==0";
      objCutMap_["zee"] = "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
      objCutMap_["zll"] = "((@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0) || (@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
      objCutMap_["ttz"] = "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "@Muons.size()==1 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0";
      objCutMap_["sle"] = "@Muons.size()==0 && @Electrons.size()==1 && isoMuonTracks==0 && isoPionTracks==0";

    } else if (ntupleVersion_ == "V15" || ntupleVersion_ == "V16" || ntupleVersion_ == "V17") {

      objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // zmm skim cuts:  NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0
      objCutMap_["zmm"] = "1";
      // zee skim cuts:  NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0
      objCutMap_["zee"] = "1";
      objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0) || (NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0))";
      // photon skim cuts:  "Sum$(Photons_fullID)==1 && Photons_hasPixelSeed==0 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0"
      objCutMap_["photon"] = "@Photons.size()==1";  // Andrew recommendation, no skim cuts
      if (isMC_) objCutMap_["photon"] += " && !Photons_nonPrompt && madMinPhotonDeltaR>=0.4";  // Andrew recommendation, no skim cuts
      // objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      objCutMap_["photonqcd"] = "@Photons.size()==1";
      if (isMC_) objCutMap_["photonqcd"] += " && Photons_nonPrompt";
      objCutMap_["ttz"] = "NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "NMuons==1 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0";
      objCutMap_["sle"] = "NMuons==0 && NElectrons==1 && isoMuonTracks==0 && isoPionTracks==0";
    }

    minDphiCutMap_["nominal"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
    minDphiCutMap_["hdp"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
    minDphiCutMap_["ldp"] = "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3 || DeltaPhi4<0.3)";
    minDphiCutMap_["ldpnominal"] = "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3 || DeltaPhi4<0.3)";

    MHTCutMap_["nominal"] = "MHT>=300";
    MHTCutMap_["hdp"] = "MHT>=250";
    MHTCutMap_["ldp"] = "MHT>=250";
    MHTCutMap_["ldpnominal"] = "MHT>=300";

  } else {  // ntuple

    if (ntupleVersion_ == "V12") {
    } else if (ntupleVersion_ == "V15" || ntupleVersion_ == "V16" || ntupleVersion_ == "V17") {
      objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      // objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoMuonTracksclean==0";
      objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0";
      // objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoElectronTracksclean==0";
      objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0) || (NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1 && Photons_hasPixelSeed==0) && NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["ttz"] = "NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "NMuons==1 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["slm"] += " && TransverseMass(METPt,METPhi,@Muons.at(0).Pt(),@Muons.at(0).Phi()) < 100";  // add'l skim cut
      objCutMap_["sle"] = "NMuons==0 && NElectrons==1 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["sle"] += " && TransverseMass(METPt,METPhi,@Electrons.at(0).Pt(),@Electrons.at(0).Phi()) < 100";  // add'l skim cut
    }

    minDphiCutMap_["nominal"] = "DeltaPhi1clean>0.5 && DeltaPhi2clean>0.5 && DeltaPhi3clean>0.3 && DeltaPhi4clean>0.3";
    minDphiCutMap_["hdp"] = "DeltaPhi1clean>0.5 && DeltaPhi2clean>0.5 && DeltaPhi3clean>0.3 && DeltaPhi4clean>0.3";
    minDphiCutMap_["ldp"] = "(DeltaPhi1clean<0.5 || DeltaPhi2clean<0.5 || DeltaPhi3clean<0.3 || DeltaPhi4clean<0.3)";
    minDphiCutMap_["ldpnominal"] = "(DeltaPhi1clean<0.5 || DeltaPhi2clean<0.5 || DeltaPhi3clean<0.3 || DeltaPhi4clean<0.3)";

    MHTCutMap_["nominal"] = "MHTclean>=300 && MHTclean<=HTclean";
    MHTCutMap_["hdp"] = "MHTclean>=250 && MHTclean<=HTclean";
    MHTCutMap_["ldp"] = "MHTclean>=250 && MHTclean<=HTclean";
    MHTCutMap_["ldpnominal"] = "MHTclean>=300 && MHTclean<=HTclean";
  }

  triggerMapByName_["zmm"] = {
    "HLT_IsoMu20_v",  // Sam+
    "HLT_IsoMu22_v",  // Sam+
    "HLT_IsoMu24_v",  // prescaled in late 2017 --Owen
    "HLT_IsoMu27_v",
    "HLT_IsoMu22_eta2p1_v",  // Sam+
    "HLT_IsoMu24_eta2p1_v",  // Sam+
    "HLT_IsoTkMu22_v",  // Sam+
    "HLT_IsoTkMu24_v",
    "HLT_Mu15_IsoVVVL_PFHT350_v",
    "HLT_Mu15_IsoVVVL_PFHT400_v",
    "HLT_Mu15_IsoVVVL_PFHT450_v",  // Sam+
    "HLT_Mu15_IsoVVVL_PFHT600_v",  // Sam+
    "HLT_Mu50_IsoVVVL_PFHT400_v",  // Sam+
    "HLT_Mu50_IsoVVVL_PFHT450_v",  // Sam+
    "HLT_Mu50_v",
    "HLT_Mu55_v"  // Sam+
    // From AN-18-271
    // HLT_IsoMuX_v* (X=20,22,24,27);
    // HLT_IsoMuX_eta2p1_v* (X=22,24) ;
    // HLT_IsoTkMuX_v* (X=22,24) ;
    // HLT_Mu15_IsoVVVL_PFHTX_v* (X=350,400,450,600) ;
    // HLT_Mu50_IsoVVVL_PFHTX_v* (400,450) ;
    // HLT_MuX_v* (X=50,55) .
  };

  triggerMapByName_["zee"] = {
    "HLT_Ele105_CaloIdVT_GsfTrkIdT_v",
    "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    "HLT_Ele135_CaloIdVT_GsfTrkIdT_v",  // Sam+
    "HLT_Ele145_CaloIdVT_GsfTrkIdT_v",  // Sam+
    "HLT_Ele15_IsoVVVL_PFHT350_v",
    "HLT_Ele15_IsoVVVL_PFHT400_v",
    "HLT_Ele15_IsoVVVL_PFHT450_v",  // Sam+
    "HLT_Ele15_IsoVVVL_PFHT600_v",  // Sam+
    "HLT_Ele27_WPTight_Gsf_v",
    "HLT_Ele35_WPTight_Gsf_v",  // Sam+
    "HLT_Ele20_WPLoose_Gsf_v",  // Sam+
    "HLT_Ele45_WPLoose_Gsf_v",  // Sam+
    "HLT_Ele25_eta2p1_WPTight_Gsf_v",  // Sam+
    "HLT_Ele20_eta2p1_WPLoose_Gsf_v",  // Sam+
    "HLT_Ele27_eta2p1_WPLoose_Gsf_v"
    // From AN-18-271
    // HLT_EleX_CaloIdVT_GsfTrkIdT_v*(X=105, 115, 135, 145);
    // HLT_Ele25_eta2p1_WPTight_Gsf_v*;
    // HLT_EleX_eta2p1_WPLoose_Gsf_v*(X=20, 27);
    // HLT_Ele15_IsoVVVL_PFHTX_v*(X=350, 400, 450, 600);
    // HLT_EleX_WPTight_Gsf_v*(X=27, 35); and
    // HLT_EleX_WPLoose_Gsf_v*(X=20,45).
  };
  triggerMapByName_["photon"] = {
    "HLT_Photon175_v",
    "HLT_Photon200_v"
  };
  triggerMapByName_["sig"] = {
    "HLT_PFMET100_PFMHT100_IDTight_v",
    "HLT_PFMET110_PFMHT110_IDTight_v",
    "HLT_PFMET120_PFMHT120_IDTight_v",
    "HLT_PFMET130_PFMHT130_IDTight_v",  // Sam+
    "HLT_PFMET140_PFMHT140_IDTight_v",  // Sam+
    "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v",
    "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v",  // Sam+
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v",  // Sam+
    "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET110_PFMHT110_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET130_PFMHT130_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET140_PFMHT140_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60_v"  // Sam+
    // From AN-18-271
    // HLT_PFMETX_PFMHTX_IDTight_v* (X=100,110,120,130,140).
    // HLT_PFMETNoMuX_PFMHTNoMuX_IDTight_v* (X=100,110,120,130,140)
    // HLT_PFMETX_PFMHTX_IDTight_PFHT60_v* (X=100,110,120,130,140)
    // HLT_PFMETNoMuX_PFMHTNoMuX_IDTight_PFHT60_v* (X=100,110,120,130)
  };
  triggerMapByName_["zll"].reserve(triggerMapByName_["zmm"].size() + triggerMapByName_["zee"].size());
  triggerMapByName_["zll"] = triggerMapByName_["zmm"];
  triggerMapByName_["zll"].insert(triggerMapByName_["zll"].end(), triggerMapByName_["zee"].begin(), triggerMapByName_["zee"].end());
  triggerMapByName_["sle"] = triggerMapByName_["sig"];
  triggerMapByName_["slm"] = triggerMapByName_["sig"];

}  // ======================================================================================

void
CutManager::setTriggerIndexList(const char* sample, vector<unsigned>* triggerIndexList,
				vector<string>* TriggerNames, vector<int>* TriggerPrescales) {

  triggerIndexList->clear();
  std::vector<TString> triggers;
  vstring_map trigMap = triggerMapByName();
  if (trigMap.count(sample) > 0) {
    triggers = trigMap.at(sample);
  } else {
    cout << "No matches in triggerMapByName for sample " << sample << endl;
    return;
  }
  for (auto myTrigName : triggers) {
    for (unsigned int ti = 0; ti < TriggerNames->size(); ++ti) {
      if (TString(TriggerNames->at(ti)).Contains(myTrigName)) {
	triggerIndexList->push_back(ti);
	Int_t prescale = TriggerPrescales->at(ti);
	if (verbosity_ >= 2 || (verbosity_ >= 1 && prescale != 1))
	  cout << "Trigger " << TriggerNames->at(ti) << " (" << ti << ") prescaled by " << prescale << endl;
	continue;
      }
    }
  }

}  // ======================================================================================

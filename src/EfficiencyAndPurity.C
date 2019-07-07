//
//  Weights to correct MC for efficiencies, data for purity
//

#include "EfficiencyAndPurity.h"
#include <TF1.h>
#include <TRegexp.h>
// #include <algorithm>  // for std::find()

// ClassImp(EfficiencyAndPurity)

EfficiencyAndPurity::EfficiencyAndPurity(string deltaPhi) : deltaPhi_(deltaPhi) {
  openFiles();
}

void
EfficiencyAndPurity::openFiles() {
  TString plotDir("../plots/histograms/");

  prefiringWeightFile_ = new TFile((plotDir+"L1PrefiringMaps_new.root").Data(), "read");

  purityTrigEffFile_.push_back(new TFile((plotDir+"effHists.root").Data(), "read"));
  purityTrigEffFile_.push_back(new TFile((plotDir+"effHists.root").Data(), "read"));
  purityTrigEffFile_.push_back(new TFile((plotDir+"effHists.root").Data(), "read"));

  photonTrigEffFile_.push_back(new TFile((plotDir+"triggersRa2bRun2_v4_withTEffs.root").Data(), "read"));
  photonTrigEffFile_.push_back(new TFile((plotDir+"triggersRa2bRun2_v4_withTEffs.root").Data(), "read"));
  photonTrigEffFile_.push_back(new TFile((plotDir+"triggersRa2bRun2_v4_withTEffs.root").Data(), "read"));

  photonSFFile_.push_back(new TFile((plotDir+"Photons2016_SF_all.root").Data(), "read"));
  photonSFFile_.push_back(new TFile((plotDir+"Photons2017_SF_all.root").Data(), "read"));
  photonSFFile_.push_back(new TFile((plotDir+"Photons2018_SF_all.root").Data(), "read"));

  // elecSFFile_.push_back(new TFile((plotDir+"ElectronScaleFactors_Run2017.root").Data(), "read"));
  // elecSFFile_.push_back(new TFile((plotDir+"ElectronScaleFactors_Run2017.root").Data(), "read"));
  // elecSFFile_.push_back(new TFile((plotDir+"ElectronScaleFactors_Run2017.root").Data(), "read"));

  elecIDandIsoSFFile_.push_back(new TFile((plotDir+"Electrons2016_SF_ISOandID.root").Data(), "read"));
  elecIDandIsoSFFile_.push_back(new TFile((plotDir+"Electrons2017_SF_ISOandID.root").Data(), "read"));
  elecIDandIsoSFFile_.push_back(new TFile((plotDir+"Electrons2018_SF_ISOandID.root").Data(), "read"));

  elecRecoLowSFFile_.push_back(new TFile((plotDir+"Electrons2016_SF_RECOlow.root").Data(), "read"));
  elecRecoLowSFFile_.push_back(new TFile((plotDir+"Electrons2017_SF_RECOlow.root").Data(), "read"));
  //elecRecoLowSFFile_.push_back(new TFile((plotDir+"Electrons2018_SF_RECOlow.root").Data(), "read"));//not available yet for 2018

  elecRecoHighSFFile_.push_back(new TFile((plotDir+"Electrons2016_SF_RECOhigh.root").Data(), "read"));//pt>20
  elecRecoHighSFFile_.push_back(new TFile((plotDir+"Electrons2017_SF_RECOhigh.root").Data(), "read"));//pt>20
  elecRecoHighSFFile_.push_back(new TFile((plotDir+"Electrons2018_SF_RECOhigh.root").Data(), "read"));//pt>10

  muonIDSFFile_.push_back(new TFile((plotDir+"Muons2016_SF_ID.root").Data(), "read"));
  muonIDSFFile_.push_back(new TFile((plotDir+"Muons2017_SF_ID.root").Data(), "read"));
  muonIDSFFile_.push_back(new TFile((plotDir+"Muons2018_SF_ID.root").Data(), "read"));

  muonIsoSFFile_.push_back(new TFile((plotDir+"Muons2016_SF_ISO.root").Data(), "read"));
  muonIsoSFFile_.push_back(new TFile((plotDir+"Muons2017_SF_ISO.root").Data(), "read"));
  muonIsoSFFile_.push_back(new TFile((plotDir+"Muons2018_SF_ISO.root").Data(), "read"));

  DRfun_ = new TF1("DRfun", "[0] + [1]*min([3], x)");
  // DRfun_ = new TF1("DRfun", "[2] + (1/3)*[1]*(min([3], x) - ([2] - [0]) / [1])");
  DRpars_.push_back({0.8566, 0.0001950, 0.9530, 900});  // v17 fragmentation and Zll purity
  DRpars_.push_back({0.8566, 0.0001950, 0.9530, 900});
  DRpars_.push_back({0.8566, 0.0001950, 0.9530, 900});
  // Graph_from_hHT_DR_zmm  -- v17 fragmentation and Zll purity
  // Function parameter 0:  0.856646880781 +/- 0.0197955271053  // 0.857 +/- 0.020
  // Function parameter 1:  0.000194995307077 +/- 3.78488782632e-05  // (0.195 +/- 0.038) e-3
  // Average y0 = 0.9530 +/- 0.0069;  x0 = (y0 -p0) / p1  // 0.953 +/- 0.007
  // DRpars_.push_back({0.8547, 0.0001983, 0.9526, 900});  // v17
  // DRpars_.push_back({0.8547, 0.0001983, 0.9526, 900});
  // DRpars_.push_back({0.8547, 0.0001983, 0.9526, 900});
  // Graph_from_hHT_DR_zmm  -- with v17 for 2016, 17 also
  // Function parameter 0:  0.854716822646 +/- 0.0198007015338
  // Function parameter 1:  0.000198285579166 +/- 3.78772029776e-05
  // Average y0 = 0.9526 +/- 0.0069;  x0 = (y0 -p0) / p1
  // Graph_from_hHT_DR_zmm  for ARC freeze
  // Function parameter 0:  0.847711708842 +/- 0.0195782643702
  // Function parameter 1:  0.00019064466223 +/- 3.74122384363e-05
  // Average y0 = 0.9419;  x0 = (y0 -p0) / p1
  // DRpars_.push_back({0.8386, 0.0001812, 0.9290, 900});
  // Graph_from_hHT_DR_zmm  Run 2 19 Mar 2019
  // Function parameter 0:  0.838552899836 +/- 0.0192480033502
  // Function parameter 1:  0.000181188022872 +/- 3.63880256987e-05
  // Average y0 = 0.9290;  x0 = (y0 -p0) / p1
  // DRpars_.push_back({0.8229, 0.0001665, 0.9061, 900});
  // Graph_from_hHT_DR_zmm  Run 2 21 Feb 2019
  // Function parameter 0:  0.822859680122 +/- 0.0188240089302
  // Function parameter 1:  0.000166574459092 +/- 3.55474346141e-05
  // DRpars_.push_back({0.8378, 0.0001363, 0.9054, 900});
  // Graph_from_hHT_DR_zmm  2016 noPU
  // Function parameter 0:  0.837859796025 +/- 0.0352792948296
  // Function parameter 1:  0.000136284149882 +/- 6.6930863488e-05
  // Average y0 = 1.0871.  x0 = (y0 -p0) / p1

  // effWt /= (min(HT, 900.0) - 497.4)*(0.0002288) + 1.0395;  // Run2 Z Pt weighted
  // Graph_from_hHT_DR_zmm
  // Function parameter 0:  0.92567079283 +/- 0.021512658121
  // Function parameter 1:  0.000228788345936 +/- 4.08070918994e-05
  // Average y0 = 1.0395.  x0 = (y0 -p0) / p1

}  // ======================================================================================

void
EfficiencyAndPurity::getHistos(const char* sample, int currentYear) {

  std::vector<const char*> year;  year.push_back("2016");  year.push_back("2017");  year.push_back("2018");

  // For L1 prefiring weight
  if (currentYear == Year2016) {
    hPrefiring_photon_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_photonptvseta_2016BtoH");
    if (hPrefiring_photon_ == nullptr) cout << "***** Histogram L1prefiring_photonptvseta_2016BtoH not found *****" << endl;
    hPrefiring_jet_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_jetptvseta_2016BtoH");
    if (hPrefiring_jet_ == nullptr) cout << "***** Histogram L1prefiring_jetptvseta_2016BtoH not found *****" << endl;
  } else if (currentYear == Year2017) {
    hPrefiring_photon_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_photonptvseta_2017BtoF");
    if (hPrefiring_photon_ == nullptr) cout << "***** Histogram L1prefiring_photonptvseta_2017BtoF not found *****" << endl;
    hPrefiring_jet_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_jetptvseta_2017BtoF");
    if (hPrefiring_jet_ == nullptr) cout << "***** Histogram L1prefiring_jetptvseta_2017BtoF not found *****" << endl;
  }

  // For purity, Fdir, trigger eff, reco eff
  theSample_ = TString(sample);
  // hSFeff_ = nullptr;
  hSFeff_.clear();
  FdirGraph_ = nullptr;
  hPurity_.clear();
  hTrigEff_.clear();
  if (theSample_.Contains("zmm") && !theSample_.Contains("tt")) {
    TString hn("h_pur_m");
    if (deltaPhi_.find("ldp") != string::npos) hn += "_ldp";
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get(hn));
    if (hPurity_.back() == nullptr) cout << "***** Histogram " << hn << " not found *****" << endl;
  } else if (theSample_.Contains("zee") && !theSample_.Contains("tt")) {
    TString hn("h_pur_e");
    if (deltaPhi_.find("ldp") != string::npos) hn += "_ldp";
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get(hn));
    if (hPurity_.back() == nullptr) cout << "***** Histogram " << hn << " not found *****" << endl;
  } else if (theSample_.Contains("photon")) {
    TString hn("h_pur_eb_");
    if (deltaPhi_.find("ldp") != string::npos) hn += "ldp_";
    hn += year.at(currentYear);
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get(hn));
    if (hPurity_.back() == nullptr) cout << "***** Histogram " << hn << " not found *****" << endl;
    hn("eb") = "ec";
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get(hn));
    if (hPurity_.back() == nullptr) cout << "***** Histogram " << hn << " not found *****" << endl;
    TString gFdirName;
    if (deltaPhi_.find("ldpnominal") != string::npos) gFdirName = "h_bin46_NJets8910_ldpnominal";
    else if (deltaPhi_.find("ldp") != string::npos &&
	     deltaPhi_.find("nominal") == string::npos) gFdirName = "h_bin59_NJets8910_ldp";
    else gFdirName = "h_bin46_NJets8910";
    FdirGraph_ = (TGraphErrors*) purityTrigEffFile_.at(currentYear)->Get(gFdirName);
    if (FdirGraph_ == nullptr) cout << "***** Histogram " << gFdirName << " not found *****" << endl;
  } else if (theSample_.Contains("dymm") || theSample_.Contains("ttmm") ||
	     theSample_.Contains("ttzmm") || theSample_.Contains("VVmm")) {
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_m"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_m not found *****" << endl;
    // hSFeff_ = (TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_SFm_MHT");
    // if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;

    //These change based on year
    if (currentYear == Year2016) {
      hSFeff_.push_back((TH2F*) muonIDSFFile_.at(currentYear)->Get("NUM_MediumID_DEN_genTracks_eta_pt"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for 2016 muon ID SFs not found *****"
					  << " file " << muonIDSFFile_.at(currentYear)->GetName() << endl;
      hSFeff_.push_back((TH2F*) muonIsoSFFile_.at(currentYear)->Get("NUM_TightRelIso_DEN_MediumID_eta_pt"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for 2016 muon iso SFs not found *****" << endl;
    } else if (currentYear == Year2017) {
      hSFeff_.push_back((TH2F*) muonIDSFFile_.at(currentYear)->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for 2017 muon ID SFs not found *****"
					  << " file " << muonIDSFFile_.at(currentYear)->GetName() << endl;
      hSFeff_.push_back((TH2F*) muonIsoSFFile_.at(currentYear)->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for 2017 muon iso SFs not found *****" << endl;
    }
    else {
      hSFeff_.push_back((TH2F*) muonIDSFFile_.at(currentYear)->Get("NUM_MediumID_DEN_TrackerMuons_pt_abseta"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for 2018 muon ID SFs not found *****"
					  << " file " << muonIDSFFile_.at(currentYear)->GetName() << endl;
      hSFeff_.push_back((TH2F*) muonIsoSFFile_.at(currentYear)->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for 2018 muon iso SFs not found *****" << endl;
    }

  } else if (theSample_.Contains("dyee") || theSample_.Contains("ttee") ||
	     theSample_.Contains("ttzee") || theSample_.Contains("VVee")) {
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_e"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_e not found *****" << endl;
    // hSFeff_ = (TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_SFe_MHT");  // Maybe this should be h_NJets
    // if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;

    //These change based on year
    if (currentYear == Year2016) {
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2016_CutBasedVetoNoIso94XV2"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2016_Mini"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron iso SFs not found *****" << endl;
    }
    else if (currentYear == Year2017) {
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2017_CutBasedVetoNoIso94XV2"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2017_MVAVLooseTightIP2DMini"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron iso SFs not found *****" << endl;
    }
    else if (currentYear == Year2018) {
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2018_CutBasedVetoNoIso94XV2"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2018_Mini"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron iso SFs not found *****" << endl;
    }
    hSFeff_.push_back((TH2F*) elecRecoHighSFFile_.at(currentYear)->Get("EGamma_SF2D"));
    if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron reco (high pT) SFs not found *****" << endl;

    //Electron reco (pt<10 GeV) SFs for 2018 not (yet) available
    if (currentYear != Year2018) {
      hSFeff_.push_back((TH2F*) elecRecoLowSFFile_.at(currentYear)->Get("EGamma_SF2D"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron reco (low pT) SFs not found *****" << endl;
    }


  } else if (theSample_.Contains("gjets")) {
    // // Troy-style trigger efficiency histograms:
    // hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_eb"));
    // if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_eb not found *****" << endl;
    // hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_ec"));
    // if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_ec not found *****" << endl;
    // // Sam-style trigger efficiency functions:
    // fTrigEff_.push_back((TF1*) purityTrigEffFile_.at(currentYear)->Get("f_trig_eb"));
    // if (fTrigEff_.back() == nullptr) {
    //   cout << "***** Histogram f_trig_eb not found *****" << endl;
    // } else {
    //   cout << "Asymptotic EB trigger efficiency = " << 1 + fTrigEff_[0]->GetParameter(0)
    // 	   << "+/-" << fTrigEff_[0]->GetParError(0) << endl;
    // }
    // fTrigEff_.push_back((TF1*) purityTrigEffFile_.at(currentYear)->Get("f_trig_ec"));
    // if (fTrigEff_.back() == nullptr) {
    //   cout << "***** Histogram f_trig_ec not found *****" << endl;
    // } else {
    //   cout << "Asymptotic EC trigger efficiency = " << 1 + fTrigEff_[1]->GetParameter(0)
    // 	   << "+/-" << fTrigEff_[1]->GetParError(0) << endl;
    // }
    // Sam-style trigger efficiency TEfficiency's
    TString te("teff_SinglePhotonBarrelLooseHdp_hists_Run");  te += year.at(currentYear);  te += "_JetHT";
    if (deltaPhi_.find("ldp") != string::npos) te("Hdp") = "Ldp";
    eTrigEff_.push_back((TEfficiency*) photonTrigEffFile_.at(currentYear)->Get(te.Data()));
    if (eTrigEff_.back() == nullptr)  cout << "***** Histogram " << te << " not found *****" << endl;
    te("Barrel") = "Endcap";
    eTrigEff_.push_back((TEfficiency*) photonTrigEffFile_.at(currentYear)->Get(te.Data()));
    if (eTrigEff_.back() == nullptr)  cout << "***** Histogram " << te << " not found *****" << endl;
    // hSFeff_ = (TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_SFg_MHT");  // Maybe this should be h_NJets
    // if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;

    hSFeff_.push_back((TH2F*) photonSFFile_.at(currentYear)->Get("EGamma_SF2D"));
    if (hSFeff_.back() == nullptr) cout << "***** Histogram for photon SFs not found *****" << endl;

    DRfun_->SetParameters(DRpars_.at(currentYear).data());
    cout << "DR intercept = " << DRpars_.at(currentYear)[0] << ", slope = " << DRpars_.at(currentYear)[1]
	 << ", Ave. value = " << DRpars_.at(currentYear)[2] << ", saturation point = " << DRpars_.at(currentYear)[3] << endl;
  }

}  // ======================================================================================

pair<double, double>
EfficiencyAndPurity::weight(CCbinning* CCbins,
					      Int_t NJets, Int_t BTags,
					      Double_t MHT, Double_t HT,
					      vector<TLorentzVector> ZCandidates,
					      vector<TLorentzVector> Photons,
					      vector<TLorentzVector> Electrons,
					      vector<TLorentzVector> Muons,
					      vector<double> EBphoton,
					      bool applyDRfitWt,
					      int currentYear) {
  // For double ratio, apply weights for purity, Fdir, trigger eff, reco eff.
  // Systematic error is relative.
  double effWt = 1, effSys = 0;

  if ((theSample_.Contains("zmm") || theSample_.Contains("zee")) && !theSample_.Contains("tt")) {
    int binCCjb = CCbins->jb(CCbins->jbin(NJets), CCbins->bbin(NJets, BTags));
    if (hPurity_[0] != nullptr) {
      double eff = hPurity_[0]->GetBinContent(binCCjb);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[0]->GetBinError(binCCjb) / eff : 0;
    }

  } else if (theSample_.Contains("photon")) {
    if(EBphoton.at(0) == 1 && hPurity_[0] != nullptr) {
      int bin = hPurity_[0]->GetNbinsX();  while (MHT < hPurity_[0]->GetBinLowEdge(bin)) bin--;
      double eff = hPurity_[0]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[0]->GetBinError(bin) / eff : 0;
    }
    if(EBphoton.at(0) == 0 && hPurity_[1] != nullptr) {
      int bin = hPurity_[1]->GetNbinsX();  while (MHT < hPurity_[1]->GetBinLowEdge(bin)) bin--;
      double eff = hPurity_[1]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[1]->GetBinError(bin) / eff : 0;
    }
    int CCbinjk = CCbins->jk(CCbins->jbin(NJets), CCbins->kinBin(HT, MHT));
    // if (CCbinjk > 0 && FdirHist_ != nullptr && CCbinjk <= FdirHist_->GetNbinsX()) effWt *= FdirHist_->GetBinContent(CCbinjk);
    if (CCbinjk > 0 && FdirGraph_ != nullptr && CCbinjk < FdirGraph_->GetN()) {
      double eff = FdirGraph_->GetY()[CCbinjk-1];
      effWt *= eff;
      effSys = eff > 0 ? quadSum(effSys, FdirGraph_->GetErrorY(CCbinjk-1)) / eff : 0;
    }
    // FIXME:  For LDP, fill extra bins with value 0.825

  } else if (theSample_.Contains("dymm") || theSample_.Contains("dyee") ||
	     theSample_.Contains("ttmm") || theSample_.Contains("ttee") ||
	     theSample_.Contains("ttzmm") || theSample_.Contains("ttzee") ||
	     theSample_.Contains("VVmm") || theSample_.Contains("VVee")) {
    if (hTrigEff_[0] != nullptr) {
      Double_t zpt = ZCandidates.at(0).Pt();  zpt = max(hTrigEff_[0]->GetBinLowEdge(1), zpt);
      int bin = hTrigEff_[0]->GetNbinsX();  while (zpt < hTrigEff_[0]->GetBinLowEdge(bin)) bin--;
      double eff = hTrigEff_[0]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hTrigEff_[0]->GetBinError(bin) / eff : 0;
    }
    if (hSFeff_[0] != nullptr && hSFeff_[1] != nullptr) {
      float pt = 0; float eta = 0;
      int globalbin_ID = 0; int globalbin_ISO = 0;
      
      double sysCorr = 0;  // Errors are correlated between Z daughter leptons
      if (theSample_.Contains("ee") && hSFeff_[2]!=nullptr) {
	int globalbin_RECO = 0;
	int numElectrons = Electrons.size();
	for (int i=0; i<numElectrons; i++){
	  pt  = Electrons.at(i).Pt(); if (pt>500) pt=499.9;
	  eta = Electrons.at(i).Eta();
	  globalbin_ID  = hSFeff_[0]->FindBin(eta,pt); globalbin_ISO = hSFeff_[1]->FindBin(eta,pt);
	  double effID = hSFeff_[0]->GetBinContent(globalbin_ID);
	  double effISO = hSFeff_[1]->GetBinContent(globalbin_ISO);
	  effWt *= effID * effISO;
	  double sysOne = quadSum(effID > 0 ? hSFeff_[0]->GetBinError(globalbin_ID)/effID : 0,
				  effISO > 0 ? hSFeff_[1]->GetBinError(globalbin_ISO)/effISO : 0);
	  if (currentYear == Year2018) { //2018 only has reco SFs for pt>10
	    if (pt<=10.0) pt = 10.1;
	    globalbin_RECO = hSFeff_[2]->FindBin(eta,pt);
	    double eff = hSFeff_[2]->GetBinContent(globalbin_RECO);
	    effWt *= eff;
	    sysOne = eff > 0 ? quadSum(sysOne, hSFeff_[2]->GetBinError(globalbin_RECO)/eff) : sysOne;
	  }
	  else { //2016 and 2017 have reco SFs for 0-20 GeV and >20 GeV
	    if (pt>20) {
	      globalbin_RECO = hSFeff_[2]->FindBin(eta,pt);
	      double eff = hSFeff_[2]->GetBinContent(globalbin_RECO);
	      effWt *= eff;
	      sysOne = eff > 0 ? quadSum(sysOne, hSFeff_[2]->GetBinError(globalbin_RECO)/eff) : sysOne;
	    }
	    else {
	      globalbin_RECO = hSFeff_[3]->FindBin(eta,pt);
	      double eff = hSFeff_[3]->GetBinContent(globalbin_RECO);
	      effWt *= eff;
	      sysOne = eff > 0 ? quadSum(sysOne, hSFeff_[3]->GetBinError(globalbin_RECO)/eff) : 0;
	    }
	  } //2016 or 2017
	  sysCorr += sysOne;
	} //loop over all electrons (should be two)
      }

      else if (theSample_.Contains("mm")) {
	double muIDISOsys = 0.003;  //  Hard-wired based on _syst histo for 2017 from POG
	int numMuons = Muons.size();
	for (int i=0; i<numMuons; i++){
	  pt = Muons.at(i).Pt(); if (pt>120) pt=119.9;
	  if (currentYear == Year2016) {
	    eta = Muons.at(i).Eta();
	    globalbin_ID = hSFeff_[0]->FindBin(eta,pt); globalbin_ISO = hSFeff_[1]->FindBin(eta,pt);
	  }
	  else { //2017 and 2018 have abs(eta), and the x- and y-axis are swapped compared to 2016
	    eta = abs(Muons.at(i).Eta());
	    globalbin_ID = hSFeff_[0]->FindBin(pt,eta); globalbin_ISO = hSFeff_[1]->FindBin(pt,eta);
	  }
	  double effSF = hSFeff_[0]->GetBinContent(globalbin_ID);
	  double effISO = hSFeff_[1]->GetBinContent(globalbin_ISO);
	  effWt *= effSF * effISO;
	  sysCorr += quadSum(effSF > 0 ? hSFeff_[0]->GetBinError(globalbin_ID)/effSF : 0,
			     effISO > 0 ? hSFeff_[1]->GetBinError(globalbin_ISO)/effISO : 0,
			     muIDISOsys);
	}
      }
      effSys = quadSum(effSys, sysCorr);

    }

  } else if (theSample_.Contains("gjets")) {
    // if(EBphoton.at(0) == 1 && hTrigEff_[0] != nullptr) {
    //   int bin = hTrigEff_[0]->GetNbinsX();  while (MHT < hTrigEff_[0]->GetBinLowEdge(bin)) bin--;
    //   effWt *= hTrigEff_[0]->GetBinContent(bin);
    // }
    // if(EBphoton.at(0) == 0 && hTrigEff_[1] != nullptr) {
    //   int bin = hTrigEff_[1]->GetNbinsX();  while (MHT < hTrigEff_[1]->GetBinLowEdge(bin)) bin--;
    //   effWt *= hTrigEff_[1]->GetBinContent(bin);
    // }
    // if(EBphoton.at(0) == 1 && fTrigEff_[0] != nullptr) {
    //   effWt *= fTrigEff_[0]->Eval(max(double(205), Photons.at(0).Pt()));  // hard-wired cutoff
    // }
    // if(EBphoton.at(0) == 0 && fTrigEff_[1] != nullptr) {
    //   effWt *= fTrigEff_[1]->Eval(max(double(205), Photons.at(0).Pt()));  // hard-wired cutoff
    // }
    if(EBphoton.at(0) == 1 && eTrigEff_[0] != nullptr) {
      TH1F* htot = (TH1F*) eTrigEff_[0]->GetTotalHistogram();
      Int_t bin = min(htot->GetNbinsX(), htot->FindBin(Photons.at(0).Pt()));
      double eff = eTrigEff_[0]->GetEfficiency(bin);
      effWt *= eff;
      effSys = eff > 0 ? eTrigEff_[0]->GetEfficiencyErrorLow(bin) / eff : 0;
    } else if(EBphoton.at(0) == 0 && eTrigEff_[1] != nullptr) {
      TH1F* htot = (TH1F*) eTrigEff_[1]->GetTotalHistogram();
      Int_t bin = min(htot->GetNbinsX(), htot->FindBin(Photons.at(0).Pt()));
      double eff = eTrigEff_[1]->GetEfficiency(bin);
      effWt *= eff;
      effSys = eff > 0 ? eTrigEff_[1]->GetEfficiencyErrorLow(bin) / eff : 0;
    }

    if (applyDRfitWt) effWt /= DRfun_->Eval(HT);

    if (hSFeff_[0] != nullptr) {
      float photon_pt = 0; float photon_eta = 0;
      int globalbin_photon = 0;
      int numPhotons = Photons.size();
      for (int i=0; i<numPhotons; i++){
	photon_pt = Photons.at(i).Pt(); photon_eta = Photons.at(i).Eta();
	if (photon_pt>500) photon_pt=499.9;
	globalbin_photon  = hSFeff_[0]->FindBin(photon_eta, photon_pt);
	double eff = hSFeff_[0]->GetBinContent(globalbin_photon);
	effWt *= eff;
	effSys = eff > 0 ? quadSum(effSys, hSFeff_[0]->GetBinError(globalbin_photon)/eff) : effSys;
      }
    }
  }
  return pair<double, double>(effWt, effSys);

}  // ======================================================================================

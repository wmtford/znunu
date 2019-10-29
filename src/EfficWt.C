//
//  Weights to correct MC for efficiencies, data for purity
//

#include "EfficWt.h"
#include <TRegexp.h>

// ClassImp(EfficWt)

EfficWt::EfficWt(string deltaPhi) : deltaPhi_(deltaPhi) {
  openFiles();
}

void
EfficWt::openFiles() {
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
  // Graph_from_hHT_DR_zmm  -- with v17 for 2016, 17 also
  // Function parameter 0:  0.854716822646 +/- 0.0198007015338
  // Function parameter 1:  0.000198285579166 +/- 3.78772029776e-05
  // Average y0 = 0.9526 +/- 0.0069;  x0 = (y0 -p0) / p1
  // Graph_from_hHT_DR_zmm  for ARC freeze
  // Function parameter 0:  0.847711708842 +/- 0.0195782643702
  // Function parameter 1:  0.00019064466223 +/- 3.74122384363e-05
  // Average y0 = 0.9419;  x0 = (y0 -p0) / p1

}  // ======================================================================================

void
EfficWt::getHistos(const char* sample, int currentYear) {

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

    hSFeff_.push_back((TH2F*) photonSFFile_.at(currentYear)->Get("EGamma_SF2D"));
    if (hSFeff_.back() == nullptr) cout << "***** Histogram for photon SFs not found *****" << endl;

    DRfun_->SetParameters(DRpars_.at(currentYear).data());
    cout << "DR intercept = " << DRpars_.at(currentYear)[0] << ", slope = " << DRpars_.at(currentYear)[1]
	 << ", Ave. value = " << DRpars_.at(currentYear)[2] << ", saturation point = " << DRpars_.at(currentYear)[3] << endl;
  }

}  // ======================================================================================

pair<double, double>
EfficWt::weight(Ntuple* Tupl, CCbinning* CCbins, bool applyDRfitWt, int currentYear) {
  // For double ratio, apply weights for purity, Fdir, trigger eff, reco eff.
  // Systematic error is relative.
  double effWt = 1, effSys = 0;

  if ((theSample_.Contains("zmm") || theSample_.Contains("zee")) && !theSample_.Contains("tt")) {
    int binCCjb = CCbins->jb(CCbins->jbin(Tupl->NJets), CCbins->bbin(Tupl->NJets, Tupl->BTags));
    if (hPurity_[0] != nullptr) {
      double eff = hPurity_[0]->GetBinContent(binCCjb);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[0]->GetBinError(binCCjb) / eff : 0;
    }

  } else if (theSample_.Contains("photon") && Tupl->Photons->size() == 1) {
    if(Tupl->Photons_isEB->at(0) == 1 && hPurity_[0] != nullptr) {
      int bin = hPurity_[0]->GetNbinsX();  while (Tupl->MHT < hPurity_[0]->GetBinLowEdge(bin)) bin--;
      double eff = hPurity_[0]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[0]->GetBinError(bin) / eff : 0;
    }
    if(Tupl->Photons_isEB->at(0) == 0 && hPurity_[1] != nullptr) {
      int bin = hPurity_[1]->GetNbinsX();  while (Tupl->MHT < hPurity_[1]->GetBinLowEdge(bin)) bin--;
      double eff = hPurity_[1]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[1]->GetBinError(bin) / eff : 0;
    }
    int CCbinjk = CCbins->jk(CCbins->jbin(Tupl->NJets), CCbins->kinBin(Tupl->HT, Tupl->MHT));
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
    if (hTrigEff_[0] != nullptr && Tupl->ZCandidates->size() == 1) {
      Double_t zpt = Tupl->ZCandidates->at(0).Pt();  zpt = max(hTrigEff_[0]->GetBinLowEdge(1), zpt);
      int bin = hTrigEff_[0]->GetNbinsX();  while (zpt < hTrigEff_[0]->GetBinLowEdge(bin)) bin--;
      double eff = hTrigEff_[0]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hTrigEff_[0]->GetBinError(bin) / eff : 0;
    }
    if (hSFeff_[0] != nullptr && hSFeff_[1] != nullptr) {
      float pt = 0; float eta = 0;
      int globalbin_ID = 0; int globalbin_ISO = 0;
      
      double sysCorr = 0;  // Errors are correlated between Z daughter leptons
      if (theSample_.Contains("ee") && hSFeff_[2]!=nullptr && Tupl->NElectrons == 2) {
	int globalbin_RECO = 0;
	int numElectrons = Tupl->Electrons->size();
	for (int i=0; i<numElectrons; i++){
	  // The following line was absent for SUS-19-006
	  if (!(Tupl->Electrons_mediumID->at(i) && Tupl->Electrons_passIso->at(i))) continue;

	  pt  = Tupl->Electrons->at(i).Pt(); if (pt>500) pt=499.9;
	  eta = Tupl->Electrons->at(i).Eta();
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
	int numMuons = Tupl->Muons->size();
	for (int i=0; i<numMuons; i++){
	  // The following line was absent for SUS-19-006
	  if (!(Tupl->Muons_mediumID->at(i) && Tupl->Muons_passIso->at(i))) continue;

	  pt = Tupl->Muons->at(i).Pt(); if (pt>120) pt=119.9;
	  if (currentYear == Year2016) {
	    eta = Tupl->Muons->at(i).Eta();
	    globalbin_ID = hSFeff_[0]->FindBin(eta,pt); globalbin_ISO = hSFeff_[1]->FindBin(eta,pt);
	  }
	  else { //2017 and 2018 have abs(eta), and the x- and y-axis are swapped compared to 2016
	    eta = abs(Tupl->Muons->at(i).Eta());
	    globalbin_ID = hSFeff_[0]->FindBin(pt,eta); globalbin_ISO = hSFeff_[1]->FindBin(pt,eta);
	  }
	  double effSF = hSFeff_[0]->GetBinContent(globalbin_ID);
	  double effISO = hSFeff_[1]->GetBinContent(globalbin_ISO);
	  effWt *= effSF * effISO;
	  sysCorr += quadSum(effSF > 0 ? hSFeff_[0]->GetBinError(globalbin_ID)/effSF : 0,
			     effISO > 0 ? hSFeff_[1]->GetBinError(globalbin_ISO)/effISO : 0,
			     muIDISOsys);
	} //loop over all muons (should be two)
      }
      effSys = quadSum(effSys, sysCorr);

    }

  } else if (theSample_.Contains("gjets") && Tupl->Photons->size() == 1) {
    if(Tupl->Photons_isEB->at(0) == 1 && eTrigEff_[0] != nullptr) {
      TH1F* htot = (TH1F*) eTrigEff_[0]->GetTotalHistogram();
      Int_t bin = min(htot->GetNbinsX(), htot->FindBin(Tupl->Photons->at(0).Pt()));
      double eff = eTrigEff_[0]->GetEfficiency(bin);
      effWt *= eff;
      effSys = eff > 0 ? eTrigEff_[0]->GetEfficiencyErrorLow(bin) / eff : 0;
    } else if(Tupl->Photons_isEB->at(0) == 0 && eTrigEff_[1] != nullptr) {
      TH1F* htot = (TH1F*) eTrigEff_[1]->GetTotalHistogram();
      Int_t bin = min(htot->GetNbinsX(), htot->FindBin(Tupl->Photons->at(0).Pt()));
      double eff = eTrigEff_[1]->GetEfficiency(bin);
      effWt *= eff;
      effSys = eff > 0 ? eTrigEff_[1]->GetEfficiencyErrorLow(bin) / eff : 0;
    }

    if (applyDRfitWt) effWt /= DRfun_->Eval(Tupl->HT);

    if (hSFeff_[0] != nullptr && Tupl->Photons->size() == 1) {
      float photon_pt = 0; float photon_eta = 0;
      int globalbin_photon = 0;
      photon_pt = Tupl->Photons->at(0).Pt(); photon_eta = Tupl->Photons->at(0).Eta();
      if (photon_pt>500) photon_pt=499.9;
      globalbin_photon  = hSFeff_[0]->FindBin(photon_eta, photon_pt);
      double eff = hSFeff_[0]->GetBinContent(globalbin_photon);
      effWt *= eff;
      effSys = eff > 0 ? quadSum(effSys, hSFeff_[0]->GetBinError(globalbin_photon)/eff) : effSys;
    }
  }
  return pair<double, double>(effWt, effSys);

}  // ======================================================================================

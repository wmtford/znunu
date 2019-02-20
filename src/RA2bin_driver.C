{
  gROOT->Reset();
  gROOT->ProcessLine(".L ../src/RA2bin_inputs_Zinv.C+");

  bool do2016 = false;
  bool do2017 = false;
  bool do2018 = false;
  bool do2018AB = true;
  bool do2018CD = false;
  bool doRunw = false;

  if (do2016) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair("../datFiles/gJets_2016v16_DR2016wt", 35.9)
    };
    RA2bin_inputs_Zinv(Signal, gJetsInputs, "../datFiles/DR_2016v16_DR2016wt_noPU", "../datFiles/DY_2016v16_noPU",
		       "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 35.9/35.9, true);
  }
  else if (do2017) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair("../datFiles/gJets_2017v16_ZptWt_DR2017wt", 41.5)
    };
    RA2bin_inputs_Zinv(Signal, gJetsInputs, "../datFiles/DR_2017v16_ZptWt_DR2017wt", "../datFiles/DY_2017v16_ZptWt",
		       "../plots/histograms/ZinvMCttzMC174bin_2017v16_ZptWt.root", 41.5/41.5, true);
  }
  else if (do2018AB) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair("../datFiles/gJets_2018v16_ZptWt_DR2017wt", 21.1)
    };
    RA2bin_inputs_Zinv(Signal, gJetsInputs, "../datFiles/DR_2018ABv16_ZptWt_DR2017wt", "../datFiles/DY_2018ABv16_ZptWt_noJ",
		       "../plots/histograms/ZinvMCttzMC174bin_2017v16_ZptWt.root", 21.1/41.5, true);
  }
  else if (do2018CD) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair("../datFiles/gJets_2018HEMv16_ZptWt_DR2017wt", 38.3)
    };
    RA2bin_inputs_Zinv(Signal, gJetsInputs, "../datFiles/DR_2018CDv16_ZptWt_DR2017wt", "../datFiles/DY_2018CDv16_ZptWt_noJ",
		       "../plots/histograms/ZinvMCttzMC174bin_2017v16_ZptWt.root", 38.3/41.5, true);
  }
  else if (do2018) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair("../datFiles/gJets_2018v16_ZptWt_DR2017wt", 21.1),
      std::make_pair("../datFiles/gJets_2018HEMv16_ZptWt_DR2017wt", 38.3)
    };
    cout << "No DR, DY file for 2018" << endl;
    // RA2bin_inputs_Zinv(Signal, gJetsInputs, "../datFiles/DR_Run2v16_ZptWt_DR2017wt_noPU", "../datFiles/DY_Run2v16_ZptWt_noJ",
    // 		       "../plots/histograms/ZinvMCttzMC174bin_Run2v16_ZptWt.root", 136.8/136.8, true);
  }
  else if (doRun2) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair("../datFiles/gJets_2016v16_DR2016wt", 35.9),
      std::make_pair("../datFiles/gJets_2017v16_ZptWt_DR2017wt", 41.5),
      std::make_pair("../datFiles/gJets_2018v16_ZptWt_DR2017wt", 21.1),
      std::make_pair("../datFiles/gJets_2018HEMv16_ZptWt_DR2017wt", 38.3)
    };
    RA2bin_inputs_Zinv(Signal, gJetsInputs, "../datFiles/DR_Run2v16_ZptWt_DR1617wt_noPU", "../datFiles/DY_Run2v16_ZptWt_noJ",
		       "../plots/histograms/ZinvMCttzMC174bin_Run2v16_ZptWt.root", 136.8/136.8, true);
  }


  // RA2bin_inputs_Zinv(LDP, gJetsInputs, "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 136.8/35.9, true);
  // RA2bin_inputs_Zinv(HDP, gJetsInputs, "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 136.8/35.9, true);



  
}

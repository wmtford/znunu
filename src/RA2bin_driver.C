{
  gROOT->Reset();
  gROOT->ProcessLine(".L RA2bin_inputs_Zinv.C+");

  std::vector< std::pair<TString, float> > gJetsInputs = {std::make_pair("../datFiles/gJets_2016", 35.9),
							  std::make_pair("../datFiles/gJets_2017", 41.5),
							  std::make_pair("../datFiles/gJets_2018", 59.4)};
  // RA2bin_inputs_Zinv(LDP, gJetsInputs, "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 136.8/35.9, true);
  // RA2bin_inputs_Zinv(HDP, gJetsInputs, "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 136.8/35.9, true);
  RA2bin_inputs_Zinv(Signal, gJetsInputs, "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 136.8/35.9, true);

  
}

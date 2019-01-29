{
  gROOT->Reset();
  gROOT->ProcessLine(".L RA2bin_inputs_Zinv.C+");

  // RA2bin_inputs_Zinv(LDP, "../datFiles/gJets_2017", "../datFiles/DR_2017", "../datFiles/DY_2017", "../plots/histograms/ZinvMCttzMC174bin.root", 1, true);
  // RA2bin_inputs_Zinv(HDP, "../datFiles/gJets_2017", "../datFiles/DR_2017", "../datFiles/DY_2017", "../plots/histograms/ZinvMCttzMC174 bin.root", 1, true);
  RA2bin_inputs_Zinv(Signal, "../datFiles/gJets_2017", "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin.root", 136.8/41.5, true);

  
}

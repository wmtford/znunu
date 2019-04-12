{
  gROOT->Reset();
  gROOT->ProcessLine(".L ../src/RA2bin_inputs_Zinv.C+");

  bool do2016 = false;
  bool do2017 = false;
  bool do2018AB = false;
  bool do2018CD = false;
  bool do2018 = false;
  bool doRun2 = true;

  // MC k-factor correction (Run 2):  (1.257*35.9+1.165*(136.8-35.9))/(1.23*136.8)

  const char* gJets16 = "../datFiles/gJets_2016v16_DRr2wt";
  const char* gJets17 = "../datFiles/gJets_2017v16_DRr2wt";
  const char* gJets18AB = "../datFiles/gJets_2018ABv17_DRr2wt";
  const char* gJets18CD = "../datFiles/gJets_2018CDv17_DRr2wt";
  const char* DR = "../datFiles/DR_Run2v161617_DRr2wt";
  const char* DY16 = "../datFiles/DY_2016v16_noPU";
  const char* DY17 = "../datFiles/DY_2017v16";
  const char* DY18AB = "../datFiles/DY_2018ABv16";
  const char* DY18CD = "../datFiles/DY_2018CDv16";
  const char* DYRun2 = "../datFiles/DY_Run2v161617";
  const char* zJets16 = "../plots/histograms/ZinvMCttzMC174bin_2016v16_noPU.root";
  const char* zJets17 = "../plots/histograms/ZinvMCttzMC174bin_2017v16_HT17_ZptWt.root";
  const char* zJets18AB = "../plots/histograms/ZinvMCttzMC174bin_2018v16_HT17_ZptWt.root";
  const char* zJets18CD = "../plots/histograms/ZinvMCttzMC174bin_2018HEMv16_HT17_ZptWt.root";
  const char* zJetsRun2 = "../plots/histograms/ZinvMCttzMC174bin_Run2v161617.root";

  if (do2016) {
    std::vector< std::pair<TString, float> > gJetsInputs = {std::make_pair(gJets16, 35.9)};
    RA2bin_inputs_Zinv(Signal, gJetsInputs, DR, DY16, zJets16, 35.9/35.9, true);
  }
  else if (do2017) {
    std::vector< std::pair<TString, float> > gJetsInputs = {std::make_pair(gJets17, 41.5)};
    RA2bin_inputs_Zinv(Signal, gJetsInputs, DR, DY17, zJets17, 41.5/41.5, true);
  }
  else if (do2018AB) {
    std::vector< std::pair<TString, float> > gJetsInputs = {std::make_pair(gJets18AB, 21.1)};
    RA2bin_inputs_Zinv(Signal, gJetsInputs, DR, DY18AB, zJets18AB, 21.0/59.6, true);
  }
  else if (do2018CD) {
    std::vector< std::pair<TString, float> > gJetsInputs = {std::make_pair(gJets18CD, 38.3)};
    RA2bin_inputs_Zinv(Signal, gJetsInputs, DR, DY18CD, zJets18CD, 38.6/59.6, true);
  }
  else if (do2018) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair(gJets18AB, 21.1),
      std::make_pair(gJets18CD, 38.3)
    };
    cout << "No DR, DY file for 2018" << endl;
    // RA2bin_inputs_Zinv(Signal, gJetsInputs, DR, DY18, zJets18, 59.6/59.6, true);
  }
  else if (doRun2) {
    std::vector< std::pair<TString, float> > gJetsInputs = {
      std::make_pair(gJets16, 35.9),
      std::make_pair(gJets17, 41.5),
      std::make_pair(gJets18AB, 21.0),
      std::make_pair(gJets18CD, 38.6)
    };
    RA2bin_inputs_Zinv(Signal, gJetsInputs, DR, DYRun2, zJetsRun2, 137/137, true);
  }


  // RA2bin_inputs_Zinv(LDP, gJetsInputs, "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 136.8/35.9, true);
  // RA2bin_inputs_Zinv(HDP, gJetsInputs, "../datFiles/DR_Run2", "../datFiles/DY_Run2", "../plots/histograms/ZinvMCttzMC174bin_2016v16.root", 136.8/35.9, true);




  gApplication->Terminate(0);
  
}

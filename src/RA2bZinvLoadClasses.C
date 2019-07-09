{
  gROOT->Reset();
  TString arch = gSystem->Getenv("SCRAM_ARCH");
  if (arch.Contains("slc6")) {
    gSystem->Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0/lib/libboost_program_options.so");
    gSystem->AddIncludePath(" -I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0-cms/include ");
  } else if (arch.Contains("slc7")) {
    gSystem->Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0/lib/libboost_program_options.so");
    gSystem->AddIncludePath(" -I/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0/include ");
  } else {
    cout << "Unknown SCRAM_ARCH = " << arch << "; can't load boost library.";
  }
  gROOT->ProcessLine(".L TreeConfig.C+");
  gROOT->ProcessLine(".L CCbinning.C+");
  gROOT->ProcessLine(".L CutManager.C+");
  gROOT->ProcessLine(".L EfficiencyAndPurity.C+");
  gROOT->ProcessLine(".L RA2bZinvAnalysis.C+");
}

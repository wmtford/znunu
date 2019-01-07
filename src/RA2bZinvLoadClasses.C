{
  gROOT->Reset();
  bool forSL6 = false;
  if (forSL6) {
    gSystem->Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0-cms/lib/libboost_program_options.so");
    gSystem->AddIncludePath(" -I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0-cms/include ");
  } else {  // sl7
    gSystem->Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0/lib/libboost_program_options.so");
    gSystem->AddIncludePath(" -I/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0/include ");
  }
  gROOT->ProcessLine(".L CCbinning.C+");
  gROOT->ProcessLine(".L RA2bZinvAnalysis.C+");
}

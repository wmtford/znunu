{
  gROOT->Reset();
  gSystem->Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0-cms/lib/libboost_program_options.so");
  gSystem->AddIncludePath(" -I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0-cms/include ");
  gROOT->ProcessLine(".L CCbinning.C+");
  gROOT->ProcessLine(".L RA2bZinvAnalysis.C+");
}

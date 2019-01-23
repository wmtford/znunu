/*
  Run from the command line with
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C
  or
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C\(\"2016B\"\)
*/

#include "TROOT.h"
#include "TEnv.h"

void RA2bZinvDriver(const std::string& runBlock = "") {
  cout << "runBlock passed to Driver = " << runBlock << endl;

  gEnv->SetValue("TFile.AsyncPrefetching", 1);

  bool doHzvv = false;
  bool doHttzvv = false;
  bool doHzmm = false;
  bool doHzee = false;
  bool doHphoton = true;
  bool doHdymm = false;
  bool doHdyee = false;
  bool doHttzmm = false;
  bool doHttzee = false;
  bool doHVVmm = false;
  bool doHVVee = false;
  bool doHttmm = false;
  bool doHttee = false;
  bool doHgjets = false;
  bool doHgjetsqcd = false;
  const std::string makeClassSample = "";  // Must be compatible with compiler directives
  bool doListTrigPrescales = false;
  const std::string dumpSelEvIDsample("");

  // RA2bZinvAnalysis analyzer("", runBlock);  // Default configuration, V12
  RA2bZinvAnalysis analyzer("data2016.cfg", runBlock);
  // RA2bZinvAnalysis analyzer("data2017.cfg", runBlock);
  // RA2bZinvAnalysis analyzer("data2018.cfg", runBlock);
  // RA2bZinvAnalysis analyzer("lowDphi.cfg", runBlock);

  std::string fnstr;

  if (doHzvv || doHttzvv) {
    fnstr = "histsZjets";  fnstr += runBlock + ".root";
    char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
    TFile *histoOutFile = TFile::Open(outfn, "RECREATE");
    TH1F *hCCzvv = nullptr;
    if (doHzvv) {
      std::vector<TH1*> h_zinv = analyzer.makeHistograms("zinv");
      for (auto& theHist : h_zinv) {
	theHist->Print();
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk") && !hName.Contains("Jb") && !hName.Contains("spl")) {
	  hCCzvv = (TH1F*) theHist;
	}
      }
    }
    if (doHttzvv) {
      std::vector<TH1*> h_ttzvv = analyzer.makeHistograms("ttzvv");
      for (auto& theHist : h_ttzvv) {
	theHist->Print();
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk") && !hName.Contains("Jb") && !hName.Contains("spl") && hCCzvv != nullptr) {
	  TH1F* hCCzinvAll = (TH1F*) hCCzvv->Clone();  hCCzinvAll->SetName("hCCzinvAll");  hCCzinvAll->Sumw2();
	  hCCzinvAll->Add(theHist);
	  hCCzinvAll->SetName("hCCzinvAll");
	  hCCzinvAll->Print();
	  hCCzinvAll->Draw();
	}
      }
    }
    histoOutFile->Write();
  }

  if (doHzmm || doHzee) {
    TFile *histoOutFile;
    fnstr = "histsZ";
    TH1F *hCCzmm = nullptr, *hCCjbzmm = nullptr, *hCCJbzmm = nullptr;
    if (doHzmm && doHzee) {
      fnstr += "ll" + runBlock + ".root";
      char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
      histoOutFile = TFile::Open(outfn, "RECREATE");
    }
    if (doHzmm) {
      if (!doHzee) {
	fnstr += "mm" + runBlock + ".root";
	char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
	histoOutFile = TFile::Open(outfn, "RECREATE");
      }
      std::vector<TH1*> h_zmm = analyzer.makeHistograms("zmm");
      for (auto& theHist : h_zmm) {
	theHist->Print();
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk") && !hName.Contains("Jb") && !hName.Contains("spl"))
	  hCCzmm = (TH1F*) theHist;
	if (hName.Contains("hCC") && hName.Contains("jb") && !hName.Contains("jk") && !hName.Contains("Jb") && !hName.Contains("spl"))
	  hCCjbzmm = (TH1F*) theHist;
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk") && hName.Contains("Jb") && !hName.Contains("spl"))
	  hCCJbzmm = (TH1F*) theHist;
      }
    }
    if (doHzee) {
      if (!doHzmm) {
	fnstr += "ee" + runBlock + ".root";
	char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
	histoOutFile = TFile::Open(outfn, "RECREATE");
      }
      std::vector<TH1*> h_zee = analyzer.makeHistograms("zee");
      for (auto& theHist : h_zee) {
	theHist->Print();
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk") && !hName.Contains("Jb") && !hName.Contains("spl") && hCCzmm != nullptr) {
	  TH1F* hCC_zll = (TH1F*) hCCzmm->Clone();  hCC_zll->SetName("hCC_zll");  hCC_zll->Sumw2();
	  hCC_zll->Add(theHist);
	  hCC_zll->Print();
	  hCC_zll->Draw();
	}
	if (hName.Contains("hCC") && hName.Contains("jb") && !hName.Contains("jk") && !hName.Contains("Jb") && !hName.Contains("spl") && hCCjbzmm != nullptr) {
	  TH1F* hCCjb_zll = (TH1F*) hCCjbzmm->Clone();  hCCjb_zll->SetName("hCCjb_zll");  hCCjb_zll->Sumw2();
	  hCCjb_zll->Add(theHist);
	  hCCjb_zll->Print();
	  hCCjb_zll->Draw();
	}
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk") && hName.Contains("Jb") && !hName.Contains("spl") && hCCJbzmm != nullptr) {
	  TH1F* hCCJb_zll = (TH1F*) hCCJbzmm->Clone();  hCCJb_zll->SetName("hCCJb_zll");  hCCJb_zll->Sumw2();
	  hCCJb_zll->Add(theHist);
	  hCCJb_zll->Print();
	  hCCJb_zll->Draw();
	}
      }
    }
    histoOutFile->Write();
  }

  if (doHdymm || doHdyee || doHttzmm || doHttzee || doHVVmm || doHVVee || doHttmm || doHttee) {
    fnstr = "histsDYMC";  fnstr += runBlock + ".root";
    char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
    TFile *histoOutFile = TFile::Open(outfn, "RECREATE");

    if (doHdymm) {
      std::vector<TH1*> h_dymm = analyzer.makeHistograms("dymm");
      for (auto& theHist : h_dymm) {theHist->Print();  theHist->Draw();}
    }
    if (doHdyee) {
      std::vector<TH1*> h_dyee = analyzer.makeHistograms("dyee");
      for (auto& theHist : h_dyee) {theHist->Print();  theHist->Draw();}
    }

    if (doHttzmm) {
      std::vector<TH1*> h_ttzmm = analyzer.makeHistograms("ttzmm");
      for (auto& theHist : h_ttzmm) {theHist->Print();  theHist->Draw();}
    }
    if (doHttzee) {
      std::vector<TH1*> h_ttzee = analyzer.makeHistograms("ttzee");
      for (auto& theHist : h_ttzee) {theHist->Print();  theHist->Draw();}
    }

    if (doHVVmm) {
      std::vector<TH1*> h_VVmm = analyzer.makeHistograms("VVmm");
      for (auto& theHist : h_VVmm) {theHist->Print();  theHist->Draw();}
    }
    if (doHVVee) {
      std::vector<TH1*> h_VVee = analyzer.makeHistograms("VVee");
      for (auto& theHist : h_VVee) {theHist->Print();  theHist->Draw();}
    }

    if (doHttmm) {
      std::vector<TH1*> h_ttmm = analyzer.makeHistograms("ttmm");
      for (auto& theHist : h_ttmm) {theHist->Print();  theHist->Draw();}
    }
    if (doHttee) {
      std::vector<TH1*> h_ttee = analyzer.makeHistograms("ttee");
      for (auto& theHist : h_ttee) {theHist->Print();  theHist->Draw();}
    }

    histoOutFile->Write();
  }

  if (doHphoton) {
    fnstr = "histsPhoton";  fnstr += runBlock + ".root";
    char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
    TFile *histoOutFile = TFile::Open(outfn, "RECREATE");
    std::vector<TH1*> h_photon = analyzer.makeHistograms("photon");
    for (auto& theHist : h_photon) {
      theHist->Print();
      theHist->Draw();
    }
    histoOutFile->Write();
  }

  if (doHgjets || doHgjetsqcd) {
    fnstr = "histsGjets";  fnstr += runBlock + ".root";
    char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
    TFile *histoOutFile = TFile::Open(outfn, "RECREATE");
    if (doHgjets) {
      std::vector<TH1*> h_gJets = analyzer.makeHistograms("gjets");
      for (auto& theHist : h_gJets) {theHist->Print();  theHist->Draw();}
    }
    if (doHgjetsqcd) {
      std::vector<TH1*> h_gjetsqcd = analyzer.makeHistograms("gjetsqcd");
      for (auto& theHist : h_gjetsqcd) {theHist->Print();  theHist->Draw();}
    }
    histoOutFile->Write();
  }

  if (doListTrigPrescales) {
    analyzer.checkTrigPrescales("zmm");
  }

  if (!makeClassSample.empty()) {
    analyzer.runMakeClass(makeClassSample);
  }

  if (!dumpSelEvIDsample.empty())
    analyzer.dumpSelEvIDs(dumpSelEvIDsample.data(), (std::string("evtIDs_") + runBlock + ".txt").data());

  gApplication->Terminate(0);

}

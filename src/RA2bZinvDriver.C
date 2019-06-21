/*
  Run from the command line with
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C
  or
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C\(\"2018AB\"\)
  or
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C\(\"2018AB,1\")
*/

#include "TROOT.h"
#include "TEnv.h"

void RA2bZinvDriver(const std::string& runBlock = "", int toDo = -1) {
  cout << "runBlock passed to Driver = " << runBlock << endl;

  gEnv->SetValue("TFile.AsyncPrefetching", 1);

  enum samples {ZllData, PhotonData, ZllMC, PhotonMC, ZinvMC};
  bool doHZllData = false;
  bool doHphotonData = false;
  bool doHZllMC = false;
  bool doHphotonMC = false;
  bool doHzinvMC = false;
  switch(toDo) {
  case ZllData: doHZllData = true;  cout << "Zll Data" << endl;  break;
  case PhotonData: doHphotonData = true;  cout << "Photon Data" << endl;  break;
  case ZllMC: doHZllMC = true;  cout << "Zll MC" << endl;  break;
  case PhotonMC: doHphotonMC = true;  cout << "Photon MC" << endl;;  break;
  case ZinvMC: doHzinvMC = true;  cout << "Zinv MC" << endl;  break;
  default: cout << "Default" << endl;  break;
  }

  bool doHzvv = false;
  bool doHttzvv = false;
  bool doHzmm = false;
  bool doHzee = false;
  bool doHphoton = false;
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
  if (doHzinvMC) {
    doHzvv = true;
    doHttzvv = true;
  }
  if (doHZllData) {
    doHzmm = true;
    doHzee = true;
  }
  if (doHZllMC) {
    doHdymm = true;
    doHdyee = true;
    doHttzmm = true;
    doHttzee = true;
    doHVVmm = true;
    doHVVee = true;
    doHttmm = true;
    doHttee = true;
  }
  if (doHphotonData) doHphoton = true;
  if (doHphotonMC) {
    doHgjets = true;
    doHgjetsqcd = true;
  }

  bool doListTrigPrescales = false;
  const std::string dumpSelEvIDsample("");

  std::string cfgName("data");  cfgName += runBlock;  cfgName += ".cfg";
  cfgName = regex_replace(cfgName, regex("HEM"), "");
  cfgName = regex_replace(cfgName, regex("HEP"), "");
  cfgName = regex_replace(cfgName, regex("AB"), "");
  cfgName = regex_replace(cfgName, regex("CD"), "");

  void combine(std::vector<TH1*> hl1, std::vector<TH1*> hl2);
  std::string fnstr;
  RA2bZinvAnalysis analyzer(cfgName, runBlock);
    
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
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk")
	    && !hName.Contains("Jb") && !hName.Contains("spl")) {
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
	if (hName.Contains("hCC") && !hName.Contains("jb") && !hName.Contains("jk")
	    && !hName.Contains("Jb") && !hName.Contains("spl") && hCCzvv != nullptr) {
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
    std::vector<TH1*> h_zmm, h_zee;
    TFile *histoOutFile;
    fnstr = "histsZ";
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
      h_zmm = analyzer.makeHistograms("zmm");
      for (auto& theHist : h_zmm) {
	theHist->Print();
	theHist->Draw();
      }
    }
    if (doHzee) {
      if (!doHzmm) {
	fnstr += "ee" + runBlock + ".root";
	char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
	histoOutFile = TFile::Open(outfn, "RECREATE");
      }
      h_zee = analyzer.makeHistograms("zee");
      for (auto& theHist : h_zee) {
	theHist->Print();
	theHist->Draw();
      }
      if (doHzmm && doHzee) combine(h_zmm, h_zee);
    }
    histoOutFile->Write();
  }

  if (doHdymm || doHdyee || doHttzmm || doHttzee || doHVVmm || doHVVee || doHttmm || doHttee) {
    fnstr = "histsDYMC";  fnstr += runBlock + ".root";
    char* outfn = new char[fnstr.length()+1];  std::strcpy (outfn, fnstr.c_str());
    TFile *histoOutFile = TFile::Open(outfn, "RECREATE");
    std::vector<TH1*> h_dymm, h_dyee, h_ttzmm, h_ttzee, h_VVmm, h_VVee, h_ttmm, h_ttee;

    if (doHdymm) {
      h_dymm = analyzer.makeHistograms("dymm");
      for (auto& theHist : h_dymm) {theHist->Print();  theHist->Draw();}
    }
    if (doHdyee) {
      h_dyee = analyzer.makeHistograms("dyee");
      for (auto& theHist : h_dyee) {theHist->Print();  theHist->Draw();}
    }
    if (doHdymm && doHdyee) combine(h_dymm, h_dyee);

    if (doHttzmm) {
      h_ttzmm = analyzer.makeHistograms("ttzmm");
      for (auto& theHist : h_ttzmm) {theHist->Print();  theHist->Draw();}
    }
    if (doHttzee) {
      h_ttzee = analyzer.makeHistograms("ttzee");
      for (auto& theHist : h_ttzee) {theHist->Print();  theHist->Draw();}
    }
    if (doHttzmm && doHttzee) combine(h_ttzmm, h_ttzee);

    if (doHVVmm) {
      h_VVmm = analyzer.makeHistograms("VVmm");
      for (auto& theHist : h_VVmm) {theHist->Print();  theHist->Draw();}
    }
    if (doHVVee) {
      h_VVee = analyzer.makeHistograms("VVee");
      for (auto& theHist : h_VVee) {theHist->Print();  theHist->Draw();}
    }
    if (doHVVmm && doHVVee) combine(h_VVmm, h_VVee);

    if (doHttmm) {
      h_ttmm = analyzer.makeHistograms("ttmm");
      for (auto& theHist : h_ttmm) {theHist->Print();  theHist->Draw();}
    }
    if (doHttee) {
      h_ttee = analyzer.makeHistograms("ttee");
      for (auto& theHist : h_ttee) {theHist->Print();  theHist->Draw();}
    }
    if (doHttmm && doHttee) combine(h_ttmm, h_ttee);

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

  if (!dumpSelEvIDsample.empty()) {
    analyzer.dumpSelEvIDs(dumpSelEvIDsample.data(), (std::string("evtIDs_") + runBlock + ".txt").data());
  }

  gApplication->Terminate(0);

}

void combine(std::vector<TH1*> hl1, std::vector<TH1*> hl2) {
  for (auto& h1 : hl1) {
    TString h1Name(h1->GetName());
    if (!h1Name.Contains("hCC")) continue;
    TString heeName = h1Name;  heeName("mm") = "ee";
    TString hllName = h1Name;  hllName("mm") = "ll";
    for (auto& h2 : hl2) {
      TString h2Name(h2->GetName());
      if (h2Name.EqualTo(heeName)) {
	TH1D* hll = (TH1D*) h1->Clone();  hll->SetNameTitle(hllName, hllName);
	hll->Add(h2);
	hll->Print();
	hll->Draw();
      }
    }
  }
}

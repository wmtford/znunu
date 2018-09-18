{
  /*
  Run from the command line with
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C
  */

#include "TROOT.h"
#include "TEnv.h"

  gEnv->SetValue("TFile.AsyncPrefetching", 1);

  bool doHzvv = false;
  bool doHttzvv = false;
  bool doHzmm = true;
  bool doHzee = true;
  bool doHdymm = false;
  bool doHdyee = false;
  bool doHttzmm = false;
  bool doHttzee = false;
  bool doHVVmm = false;
  bool doHVVee = false;
  bool doHttmm = false;
  bool doHttee = false;
  const std::string makeClassSample = "";  // Must be compatible with compiler directives
  bool doListTrigPrescales = false;

  RA2bZinvAnalysis analyzer;  // Default configuration
  // RA2bZinvAnalysis analyzer("data2018.cfg");
  // RA2bZinvAnalysis analyzer("lowDphi.cfg");

  if (doHzvv || doHttzvv) {
    TFile *histoOutFile = TFile::Open("histsZjets.root", "RECREATE");
    TH1F *hCCzvv = nullptr;
    if (doHzvv) {
      std::vector<TH1*> h_zinv = analyzer.makeHistograms("zinv");
      for (auto& theHist : h_zinv) {
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC")) {
	  hCCzvv = (TH1F*) theHist;
	}
      }
    }
    if (doHttzvv) {
      std::vector<TH1*> h_ttzvv = analyzer.makeHistograms("ttzvv");
      for (auto& theHist : h_ttzvv) {
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC") && hCCzvv) {
	  TH1F* hCCzinvAll = (TH1F*) hCCzvv->Clone();
	  hCCzinvAll->Add(theHist);
	  hCCzinvAll->SetName("hCCzinvAll");
	  hCCzinvAll->Draw();
	}
      }
    }
    histoOutFile->Write();
  }

  if (doHzmm || doHzee) {
    TFile *histoOutFile;
    TH1F *hCCzmm = nullptr;
    if (doHzmm && doHzee) histoOutFile = TFile::Open("histsZll.root", "RECREATE");
    if (doHzmm) {
      if (!doHzee) histoOutFile = TFile::Open("histsZmm.root", "RECREATE");
      std::vector<TH1*> h_zmm = analyzer.makeHistograms("zmm");
      for (auto& theHist : h_zmm) {
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC")) {
	  hCCzmm = (TH1F*) theHist;
	}
      }
    }
    if (doHzee) {
      if (!doHzmm) histoOutFile = TFile::Open("histsZee.root", "RECREATE");
      std::vector<TH1*> h_zee = analyzer.makeHistograms("zee");
      for (auto& theHist : h_zee) {
	theHist->Draw();
	TString hName(theHist->GetName());
	if (hName.Contains("hCC") && hCCzmm) {
	  TH1F* hCCll = (TH1F*) hCCzmm->Clone();
	  hCCll->Add(theHist);
	  hCCll->SetName("hCCll");
	  hCCll->Draw();
	}
      }
    }
    histoOutFile->Write();
  }

  if (doHdymm || doHdyee || doHttzmm || doHttzee || doHVVmm || doHVVee || doHttmm || doHttee) {
    TFile *histoOutFile = TFile::Open("histsDYMC.root", "RECREATE");

    if (doHdymm) {
      std::vector<TH1*> h_dymm = analyzer.makeHistograms("dymm");
      for (auto& theHist : h_dymm) theHist->Draw();
    }
    if (doHdyee) {
      std::vector<TH1*> h_dyee = analyzer.makeHistograms("dyee");
      for (auto& theHist : h_dyee) theHist->Draw();
    }

    if (doHttzmm) {
      std::vector<TH1*> h_ttzmm = analyzer.makeHistograms("ttzmm");
      for (auto& theHist : h_ttzmm) theHist->Draw();
    }
    if (doHttzee) {
      std::vector<TH1*> h_ttzee = analyzer.makeHistograms("ttzee");
      for (auto& theHist : h_ttzee) theHist->Draw();
    }

    if (doHVVmm) {
      std::vector<TH1*> h_VVmm = analyzer.makeHistograms("VVmm");
      for (auto& theHist : h_VVmm) theHist->Draw();
    }
    if (doHVVee) {
      std::vector<TH1*> h_VVee = analyzer.makeHistograms("VVee");
      for (auto& theHist : h_VVee) theHist->Draw();
    }

    if (doHttmm) {
      std::vector<TH1*> h_ttmm = analyzer.makeHistograms("ttmm");
      for (auto& theHist : h_ttmm) theHist->Draw();
    }
    if (doHttee) {
      std::vector<TH1*> h_ttee = analyzer.makeHistograms("ttee");
      for (auto& theHist : h_ttee) theHist->Draw();
    }

    histoOutFile->Write();
  }

  if (doListTrigPrescales) {
    analyzer.checkTrigPrescales("zmm");
  }

  if (!makeClassSample.empty()) {
    analyzer.runMakeClass(makeClassSample);
  }

  gApplication->Terminate(0);

}

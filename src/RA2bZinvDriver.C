{
  /*
  Run from the command line with
  root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C
  */

#include "TROOT.h"
#include "TEnv.h"

  gEnv->SetValue("TFile.AsyncPrefetching", 1);

  bool doCCzvv = true;
  bool doCCttzvv = true;
  bool do1Dzvv = false;
  bool do1Dttzvv = false;
  bool do1Dzmm = false;
  bool do1Dzee = false;
  bool do1Ddymm = false;
  bool do1Ddyee = false;
  bool do1Dttzmm = false;
  bool do1Dttzee = false;
  bool do1DVVmm = false;
  bool do1DVVee = false;
  bool do1Dttmm = false;
  bool do1Dttee = false;
  bool doMakeClass = false;
  bool doListTrigPrescales = false;

  RA2bZinvAnalysis analyzer;

  if (doCCzvv || doCCttzvv) {
    // Output file
    TFile *CChistos = TFile::Open("hCCzinv.root", "RECREATE");

    if (doCCzvv) {
      TH1F* hCCzvv = analyzer.makeCChist("zinv");
      hCCzvv->SetName("hCCzvv");  hCCzvv->Draw();
    }
    if (doCCttzvv) {
      TH1F* hCCttzvv = analyzer.makeCChist("ttzvv");
      hCCttzvv->SetName("hCCttzvv");  hCCttzvv->Draw();
    }
    if (doCCzvv && doCCttzvv) {
      TH1F* hCCzinv = (TH1F*) hCCzvv->Clone();
      hCCzinv->Add(hCCttzvv);
      hCCzinv->SetName("hCCzinv");  hCCzinv->Draw();
    }
    CChistos->Write();
  }

  if (do1Dzvv || do1Dttzvv) {
    TFile *histos1D = TFile::Open("histsZjets.root", "RECREATE");
    if (do1Dzvv) {
      std::vector<TH1F*> h_zinv = analyzer.makeHistograms("zinv");
      for (auto& theHist : h_zinv) theHist->Draw();
    }
    if (do1Dttzvv) {
      std::vector<TH1F*> h_ttzvv = analyzer.makeHistograms("ttzvv");
      for (auto& theHist : h_ttzvv) theHist->Draw();
    }
    histos1D->Write();
  }

  if (do1Dzmm || do1Dzee) {
    TFile *histos1D = TFile::Open("histsDY.root", "RECREATE");
    if (do1Dzmm) {
      std::vector<TH1F*> h_zmm = analyzer.makeHistograms("zmm");
      for (auto& theHist : h_zmm) theHist->Draw();
    }
    if (do1Dzee) {
      std::vector<TH1F*> h_zee = analyzer.makeHistograms("zee");
      for (auto& theHist : h_zee) theHist->Draw();
    }
    histos1D->Write();
  }

  if (do1Ddymm || do1Ddyee || do1Dttzmm || do1Dttzee || do1DVVmm || do1DVVee || do1Dttmm || do1Dttee) {
    TFile *histos1D = TFile::Open("histsDYMC.root", "RECREATE");

    if (do1Ddymm) {
      std::vector<TH1F*> h_dymm = analyzer.makeHistograms("dymm");
      for (auto& theHist : h_dymm) theHist->Draw();
    }
    if (do1Ddyee) {
      std::vector<TH1F*> h_dyee = analyzer.makeHistograms("dyee");
      for (auto& theHist : h_dyee) theHist->Draw();
    }

    if (do1Dttzmm) {
      std::vector<TH1F*> h_ttzmm = analyzer.makeHistograms("ttzmm");
      for (auto& theHist : h_ttzmm) theHist->Draw();
    }
    if (do1Dttzee) {
      std::vector<TH1F*> h_ttzee = analyzer.makeHistograms("ttzee");
      for (auto& theHist : h_ttzee) theHist->Draw();
    }

    if (do1DVVmm) {
      std::vector<TH1F*> h_VVmm = analyzer.makeHistograms("VVmm");
      for (auto& theHist : h_VVmm) theHist->Draw();
    }
    if (do1DVVee) {
      std::vector<TH1F*> h_VVee = analyzer.makeHistograms("VVee");
      for (auto& theHist : h_VVee) theHist->Draw();
    }

    if (do1Dttmm) {
      std::vector<TH1F*> h_ttmm = analyzer.makeHistograms("ttmm");
      for (auto& theHist : h_ttmm) theHist->Draw();
    }
    if (do1Dttee) {
      std::vector<TH1F*> h_ttee = analyzer.makeHistograms("ttee");
      for (auto& theHist : h_ttee) theHist->Draw();
    }

    histos1D->Write();
  }

  if (doListTrigPrescales) {
    analyzer.checkTrigPrescales("zmm");
  }

  if (doMakeClass) {
    // analyzer.runMakeClass("zinv", "MC_V12");
    // analyzer.runMakeClass("zmm", "data_V12");
    analyzer.runMakeClass("zmm", "data_V15");
  }

  gApplication->Terminate(0);

}

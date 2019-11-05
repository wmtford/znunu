
// Combine histograms from multiple files

void mergeHists() {
  gROOT->Reset();

  bool addLinearXlast = false;
  std::vector<Double_t> scale;
  TString matchHname(""), newHname("");

  // // Zinv expectation histogram
  // // std::vector<const char*> mmFiles = {"../outputs/histsZjets_2016v16_noPU.root"};
  // // std::vector<const char*> mmFiles = {"../outputs/histsZjets_2017v16_HT17wt_ZptWt.root"};
  // // std::vector<const char*> mmFiles = {"histsZjets2017.root"};
  // // std::vector<const char*> mmFiles = {"../outputs/histsZjets_2018v16_HT17wt_ZptWt.root"};  scale = {21.0/59.2};
  // // std::vector<const char*> mmFiles = {"../outputs/histsZjets_2018HEMv16_HT17wt_ZptWt.root"};scale = {38.2/59.2};
  // std::vector<const char*> mmFiles = {"../outputs/histsZjets_2016v17.root",
  //                                     "../outputs/histsZjets_2017v17.root",
  // 				      "../outputs/histsZjets_2018v17.root"};
  // // scale = {1, 21.0/59.2, 38.2/59.2, 1};
  // // addLinearXlast = true;
  // // matchHname = "hCCzinvAll";
  // // newHname += "plot_zinv_nj5_nb4_kin10_1";

  // // gJets MC
  // std::vector<const char*> mmFiles = {
  //   "../outputs/histsGjets_2016v17_DRr2wt.root",
  //   "../outputs/histsGjets_2017v17_DRr2wt.root",
  //   "../outputs/histsGjets_2018v17_DRr2wt.root"
  //   // "../outputs/histsGjetsldpnominal_2016v16_DRr2wt.root",
  //   // "../outputs/histsGjetsldpnominal_2017v16_DRr2wt.root",
  //   // "../outputs/histsGjetsldpnominal_2018v17_DRr2wt.root"
  // };
  // // scale = {1, 21.0/59.2, 38.2/59.2, 1};
  // // addLinearXlast = true;

  // // DY MC (mm and ee combined in each file)
  // std::vector<const char*> mmFiles = {
  //   "../outputs/histsDYMC_2016v17_corrSF.root",
  //   "../outputs/histsDYMC_2017v17_corrSF.root",
  //   "../outputs/histsDYMC_2018v17_corrSF.root",
  // };
  // // Consider multiplying DYMC by 2017 k factor instead of old 2016:  scale = 1.165/1.23
  // // To update 2016, scale = 1.257/1.23
  // // scale = {1, 21.0/59.2, 38.2/59.2, 1};
  // // addLinearXlast = true;

  // // DY data (mm and ee combined in each file)
  // std::vector<const char*> mmFiles = {
  //   "../outputs/histsDY_2016v17.root",
  //   "../outputs/histsDY_2017v17.root",
  //   "../outputs/histsDY_2018v17.root"
  //   // "../outputs/histsDY_2018ABv16.root",
  //   // "../outputs/histsDY_2018CDv16.root"
  // };

  // // Photon data
  // std::vector<const char*> mmFiles = {
  //   "../outputs/histsPhoton_2016v17.root",
  //   "../outputs/histsPhoton_2017v17.root",
  //   "../outputs/histsPhoton_2018v17.root"
  //   // "../outputs/histsPhoton_2018ABv16.root",
  //   // "../outputs/histsPhoton_2018CDv16.root"
  // };

  // DY data in separate mm, ee files
  std::vector<const char*> eeFiles;
  // std::vector<const char*> mmFiles = {
  //   "histsZmm2017B.root",
  //   "histsZmm2017C.root",
  //   "histsZmm2017D.root",
  //   "histsZmm2017E.root",
  //   "histsZmm2017F.root"
  // };
  // std::vector<const char*> eeFiles = {
  //   "histsZee2017B.root",
  //   "histsZee2017C.root",
  //   "histsZee2017D.root",
  //   "histsZee2017E.root",
  //   "histsZee2017F.root"
  // };

  // std::vector<const char*> mmFiles = {
  //   "histsZmm2018A.root",
  //   "histsZmm2018B.root"
  // };
  // std::vector<const char*> eeFiles = {
  //   "histsZee2018A.root",
  //   "histsZee2018B.root"
  // };

  std::vector<std::vector<const char*>> inFiles;
  inFiles.push_back(mmFiles);
  if (!eeFiles.empty()) inFiles.push_back(eeFiles);
  std::vector<TH1*> histos;
  size_t hdx, hdx0 = 0;
  for (size_t iem = 0; iem < inFiles.size(); ++iem) {
    for (size_t fdx = 0; fdx < inFiles[iem].size(); ++fdx) {
      Double_t wt = 1;
      if (!scale.empty()) wt = scale[fdx];
      cout << inFiles[iem][fdx] << endl;
      TFile *f1 = TFile::Open(inFiles[iem][fdx]);
      TIter keyList(f1->GetListOfKeys());
      TKey *key;
      hdx = hdx0;
      while ((key = (TKey*)keyList())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	// if (!cl->InheritsFrom("TH1")) continue;
	if (!(cl->InheritsFrom("TH1F") || cl->InheritsFrom("TH1D"))) continue;
	TH1 *h = (TH1*)key->ReadObj();
	if (!matchHname.IsNull() && !matchHname.EqualTo(h->GetName())) continue;
	if (!newHname.IsNull()) h->SetNameTitle(newHname, newHname);
	if (fdx == 0) {
	  histos.push_back((TH1*) h->Clone());
	  histos.back()->Scale(wt);
	} else {
	  if (addLinearXlast && fdx < inFiles[iem].size() - 1) {
	    for (int binx = 1; binx <= h->GetNbinsX(); ++binx) {
	      histos[hdx]->SetBinContent(binx, histos[hdx]->GetBinContent(binx) + wt*(h->GetBinContent(binx)));
	      histos[hdx]->SetBinError(binx, histos[hdx]->GetBinError(binx) + wt*(h->GetBinError(binx)));
	    }
	  }
	  else
	    histos[hdx]->Add(h, wt);
	}
	++hdx;
      }
    }
    hdx0 = hdx;
  }
  TFile *histoOutFile = TFile::Open("mergedHists.root", "RECREATE");
  for (auto& theHist : histos) theHist->SetDirectory(gDirectory);
  histoOutFile->Write();

  gApplication->Terminate(0);
}

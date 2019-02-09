
// Combine histograms from multiple files

void mergeHists() {
  gROOT->Reset();

  bool addLinear12 = false;
  std::vector<Double_t> scale;
  TString matchHname(""), newHname("");

  std::vector<const char*> mmFiles = {"../outputs/histsZjets_2016v16.root"};
  matchHname = "hCCzinvAll";
  newHname += "plot_zinv_nj5_nb4_kin10_1";

  // std::vector<const char*> mmFiles = {
  //   "../outputs/histsGjets_2017v16_DRwt2.root",
  //   "../outputs/histsGjets_2018v16_DRwt2.root",
  //   "../outputs/histsGjets_2016v16_DRwt2.root"
  // };
  // addLinear12 = true;

  // std::vector<const char*> mmFiles = {
  //   "../outputs/histsDYMC_2016v16.root",
  //   "../outputs/histsDYMC_2017v16.root"
  // };
  // scale = {1, (41.5+59.4)/41.5};

  // std::vector<const char*> mmFiles = {
  //   "../outputs/histsDY_2016v16.root",
  //   "../outputs/histsDY_2017v16.root",
  //   "../outputs/histsDY_2018v16.root"
  // };

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
	h->SetNameTitle(newHname, newHname);
	if (fdx == 0) {
	  histos.push_back((TH1*) h->Clone());
	  histos.back()->Scale(wt);
	} else {
	  if (addLinear12 && fdx == 1) {
	    for (int binx = 1; binx <= h->GetNbinsX(); ++binx) {
	      histos[hdx]->SetBinContent(binx, histos[hdx]->GetBinContent(binx) + h->GetBinContent(binx));
	      histos[hdx]->SetBinError(binx, histos[hdx]->GetBinError(binx) + h->GetBinError(binx));
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
}

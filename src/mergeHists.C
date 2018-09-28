
// Combine histograms from multiple files

void mergeHists() {
  gROOT->Reset();
  std::vector<const char*> mmFiles = {
    "histsZmm2017B.root",
    "histsZmm2017C.root",
    "histsZmm2017D.root",
    "histsZmm2017E.root",
    "histsZmm2017F.root"
  };
  std::vector<const char*> eeFiles = {
    "histsZee2017B.root",
    "histsZee2017C.root",
    "histsZee2017D.root",
    "histsZee2017E.root",
    "histsZee2017F.root"
  };
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
  inFiles.push_back(eeFiles);
  std::vector<TH1*> histos;
  size_t hdx, hdx0 = 0;
  for (size_t iem = 0; iem < inFiles.size(); ++iem) {
    for (size_t fdx = 0; fdx < inFiles[iem].size(); ++fdx) {
      cout << inFiles[iem][fdx] << endl;
      TFile *f1 = TFile::Open(inFiles[iem][fdx]);
      TIter keyList(f1->GetListOfKeys());
      TKey *key;
      hdx = hdx0;
      while ((key = (TKey*)keyList())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();
	if (fdx == 0) histos.push_back((TH1*) h->Clone());
	else histos[hdx]->Add(h);
	++hdx;
      }
    }
    hdx0 = hdx;
  }
  TFile *histoOutFile = TFile::Open("mergedHists.root", "RECREATE");
  for (auto& theHist : histos) theHist->SetDirectory(gDirectory);
  histoOutFile->Write();
}

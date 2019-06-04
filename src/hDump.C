void hDump(TH1* histo) {

  for (int i=1; i<= histo->GetNbinsX(); ++i) {
    double y = histo->GetBinContent(i);
    double yerr = histo->GetBinError(i);
    double ferr = y > 0 ? yerr/y : 0;
    char buf[256];
    sprintf(buf, "%s %11.3f %8.3f (%5.1f\%)", histo->GetXaxis()->GetBinLabel(i), y, yerr, 100*ferr);
    cout << buf << endl;
  }

}

#include <vector>
#include <cstdio>
#include <iostream>   // std::cout
#include <string>     // std::string, std::to_string
#include "tdrstyle.C"
#include "CMS_lumi.C"

using namespace std;

/*

root.exe -b -q 'Plot_searchBin_full.C("stacked","searchH_b","Elog408_","Elog404_")'
root.exe -b -q 'Plot_searchBin_full.C("stacked","QCD_Low",  "Elog408_","Elog404_")'
root.exe -b -q 'Plot_searchBin_full.C("stacked","QCD_Up",   "Elog408_","Elog404_")'

.L Plot_searchBin_full.C
Plot_searchBin_full("stacked","searchH_b","Elog365_");
Plot_searchBin_full("stacked","QCD_Low","Elog365_");
Plot_searchBin_full("stacked","QCD_Up","Elog365_");

root.exe -b -q 'Plot_searchBin_full.C("stacked","searchH_b","Elog365_")'
root.exe -b -q 'Plot_searchBin_full.C("stacked","QCD_Low","Elog365_")'
root.exe -b -q 'Plot_searchBin_full.C("stacked","QCD_Up","Elog365_")'

*/

void shift_bin(TH1* input, TH1* output){
  for (int ibin=1;ibin<=input->GetNbinsX();ibin++){
    output->SetBinContent(ibin,input->GetBinContent(ibin));    
    output->SetBinError(ibin,input->GetBinError(ibin));    
  }
}

void Plot_searchBin_full(string sample="signal",
			 string histname="ZinvBGpred",
			 string MChistname="hMCexp",
			 // string MChistname="plot_zinv_nj5_nb4_kin10_1",
			 string elog="",string elogExp="", int pull=0) {

  bool ZnnMCvsData = false;
  bool ZnnMCvsZllMC = false;
  bool ZllMCvsZllMC = true;
  bool ZllDataVsZllMC = false;
  bool PhoDataVsPhoMC = false;

  bool Y2016 = false;
  bool Y2017 =false;
  bool Y2018 = false;
  bool Run2 = true;

  TString CMSlabel, legendHeader, predLegend, expectLegend, ratioOrdinateTitle;
  char tempname[200], MCtempname[200];
  std::vector<Float_t> legCoord = {0.62, 0.55, 0.95, 0.77};  // (x1, y1, x2, y2)

  if (ZnnMCvsData) {
    legendHeader = "Z#rightarrow#nu#bar{#nu} background";
    MChistname="plot_zinv_nj5_nb4_kin10_1";
    CMSlabel = "#bf{CMS} #scale[0.76]{#it{Preliminary}}";
    ratioOrdinateTitle = "#scale[0.75]{#frac{Expectation}{Prediction}} ";
    expectLegend = "Expectation from simulation";
    predLegend = "Prediction from data";
  } else {
    MChistname="hMCexp";
    // CMSlabel = "#bf{CMS} #scale[0.76]{#it{Simulation Preliminary}}";
    CMSlabel = "#bf{CMS} #scale[0.76]{#it{Simulation}}";
    ratioOrdinateTitle = "#frac{Direct}{Prediction} ";
    if (ZnnMCvsZllMC) {
      legendHeader = "Z#rightarrow#nu#bar{#nu} background";
      expectLegend = "Z#rightarrow#nu#bar{#nu} from simulation";
      predLegend = "Treat Z#rightarrow ll simulation as data";
    } else if (ZllMCvsZllMC) {
      legCoord[0] = 0.67;
      legendHeader = "Z#rightarrow l^{+}l^{-}+jets yield";
      expectLegend = "Direct from simulation";
      predLegend = "Treat simulation as data";
    } else if (ZllDataVsZllMC) {
      legendHeader = "Z#rightarrow ll control sample";
      expectLegend = "Direct from Z#rightarrow ll data";
      predLegend = "Predict from Z#rightarrow ll data";
    } else if (PhoDataVsPhoMC) {
      legendHeader = "Photon control sample";
      expectLegend = "Direct from photon data";
      predLegend = "Predict from photon data";
    }
  }
  double lumi;
  if (Y2016) {
    if (ZnnMCvsData) {
      ;
    } else {
      if (ZnnMCvsZllMC) ;
      else if (ZllMCvsZllMC)
	sprintf(tempname, "/usr/users/wtford/cms/znunu/outputs/hClosure_Zll_2016v16.root");
      sprintf(MCtempname, "%s", tempname);  // (same as prediction file for closure)
    }
    lumi = 35.9;
  }
  else if (Y2017) {
    if (ZnnMCvsData) {
      ;
    } else {
      if (ZnnMCvsZllMC) ;
      else if (ZllMCvsZllMC)
	sprintf(tempname, "/usr/users/wtford/cms/znunu/outputs/hClosure_Zll_2017v16.root");
      sprintf(MCtempname, "%s", tempname);  // (same as prediction file for closure)
    }
    lumi = 41.5;
  }
  else if (Y2018) {
    if (ZnnMCvsData) {
      ;
    } else {
      if (ZnnMCvsZllMC) ;
      else if (ZllMCvsZllMC)
	sprintf(tempname, "/usr/users/wtford/cms/znunu/outputs/hClosure_Zll_2018v17.root");
      sprintf(MCtempname, "%s", tempname);  // (same as prediction file for closure)
    }
    lumi = 59.6;
  }
  else if (Run2) {
    if (ZnnMCvsData) {
      sprintf(tempname, "/usr/users/wtford/cms/znunu/zinvData_2019May7_Run2/ZinvHistos.root");
      sprintf(MCtempname, "/usr/users/wtford/cms/znunu/zinvData_2019May7_Run2/ZinvMCttzMC174bin_Run2v17.root");
    } else {
      if (ZnnMCvsZllMC)
	sprintf(tempname, "/usr/users/wtford/cms/znunu/outputs/hClosure_Zinv_Run2v161617.root");
      else if (ZllMCvsZllMC)
	sprintf(tempname, "/usr/users/wtford/cms/znunu/outputs/hClosure_Zll_Run2v17.root");
      sprintf(MCtempname, "%s", tempname);  // (same as prediction file for closure)
    }
    lumi = 137;
  }
  double lumi_ref = lumi; // Histos already normalized

  // Prediction file
  TFile * EstFile = new TFile(tempname,"R");
  printf("Opened %s\n", tempname);
  // Expectation file
  TFile * GenFile = new TFile(MCtempname,"R");
  printf("Opened %s\n", MCtempname);

  bool doBinShift = false;

  ///////////////////////////////////////////////////////////////////////////////////////////
  ////Some cosmetic work for official documents.
  //
  // Set basic style
  //
  //  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  gStyle->SetPalette(1) ; // for better color output
  //  gROOT->LoadMacro("CMS_lumi.C");

  //
  // Canvas size
  // int W = 1200;
  int W = 960;
  int H = 600;
  int H_ref = 600;
  int W_ref = 800;
  float T = 0.10*H_ref;
  float B = 0.06*H_ref;
  float L = 0.16*W_ref;
  float R = 0.04*W_ref;

  //
  // Various vertical line coordinates
  // float ymax_top = 200000.;
  float ymax_top = 2000000.;
  // float ymax_top = 40000.;
  // float ymin_top = 0.015 -0.01;  // wtf
  float ymin_top = 0.015;

  float ymax2_top = 20000.;
  // float ymax2_top = 1000.;
  float ymax3_top = 200.;
  float ymax4_top = 30.;

  // float ymax_bottom = 1.99 +2.3;
  float ymax_bottom = 1.99;
  float ymin_bottom = 0.01;

  float ymax2_bottom = 2.15;
  float ymax3_bottom = 2.15;
  float ymax4_bottom = 2.15;
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  //
  // More specific style set, opening input files etc

  gStyle->SetOptStat(0);  ///to avoid the stat. on the plots
  //gStyle->SetErrorX(0);

  //
  // Define legend
  //

  TLegend* catLeg1 = new TLegend(legCoord[0], legCoord[1], legCoord[2], legCoord[3]);
  catLeg1->SetTextSize(0.04);
  catLeg1->SetTextFont(42);
  catLeg1->SetFillColor(0);
  catLeg1->SetLineColor(1);
  catLeg1->SetBorderSize(1);

  //
  // Define canvas
  //
  TCanvas *canvas = new TCanvas("canvas","canvas",10,10,W,H);

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( R/W );
  canvas->SetTopMargin( T/H );
  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

  canvas->Divide(1, 2);
  
  //
  // Define pads
  //
  TPad* canvas_up = (TPad*) canvas->GetListOfPrimitives()->FindObject("canvas_1");
  TPad* canvas_dw = (TPad*) canvas->GetListOfPrimitives()->FindObject("canvas_2");

  //
  // define the size
  double up_height     = 0.8;  // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.30; // please tune so that the smaller canvas size will work in your environment
  double font_size_dw  = 0.1;  // please tune the font size parameter for bottom figure
  double dw_height     = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.04; // KH, added to put the bottom one closer to the top panel

  //
  // set pad size
  canvas_up->SetPad(0., 1 - up_height,    1., 1.00);
  canvas_dw->SetPad(0., 0.,               1., dw_height+dw_height_offset);
  //
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_up->SetTopMargin(0.12);
  canvas_up->SetLeftMargin(0.1);
  //
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  canvas_dw->SetBottomMargin(0.35);
  canvas_dw->SetTopMargin(0);
  canvas_dw->SetLeftMargin(0.1);
  
  //
  // draw top figure
  canvas_up->cd();

  TH1D * GenHist, * EstHist,* thist, *EstSysUpHist, *EstSysLowHist;
  TH1D * GenHist_input, * EstHist_input, *EstSysUpHist_input, *EstSysLowHist_input;
  TH1D * histTemplate;
  THStack *tempstack;

  double HT_x_max=2500.;
  double HT_x_min=400.;
  double MHT_x_max=1000.;
  double NJet_x_max=15.;
  double NBtag_x_max=4.;
  // double search_x_max=73.-0.5;
  // double search_x_max=161.-0.5;
  double search_x_max=175.-0.5;
  if(histname.find("QCD")!=string::npos)search_x_max=224.;
  double search_x_min=1.-0.5;

  sprintf(tempname,"%s",histname.c_str());
  sprintf(MCtempname,"%s",MChistname.c_str());
  cout << "histname = " << tempname << ", MChistname = " << MChistname << endl;
  if(sample.find("stacked")!=string::npos){
    tempstack=(THStack*)EstFile->Get(tempname)->Clone();
    EstHist_input=(TH1D*) tempstack->GetStack()->Last();
    tempstack=(THStack*)GenFile->Get(MCtempname)->Clone();   
    GenHist_input=(TH1D*) tempstack->GetStack()->Last();
    /*
    tempstack=(THStack*)EstFile->Get(tempname)->Clone();
    EstHistD=(TH1D*) tempstack->GetStack()->Last();
    tempstack=(THStack*)GenFile->Get(tempname)->Clone();
    GenHistD=(TH1D*) tempstack->GetStack()->Last();    
    */
  }
  else{
    if (!EstFile->Get(tempname)) {
      cout << "No prediction histogram found in this file and directory" << endl;
      return;
    }
    EstHist_input=(TH1D*) EstFile->Get(tempname)->Clone();
    string hsysUp = tempname;  hsysUp = regex_replace(hsysUp, regex("pred"), "sysUp");
    if (!EstFile->Get(hsysUp.c_str())) {
    // If (!EstFile->Get("ZinvBGsysUp")) {
      cout << "No prediction upper error histogram found in this file and directory" << endl;
      return;
    }
    EstSysUpHist_input=(TH1D*) EstFile->Get(hsysUp.c_str())->Clone();
    string hsysLow = tempname;  hsysLow = regex_replace(hsysLow, regex("pred"), "sysLow");
    if (!EstFile->Get(hsysLow.c_str())) {
    // if (!EstFile->Get("ZinvBGsysLow")) {
      cout << "No prediction lower error histogram found in this file and directory" << endl;
      return;
    }
    EstSysLowHist_input=(TH1D*) EstFile->Get(hsysLow.c_str())->Clone();
    if (!GenFile->Get(MCtempname)) {
      cout << "No MC histogram found in this file and directory" << endl;
      return;
    }
    GenHist_input=(TH1D*) GenFile->Get(MCtempname)->Clone();
    /*
    EstHistD_input=(TH1D*) EstFile->Get(tempname)->Clone();
    GenHistD_input=(TH1D*) GenFile->Get(tempname)->Clone();
    */
  }

  // GenHist_input->Print("all");
  //  TH1D * GenHist = static_cast<TH1D*>(GenHist_input->Clone("GenHist"));
  GenHist = static_cast<TH1D*>(GenHist_input->Clone("GenHist"));
  //  TH1D * EstHist = static_cast<TH1D*>(EstHist_input->Clone("EstHist"));
  EstHist = static_cast<TH1D*>(EstHist_input->Clone("EstHist"));
  EstSysUpHist = static_cast<TH1D*>(EstSysUpHist_input->Clone("EstSysUpHist"));
  EstSysLowHist = static_cast<TH1D*>(EstSysLowHist_input->Clone("EstSysLowHist"));
  float ls = 0.5;
  if (doBinShift) {
    shift_bin(GenHist_input,GenHist);
    shift_bin(EstHist_input,EstHist);
    ls = 0.0;
  }

  // GenHist->Print("all");
  // return;

  // TGraphAsymmErrors* ZinvBGsyst = new TGraphAsymmErrors(ZinvBGpred);
  TGraphAsymmErrors* ZinvBGsyst = new TGraphAsymmErrors(EstHist);
  ZinvBGsyst->SetMarkerSize(0);
  ZinvBGsyst->SetLineColorAlpha(kRed, 0);
  // ZinvBGsyst->SetLineWidth(-1);
  ZinvBGsyst->SetLineWidth(2);
  ZinvBGsyst->SetFillColor(kRed-9);
  ZinvBGsyst->SetFillStyle(1001);
  TGraphAsymmErrors* ZinvBGsystFrac = new TGraphAsymmErrors(*ZinvBGsyst);
  for (Int_t bin0=0; bin0<EstSysLowHist->GetNbinsX(); ++bin0) {
    if (EstHist->GetBinContent(bin0+1) > 0) {
      ZinvBGsyst->SetPointError(bin0, 0.5, 0.5, EstSysLowHist->GetBinContent(bin0+1), EstSysUpHist->GetBinContent(bin0+1));
      ZinvBGsystFrac->SetPointError(bin0, 0.5, 0.5, EstSysLowHist->GetBinContent(bin0+1)/EstHist->GetBinContent(bin0+1), EstSysUpHist->GetBinContent(bin0+1)/EstHist->GetBinContent(bin0+1));
      Double_t xf, yf;
      ZinvBGsystFrac->GetPoint(bin0, xf, yf);
      ZinvBGsystFrac->SetPoint(bin0, xf, 1.0);
    }
  }
  cout << "ZinvBGsystFrac" << endl;
  // ZinvBGsystFrac->Print("all");

  GenHist->SetLineColor(4);
  EstHist->SetLineColor(4);
  //GenHist->GetXaxis()->SetLabelFont(42);
  //GenHist->GetXaxis()->SetLabelOffset(0.007);
  //GenHist->GetXaxis()->SetLabelSize(0.04);
  //GenHist->GetXaxis()->SetTitleSize(0.05);
  //GenHist->GetXaxis()->SetTitleOffset(0.9);
  //GenHist->GetXaxis()->SetTitleOffset(0.5);
  //GenHist->GetXaxis()->SetTitleFont(42);
  //GenHist->GetYaxis()->SetLabelFont(42);
  //GenHist->GetYaxis()->SetLabelOffset(0.007);
  //GenHist->GetYaxis()->SetLabelSize(0.04);
  GenHist->GetYaxis()->SetLabelSize(0.045*1.15);
  GenHist->GetYaxis()->SetTitleSize(0.06*1.15);
  GenHist->GetYaxis()->SetTitleOffset(0.6);
  GenHist->GetYaxis()->SetTitleFont(42);


  //EstHist->GetXaxis()->SetLabelFont(42);
  //EstHist->GetXaxis()->SetLabelOffset(0.007);
  //EstHist->GetXaxis()->SetLabelSize(0.04);
  //EstHist->GetXaxis()->SetTitleSize(0.05);
  //EstHist->GetXaxis()->SetTitleOffset(0.9);
  //EstHist->GetXaxis()->SetTitleFont(42);
  //EstHist->GetYaxis()->SetLabelFont(42);
  //EstHist->GetYaxis()->SetLabelOffset(0.007);
  //EstHist->GetYaxis()->SetLabelSize(0.04);
  //EstHist->GetYaxis()->SetTitleSize(0.08);
  //EstHist->GetYaxis()->SetTitleOffset(2.0);
  //EstHist->GetYaxis()->SetTitleFont(42);
  char xtitlename[200];
  char ytitlename[200];
  sprintf(xtitlename,"Search region bin number");
  sprintf(ytitlename,"Events");
  gPad->SetLogy();
  GenHist->SetMaximum(ymax_top);
  GenHist->SetMinimum(ymin_top);
  GenHist->GetXaxis()->SetRangeUser(search_x_min,search_x_max);

  //GenHist->GetYaxis()->SetTickLength(0.015);
  //GenHist->GetXaxis()->SetTickLength(0.02);

  //gPad->SetGridx(1);
  TExec *ex1 = new TExec("ex1","gStyle->SetErrorX(0);");
  TExec *ex2 = new TExec("ex2","gStyle->SetErrorX(0.5);");

  GenHist->SetTitle("");
  GenHist->SetMarkerStyle(20);
  GenHist->SetMarkerSize(0.7);
  // GenHist->SetMarkerSize(1.2);
  GenHist->SetLineColor(1);
  GenHist->GetXaxis()->SetTitle(xtitlename);
  GenHist->GetYaxis()->SetTitle(ytitlename);
  GenHist->Scale(lumi/lumi_ref);
  cout << lumi << " / " << lumi_ref << endl;
  EstHist->Scale(lumi/lumi_ref);
  TH1D * GenHist_Normalize = static_cast<TH1D*>(GenHist->Clone("GenHist_Normalize"));
  //  GenHist_Normalize = static_cast<TH1D*>(GenHist->Clone("GenHist_Normalize"));
  GenHist_Normalize->SetMaximum(ymax_top);
  GenHist_Normalize->SetMinimum(ymin_top);
  ex1->Draw();
  //GenHist_Normalize->GetListOfFunctions()->Add(ex1);
  GenHist_Normalize->DrawCopy("e");

  // EstHist->SetFillStyle(3144);
  // EstHist->SetFillColor(kRed-10);
  EstHist->SetMarkerStyle(0);
  EstHist->SetMarkerSize(0.7);
  TH1D * EstHist_Normalize = static_cast<TH1D*>(EstHist->Clone("EstHist_Normalize"));
  //  EstHist_Normalize = static_cast<TH1D*>(EstHist->Clone("EstHist_Normalize"));
  ex2->Draw();
  // ex1->Draw();
  //EstHist_Normalize->GetListOfFunctions()->Add(ex2);
  ZinvBGsyst->Draw("2 same");
  // ZinvBGsyst->Draw("e2 same");
  EstHist_Normalize->DrawCopy("e2same");
  EstHist_Normalize->DrawCopy("esame");

  // GenHist->Print("all");
  // EstHist->Print("all");
  
  //
  // Re-draw to have "expectation" on top of "prediction"
  ex1->Draw();
  GenHist_Normalize->DrawCopy("esame");
  //

  TString line = "";
  if (Run2) sprintf(tempname,"%8.0f",lumi);
  else sprintf(tempname,"%8.1f",lumi);
  line+=tempname;
  line+=" fb^{-1} (13 TeV)";
  
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=0;
    
  // writeExtraText = true;
  writeExtraText = false;
  extraText   = "      Simulation";
  //float extraTextFont = 52;  // default is helvetica-italics

  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  //lumiTextSize     = 0.5;
  //float lumiTextOffset   = 0.2;
  //cmsTextSize      = 0.65;
  //float cmsTextOffset    = 0.1;  // only used in outOfFrame version
  
  //relPosX    = 0.045;
  //relPosY    = 0.035;
  //relExtraDY = 1.2;
  
  // ratio of "CMS" and extra text size
  //float extraOverCmsTextSize  = 0.76;
    
  //TString lumi_13TeV = "20.1 fb^{-1}";
  //TString lumi_8TeV  = "19.7 fb^{-1}";
  //TString lumi_7TeV  = "5.1 fb^{-1}";
  TString lumi_sqrtS = line;

  //
  if(histname.find("QCD")==string::npos ){
    // Signal binning
    
    //-----------------------------------------------------------
    // Putting lines and labels explaining search region definitions
    //-----------------------------------------------------------

    // TString CMSlabel = "";
    // cmsText = "#bf{CMS} #it{Simulation}";
    // CMSlabel += "#splitline{#bf{CMS}}{#scale[0.6]{#it{Simulation}}}";

    
    double x0 = gStyle->GetPadLeftMargin();
    double x1 = 1.-gStyle->GetPadRightMargin();
    // double y0 = 1.005-gStyle->GetPadTopMargin();
    // double y1 = 0.96;
    double y0 = 0.930-gStyle->GetPadTopMargin();
    // y1 = y0 + 0.045;
    double y1 = y0 + 0.045;
    TPaveText *Lumitxt = new TPaveText(x0,y0,x1,y1,"NDC");
    Lumitxt->SetBorderSize(0);
    Lumitxt->SetFillColor(0);
    Lumitxt->SetTextFont(42);
    Lumitxt->SetTextAlign(31);
    // Lumitxt->SetTextSize(1.7*gStyle->GetPadTopMargin());
    Lumitxt->SetTextSize(1.5*gStyle->GetPadTopMargin());
    Lumitxt->SetMargin(0.);
    Lumitxt->AddText(line);
    Lumitxt->Draw("same");

    x0 = gStyle->GetPadLeftMargin()-0.062;
    x1 = x0 + 0.10;
    // x0 = gStyle->GetPadLeftMargin()+0.03;
    // x1 = gStyle->GetPadLeftMargin()+0.13;
    // y0 = 0.905-gStyle->GetPadTopMargin();
    // y1 = 0.88;
    TPaveText *CMStxt = new TPaveText(x0,y0,x1,y1,"NDC");
    CMStxt->SetBorderSize(0);
    CMStxt->SetFillColor(0);
    CMStxt->SetTextFont(42);
    CMStxt->SetTextAlign(11);
    CMStxt->SetTextSize(1.9*gStyle->GetPadTopMargin());
    // CMStxt->SetTextSize(2.0*gStyle->GetPadTopMargin());
    CMStxt->SetMargin(0.);
    CMStxt->AddText(CMSlabel);
    CMStxt->Draw("same");

    // Njet separation lines
    TLine *tl_njet = new TLine();
    tl_njet->SetLineStyle(2);
    tl_njet->DrawLine(30.+ls,ymin_top,30.+ls,ymax_top);
    tl_njet->DrawLine(70.+ls,ymin_top,70.+ls,ymax_top); 
    tl_njet->DrawLine(110.+ls,ymin_top,110.+ls,ymax_top);
    tl_njet->DrawLine(142.+ls,ymin_top,142.+ls,ymax_top);
    // tl_njet->DrawLine(40.,ymin_top,40.,ymax_top);
    // tl_njet->DrawLine(80.,ymin_top,80.,ymax_top); 
    // tl_njet->DrawLine(120.,ymin_top,120.,ymax_top);
    // tl_njet->DrawLine(25.-0.5,ymin_top,25.-0.5,ymax_top); 
    // tl_njet->DrawLine(49.-0.5,ymin_top,49.-0.5,ymax_top); 

    // Njet labels
    TLatex * ttext_njet = new TLatex();
    ttext_njet->SetTextFont(42);
    // ttext_njet->SetTextSize(0.060);
    ttext_njet->SetTextSize(0.040);
    ttext_njet->SetTextAlign(22);
    ttext_njet->DrawLatex(15. , ymax_top/4. , "2 #leq N_{#scale[0.2]{ }jet} #leq 3");
    ttext_njet->DrawLatex(50. , ymax_top/4. , "4 #leq N_{#scale[0.2]{ }jet} #leq 5");
    ttext_njet->DrawLatex(90. , ymax_top/4. , "6 #leq N_{#scale[0.2]{ }jet} #leq 7");
    ttext_njet->DrawLatex(126. , ymax_top/4. , "8 #leq N_{#scale[0.2]{ }jet} #leq 9");
    ttext_njet->DrawLatex(158. , ymax_top/4. , "N_{#scale[0.2]{ }jet} #geq 10");
    // ttext_njet->DrawLatex(15. , ymax_top/4. , "N_{#scale[0.2]{ }jet} = 2");
    // ttext_njet->DrawLatex(50. , ymax_top/4. , "3 #leq N_{#scale[0.2]{ }jet} #leq 4");
    // ttext_njet->DrawLatex(90. , ymax_top/4. , "5 #leq N_{#scale[0.2]{ }jet} #leq 6");
    // ttext_njet->DrawLatex(126. , ymax_top/4. , "7 #leq N_{#scale[0.2]{ }jet} #leq 8");
    // ttext_njet->DrawLatex(158. , ymax_top/4. , "N_{#scale[0.2]{ }jet} #geq 9");
    // ttext_njet->DrawLatex(13.-0.5 , ymax_top/4. , "4 #leq N_{#scale[0.2]{ }jet} #leq 6");
    // ttext_njet->DrawLatex(37.-0.5 , ymax_top/4. , "7 #leq N_{#scale[0.2]{ }jet} #leq 8");
    // ttext_njet->DrawLatex(61.-0.5 , ymax_top/4. , "N_{#scale[0.2]{ }jet} #geq 9");

    // Nb separation lines
    TLine *tl_nb = new TLine();
    tl_nb->SetLineStyle(3);
    tl_nb->DrawLine(10.+ls,ymin_top,10.+ls,ymax2_top); 
    tl_nb->DrawLine(20.+ls,ymin_top,20.+ls,ymax2_top); 

    tl_nb->DrawLine(40.+ls,ymin_top,40.+ls,ymax2_top); 
    tl_nb->DrawLine(50.+ls,ymin_top,50.+ls,ymax2_top); 
    tl_nb->DrawLine(60.+ls,ymin_top,60.+ls,ymax2_top); 

    tl_nb->DrawLine(80.+ls,ymin_top,80.+ls,ymax2_top); 
    tl_nb->DrawLine(90.+ls,ymin_top,90.+ls,ymax2_top); 
    tl_nb->DrawLine(100.+ls,ymin_top,100.+ls,ymax2_top); 

    tl_nb->DrawLine(118.+ls,ymin_top,118.+ls,ymax2_top); 
    tl_nb->DrawLine(126.+ls,ymin_top,126.+ls,ymax2_top); 
    tl_nb->DrawLine(134.+ls,ymin_top,134.+ls,ymax2_top); 

    tl_nb->DrawLine(150.+ls,ymin_top,150.+ls,ymax2_top); 
    tl_nb->DrawLine(158.+ls,ymin_top,158.+ls,ymax2_top); 
    tl_nb->DrawLine(166.+ls,ymin_top,166.+ls,ymax2_top); 

    // tl_nb->DrawLine(10.,ymin_top,10.,ymax2_top); 
    // tl_nb->DrawLine(20.,ymin_top,20.,ymax2_top); 
    // tl_nb->DrawLine(30.,ymin_top,30.,ymax2_top); 

    // tl_nb->DrawLine(50.,ymin_top,50.,ymax2_top); 
    // tl_nb->DrawLine(60.,ymin_top,60.,ymax2_top); 
    // tl_nb->DrawLine(70.,ymin_top,70.,ymax2_top); 

    // tl_nb->DrawLine(90.,ymin_top,90.,ymax2_top); 
    // tl_nb->DrawLine(100.,ymin_top,100.,ymax2_top); 
    // tl_nb->DrawLine(110.,ymin_top,110.,ymax2_top); 

    // tl_nb->DrawLine(130.,ymin_top,130.,ymax2_top); 
    // tl_nb->DrawLine(140.,ymin_top,140.,ymax2_top); 
    // tl_nb->DrawLine(150.,ymin_top,150.,ymax2_top); 
    // tl_nb->DrawLine( 7.-0.5,ymin_top, 7.-0.5,ymax2_top); 
    // tl_nb->DrawLine(13.-0.5,ymin_top,13.-0.5,ymax2_top); 
    // tl_nb->DrawLine(19.-0.5,ymin_top,19.-0.5,ymax2_top); 
    // tl_nb->DrawLine(31.-0.5,ymin_top,31.-0.5,ymax3_top); 
    // tl_nb->DrawLine(37.-0.5,ymin_top,37.-0.5,ymax3_top); 
    // tl_nb->DrawLine(43.-0.5,ymin_top,43.-0.5,ymax3_top); 
    // tl_nb->DrawLine(55.-0.5,ymin_top,55.-0.5,ymax4_top); 
    // tl_nb->DrawLine(61.-0.5,ymin_top,61.-0.5,ymax4_top); 
    // tl_nb->DrawLine(67.-0.5,ymin_top,67.-0.5,ymax4_top); 
    
    // Nb labels
    TLatex * ttext_nb = new TLatex();
    ttext_nb->SetTextFont(42);
    // ttext_nb->SetTextSize(0.060);
    ttext_nb->SetTextSize(0.040);
    ttext_nb->SetTextAlign(22);
    
    ttext_nb->DrawLatex( 9.+0., ymax_top/12. , "N_{#scale[0.2]{ }b-jet}");
    ttext_nb->DrawLatex( 6.+0., ymax_top/40. , "0");
    ttext_nb->DrawLatex(15.+0., ymax_top/40. , "1");
    ttext_nb->DrawLatex(25.+0., ymax_top/40. , "#geq 2");
    ttext_nb->DrawLatex(5.+30., ymax_top/40. , "0");
    ttext_nb->DrawLatex(15.+30., ymax_top/40. , "1");
    ttext_nb->DrawLatex(25.+30., ymax_top/40. , "2");
    ttext_nb->DrawLatex(35.+30., ymax_top/40. , "#geq 3");

    // ttext_nb->DrawLatex( 9.+40., ymax_top/12. , "N_{#scale[0.2]{ }b-jet}");
    // ttext_nb->DrawLatex( 6.+40., ymax_top/40. , "0");
    // ttext_nb->DrawLatex(16.+40., ymax_top/40. , "1");
    // ttext_nb->DrawLatex(26.+40., ymax_top/40. , "2");
    // ttext_nb->DrawLatex(36.+40., ymax_top/40. , "#geq 3");

    //
  } else {
    // LDP, HDP binning
    
    //-----------------------------------------------------------
    // Putting lines and labels explaining search region definitions
    //-----------------------------------------------------------

    // Njet separation lines
    TLine *tl_njet = new TLine();
    tl_njet->SetLineStyle(2);
    cout << "setting Njet separaton lines" << endl;
    tl_njet->DrawLine(53.-0.5,ymin_top,53.-0.5,ymax_top); 
    tl_njet->DrawLine(105.-0.5,ymin_top,105.-0.5,ymax_top); 
    tl_njet->DrawLine(157.-0.5,ymin_top,157.-0.5,ymax_top); 
    // tl_njet->DrawLine( 45.,ymin_top, 45.,ymax_top); 
    // tl_njet->DrawLine( 89.,ymin_top, 89.,ymax_top); 
    // tl_njet->DrawLine(133.,ymin_top,133.,ymax_top); 
    // tl_njet->DrawLine(177.,ymin_top,177.,ymax_top); 

    // Njet labels
    TLatex * ttext_njet = new TLatex();
    ttext_njet->SetTextFont(42);
    ttext_njet->SetTextSize(0.04);
    ttext_njet->SetTextAlign(22);
    ttext_njet->DrawLatex(26.-0.5 , ymax_top/4. , "3 #leq N_{#scale[0.2]{ }jet} #leq 4");
    ttext_njet->DrawLatex(78.-0.5 , ymax_top/4. , "5 #leq N_{#scale[0.2]{ }jet} #leq 6");
    ttext_njet->DrawLatex(130.-0.5 , ymax_top/4. , "7 #leq N_{#scale[0.2]{ }jet} #leq 8");
    ttext_njet->DrawLatex(182.-0.5 , ymax_top/4. , "N_{#scale[0.2]{ }jet} #geq 9");
    // ttext_njet->DrawLatex(23. , ymax_top/4. , "N_{jets} = 4");
    // ttext_njet->DrawLatex(67. , ymax_top/4. , "N_{jets} = 5");
    // ttext_njet->DrawLatex(111., ymax_top/4. , "N_{jets} = 6");
    // ttext_njet->DrawLatex(155., ymax_top/4. , "7 #leq N_{jets} #leq 8");
    // ttext_njet->DrawLatex(199., ymax_top/4. , "N_{jets} #geq 9");

    // Nb separation lines
    TLine *tl_nb = new TLine();
    tl_nb->SetLineStyle(3);
    tl_nb->DrawLine(12.,ymin_top,12.,ymax2_top); 
    tl_nb->DrawLine(23.,ymin_top,23.,ymax2_top); 
    tl_nb->DrawLine(34.,ymin_top,34.,ymax2_top); 

    tl_nb->DrawLine(56.,ymin_top,56.,ymax2_top); 
    tl_nb->DrawLine(67.,ymin_top,67.,ymax2_top); 
    tl_nb->DrawLine(78.,ymin_top,78.,ymax2_top); 

    tl_nb->DrawLine(100.,ymin_top,100.,ymax2_top); 
    tl_nb->DrawLine(111.,ymin_top,111.,ymax2_top); 
    tl_nb->DrawLine(122.,ymin_top,122.,ymax2_top); 

    tl_nb->DrawLine(144.,ymin_top,144.,ymax3_top); 
    tl_nb->DrawLine(155.,ymin_top,155.,ymax3_top); 
    tl_nb->DrawLine(166.,ymin_top,166.,ymax3_top); 

    tl_nb->DrawLine(188.,ymin_top,188.,ymax4_top); 
    tl_nb->DrawLine(199.,ymin_top,199.,ymax4_top); 
    tl_nb->DrawLine(210.,ymin_top,210.,ymax4_top); 

    // Nb labels
    TLatex * ttext_nb = new TLatex();
    ttext_nb->SetTextFont(42);
    ttext_nb->SetTextSize(0.04);
    ttext_nb->SetTextAlign(22);
    ttext_nb->SetTextAngle(90);

    ttext_nb->DrawLatex( 6. , ymax_top/50. , "N_{b} = 0");
    ttext_nb->DrawLatex(17. , ymax_top/50. , "N_{b} = 1");
    ttext_nb->DrawLatex(28. , ymax_top/50. , "N_{b} = 2");
    ttext_nb->DrawLatex(39. , ymax_top/50. , "N_{b} #geq 3");
    
    TText * ttext = new TLatex(160. , ymax_top/50. , "Normalized to 10 fb^{-1}");
    ttext->SetTextFont(42);
    ttext->SetTextSize(0.045);
    ttext->SetTextAlign(22);
    ttext->Draw();

  }

  // Legend & texts
  sprintf(tempname, "%s", legendHeader.Data());
  catLeg1->SetHeader(tempname);
  sprintf(tempname, "%s", expectLegend.Data());
  catLeg1->AddEntry(GenHist,tempname,"ep");
  sprintf(tempname, "%s", predLegend.Data());
  EstHist->SetFillColor(kRed-9);
  EstHist->SetLineColor(kBlue);
  EstHist->SetFillStyle(1001);
  catLeg1->AddEntry(EstHist,tempname);
  catLeg1->Draw();

  gPad->RedrawAxis();

  //
  // Bottom ratio plot
  //
  // ----------

    //
    // Preparing ratio histograms
      TH1D * numerator   = static_cast<TH1D*>(GenHist->Clone("numerator"));
      TH1D * numerator_fullstaterr   = static_cast<TH1D*>(GenHist->Clone("numerator_fullstaterr"));
      TH1D * denominator = static_cast<TH1D*>(EstHist->Clone("denominator"));

      TH1D * GenHist_Clone = static_cast<TH1D*>(GenHist->Clone("GenHist_Clone"));
      TH1D * EstHist_Clone = static_cast<TH1D*>(EstHist->Clone("EstHist_Clone"));
      TH1D * EstHist_NoError = static_cast<TH1D*>(EstHist->Clone("EstHist_NoError"));
      TH1D * One_NoError = static_cast<TH1D*>(EstHist->Clone("One_NoError"));
      for (int ibin=1; ibin<EstHist_NoError->GetNbinsX(); ibin++){ // scan excluding underflow and overflow bins
	EstHist_NoError->SetBinError(ibin,0.);
	One_NoError->SetBinContent(ibin,1.);
	One_NoError->SetBinError(ibin,0.);
      }

      //EstHistD->Add(GenHistD,-1);
      numerator->Divide(GenHist_Clone,EstHist_NoError,1,1,"");
      denominator->Divide(EstHist_Clone,EstHist_NoError,1,1,"");

      numerator_fullstaterr->Divide(GenHist_Clone,EstHist_Clone,1,1,"");  // Expectation/Prediction
      numerator_fullstaterr->Add(One_NoError,-1.);                        // Expectation/Prediction-1

      // draw bottom figure
      canvas_dw->cd();
      numerator->SetOption("ex0");  // wtf
      // font size
      numerator->GetXaxis()->SetLabelSize(font_size_dw);
      numerator->GetXaxis()->SetTitleSize(font_size_dw);
      numerator->GetYaxis()->SetLabelSize(font_size_dw);
      numerator->GetYaxis()->SetTitleSize(font_size_dw);

      //
      // Horizontal Lines
      TLine *tline  = new TLine(search_x_min,1.,search_x_max,1.);
      TLine *tline0 = new TLine(search_x_min,0.,search_x_max,0.);

      //
      // Common to all bottom plots
      //
      //sprintf(ytitlename,"#frac{Estimate - #tau_{had} BG}{#tau_{had} BG} ");
      sprintf(ytitlename, "%s", ratioOrdinateTitle.Data());
      numerator->SetMaximum(ymax_bottom);
      numerator->SetMinimum(ymin_bottom);

      //
      // Specific to each bottom plot
      //
      // Setting style
      //numerator->SetMaximum(1.4);
      //numerator->GetXaxis()->SetLabelFont(42);
      //numerator->GetXaxis()->SetLabelOffset(0.007);
      numerator->GetXaxis()->SetLabelSize(0.18*0.045/0.06);
      numerator->GetXaxis()->SetTitleSize(0.18);
      numerator->GetXaxis()->SetTitleOffset(0.9);
      numerator->GetXaxis()->SetTitleFont(42);
      //numerator->GetYaxis()->SetLabelFont(42);
      //numerator->GetYaxis()->SetLabelOffset(0.007);
      numerator->GetYaxis()->SetLabelSize(0.18*0.045/0.06);
      numerator->GetYaxis()->SetTitleSize(0.18);
      //numerator->GetYaxis()->SetTitleOffset(0.5);
      numerator->GetYaxis()->SetTitleOffset(0.25);
      numerator->GetYaxis()->SetTitleFont(42);

      numerator->GetXaxis()->SetTitle(xtitlename);
      numerator->GetYaxis()->SetTitle(ytitlename);

      //gPad->SetGridx(1);


      if (pull==1){

	sprintf(ytitlename,"#frac{Exp - Pre}{Stat Error} ");
	numerator->SetMaximum(8.);
	numerator->SetMinimum(-8.);
	
	//
	// Specific to each bottom plot
	//
	// Setting style

	for (int ibin=0; ibin<numerator_fullstaterr->GetNbinsX()+2; ibin++){ // scan including underflow and overflow bins
	  numerator_fullstaterr->SetBinContent(ibin,numerator_fullstaterr->GetBinContent(ibin)/numerator_fullstaterr->GetBinError(ibin));
	  numerator_fullstaterr->SetBinError(ibin,0.);
	}
	numerator_fullstaterr->Print("all");
	
	numerator_fullstaterr->GetXaxis()->SetLabelSize(font_size_dw);
	numerator_fullstaterr->GetXaxis()->SetTitleSize(font_size_dw);
	numerator_fullstaterr->GetYaxis()->SetLabelSize(font_size_dw);
	numerator_fullstaterr->GetYaxis()->SetTitleSize(font_size_dw);

	numerator_fullstaterr->GetXaxis()->SetTitleSize(0.12);
	numerator_fullstaterr->GetXaxis()->SetTitleOffset(0.9);
	numerator_fullstaterr->GetXaxis()->SetTitleFont(42);
	numerator_fullstaterr->GetYaxis()->SetTitleSize(0.13);
	numerator_fullstaterr->GetYaxis()->SetTitleOffset(0.5);
	numerator_fullstaterr->GetYaxis()->SetTitleFont(42);
	
	numerator_fullstaterr->GetXaxis()->SetTitle(xtitlename);
	numerator_fullstaterr->GetYaxis()->SetTitle(ytitlename);
	//numerator_fullstaterr->SetFillColor(kGreen-3);
	numerator_fullstaterr->SetFillColor(kRed-10);
	numerator_fullstaterr->DrawCopy();

	//
	// Drawing lines
	tline0->SetLineStyle(2);
	//tline0->Draw();

      }
      else {

      //
      // Plotting
      tline->SetLineStyle(1);
      tline->SetLineColor(kBlue);
      tline->Draw();
      numerator->GetYaxis()->SetNdivisions(505);
      numerator->GetYaxis()->SetTickLength(0.015);
      numerator->GetXaxis()->SetTickLength(0.08);
      numerator->SetTitle("");
      ex1->Draw();
      numerator->DrawCopy();

      ex1->Draw();
      ZinvBGsystFrac->Draw("2 same");
      denominator->DrawCopy("e2same");
      denominator->DrawCopy("same");
      ex1->Draw();
      tline->SetLineStyle(1);
      tline->SetLineColor(kBlue);
      tline->Draw();
      // numerator->DrawCopy("same");
      numerator->DrawCopy("e0same");  // Show errors for points outside plot range

      // numerator->Print("all");
      // denominator->Print("all");
      // numerator_fullstaterr->Print("all");

      }
      
      //
      // Drawing lines

      //
      if(histname.find("QCD")==string::npos ){
	// Signal binning

        // Njet separation lines
	TLine *tl_njet = new TLine();
	tl_njet->SetLineStyle(2);
        tl_njet->DrawLine(30.+ls, ymin_bottom,30.+ls, ymax_bottom); 
        tl_njet->DrawLine(70.+ls, ymin_bottom,70.+ls, ymax_bottom); 
        tl_njet->DrawLine(110.+ls, ymin_bottom,110.+ls, ymax_bottom); 
        tl_njet->DrawLine(142.+ls, ymin_bottom,142.+ls, ymax_bottom); 
        // tl_njet->DrawLine(40.+ls, ymin_bottom,40.+ls, ymax_bottom); 
        // tl_njet->DrawLine(80.+ls, ymin_bottom,80.+ls, ymax_bottom); 
        // tl_njet->DrawLine(120.+ls, ymin_bottom,120.+ls, ymax_bottom); 
        // tl_njet->DrawLine( 25.-0.5,ymin_bottom, 25.-0.5,ymax_bottom); 
        // tl_njet->DrawLine( 49.-0.5,ymin_bottom, 49.-0.5,ymax_bottom); 

        // Nb separation lines
	TLine *tl_nb = new TLine();
	tl_nb->SetLineStyle(3);
	tl_nb->DrawLine(10.+ls, ymin_bottom,10.+ls, ymax2_bottom); 
        tl_nb->DrawLine(20.+ls, ymin_bottom,20.+ls, ymax2_bottom); 	

    	tl_nb->DrawLine(40.+ls, ymin_bottom,40.+ls, ymax2_bottom); 	
    	tl_nb->DrawLine(50.+ls, ymin_bottom,50.+ls, ymax2_bottom); 	
    	tl_nb->DrawLine(60.+ls, ymin_bottom,60.+ls, ymax2_bottom); 	

    	tl_nb->DrawLine(80.+ls, ymin_bottom,80.+ls, ymax2_bottom); 	
     	tl_nb->DrawLine(90.+ls, ymin_bottom,90.+ls, ymax2_bottom);
     	tl_nb->DrawLine(100.+ls, ymin_bottom,100.+ls, ymax2_bottom);

     	tl_nb->DrawLine(118.+ls, ymin_bottom,118.+ls, ymax2_bottom);
     	tl_nb->DrawLine(126.+ls, ymin_bottom,126.+ls, ymax2_bottom);
     	tl_nb->DrawLine(134.+ls, ymin_bottom,134.+ls, ymax2_bottom);

     	tl_nb->DrawLine(150.+ls, ymin_bottom,150.+ls, ymax2_bottom);
     	tl_nb->DrawLine(158.+ls, ymin_bottom,158.+ls, ymax2_bottom);
     	tl_nb->DrawLine(166.+ls, ymin_bottom,166.+ls, ymax2_bottom);

	// tl_nb->DrawLine(10.,ymin_bottom,10.,ymax2_bottom); 
        // tl_nb->DrawLine(20.,ymin_bottom,20.,ymax2_bottom); 	
    	// tl_nb->DrawLine(30.,ymin_bottom,30.,ymax2_bottom); 	

    	// tl_nb->DrawLine(50.,ymin_bottom,50.,ymax2_bottom); 	
    	// tl_nb->DrawLine(60.,ymin_bottom,60.,ymax2_bottom); 	
    	// tl_nb->DrawLine(70.,ymin_bottom,70.,ymax2_bottom); 	

    	// tl_nb->DrawLine(90.,ymin_bottom,90.,ymax2_bottom); 	
     	// tl_nb->DrawLine(100.,ymin_bottom,100.,ymax2_bottom);
     	// tl_nb->DrawLine(110.,ymin_bottom,110.,ymax2_bottom);

     	// tl_nb->DrawLine(130.,ymin_bottom,130.,ymax2_bottom);
     	// tl_nb->DrawLine(140.,ymin_bottom,140.,ymax2_bottom);
     	// tl_nb->DrawLine(150.,ymin_bottom,150.,ymax2_bottom);
      // tl_nb->DrawLine( 7.-0.5,ymin_bottom, 7.-0.5,ymax2_bottom); 
      // tl_nb->DrawLine(13.-0.5,ymin_bottom,13.-0.5,ymax2_bottom); 
      // tl_nb->DrawLine(19.-0.5,ymin_bottom,19.-0.5,ymax2_bottom); 
      
      // tl_nb->DrawLine(31.-0.5,ymin_bottom,31.-0.5,ymax2_bottom); 
      // tl_nb->DrawLine(37.-0.5,ymin_bottom,37.-0.5,ymax2_bottom); 
      // tl_nb->DrawLine(43.-0.5,ymin_bottom,43.-0.5,ymax2_bottom); 
      
      // tl_nb->DrawLine(55.-0.5,ymin_bottom,55.-0.5,ymax2_bottom); 
      // tl_nb->DrawLine(61.-0.5,ymin_bottom,61.-0.5,ymax2_bottom); 
      // tl_nb->DrawLine(67.-0.5,ymin_bottom,67.-0.5,ymax2_bottom); 

      } else {
	// LDP, HDP binning
	
      // Njet separation lines
      TLine *tl_njet = new TLine();
      tl_njet->SetLineStyle(2);
      tl_njet->DrawLine(53.-0.5,ymin_bottom,53.-0.5,ymax_bottom); 
      tl_njet->DrawLine(105.-0.5,ymin_bottom,105.-0.5,ymax_bottom); 
      tl_njet->DrawLine(157.-0.5,ymin_bottom,157.-0.5,ymax_bottom); 
      // tl_njet->DrawLine( 45.,ymin_bottom, 45.,ymax_bottom); 
      // tl_njet->DrawLine( 89.,ymin_bottom, 89.,ymax_bottom); 
      // tl_njet->DrawLine(133.,ymin_bottom,133.,ymax_bottom); 
      // tl_njet->DrawLine(177.,ymin_bottom,177.,ymax_bottom); 


      // Nb separation lines
      TLine *tl_nb = new TLine();
      tl_nb->SetLineStyle(3);
      tl_nb->DrawLine(14.-0.5,ymin_bottom,14.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(27.-0.5,ymin_bottom,27.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(40.-0.5,ymin_bottom,40.-0.5,ymax2_bottom);
      tl_nb->DrawLine(53.-0.5,ymin_bottom,53.-0.5,ymax2_bottom);
      tl_nb->DrawLine(66.-0.5,ymin_bottom,66.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(79.-0.5,ymin_bottom,79.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(92.-0.5,ymin_bottom,92.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(105.-0.5,ymin_bottom,105.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(118.-0.5,ymin_bottom,118.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(131.-0.5,ymin_bottom,131.-0.5,ymax2_bottom); 
      tl_nb->DrawLine(144.-0.5,ymin_bottom,144.-0.5,ymax2_bottom);
      tl_nb->DrawLine(157.-0.5,ymin_bottom,157.-0.5,ymax2_bottom);
      tl_nb->DrawLine(170.-0.5,ymin_bottom,170.-0.5,ymax2_bottom);
      tl_nb->DrawLine(183.-0.5,ymin_bottom,183.-0.5,ymax2_bottom);
      tl_nb->DrawLine(196.-0.5,ymin_bottom,196.-0.5,ymax2_bottom);
     // tl_nb->DrawLine(12.,ymin_bottom,12.,ymax2_bottom); 
      // tl_nb->DrawLine(23.,ymin_bottom,23.,ymax2_bottom); 
      // tl_nb->DrawLine(34.,ymin_bottom,34.,ymax2_bottom); 
      
      // tl_nb->DrawLine(56.,ymin_bottom,56.,ymax2_bottom); 
      // tl_nb->DrawLine(67.,ymin_bottom,67.,ymax2_bottom); 
      // tl_nb->DrawLine(78.,ymin_bottom,78.,ymax2_bottom); 
      
      // tl_nb->DrawLine(100.,ymin_bottom,100.,ymax2_bottom); 
      // tl_nb->DrawLine(111.,ymin_bottom,111.,ymax2_bottom); 
      // tl_nb->DrawLine(122.,ymin_bottom,122.,ymax2_bottom); 

      // tl_nb->DrawLine(144.,ymin_bottom,144.,ymax2_bottom); 
      // tl_nb->DrawLine(155.,ymin_bottom,155.,ymax2_bottom); 
      // tl_nb->DrawLine(166.,ymin_bottom,166.,ymax2_bottom); 
      
      // tl_nb->DrawLine(188.,ymin_bottom,188.,ymax2_bottom); 
      // tl_nb->DrawLine(199.,ymin_bottom,199.,ymax2_bottom); 
      // tl_nb->DrawLine(210.,ymin_bottom,210.,ymax2_bottom); 
    
      }

      gPad->RedrawAxis();

      //
      //

  // CMS_lumi( canvas, iPeriod, iPos );

  // sprintf(tempname,"Closure_%s_%s_Full_%s%sPlot.png",histname.c_str(),sample.c_str(),elog.c_str(),elogExp.c_str());
  // if (pull==1) 
  //   sprintf(tempname,"ClosurePull_%s_%s_Full_%s%sPlot.png",histname.c_str(),sample.c_str(),elog.c_str(),elogExp.c_str());
  // canvas->Print(tempname);

  sprintf(tempname,"Closure_%s_%s%s%s.pdf",histname.c_str(),sample.c_str(),elog.c_str(),elogExp.c_str());
  if (pull==1)
    sprintf(tempname,"ClosurePull_%s_%s_Full_%s%sPlot.pdf",histname.c_str(),sample.c_str(),elog.c_str(),elogExp.c_str());
  canvas->Print(tempname);

}

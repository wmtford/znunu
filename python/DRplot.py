#!/usr/bin/python
"""
getDoubleRatioDatFiles.py generates the double ratio
dat files for integration. Can be run from
command line, interactively or with the driver script:
scripts/extrapDatFileDriver.csh
"""
import sys
import ROOT
import RA2b

ROOT.gROOT.SetBatch(1)

## doSample is a list of dat files to generate ##
doSample = []
## Check if samples are given at runtime (e.g. sig hdp ldp)
## can give more than one and will run over each sequentially
if len(sys.argv) > 1:
    for i in range(1,len(sys.argv)):
        doSample.append(sys.argv[i])
else:
    ## default runs over all 
    # doSample = ['sig','hdp','ldp']
    doSample = ['sig']
#################################################

haveHistograms = True
runYear = "Run2"
postDRfit = True
plotRefLine = True
Ymin = 0.5
Ymax = 1.5
extraText = '    Supplementary'

DYkfactor = 1.23
yTitle = "R^{data}_{Z#rightarrow l^{+}l^{-}/#gamma}/R^{sim}_{Z#rightarrow l^{+}l^{-}/#gamma}"

for sample in doSample:

    histoNJets = {}
    histoHT = {}
    histoMHT = {}

    if (not haveHistograms):
        histoNJets = None
        histoHT = None
        histoMHT = None
        ## no 9+ jets events in ldp sample
        if(sample=='ldp'):
            nj_binning = [1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]
        else:
            nj_binning = [1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]

        ## include sideband if ldp or hdp
            if(sample=='ldp' or sample=='hdp'):
            #mht_binning = [250.,300.,400.,600.,900.]
                mht_binning = [250.,900.]
                mht_binning = [300.,350.,400.,450.,600.,750,900.]
                mhtCut = True
            else:
                mhtCut = True
                mht_binning = [300.,350.,400.,450.,600.,750,900.]

        ## get the double ratio graphs
        nj_dr = RA2b.getDoubleRatioGraph('NJets',applyPuWeight=True,dphiCut=sample,binning=nj_binning,applyMHTCut=mhtCut)
        mht_dr = RA2b.getDoubleRatioGraph('MHT',applyPuWeight=True,dphiCut=sample,binning=mht_binning,applyMHTCut=mhtCut)
        ht_dr = RA2b.getDoubleRatioGraph('HT',applyPuWeight=True,dphiCut=sample,applyMHTCut=mhtCut)
        
        ## get the double ratio plots with values and uncertainties
        dr_out = RA2b.getDoubleRatioPlot([nj_dr,mht_dr,ht_dr])

    else:
        # Histograms from RA2bZinvAnalysis
        mcLumiRatio_mm = 1
        mcLumiRatio_ee = 1
        mcLumiRatio_photon = 1
        if (runYear is "2016"):
            iPeriod = 5
            dataZllFile = ROOT.TFile('../outputs/histsDY_2016v16.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_2016v16.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMC_2016v16_noPU.root')
            if (postDRfit):
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2016v16_noPU_DR2016wt.root')
            else:
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2016v16_noPU.root')
        elif (runYear is "2017"):
            iPeriod = 6
            dataZllFile = ROOT.TFile('../outputs/histsDY_2017v16.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_2017v16.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMC_2017v16_HT17wt_ZptWt.root')
            if (postDRfit):
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2017v16_DR2017wt.root')
            else:
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2017v16.root')
        elif (runYear is "2018AB"):
            iPeriod = 8
            dataZllFile = ROOT.TFile('../outputs/histsDY_2018ABv16.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_2018ABv16.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMC_2018v16_HT17wt_ZptWt.root')
            if (postDRfit):
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2018v16_DR2017wt.root')
            else:
                 mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2018v16.root')
            mcLumiRatio_mm = 21.0/59.6
            mcLumiRatio_ee =  21.0/59.6
            mcLumiRatio_photon =  21.0/59.6
        elif (runYear is "2018CD"):
            iPeriod = 8
            dataZllFile = ROOT.TFile('../outputs/histsDY_2018CDv16.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_2018CDv16.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMC_2018HEMv16_HT17wt_ZptWt.root')
            if (postDRfit):
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2018HEMv16_DR2017wt.root')
            else:
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2018HEMv16.root')
            mcLumiRatio_mm = 38.6/59.6
            mcLumiRatio_ee =  38.6/59.6
            mcLumiRatio_photon =  38.6/59.6
        elif (runYear is "Run2"):
            iPeriod = 9
            dataZllFile = ROOT.TFile('../outputs/histsDY_Run2v17.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_Run2v17.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMC_Run2v17.root')
            if (postDRfit):
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_Run2v17_DRr2wt.root')
            else:
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_Run2v17.root')
        elif (runYear is "Run2ldpnominal"):
            iPeriod = 9
            dataZllFile = ROOT.TFile('../outputs/histsDYldpnominal_Run2v161617.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhotonldpnominal_Run2v161617.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMCldpnominal_Run2v161617.root')
            if (postDRfit):
                mcPhotonFile = ROOT.TFile('../outputs/histsGjetsldpnominal_Run2v161617_DRr2wt.root')
            else:
                mcPhotonFile = ROOT.TFile('../outputs/histsGjetsldpnominal_Run2v161617.root')
        histoNJets['pho_da'] = dataPhotonFile.Get("hNJets_DR_photon")
        histoNJets['pho_cg'] = histoNJets['pho_da'].Clone()
        for i in range(1, histoNJets['pho_cg'].GetNbinsX()+1):
            histoNJets['pho_cg'].SetBinContent(i, histoNJets['pho_cg'].GetBinCenter(i))
        histoNJets['pho_mc'] = mcPhotonFile.Get("hNJets_DR_gjets")
        histoNJets['pho_mc'].Scale(mcLumiRatio_photon)
        histoNJets['zmm_da'] = dataZllFile.Get("hNJets_DR_zmm")
        histoNJets['zmm_mc'] = mcZllFile.Get("hNJets_DR_dymm")
        histoNJets['zmm_mc'].Scale(mcLumiRatio_mm/DYkfactor)
        histoNJets['zee_da'] = dataZllFile.Get("hNJets_DR_zee")
        histoNJets['zee_mc'] = mcZllFile.Get("hNJets_DR_dyee")
        histoNJets['zee_mc'].Scale(mcLumiRatio_ee/DYkfactor)
        histoMHT['pho_da'] = dataPhotonFile.Get("hMHT_DR_photon")
        histoMHT['pho_cg'] = dataPhotonFile.Get("hMHT_DR_xWt_photon")
        histoMHT['pho_cg'].Divide(histoMHT['pho_da'])
        histoMHT['pho_mc'] = mcPhotonFile.Get("hMHT_DR_gjets")
        histoMHT['pho_mc'].Scale(mcLumiRatio_photon)
        histoMHT['zmm_da'] = dataZllFile.Get("hMHT_DR_zmm")
        histoMHT['zmm_mc'] = mcZllFile.Get("hMHT_DR_dymm")
        histoMHT['zmm_mc'].Scale(mcLumiRatio_mm/DYkfactor)
        histoMHT['zee_da'] = dataZllFile.Get("hMHT_DR_zee")
        histoMHT['zee_mc'] = mcZllFile.Get("hMHT_DR_dyee")
        histoMHT['zee_mc'].Scale(mcLumiRatio_ee/DYkfactor)
        histoHT['pho_da'] = dataPhotonFile.Get("hHT_DR_photon")
        histoHT['pho_cg'] = dataPhotonFile.Get("hHT_DR_xWt_photon")
        histoHT['pho_cg'].Divide(histoHT['pho_da'])
        histoHT['pho_mc'] = mcPhotonFile.Get("hHT_DR_gjets")
        histoHT['pho_mc'].Scale(mcLumiRatio_photon)
        histoHT['zmm_da'] = dataZllFile.Get("hHT_DR_zmm")
        histoHT['zmm_mc'] = mcZllFile.Get("hHT_DR_dymm")
        histoHT['zmm_mc'].Scale(mcLumiRatio_mm/DYkfactor)
        histoHT['zee_da'] = dataZllFile.Get("hHT_DR_zee")
        histoHT['zee_mc'] = mcZllFile.Get("hHT_DR_dyee")
        histoHT['zee_mc'].Scale(mcLumiRatio_ee/DYkfactor)

        ## get the double ratio graphs
        nj_dr = RA2b.getDoubleRatioGraph('NJets', histos=histoNJets)
        mht_dr = RA2b.getDoubleRatioGraph('MHT', histos=histoMHT)
        ht_dr = RA2b.getDoubleRatioGraph('HT', histos=histoHT)
        
        binMeanNJets = []
        hNJetsCC_pho_da = dataPhotonFile.Get("hNJets_DRCC_photon")
        hNJetsCC_pho_cg = dataPhotonFile.Get("hNJets_DRCC_xWt_photon")
        hNJetsCC_pho_cg.Divide(hNJetsCC_pho_da)
        for bin in range(1, hNJetsCC_pho_cg.GetNbinsX()+1):
            binMeanNJets.append(hNJetsCC_pho_cg.GetBinContent(bin))
        binMeanMHT = []
        hMHTCC_pho_da = dataPhotonFile.Get("hMHT_DRCC_photon")
        hMHTCC_pho_cg = dataPhotonFile.Get("hMHT_DRCC_xWt_photon")
        hMHTCC_pho_cg.Divide(hMHTCC_pho_da)
        for bin in range(1, hMHTCC_pho_cg.GetNbinsX()+1):
            binMeanMHT.append(hMHTCC_pho_cg.GetBinContent(bin))
        binMeanHT = []
        hHT1CC_pho_da = dataPhotonFile.Get("hHT1_DRCC_photon")
        hHT1CC_pho_cg = dataPhotonFile.Get("hHT1_DRCC_xWt_photon")
        hHT1CC_pho_cg.Divide(hHT1CC_pho_da)
        for bin in range(1, hHT1CC_pho_cg.GetNbinsX()+1):
            binMeanHT.append(hHT1CC_pho_cg.GetBinContent(bin))
        hHT2CC_pho_da = dataPhotonFile.Get("hHT2_DRCC_photon")
        hHT2CC_pho_cg = dataPhotonFile.Get("hHT2_DRCC_xWt_photon")
        hHT2CC_pho_cg.Divide(hHT2CC_pho_da)
        for bin in range(1, hHT2CC_pho_cg.GetNbinsX()+1):
            binMeanHT.append(hHT2CC_pho_cg.GetBinContent(bin))
        binMeans = {}
        binMeans['NJets'] = binMeanNJets
        binMeans['MHT'] = binMeanMHT
        binMeans['HT'] = binMeanHT

        ## get the double ratio plots with values and uncertainties
        dr_out = RA2b.getDoubleRatioPlot([nj_dr, mht_dr, ht_dr], binMeans = binMeans,
                                         iPeriod = iPeriod, Ymin = Ymin, Ymax = Ymax,
                                         plotRefLine = plotRefLine, yTitle = yTitle)

    axisTitleSize = 0.033
    def drawPad(canv, pad, margins, range, fontScale, Title):  # Put one graph on a multi-graph canvas
        pad.SetTopMargin(margins[0])
        pad.SetRightMargin(margins[1])
        pad.SetBottomMargin(margins[2])
        pad.SetLeftMargin(margins[3])
        for obj in canv.GetListOfPrimitives():
            print "Object in pad:  "+obj.GetName()+" is a "+str(type(obj))
            if (type(obj) == ROOT.TGraphAsymmErrors):
                obj.SetMinimum(0.5)
                obj.SetMaximum(1.5)
                obj.GetXaxis().SetRangeUser(range[0], range[1])
                obj.GetXaxis().SetTitleSize(2*axisTitleSize)
                obj.GetXaxis().SetTitleOffset(0.030/(1.0*axisTitleSize))
                obj.GetYaxis().SetTitleSize(2*axisTitleSize)
                obj.GetYaxis().SetTitleOffset(0.040/(1.0*axisTitleSize))
                obj.GetXaxis().SetLabelSize(1.5*axisTitleSize)
                obj.GetXaxis().SetLabelOffset(fontScale*0.000165/(1.0*axisTitleSize))
                obj.GetYaxis().SetLabelSize(fontScale*1.5*axisTitleSize)
                obj.GetXaxis().SetTitle(Title[0])
                obj.GetYaxis().SetTitle(Title[1])
                obj.SetLineWidth(1)
                obj.SetMarkerStyle(20)
                obj.SetMarkerSize(1)
                obj.Draw("ap")
        for obj in canv.GetListOfPrimitives():
            if (type(obj) == ROOT.TGraph):
                obj.SetLineWidth(1)
                obj.Draw("same")
            if (type(obj) == ROOT.TF1):
                if (obj.GetName() == "refLine"):
                    obj.SetLineWidth(2)
                    obj.DrawF1(range[0], 900, "same")
                    satVal = obj.Eval(900)
                    print (obj.GetName()+" has "+str(obj.GetNumberFreeParameters())+" parameters:  "
                           +"p0 = "+str(obj.GetParameter(0))+", p1 = "
                           +str(obj.GetParameter(1))+", satVal = "+str(satVal))
                    satLine = ROOT.TLine()
                    satLine.SetLineWidth(obj.GetLineWidth())
                    satLine.SetLineStyle(obj.GetLineStyle())
                    satLine.SetLineColor(obj.GetLineColor())
                    satLine.DrawLine(900, satVal, range[1], satVal)
                else:
                    obj.SetLineWidth(1)
                    obj.Draw("same")

    def textCMS(margin):
        CMStxt = ROOT.TPaveText(margin, 0.94, 0.23, 1.0, "NDC")
        CMStxt.SetBorderSize(0)
        CMStxt.SetFillColor(0)
        CMStxt.SetTextFont(42)
        CMStxt.SetTextAlign(13)
        CMStxt.SetTextSize(0.09)
        CMStxt.SetMargin(0.)
        CMStxt.AddText("#bf{CMS}")
        return CMStxt

    def textLumi(margin):
        Lumitxt = ROOT.TPaveText(0.80, 0.92, 1-margin, 1.0, "NDC")
        Lumitxt.SetBorderSize(0)
        Lumitxt.SetFillColor(0)
        Lumitxt.SetTextFont(42)
        Lumitxt.SetTextAlign(33)
        Lumitxt.SetTextSize(0.07)
        Lumitxt.SetMargin(0.)
        Lumitxt.AddText("137 fb^{-1} (13 TeV)")
        return Lumitxt

    if (sample == "sig"):
        # Set up canvas dimensions
        Lcanv = 600
        mlpad1 = 0.19
        mrpad3 = 0.03
        mtpad = 0.10
        wscanv = 1 / (3 - 2*mlpad1 + mrpad3)
        newCanv = ROOT.TCanvas("newCanv", "newCanv", 0, 0, int(Lcanv/wscanv), int(1.2*Lcanv))
        newCanv.Range(0, 0, 1, 1)
        c3 = ROOT.gROOT.FindObject("c3")
        if (type(c3)==ROOT.TCanvas):
            print "margins T, R, B, L = "+str(c3.GetTopMargin())+", "+str(c3.GetRightMargin())+", "+str(c3.GetBottomMargin())+", "+str(c3.GetLeftMargin())+", "
            newCanv.cd()
            pad1 = ROOT.TPad("pad1", "pad1", 0, 0, wscanv, 1)
            print "pad1 width = "+str(pad1.GetWNDC())
            pad1.Draw()
            pad1.cd()
            drawPad(c3, pad1, [mtpad, 0, 0.15, mlpad1], [300., 1545.], 1,
                    ["H_{T} [GeV]", yTitle])
            # CMStxt = textCMS(mlpad1)
            # CMStxt.Draw()
            RA2b.cmsLumi(pad = pad1, iPeriod = -1, extraText = extraText)
        axisTitleSize *= 1.0/(1 - 1*mlpad1)
        c2 = ROOT.gROOT.FindObject("c2")
        if (type(c2)==ROOT.TCanvas):
            newCanv.cd()
            pad2 = ROOT.TPad("pad2", "pad2", wscanv, 0, wscanv*(2-mlpad1), 1)
            print "pad2 width = "+str(pad2.GetWNDC())
            pad2.Draw()
            pad2.cd()
            drawPad(c2, pad2, [mtpad, 0, 0.15, 0], [310., 880.], -1.5,  # This -1.5 is an emprical fudge
                    ["H^{miss}_{T} [GeV]", ""])
        c1 = ROOT.gROOT.FindObject("c1")
        if (type(c1)==ROOT.TCanvas):
            newCanv.cd()
            pad3 = ROOT.TPad("pad3", "pad3", wscanv*(2-mlpad1), 0, wscanv*(3-2*mlpad1+mrpad3), 1)
            print "pad3 width = "+str(pad3.GetWNDC())
            pad3.Draw()
            pad3.cd()
            drawPad(c1, pad3, [mtpad, mrpad3, 0.15, 0], [1.01, 9.99], -1.5,
                    ["N_{jet}", ""])
            Lumitxt = textLumi(mrpad3)
            Lumitxt.Draw()

        newCanv.SaveAs("DRtriple.pdf")
        newCanv.SaveAs("DRtriple.png")

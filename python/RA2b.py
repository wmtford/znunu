"""Module containing tools for Z->nunu estimation studies.
Designed to run on LPC produced skims
T. Mulholland: troy.mulholland@cern.ch"""

import time
import re
import math
import random
import numpy
from array import array
import ROOT

def cmsLumi(pad,  iPeriod=None,  iPosX=None, extraText=None ):

    """ modifed by TM 
    modified form of CMS_lumi
    Initiated by: Gautier Hamel de Monchenault (Saclay)
    Translated in Python by: Joshua Hardenbrook (Princeton)
    Updated by:   Dinko Ferencek (Rutgers) """

    cmsText     = "CMS";
    cmsTextFont   = 61  

    writeExtraText = True
    if(extraText==None):
        extraText   = "Preliminary"
    if(iPeriod==None):
        iPeriod = 6
    if(iPosX==None):
        iPosX=0

    extraTextFont = 52 

    lumiTextSize     = 0.6
    lumiTextOffset   = 0.2

    cmsTextSize      = 0.75
    cmsTextOffset    = 0.1

    relPosX    = 0.045
    relPosY    = 0.035
    relExtraDY = 1.2

    extraOverCmsTextSize  = 0.76
    
    lumi_13TeV = "2.3 fb^{-1}"
    lumi_13TeV_2016 = "20 fb^{-1}"
    lumi_13TeV_V8 = "4.0 fb^{-1}"
    lumi_13TeV_V9 = "24.5 fb^{-1}"
    lumi_8TeV  = "19.7 fb^{-1}" 
    lumi_7TeV  = "5.1 fb^{-1}"
    lumi_sqrtS = ""

    drawLogo      = False
    outOfFrame    = False
    if(iPosX/10==0 ): outOfFrame = True

    alignY_=3
    alignX_=2
    if( iPosX/10==0 ): alignX_=1
    if( iPosX==0    ): alignY_=1
    if( iPosX/10==1 ): alignX_=1
    if( iPosX/10==2 ): alignX_=2
    if( iPosX/10==3 ): alignX_=3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025

    pad.cd()

    lumiText = ""
    if( iPeriod==1 ):
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif ( iPeriod==2 ):
        lumiText += lumi_8TeV
        lumiText += " (8 TeV)"

    elif( iPeriod==3 ):      
        lumiText = lumi_8TeV 
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif ( iPeriod==4 ):
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
    elif ( iPeriod==5 ):
        lumiText += lumi_13TeV_2016
        lumiText += " (13 TeV)"
    elif ( iPeriod==6 ):
        lumiText += lumi_13TeV_V9
        lumiText += " (13 TeV)"
    elif ( iPeriod==7 ):
        if( outOfFrame ):lumiText += "#scale[0.85]{"
        lumiText += lumi_13TeV 
        lumiText += " (13 TeV)"
        lumiText += " + "
        lumiText += lumi_8TeV 
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
        if( outOfFrame): lumiText += "}"
    elif ( iPeriod==12 ):
        lumiText += "8 TeV"
    elif ( iPeriod==0 ):
        lumiText += lumi_sqrtS
            
    #print lumiText

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)    
    
    extraTextSize = extraOverCmsTextSize*cmsTextSize
    
    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)    

    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)

    if( outOfFrame ):
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11) 
        latex.SetTextSize(cmsTextSize*t)    
        latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText)
  
    pad.cd()

    posX_ = 0
    if( iPosX%10<=1 ):
        posX_ =   l + relPosX*(1-l-r)
    elif( iPosX%10==2 ):
        posX_ =  l + 0.5*(1-l-r)
    elif( iPosX%10==3 ):
        posX_ =  1-r - relPosX*(1-l-r)

    posY_ = 1-t - relPosY*(1-t-b)

    if( not outOfFrame ):
        if( drawLogo ):
            posX_ =   l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            CMS_logo = ROOT.TASImage("CMS-BW-label.png")
            pad_logo =  ROOT.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()          
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if( writeExtraText ) :
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize*t)
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText)
    elif( writeExtraText ):
        if( iPosX==0):
            posX_ =   0.075+l +  relPosX*(1-l-r)
            posY_ =   1-t+lumiTextOffset*t

        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)      

    pad.Update()


def tdrGrid( gridOn):
    tdrStyle.SetPadGridX(gridOn)
    tdrStyle.SetPadGridY(gridOn)

#fixOverlay: Redraws the axis
def fixOverlay(): gPad.RedrawAxis()

def setTdrStyle():
    tdrStyle =  ROOT.TStyle("tdrStyle","Style for P-TDR")

    #for the canvas:
    tdrStyle.SetCanvasBorderMode(0)
    tdrStyle.SetCanvasColor(ROOT.kWhite)
    tdrStyle.SetCanvasDefH(600) #Height of canvas
    tdrStyle.SetCanvasDefW(600) #Width of canvas
    tdrStyle.SetCanvasDefX(0)   #POsition on screen
    tdrStyle.SetCanvasDefY(0)
    
    tdrStyle.SetPadBorderMode(0)
    tdrStyle.SetPadColor(ROOT.kWhite)
    tdrStyle.SetPadGridX(False)
    tdrStyle.SetPadGridY(False)
    tdrStyle.SetGridColor(0)
    tdrStyle.SetGridStyle(3)
    tdrStyle.SetGridWidth(1)
    
    #For the frame:
    tdrStyle.SetFrameBorderMode(0)
    tdrStyle.SetFrameBorderSize(1)
    tdrStyle.SetFrameFillColor(0)
    tdrStyle.SetFrameFillStyle(0)
    tdrStyle.SetFrameLineColor(1)
    tdrStyle.SetFrameLineStyle(1)
    tdrStyle.SetFrameLineWidth(1)
    
    #For the histo:
    tdrStyle.SetHistLineColor(1)
    tdrStyle.SetHistLineStyle(0)
    tdrStyle.SetHistLineWidth(1)
    
    tdrStyle.SetEndErrorSize(2)
    #tdrStyle.SetErrorMarker(20)
    #tdrStyle.SetErrorX(0.)
    
    tdrStyle.SetMarkerStyle(20)
    
    #For the fit/function:
    tdrStyle.SetOptFit(1)
    tdrStyle.SetFitFormat("5.4g")
    tdrStyle.SetFuncColor(2)
    tdrStyle.SetFuncStyle(1)
    tdrStyle.SetFuncWidth(1)

    #For the date:
    tdrStyle.SetOptDate(0)
    
    # For the statistics box:
    tdrStyle.SetOptFile(0)
    tdrStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
    tdrStyle.SetStatColor(ROOT.kWhite)
    tdrStyle.SetStatFont(42)
    tdrStyle.SetStatFontSize(0.025)
    tdrStyle.SetStatTextColor(1)
    tdrStyle.SetStatFormat("6.4g")
    tdrStyle.SetStatBorderSize(1)
    tdrStyle.SetStatH(0.1)
    tdrStyle.SetStatW(0.15)
    
    # Margins:
    tdrStyle.SetPadTopMargin(0.05)
    tdrStyle.SetPadBottomMargin(0.13)
    tdrStyle.SetPadLeftMargin(0.16)
    tdrStyle.SetPadRightMargin(0.02)
    
    # For the Global title:
    
    tdrStyle.SetOptTitle(0)
    tdrStyle.SetTitleFont(42)
    tdrStyle.SetTitleColor(1)
    tdrStyle.SetTitleTextColor(1)
    tdrStyle.SetTitleFillColor(10)
    tdrStyle.SetTitleFontSize(0.05)
    
    # For the axis titles:
    
    tdrStyle.SetTitleColor(1, "XYZ")
    tdrStyle.SetTitleFont(42, "XYZ")
    tdrStyle.SetTitleSize(0.06, "XYZ")
    tdrStyle.SetTitleXOffset(0.9)
    tdrStyle.SetTitleYOffset(1.25)
  
    # For the axis labels:
    
    tdrStyle.SetLabelColor(1, "XYZ")
    tdrStyle.SetLabelFont(42, "XYZ")
    tdrStyle.SetLabelOffset(0.007, "XYZ")
    tdrStyle.SetLabelSize(0.05, "XYZ")
    
    # For the axis:
    
    tdrStyle.SetAxisColor(1, "XYZ")
    tdrStyle.SetStripDecimals(True)
    tdrStyle.SetTickLength(0.03, "XYZ")
    tdrStyle.SetNdivisions(510, "XYZ")
    tdrStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
    tdrStyle.SetPadTickY(1)
    
    # Change for log plots:
    tdrStyle.SetOptLogx(0)
    tdrStyle.SetOptLogy(0)
    tdrStyle.SetOptLogz(0)

    # Postscript options:
    tdrStyle.SetPaperSize(20.,20.)
    tdrStyle.SetHatchesLineWidth(5)
    tdrStyle.SetHatchesSpacing(0.05)
    
    tdrStyle.cd()


def getTreeWeight(File,  doLumi=None,  treeName=None, removeZkfactor=None, getWeight=None, setWeight=None, doNoisy=None):

    if(treeName==None):
        treeName="tree"
    if(doLumi==None):
        doLumi = 24.5
    if(doNoisy==None):
        doNoisy=False
    if(removeZkfactor==None):
        removeZkfactor=False
    if(getWeight==None): 
        getWeight=True
    if(setWeight==None): # for manually setting the weight
        setWeight=False

    Zkapa = 1.0
    if(removeZkfactor==True):
        Zkapa=1.23

    f = ROOT.TFile(File,"update")
    tree = f.Get(treeName)

    if(type(setWeight) is int or type(setWeight) is float):
        weight = setWeight
    else:
        WeightHist = ROOT.TH1F("WeightHist","WeightHist",1,-10000,10000)
        tree.Project("WeightHist","Weight");
    
        weight = WeightHist.GetMean()

        if(getWeight):
            xsecHist = ROOT.TH1F("xsecHist","xsecHist",1,0,1000000000)
            tree.Project("xsecHist","CrossSection")
            xsection = float(xsecHist.GetMean())

            nEventProc = f.Get("nEventProc")
            nEventNeg = f.Get("nEventNeg")

            nEvents = nEventProc.GetBinContent(1) - 2*nEventNeg.GetBinContent(1)

            weight = xsection/nEvents;
            if(weight < 0):
                weight *= -1;
        sampleLumi = Zkapa/(weight*1000.)

        weight *= 1000*doLumi/Zkapa

    return weight

def setTreeWeight(File,  doLumi=None,  treeName=None, removeZkfactor=None, getWeight=None, setWeight=None, doNoisy=None):
    """Sets the weight of a tree to a corresponding lumi. Default lumi is 24.5 ifb"""

    if(treeName==None):
        treeName="tree"

    f = ROOT.TFile(File,"update")
    tree = f.Get(treeName)

    weight = getTreeWeight(File, doLumi, treeName, removeZkfactor, getWeight, setWeight, doNoisy)
  
    tree.SetWeight(weight)
    tree.AutoSave() 

def getFileList(sample, doLumi=None, treeLoc=None, treeName=None, removeZkfactor=None):

    if(doLumi==None):
        doLumi = 24.5
    if(treeLoc==None):
        treeLoc = "/home/ww/work/data/lpcTrees/Skims/Run2ProductionV10"
    if(treeName==None):
        treeName = "tree"
    if(removeZkfactor==None and (sample=='zinv' or sample=='dymm' or sample=='dyee' or sample=='dyll')):
        removeZkfactor=True
    elif(removeZkfactor==None):
        removeZkfactor=False

    fileList = []

    if('sig' in sample):
        if('LDP' not in sample):
            fileList.append(treeLoc+"/tree_signal/tree_HTMHT_2016B.root")
            fileList.append(treeLoc+"/tree_signal/tree_HTMHT_2016C.root")
            fileList.append(treeLoc+"/tree_signal/tree_HTMHT_2016D.root")
            fileList.append(treeLoc+"/tree_signal/tree_HTMHT_2016E.root")
            fileList.append(treeLoc+"/tree_signal/tree_HTMHT_2016F.root")
            fileList.append(treeLoc+"/tree_signal/tree_HTMHT_2016G.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            fileList.append(treeLoc+"/tree_LDP/tree_HTMHT_2016B.root")
            fileList.append(treeLoc+"/tree_LDP/tree_HTMHT_2016C.root")
            fileList.append(treeLoc+"/tree_LDP/tree_HTMHT_2016D.root")
            fileList.append(treeLoc+"/tree_LDP/tree_HTMHT_2016E.root")
            fileList.append(treeLoc+"/tree_LDP/tree_HTMHT_2016F.root")
            fileList.append(treeLoc+"/tree_LDP/tree_HTMHT_2016G.root")
    elif(sample.startswith('zmm')):
        if('LDP' not in sample):
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016B.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016C.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016D.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016E.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016F.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016G.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016B.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016C.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016D.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016E.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016F.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016G.root") 
    elif(sample.startswith('zee')):
        if('LDP' not in sample):
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016B.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016C.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016D.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016E.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016F.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016G.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016B.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016C.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016D.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016E.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016F.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016G.root") 
    elif(sample.startswith('zll')):
        if('LDP' not in sample):
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016B.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016B.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016C.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016C.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016D.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016D.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016E.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016E.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_SingleMuon_2016F.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016F.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016G.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_SingleElectron_2016G.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016B.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016B.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016C.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016C.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016D.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016D.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016E.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016E.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_SingleMuon_2016F.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016F.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016G.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_SingleElectron_2016G.root") 
    elif(sample.startswith('photon')):
        if('LDP' not in sample):
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_SinglePhoton_2016B.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_SinglePhoton_2016C.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_SinglePhoton_2016D.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_SinglePhoton_2016E.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_SinglePhoton_2016F.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_SinglePhoton_2016G.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_SinglePhoton_2016B.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_SinglePhoton_2016C.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_SinglePhoton_2016D.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_SinglePhoton_2016E.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_SinglePhoton_2016F.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_SinglePhoton_2016G.root") 
    elif('zinv' in sample):
        if('LDP' not in sample):
            #setTreeWeight(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-600to800.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-800to1200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-1200to2500.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-2500toInf.root", doLumi, treeName, removeZkfactor)
            #fileList.append(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-600to800.root") 
            fileList.append(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-800to1200.root") 
            fileList.append(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-1200to2500.root") 
            fileList.append(treeLoc+"/tree_signal/tree_ZJetsToNuNu_HT-2500toInf.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-600to800.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-800to1200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-1200to2500.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-2500toInf.root", doLumi, treeName, removeZkfactor)
            fileList.append(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-600to800.root") 
            fileList.append(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-800to1200.root") 
            fileList.append(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-1200to2500.root") 
            fileList.append(treeLoc+"/tree_LDP/tree_ZJetsToNuNu_HT-2500toInf.root") 
    
    elif('dymm' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
    elif('dyee' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
    elif('dyll' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName, removeZkfactor)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName, removeZkfactor)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_DYJetsToLL_M-50_HT-600toInf.root") 
    elif('gjets' in sample):
        if('new' in sample):
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-100to200.root", doLumi, treeName, setWeight=2.48331699336*doLumi/7.63)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-200to400.root", doLumi, treeName, setWeight=0.139931580097*doLumi/7.63)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-400to600.root", doLumi, treeName, setWeight=0.0536453416901*doLumi/7.63)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-600toInf.root", doLumi, treeName, setWeight=0.015346951535*doLumi/7.63) 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_DR-0p4_HT-600toInf.root") 
        elif('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-600toInf.root", doLumi, treeName) 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_GJets_HT-600toInf.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-600toInf.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-100to200.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-200to400.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-400to600.root") 
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_GJets_HT-600toInf.root") 
    elif('ttgjets' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_TTGJets.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_TTGJets.root") 
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_TTGJets.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_TTGJets.root") 
    elif('gjetsqcd' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-200to300.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-300to500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-500to700.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-700to1000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-1000to1500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-1500to2000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-2000toInf.root", doLumi, treeName)
            #fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-200to300.root")
            #fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-300to500.root")
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-500to700.root")
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-700to1000.root")
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-1000to1500.root")
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-1500to2000.root")
            fileList.append(treeLoc+"/tree_GJet_CleanVars/tree_QCD_HT-2000toInf.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-200to300.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-300to500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-500to700.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-700to1000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-1000to1500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-1500to2000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-2000toInf.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-200to300.root")
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-300to500.root")
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-500to700.root")
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-700to1000.root")
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-1000to1500.root")
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-1500to2000.root")
            fileList.append(treeLoc+"/tree_GJetLDP_CleanVars/tree_QCD_HT-2000toInf.root")
    elif('ttzvv' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_signal/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_signal/tree_TTZToLLNuNu.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_LDP/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_LDP/tree_TTZToLLNuNu.root")
    elif(sample.startswith('ttzmm')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTZToLLNuNu.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTZToLLNuNu.root")
    elif(sample.startswith('ttzee')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTZToLLNuNu.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTZToLLNuNu.root")
    elif(sample.startswith('ttzll')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTZToLLNuNu.root")
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTZToLLNuNu.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTZToLLNuNu.root")
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTZToLLNuNu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTZToLLNuNu.root")
    elif('ttmm' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromT.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromT.root")
    elif('ttee' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromT.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromT.root")
    elif(sample=='ttll'):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_TTJets_SingleLeptFromT.root")
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_TTJets_SingleLeptFromT.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_DiLept.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_DiLept.root")
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_TTJets_SingleLeptFromT.root")
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromTbar.root")
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_TTJets_SingleLeptFromT.root")

    elif(sample.startswith('dibosonmm')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_WWTo2L2Nu.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_WWTo2L2Nu.root")
    elif(sample.startswith('dibosonee')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_WWTo2L2Nu.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_WWTo2L2Nu.root")
    elif(sample.startswith('dibosonll')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_WWTo2L2Nu.root")
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_WWTo2L2Nu.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_WWTo2L2Nu.root")
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_ZZTo2L2Q.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_WWTo2L2Nu.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_ZZTo2L2Q.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_WWTo2L2Nu.root")
    elif(sample.startswith('tribosonmm')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_ZZZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_WWZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYm_CleanVars/tree_WZZ.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_ZZZ.root")
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_WWZ.root")
            fileList.append(treeLoc+"/tree_DYm_CleanVars/tree_WZZ.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_ZZZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_WWZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYmLDP_CleanVars/tree_WWZ.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_ZZZ.root")
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_WWZ.root")
            fileList.append(treeLoc+"/tree_DYmLDP_CleanVars/tree_WZZ.root")
    elif(sample.startswith('tribosonee')):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_ZZZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_WWZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYe_CleanVars/tree_WZZ.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_ZZZ.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_WWZ.root")
            fileList.append(treeLoc+"/tree_DYe_CleanVars/tree_WZZ.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_ZZZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_WWZ.root", doLumi, treeName, False, True)
            setTreeWeight(treeLoc+"/tree_DYeLDP_CleanVars/tree_WZZ.root", doLumi, treeName, False, True)
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_ZZZ.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_WWZ.root")
            fileList.append(treeLoc+"/tree_DYeLDP_CleanVars/tree_WZZ.root")

    elif('top' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_signal/tree_TTJets.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_ST_s-channel.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_ST_t-channel_antitop.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_ST_t-channel_top.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_ST_tW_antitop.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_ST_tW_top.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_signal/tree_TTJets.root")
            fileList.append(treeLoc+"/tree_signal/tree_TTJets_SingleLeptFromTbar.root")
            fileList.append(treeLoc+"/tree_signal/tree_TTJets_SingleLeptFromT.root")
            fileList.append(treeLoc+"/tree_signal/tree_ST_s-channel.root")
            fileList.append(treeLoc+"/tree_signal/tree_ST_t-channel_antitop.root")
            fileList.append(treeLoc+"/tree_signal/tree_ST_t-channel_top.root")
            fileList.append(treeLoc+"/tree_signal/tree_ST_tW_antitop.root")
            fileList.append(treeLoc+"/tree_signal/tree_ST_tW_top.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_LDP/tree_TTJets.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_TTJets_SingleLeptFromTbar.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_TTJets_SingleLeptFromT.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ST_s-channel.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ST_t-channel_antitop.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ST_t-channel_top.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ST_tW_antitop.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_ST_tW_top.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_LDP/tree_TTJets.root")
            fileList.append(treeLoc+"/tree_LDP/tree_TTJets_SingleLeptFromTbar.root")
            fileList.append(treeLoc+"/tree_LDP/tree_TTJets_SingleLeptFromT.root")
            fileList.append(treeLoc+"/tree_LDP/tree_ST_s-channel.root")
            fileList.append(treeLoc+"/tree_LDP/tree_ST_t-channel_antitop.root")
            fileList.append(treeLoc+"/tree_LDP/tree_ST_t-channel_top.root")
            fileList.append(treeLoc+"/tree_LDP/tree_ST_tW_antitop.root")
            fileList.append(treeLoc+"/tree_LDP/tree_ST_tW_top.root")
    elif('wjets' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-600to800.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-800to1200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-1200to2500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-2500toInf.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-100to200.root")
            fileList.append(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-200to400.root")
            fileList.append(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-400to600.root")
            fileList.append(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-600to800.root")
            fileList.append(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-800to1200.root")
            fileList.append(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-1200to2500.root")
            fileList.append(treeLoc+"/tree_signal/tree_WJetsToLNu_HT-2500toInf.root")
        if( ('LDP' in sample) or ('IDP' in sample)):        
            setTreeWeight(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-600to800.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-800to1200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-1200to2500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-2500toInf.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-100to200.root")
            fileList.append(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-200to400.root")
            fileList.append(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-400to600.root")
            fileList.append(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-600to800.root")
            fileList.append(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-800to1200.root")
            fileList.append(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-1200to2500.root")
            fileList.append(treeLoc+"/tree_LDP/tree_WJetsToLNu_HT-2500toInf.root")
    elif('qcd' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_signal/tree_QCD_HT-200to300.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_QCD_HT-300to500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_QCD_HT-500to700.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_QCD_HT-700to1000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_QCD_HT-1000to1500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_QCD_HT-1500to2000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_QCD_HT-2000toInf.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_signal/tree_QCD_HT-200to300.root")
            fileList.append(treeLoc+"/tree_signal/tree_QCD_HT-300to500.root")
            fileList.append(treeLoc+"/tree_signal/tree_QCD_HT-500to700.root")
            fileList.append(treeLoc+"/tree_signal/tree_QCD_HT-700to1000.root")
            fileList.append(treeLoc+"/tree_signal/tree_QCD_HT-1000to1500.root")
            fileList.append(treeLoc+"/tree_signal/tree_QCD_HT-1500to2000.root")
            fileList.append(treeLoc+"/tree_signal/tree_QCD_HT-2000toInf.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_LDP/tree_QCD_HT-200to300.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_QCD_HT-300to500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_QCD_HT-500to700.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_QCD_HT-700to1000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_QCD_HT-1000to1500.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_QCD_HT-1500to2000.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_LDP/tree_QCD_HT-2000toInf.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_LDP/tree_QCD_HT-200to300.root")
            fileList.append(treeLoc+"/tree_LDP/tree_QCD_HT-300to500.root")
            fileList.append(treeLoc+"/tree_LDP/tree_QCD_HT-500to700.root")
            fileList.append(treeLoc+"/tree_LDP/tree_QCD_HT-700to1000.root")
            fileList.append(treeLoc+"/tree_LDP/tree_QCD_HT-1000to1500.root")
            fileList.append(treeLoc+"/tree_LDP/tree_QCD_HT-1500to2000.root")
            fileList.append(treeLoc+"/tree_LDP/tree_QCD_HT-2000toInf.root")
    elif('other' in sample):
        if('LDP' not in sample):
            setTreeWeight(treeLoc+"/tree_signal/tree_TTZToLLNuNu.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_TTZToQQ.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-600toInf.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_TTGJets.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_signal/tree_TTZToLLNuNu.root")
            fileList.append(treeLoc+"/tree_signal/tree_TTZToQQ.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-100to200.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-200to400.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-400to600.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-600toInf.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-100to200.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-200to400.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-400to600.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-600toInf.root")
            fileList.append(treeLoc+"/tree_signal/tree_TTGJets.root")
        if( ('LDP' in sample) or ('IDP' in sample)):
            setTreeWeight(treeLoc+"/tree_signal/tree_TTZToLLNuNu.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_TTZToQQ.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-600toInf.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-100to200.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-200to400.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-400to600.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_GJets_HT-600toInf.root", doLumi, treeName)
            setTreeWeight(treeLoc+"/tree_signal/tree_TTGJets.root", doLumi, treeName)
            fileList.append(treeLoc+"/tree_signal/tree_TTZToLLNuNu.root")
            fileList.append(treeLoc+"/tree_signal/tree_TTZToQQ.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-100to200.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-200to400.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-400to600.root")
            fileList.append(treeLoc+"/tree_signal/tree_DYJetsToLL_M-50_HT-600toInf.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-100to200.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-200to400.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-400to600.root")
            fileList.append(treeLoc+"/tree_signal/tree_GJets_HT-600toInf.root")
            fileList.append(treeLoc+"/tree_signal/tree_TTGJets.root")
    else:
        print "***** invalid/unimplemented sample code *****"

    return fileList

def getChain(sample, doLumi=None, treeName=None, doProof=None, treeLoc=None, removeZkfactor=None):
    """ returns sample TChain
    sample options:
    zmm: Z->mumu data
    zee: Z->ee data
    zll: Z->mumu + Z->ee data
    photon: photon data
    sig: signal data
    dymm: drell-yan (mumu) MC
    dyee: drell-yan (ee) MC
    dyll: drell-yan (mumu) + drell-yan (ee) MC
    ttzmm: ttZ->mumu MC
    ttzee: ttZ->ee MC
    ttzll: ttZ->mumu + ttZ->ee MC
    ttmm: ttbar->mumu MC
    ttee: ttbar->ee MC
    ttll: ttbar->mumu + ttbar->ee MC
    gjets: photon MC
    zinv: Z->nunu MC
    """    
    if(doLumi==None):
        doLumi = 24.5
    if(doProof==None):
        doProof = False
    if(treeLoc==None):
        treeLoc = "/home/ww/work/data/lpcTrees/Skims/Run2ProductionV10"
    if(treeName==None):
        treeName = "tree"
    if(removeZkfactor==None and (sample=='zinv' or sample=='dymm' or sample=='dyee' or sample=='dyll')):
        removeZkfactor=True
    elif(removeZkfactor==None):
        removeZkfactor=False

    #################################################################################
    # chain samples according to code
    #################################################################################

    chain = ROOT.TChain("tree")

    files = getFileList(sample, doLumi, treeLoc, treeName, removeZkfactor)

    for f in files:
        chain.Add(f)

    if(doProof):
        p = ROOT.TProof.Open()
        chain.SetProof()

    #################################################################################
    # end chain samples according to code
    #################################################################################

    return chain
  
def getCuts(code, nJetBin=None, bJetBin=None, kinBin=None, trig=None, dphiCut=None, extraCuts=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, applySF=None, njSplit=None, applyPuWeight=None):

    cuts = ROOT.TCut("")

    #####################################################################
    # set defaults and initialize 
    #####################################################################
    if(nJetBin==None):
        nJetBin=-1
    if(bJetBin==None):
        bJetBin=-1
    if(kinBin==None):
        kinBin=-1
    if(dphiCut==None):
        dphiCut = 'sig'
    if(applyMassCut==None):
        applyMassCut=True
    if(applyPtCut==None):
        applyPtCut=True
    if(applyHTCut==None):
        applyHTCut=True
    if(kinBin>10):
        if(applyMHTCut==None):
            applyMHTCut=False
    if(applyMHTCut==None):
        applyMHTCut=True
    if(applyNJetsCut==None):
        applyNJetsCut=False
    if(trig==None):
        trig = -1;
    if(applySF==None):
        applySF=False
    if(njSplit==None):
        njSplit=False
    if(applyPuWeight==None):
        applyPuWeight=False

    #####################################################################
    # end set defaults and initialize 
    #####################################################################

    #####################################################################
    # cut dictionaries  
    #####################################################################
    kinCut = {
        1: "HT>=300&&HT<500&&MHT>=300&&MHT<350",
        2: "HT>=500&&HT<1000&&MHT>=300&&MHT<350",
        3: "HT>=1000&&MHT>=300&&MHT<350",
        4: "HT>=350&&HT<500&&MHT>=350&&MHT<500",
        5: "HT>=500&&HT<1000&&MHT>=350&&MHT<500",
        6: "HT>=1000&&MHT>=350&&MHT<500",
        7: "HT>=500&&HT<1000&&MHT>=500&&MHT<750",
        8: "HT>=1000&&MHT>=500&&MHT<750",
        9: "HT>=750&&HT<1500&&MHT>=750",
        10: "HT>=1500&&MHT>=750",
        -1: "",
        11: "HT>=300&&HT<500&&MHT>=250&&MHT<=300",
        12: "HT>=500&&HT<1000&&MHT>=250&&MHT<=300",
        13: "HT>=1000&&MHT>=250&&MHT<=300",
    }
    nJetCuts = {
        1: "NJets==2",
        2: "NJets>=3&&NJets<=4",
        3: "NJets>=5&&NJets<=6",
        4: "NJets>=7&&NJets<=8",
        5: "NJets>=9",
        -1:"NJets>=2",
    }
    if(njSplit==True):
        nJetCuts = {
            1: "NJets==2",
            2: "NJets==3",
            3: "NJets==4",
            4: "NJets==5",
            5: "NJets==6",
            6: "NJets==7",
            7: "NJets==8",
            8: "NJets>=9",
            -1:"NJets>=2",
        }
        
    bJetCuts = {
        0: "BTags==0",
        1: "BTags==1",
        2: "BTags==2",
        3: "BTags>=3",
        -1:"BTags>=0",
        11:"BTags>=1",
        12:"BTags>=2",
    }
    bJetCutsSF = {
        0: "BTagsSF[0]*(1)",
        1: "BTagsSF[1]*(1)",
        2: "BTagsSF[2]*(1)",
        3: "BTagsSF[3]*(1)",
        -1:"(1)",
    }
    ## sig and hdp are identical
    dphiCuts = {
        "sig" : "DeltaPhi1>0.5&&DeltaPhi2>0.5&&DeltaPhi3>0.3&&DeltaPhi4>0.3",
        "hdp" : "DeltaPhi1>0.5&&DeltaPhi2>0.5&&DeltaPhi3>0.3&&DeltaPhi4>0.3",
        "ldp" : "(DeltaPhi1<0.5||DeltaPhi2<0.5||DeltaPhi3<0.3||DeltaPhi4<0.3)",
        "none" : ""
    }

    # 0 = sig, 1 = dimuon, 2 = dielectron, 3 = dilepton, 4 = photon, 5 = ttZ->nunu
    codeList = [0, 1, 2, 3, 4, 5, 6]

    if(type(code) is str):
        codeList = ["sig","zmm","zee","zll","photon","photonqcd","ttz"]

    Zcuts = {
        codeList[0]: "@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0",
        codeList[1]: "@Muons.size()==2&&@Electrons.size()==0&&isoElectronTracks==0&&isoPionTracks==0",
        codeList[2]: "@Muons.size()==0&&@Electrons.size()==2&&isoMuonTracks==0&&isoPionTracks==0",
        codeList[3]: "((@Muons.size()==2&&@Electrons.size()==0&&isoElectronTracks==0&&isoPionTracks==0)||(@Muons.size()==0&&@Electrons.size()==2&&isoMuonTracks==0&&isoPionTracks==0))",
        #codeList[4]: "(Sum$(Photons.Pt())>=200)&&(@Photons.size()==1)&&@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0",
        #codeList[4]: "Sum$(photon_nonPrompt)==0&&Photons[0].Pt()>=200&&Sum$(photon_fullID)==1&&@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0",
        # codeList[4]: "Sum$(photon_nonPrompt)==0&&Photons[0].Pt()>=200&&Sum$(photon_fullID)==1&&@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0",
        # codeList[5]: "Sum$(photon_nonPrompt)==1&&Photons[0].Pt()>=200&&Sum$(photon_fullID)==1&&@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0",
        codeList[4]: "Photons[0].Pt()>=200&Sum$(Photons_fullID)==1&&(@Photons.size()==1)&&@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0",
        codeList[5]: "Photons[0].Pt()>=200&&@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0",
        codeList[6]: "@Muons.size()==0&&@Electrons.size()==0&&isoElectronTracks==0&&isoMuonTracks==0&&isoPionTracks==0&&(@GenMuons.size()==0&&@GenElectrons.size()==0&&@GenTaus.size()==0)",
    }

    # define other cuts
    massCut = ""
    if((code == 'zmm' or code == 'zee' or code == 'zll') and applyMassCut):
        massCut = "ZCandidates.M()>=76.188&&ZCandidates.M()<=106.188"

    ptCut = ""
    if((code == 'zmm' or code == 'zee' or code == 'zll') and applyPtCut):
        ptCut = "ZCandidates.Pt()>=200."

    #commonCuts = "(JetID==1&&CSCTightHaloFilter==1&&HBHENoiseFilter==1&&HBHEIsoNoiseFilter==1&&eeBadScFilter==1&&EcalDeadCellTriggerPrimitiveFilter==1&&NVtx>0)"
    #commonCuts = "(JetID==1&&CSCTightHaloFilter==1&&globalTightHalo2016Filter && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && NVtx > 0)"
    #commonCuts = "(JetID==1)"

#(Data.globalTightHalo2016Filter==1 && Data.HBHENoiseFilter==1 && Data.HBHEIsoNoiseFilter==1 && Data.eeBadScFilter==1 && Data.EcalDeadCellTriggerPrimitiveFilter==1 && Data.BadChargedCandidateFilter && Data.BadPFMuonFilter && Data.JetID && Data.NVtx > 0

    if(trig==-1):
        commonCuts = "JetID==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0"
    else:
        commonCuts = "JetID==1 && globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0"

    trigCuts = ""
    if(type(trig) is tuple):
        trigCuts += "("
        for i in range(len(trig)):
            if(i==(len(trig)-1)):
                trigCuts+="(TriggerPass["+str(trig[i])+"]==1))"
            else:
                trigCuts+="(TriggerPass["+str(trig[i])+"]==1)||"
    elif(trig>=0):
        trigCuts = "TriggerPass["+str(trig)+"]==1"

    #####################################################################
    # end cut dictionaries  
    #####################################################################

    #####################################################################
    # add up cuts
    #####################################################################
    if(extraCuts!=None):
        cuts+=extraCuts
    if(applyHTCut==True):
        cuts+="HT>=300"
    else:
        cuts+="HT>=250"
    if(applyMHTCut==True):
        cuts+="MHT>=300"
    else:
        cuts+="MHT>=200"
    if(applyNJetsCut==True):
        cuts+="NJets>=3"
    else:
        cuts+="NJets>=2"

    cuts+=kinCut[kinBin]
    cuts+=nJetCuts[nJetBin]
    cuts+=dphiCuts[dphiCut]
    cuts+=Zcuts[code]
    cuts+=massCut
    cuts+=ptCut
    cuts+=commonCuts
    cuts+=trigCuts
    if(applySF):
        cuts*=bJetCutsSF[bJetBin]
    else:
        cuts+=bJetCuts[bJetBin]
    if(applyPuWeight):
        cuts*="puWeight*(1)"
    #####################################################################
    # end add up cuts
    #####################################################################

    #cuts.Print()

    return cuts

# use for double ratio
# scales yields by efficiency 
# but does not add in the uncertainty
def getEffWeights(sample):
    effWeights = ""
    
    effFile = ROOT.TFile.Open("../plots/histograms/effHists.root", "read")
    effFileSFm = ROOT.TFile.Open("../plots/histograms/SFcorrections.Muons.root", "read")
    effFileSFe = ROOT.TFile.Open("../plots/histograms/SFcorrections.Electrons.root", "read")

    if(sample=='zmm'):
        h_pur_m = effFile.Get("h_pur_m4")
        effWeights="((BTags==0)*"+str(h_pur_m.GetBinContent(1))+"+(BTags==1)*"+str(h_pur_m.GetBinContent(2))+"+(BTags==2)*"+str(h_pur_m.GetBinContent(3))+"+(BTags>=3)*"+str(h_pur_m.GetBinContent(4))+")"
    if(sample=='zee'):
        h_pur_e = effFile.Get("h_pur_e4")
        effWeights="((BTags==0)*"+str(h_pur_e.GetBinContent(1))+"+(BTags==1)*"+str(h_pur_e.GetBinContent(2))+"+(BTags==2)*"+str(h_pur_e.GetBinContent(3))+"+(BTags>=3)*"+str(h_pur_e.GetBinContent(4))+")"
    if(sample=='zll'):
        h_pur_m = effFile.Get("h_pur_m4")
        effWeights="(((BTags==0)*"+str(h_pur_m.GetBinContent(1))+"+(BTags==1)*"+str(h_pur_m.GetBinContent(2))+"+(BTags==2)*"+str(h_pur_m.GetBinContent(3))+"+(BTags>=3)*"+str(h_pur_m.GetBinContent(4))+")+((BTags==0)*"+str(h_pur_e.GetBinContent(1))+"+(BTags==1)*"+str(h_pur_e.GetBinContent(2))+"+(BTags==2)*"+str(h_pur_e.GetBinContent(3))+"+(BTags>=3)*"+str(h_pur_e.GetBinContent(4))+"))"
    if(sample=='photon'):
        h_pur_eb = effFile.Get("h_pur_eb")
        h_pur_ec = effFile.Get("h_pur_ec")
        h_frag = effFile.Get("h_frag1")
        effWeights = "("+str(h_frag.GetBinContent(1))+")*( (Photons_isEB==1)*((MHT<225)*"+str(h_pur_eb.GetBinContent(1))+"+(MHT>=225&&MHT<250)*"+str(h_pur_eb.GetBinContent(2))+"+(MHT>=250&&MHT<300)*"+str(h_pur_eb.GetBinContent(3))+"+(MHT>=300&&MHT<350)*"+str(h_pur_eb.GetBinContent(4))+"+(MHT>=350&&MHT<=500)*"+str(h_pur_eb.GetBinContent(5))+"+(MHT>=500)*"+str(h_pur_eb.GetBinContent(6))+") + (Photons_isEB==0)*((MHT<225)*"+str(h_pur_ec.GetBinContent(1))+"+(MHT>=225&&MHT<250)*"+str(h_pur_ec.GetBinContent(2))+"+(MHT>=250&&MHT<300)*"+str(h_pur_ec.GetBinContent(3))+"+(MHT>=300&&MHT<350)*"+str(h_pur_ec.GetBinContent(4))+"+(MHT>=350&&MHT<=500)*"+str(h_pur_ec.GetBinContent(5))+"+(MHT>=500)*"+str(h_pur_ec.GetBinContent(6))+"))"
    if(sample=='dymm'):
        h_trig_m = effFile.Get("h_trig_m1")
        h_SF_m = effFileSFm.Get("h_MHT")
        effWeights = "("+str(h_trig_m.GetBinContent(1))+")*((MHT<350)*"+str(1/h_SF_m.GetBinContent(1))+"+(MHT>=350&&MHT<450)*"+str(1/h_SF_m.GetBinContent(2))+"+(MHT>=450&&MHT<550)*"+str(1/h_SF_m.GetBinContent(3))+"+(MHT>=550&&MHT<650)*"+str(1/h_SF_m.GetBinContent(4))+"+(MHT>=650&&MHT<750)*"+str(1/h_SF_m.GetBinContent(5))+"+(MHT>=750&&MHT<850)*"+str(1/h_SF_m.GetBinContent(6))+"+(MHT>=850&&MHT<950)*"+str(1/h_SF_m.GetBinContent(7))+"+(MHT>=950&&MHT<1050)*"+str(1/h_SF_m.GetBinContent(8))+"+(MHT>=1050&&MHT<1150)*"+str(1/h_SF_m.GetBinContent(9))+"+(MHT>=1150)*"+str(1/h_SF_m.GetBinContent(10))+")"
    if(sample=='dyee'):
        h_trig_e = effFile.Get("h_trig_e1")
        h_SF_e = effFileSFe.Get("h_MHT")
        effWeights = "("+str(h_trig_e.GetBinContent(1))+")*((MHT<350)*"+str(1/h_SF_e.GetBinContent(1))+"+(MHT>=350&&MHT<450)*"+str(1/h_SF_e.GetBinContent(2))+"+(MHT>=450&&MHT<550)*"+str(1/h_SF_e.GetBinContent(3))+"+(MHT>=550&&MHT<650)*"+str(1/h_SF_e.GetBinContent(4))+"+(MHT>=650&&MHT<750)*"+str(1/h_SF_e.GetBinContent(5))+"+(MHT>=750&&MHT<850)*"+str(1/h_SF_e.GetBinContent(6))+"+(MHT>=850&&MHT<950)*"+str(1/h_SF_e.GetBinContent(7))+"+(MHT>=950&&MHT<1050)*"+str(1/h_SF_e.GetBinContent(8))+"+(MHT>=1050&&MHT<1150)*"+str(1/h_SF_e.GetBinContent(9))+"+(MHT>=1150)*"+str(1/h_SF_e.GetBinContent(10))+")"
    if(sample=='dyll'):
        h_trig_m = effFile.Get("h_trig_m1")
        h_SF_m = effFileSFm.Get("h_MHT")
        h_trig_e = effFile.Get("h_trig_e1")
        h_SF_e = effFileSFe.Get("h_MHT")
        effWeights = "(("+str(h_trig_m.GetBinContent(1))+")*((MHT<350)*"+str(1/h_SF_m.GetBinContent(1))+"+(MHT>=350&&MHT<450)*"+str(1/h_SF_m.GetBinContent(2))+"+(MHT>=450&&MHT<550)*"+str(1/h_SF_m.GetBinContent(3))+"+(MHT>=550&&MHT<650)*"+str(1/h_SF_m.GetBinContent(4))+"+(MHT>=650&&MHT<750)*"+str(1/h_SF_m.GetBinContent(5))+"+(MHT>=750&&MHT<850)*"+str(1/h_SF_m.GetBinContent(6))+"+(MHT>=850&&MHT<950)*"+str(1/h_SF_m.GetBinContent(7))+"+(MHT>=950&&MHT<1050)*"+str(1/h_SF_m.GetBinContent(8))+"+(MHT>=1050&&MHT<1150)*"+str(1/h_SF_m.GetBinContent(9))+"+(MHT>=1150)*"+str(1/h_SF_m.GetBinContent(10))+")+("+str(h_trig_e.GetBinContent(1))+")*((MHT<350)*"+str(1/h_SF_e.GetBinContent(1))+"+(MHT>=350&&MHT<450)*"+str(1/h_SF_e.GetBinContent(2))+"+(MHT>=450&&MHT<550)*"+str(1/h_SF_e.GetBinContent(3))+"+(MHT>=550&&MHT<650)*"+str(1/h_SF_e.GetBinContent(4))+"+(MHT>=650&&MHT<750)*"+str(1/h_SF_e.GetBinContent(5))+"+(MHT>=750&&MHT<850)*"+str(1/h_SF_e.GetBinContent(6))+"+(MHT>=850&&MHT<950)*"+str(1/h_SF_e.GetBinContent(7))+"+(MHT>=950&&MHT<1050)*"+str(1/h_SF_e.GetBinContent(8))+"+(MHT>=1050&&MHT<1150)*"+str(1/h_SF_e.GetBinContent(9))+"+(MHT>=1150)*"+str(1/h_SF_e.GetBinContent(10))+"))"
    if(sample=='gjets' or sample=='gjetsqcd' or sample=='ttgjets'):
        h_SF_g = effFile.Get("h_SF_g1")
        h_trig_eb = effFile.Get("h_trig_eb")
        h_trig_ec = effFile.Get("h_trig_ec")
        effWeights = "("+str(h_SF_g.GetBinContent(1))+")*( (Photons_isEB==1)*((MHT<300)*"+str(h_trig_eb.GetBinContent(1))+"+(MHT>=300&&MHT<350)*"+str(h_trig_eb.GetBinContent(2))+"+(MHT>=350&&MHT<500)*"+str(h_trig_eb.GetBinContent(3))+"+(MHT>=500&&MHT<750)*"+str(h_trig_eb.GetBinContent(4))+"+(MHT>=750)*"+str(h_trig_eb.GetBinContent(5))+")+(Photons_isEB==0)*((MHT<300)*"+str(h_trig_ec.GetBinContent(1))+"+(MHT>=300&&MHT<350)*"+str(h_trig_ec.GetBinContent(2))+"+(MHT>=350&&MHT<500)*"+str(h_trig_ec.GetBinContent(3))+"+(MHT>=500&&MHT<750)*"+str(h_trig_ec.GetBinContent(4))+"+(MHT>=750)*"+str(h_trig_ec.GetBinContent(5))+"))"

    return effWeights

def getDist(sample, dist, distRange=None, nBins=None, doVarBinning=None, binning=None, extraCuts=None, nJetBin=None, bJetBin=None, kinBin=None, dphiCut=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, doLumi=None, removeZkfactor=None, doNoisy=None, doProof=None, treeLoc=None, applyEffs=None, treeName=None, applySF=None, applyPuWeight=None):

    ###############################################################################        
    # set default parameters
    ###############################################################################        

    ###############################################################################        
    ## distribution binning, try to set as many defaults as possible ###
    ## if plotting something not below, you must give binning ###
    if(dist=='ZCandidates.M()'  or "min(ZCandidates.M()" in dist):
        if(applyMassCut==None):
            applyMassCut=False
        if(applyMassCut==False):
            if(distRange==None):
                distRange=[60.,120.]
            if(binning==None):
                binning=[60.,65.,70.,75.,80.,85.,90.,95.,100.,105.,110.,115.,120.]
        if(applyMassCut==True):
            if(distRange==None):
                distRange=[76.188,106.188]
            if(binning==None):
                binning=[75.,80.,85.,90.,95.,100.,105.]
        if(nBins==None):
            nBins=(distRange[1]-distRange[0])/2
        if(doVarBinning==None):
            doVarBinning=False
    if(applyMassCut==None):
        applyMassCut=True
    if(applyPtCut==None):
        applyPtCut=True
    if(applyHTCut==None):
        applyHTCut=True
    if(applyMHTCut==None):
        applyMHTCut=True
    if(applyNJetsCut==None):
        applyNJetsCut=False
    if(dist=='HT' or "min(HT" in dist):
        if(doVarBinning==None):
            doVarBinning=True
        if(applyHTCut==True):
            if(distRange==None):
                distRange=[250,2500]
            if(nBins==None):
                nBins=14
            if(binning==None):
                binning=[300.,350.,400.,450,500.,550.,650.,750.,850.,950.,1050,1150.,1250,1750.]
        elif(applyHTCut==False):
            if(distRange==None):
                distRange=[200,2500]
            if(nBins==None):
                nBins=14
            if(binning==None):
                binning=[250.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200.,1300.,1500.,2500.]
    if(dist=='MHT' or "min(MHT" in dist):
        if(doVarBinning==None):
            doVarBinning=True
        if(applyMHTCut==True):
            if(distRange==None):
                distRange=[250,1000]
            if(nBins==None):
                nBins=7
            if(binning==None):
                binning=[300.,350.,400.,450.,500.,600.,775.,1000.]
            else:
                nBins = len(binning)
        if(applyMHTCut==False):
            if(distRange==None):
                distRange=[100,1000]
            if(nBins==None):
                nBins=8
            if(binning==None):
                binning=[250.,300.,400.,500.,600.,775.,1000.]
            else:
                nBins = len(binning)
    if(dist=='ZCandidates.Pt()'  or "min(ZCandidates.Pt()" in dist):
        if(doVarBinning==None):
            doVarBinning=True
        if(distRange==None):
            distRange=[100,1200]
        if(nBins==None):
            nBins=14
        if(binning==None):
            binning=[0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,600.,800.,1200.]
    if(dist=='Photons.Pt()' or "min(Photons.Pt()" in dist):
        if(doVarBinning==None):
            doVarBinning=True
        if(distRange==None):
            distRange=[0,1000]
        if(nBins==None):
            nBins=13
        if(binning==None):
            binning=[0,50.,100.,150.,200.,250.,300.,350.,400.,500.,600.,775.,1000.]
    if(dist=='NJets' or "min(NJets" in dist):
        if(applyNJetsCut==True):
            if(distRange==None):
                distRange=[2.5,9.5]
            if(binning==None):
                binning=[1.5,6.5,8.5,12.5]
        elif(applyNJetsCut==False):
            if(distRange==None):
                distRange=[1.5,9.5]
            if(binning==None):
                binning=[1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]
        if(nBins==None):
            nBins=(distRange[1]-distRange[0])
        if(doVarBinning==None):
            doVarBinning=False
    if(dist=='BTags' or "min(BTags" in dist):
        if(distRange==None):
            distRange=[-0.5,3.5]
        if(binning==None):
            binning=[-0.5,0.5,1.5,2.5,3.5]
        if(nBins==None):
            nBins=(distRange[1]-distRange[0])
        if(doVarBinning==None):
            doVarBinning=False
    if('DeltaPhi' in dist):
        if(distRange==None):
            distRange=[0,3.14]
        if(binning==None):
            binning=[0.,0.5,1.0,1.5,2.0,2.5,3.0,3.14]
        if(nBins==None):
            nBins=50
        if(doVarBinning==None):
            doVarBinning=False
    if(dist=='Muons[0].Pt()' or dist=='Muons[1].Pt()' or ("min(Muons[" in dist and "].Pt()" in dist)):
        if(doVarBinning==None):
            doVarBinning=True
        if(distRange==None):
            distRange=[100,1200]
        if(nBins==None):
            nBins=14
        if(binning==None):
            binning=[0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,600.,800.,1200.]

    ## end distribution binning 
    ###############################################################################        
    if(doLumi==None):
        doLumi = 24.5
    if(removeZkfactor==None and (sample=='zinv' or sample=='dymm' or sample=='dyee' or sample=='dyll')):
       removeZkfactor=True
    elif(removeZkfactor==None):
       removeZkfactor=False
    if(doNoisy==None):
        doNoisy = False
    if(doProof==None):
        doProof = False
    if(extraCuts==None):
        extraCuts=""
    if(nJetBin==None):
        nJetBin=-1
    if(bJetBin==None):
        bJetBin=-1
    if(kinBin==None):
        kinBin=-1
    if(dphiCut==None):
        if('LDP' in sample):
            dphiCut = 'ldp'
        elif('IDP' in sample):
            dphiCut = 'none'
        else:
            dphiCut='sig'
    if(dphiCut=='ldp'):
        if('LDP' not in sample):
            sample+='LDP'
    if(treeLoc==None):
        treeLoc = "/home/ww/work/data/lpcTrees/Skims/Run2ProductionV10"
    if(applyEffs==None):
        applyEffs=False
    if(treeName==None):
        treeName = "tree"
    if(applySF==None):
        applySF=False

    nBins=int(nBins)

    ###############################################################################        
    # end set default parameters
    ###############################################################################        

    chain = getChain(sample, doLumi, treeName, doProof, treeLoc, removeZkfactor)

    #############################################################################
    # define unique root TH1F 
    #############################################################################

    plotlabel = "plot_"+str(sample)+"_"+str(re.sub('[^a-zA-Z0-9-_*.\.]','', dist))
    
    hIter = 1
    while(type(ROOT.gROOT.FindObject(plotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1
        
    plotlabel+="_"+str(hIter)

    plot = ROOT.TH1F(plotlabel,plotlabel,nBins,distRange[0],distRange[1])  
    if(doVarBinning):
        plot.Delete()
        binning_f = array('f',binning)
        plot = ROOT.TH1F(plotlabel,plotlabel,len(binning)-1,binning_f)

    #################################################################################
    # end define unique root TH1F 
    #################################################################################

    #################################################################################
    # determine sample cut code for getCuts, trigger, and additional filters/cuts
    ###############################################################################

    cutCode = {
        'sig' : 'sig',
        'zinv' : 'sig',
        'top' : 'sig',
        'qcd' : 'sig',
        'wjets' : 'sig',
        'other' : 'sig',
        'ttzvv' : 'ttz',
        'zmm' : 'zmm',
        'zmm20' : 'zmm',
        'dymm' : 'zmm',
        'ttmm' : 'zmm',
        'ttzmm' : 'zmm',
        'dibosonmm' : 'zmm',
        'tribosonmm' : 'zmm',
        'zee' : 'zee',
        'zee20' : 'zee',
        'dyee' : 'zee',
        'ttzee' : 'zee',
        'ttee' : 'zee',
        'dibosonee' : 'zee',
        'tribosonee' : 'zee',
        'zll' : 'zll',
        'zll20' : 'zll',
        'dyll' : 'zll',
        'ttzll' : 'zll',
        'ttll' : 'zll',
        'dibosonll' : 'zll',
        'photon' : 'photon',
        'photon20' : 'photon',
        'gjets' : 'photon',
        'gjetsold' : 'photon',
        'ttgjets' : 'photon',
        'gjetsqcd' : 'photonqcd',
    }
    if(sample == 'zmm'):
        trig = (20, 21, 27, 16)
    elif(sample == 'zee'):
        trig = (5, 6, 11, 3)
    elif(sample == 'zll'):
        trig = (20, 21, 27, 16, 5, 6, 11, 3)
    elif(sample == 'photon'):
        trig = 49
    elif(sample == 'sig'):
        trig = (29,33)
    else:
        trig = -1

    ExtraFilters = ROOT.TCut("")
    extraCuts = ROOT.TCut(extraCuts)
    # data filters updates came after prod V5
    if(sample=='sig' or sample=='zmm' or sample=='zee' or sample=='photon'):
        #ExtraFilters = "(CSCTightHaloFilterUpdate==1)&&(BadTrackFilterUpdate==1)&&(BadMuonFilterUpdate==1)"
        applySF = False

    extraCuts+=ExtraFilters

    if(('IDP' in sample) or ('LDP' in sample)):
        sample = sample[:-3]

    cuts = getCuts(cutCode[sample], nJetBin, bJetBin, kinBin, trig, dphiCut, extraCuts, applyMassCut, applyPtCut, applyHTCut, applyMHTCut, applyNJetsCut, applySF, applyPuWeight=applyPuWeight)

    #################################################################################
    # determine sample cut code for getCuts, trigger, and additional filters/cuts
    ###############################################################################

    if(applyEffs):
        cuts*=getEffWeights(sample)

    chain.Project(str(plotlabel),dist,str(cuts))

    Min = plot.GetBinContent(plot.FindLastBinAbove(0))

    #plot.SetMinimum(Min*0.1)

    if('20' in sample):
        plot.Sumw2(False)

    return plot


def getRemovedBins(njRange=None, nbRange=None, kinRange=None, njSplit=None):

    if(njSplit==None):
        njSplit=False
    if(njSplit==True):
        if(njRange==None):
            njRange = range(1,9)
    if(njRange==None):
        njRange = range(1,6)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    if(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)

    subtractBins = 0
    avoidBins = []

    if(1 in njRange and 3 in nbRange):
        subtractBins+=len(kinRange)
        avoidBins += range(len(kinRange)*3+1,len(kinRange)*3+1*len(kinRange)+1)
    
    kinSkip = [1,4,11]
    njSkip = [4,5]
    if(njSplit==True):
        njSkip = range(6,9)
    for avoidNJ in njSkip:
        for nb in nbRange:
            for avoidKin in kinSkip:
                if(avoidKin in kinRange and avoidNJ in njRange):
                    avoidBins+=[len(kinRange)*len(nbRange)*(avoidNJ-1)+len(kinRange)*nb+avoidKin]
                    subtractBins+=1

    return (subtractBins,avoidBins)

# code map
# 0 sig data 
# 1 Zmm data
# 2 Zee data
# 3 Zll data
# 4 photon data
# 10 Zinv MC
# 11 DY(mm) MC
# 12 DY(ee) MC
# 13 DY(ll) MC
# 14 photon MC
# 20 Zinv MC + ttZ MC 
# 21 DY(mm) + tt(mm) + ttZ(mm) MC
# 22 DY(ee) + tt(ee) + ttZ(ee) MC
# 23 DY(ll) + tt(ll) + ttZ(ll) MC
# 24 photon + QCD MC
# 25 ttZ->nunu
# 30 All bg MC

def getHist(sample, njRange=None, nbRange=None, kinRange=None, doLumi=None, removeZkfactor=None, doNoisy=None, doProof=None, dphiCut=None, treeLoc=None, applyEffs=None, treeName=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, applySF=None, extraCuts=None, doStandalone=None, njSplit=None, applyPuWeight=None):

    # timer
    start_time = time.time()

    ###############################################################################        
    # set default parameters
    ###############################################################################        
    
    if(njSplit==None):
        njSplit=False
    if(njSplit==True):
        if(njRange==None):
            njRange = range(1,9)
    if(njRange==None):
        njRange = range(1,6)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    if(nbRange==range(4)):
        if(doStandalone==None):
            doStandalone=True
    elif(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(doStandalone==None):
        doStandalone=False

    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)
    if(doLumi==None):
        doLumi = 24.5
    if(removeZkfactor==None and (sample=='zinv' or sample=='dymm' or sample=='dyee' or sample=='dyll')):
       removeZkfactor=True
    elif(removeZkfactor==None):
       removeZkfactor=False
    if(doNoisy==None):
        doNoisy = False
    if(doProof==None):
        doProof = False
    if(dphiCut==None):
        if('LDP' in sample):
            dphiCut = 'ldp'
        elif('IDP' in sample):
            dphiCut = 'none'
        else:
            dphiCut='sig'
    if(dphiCut=='ldp'):
        if('LDP' not in sample):
            sample+='LDP'
    if(treeLoc==None):
        treeLoc = "/home/ww/work/data/lpcTrees/Skims/Run2ProductionV10"
    if(applyEffs==None):
        applEffs=False
    if(treeName==None):
        treeName = "tree"
    if(applySF==None):
        applySF=False
    if(applyPuWeight==None):
        applyPuWeight=False
    if(extraCuts==None):
        extraCuts=""
    ###############################################################################        
    # end set default parameters
    ###############################################################################        

    #################################################################################
    # chain samples according to code
    #################################################################################

    chain = getChain(sample, doLumi, treeName, doProof, treeLoc, removeZkfactor)

    #################################################################################
    # end chain samples according to code
    #################################################################################

    #############################################################################
    # define unique root TH1F 
    #############################################################################

    removedBins = getRemovedBins(njRange, nbRange, kinRange, njSplit)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]
    
    nbins = len(kinRange)*len(njRange)*len(nbRange)-subtractBins
    
    plotlabel = "plot_"+str(sample)+"_nj"+str(len(njRange))+"_nb"+str(len(nbRange))+"_kin"+str(len(kinRange))
    
    hIter = 1
    while(type(ROOT.gROOT.FindObject(plotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1
        
    plotlabel+="_"+str(hIter)
  
    plot = ROOT.TH1F(plotlabel,plotlabel,nbins,0.5,nbins+0.5)

    #################################################################################
    # end define unique root TH1F 
    #################################################################################

    #################################################################################
    # determine sample cut code for getCuts, trigger, and additional filters/cuts
    ###############################################################################

    cutCode = {
        'sig' : 'sig',
        'zinv' : 'sig',
        'top' : 'sig',
        'qcd' : 'sig',
        'qcdIDP' : 'sig',
        'zinvIDP' : 'sig',
        'wjets' : 'sig',
        'other' : 'sig',
        'ttzvv' : 'ttz',
        'zmm' : 'zmm',
        'zmm20' : 'zmm',
        'dymm' : 'zmm',
        'ttmm' : 'zmm',
        'ttzmm' : 'zmm',
        'dibosonmm' : 'zmm',
        'tribosonmm' : 'zmm',
        'zee' : 'zee',
        'zee20' : 'zee',
        'dyee' : 'zee',
        'ttzee' : 'zee',
        'ttee' : 'zee',
        'dibosonee' : 'zee',
        'tribosonee' : 'zee',
        'zll' : 'zll',
        'zll20' : 'zll',
        'dyll' : 'zll',
        'ttzll' : 'zll',
        'ttll' : 'zll',
        'dibosonll' : 'zll',
        'photon' : 'photon',
        'photon20' : 'photon',
        'gjets' : 'photon',
        'gjetsold' : 'photon',
        'ttgjets' : 'photon',
        'gjetsqcd' : 'photonqcd',
    }

    if(sample == 'zmm'):
        trig = (20, 21, 27, 16)
    elif(sample == 'zee'):
        trig = (5, 6, 11, 3)
    elif(sample == 'zll'):
        trig = (20, 21, 27, 16, 5, 6, 11, 3)
    elif(sample == 'photon'):
        trig = 49
    elif(sample == 'sig'):
        trig = (29,33)
    else:
        trig = -1

    ExtraFilters = extraCuts

    # data filters updates came after prod V5
    if(sample=='sig' or sample=='zmm' or sample=='zee' or sample=='photon'):
        #ExtraFilters = "(CSCTightHaloFilterUpdate==1)&&(BadTrackFilterUpdate==1)&&(BadMuonFilterUpdate==1)"
        applySF = False

    #################################################################################
    # end determine sample cut code for getCuts, trigger, and additional filters/cuts
    #################################################################################

    #################################################################################
    # loop and fill bins with observed number of events
    #################################################################################

    # define histogram used to count events

    ecplotlabel = "eventCount"
    hIter = 1
    while(type(ROOT.gROOT.FindObject(ecplotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1

    ecplotlabel+="_"+str(hIter)

    EventCount = ROOT.TH1F(ecplotlabel,ecplotlabel,1,0,1000)

    if(('IDP' in sample) or ('LDP' in sample)):
        sample = sample[:-3]

    Bin = 1
    uncutBin = 0
    for nj in njRange:
        for nb in nbRange:
            for kin in kinRange:
                uncutBin+=1
                if(uncutBin in avoidBins):
                    continue

                cuts = getCuts(cutCode[sample], nj, nb, kin, trig, dphiCut, ExtraFilters, applyMassCut, applyPtCut, applyHTCut, applyMHTCut, applyNJetsCut, applySF,njSplit=njSplit,applyPuWeight=applyPuWeight)
                if(applyEffs):
                    cuts*=getEffWeights(sample)

                chain.Project(ecplotlabel,"nAllVertices",str(cuts))

                if(doNoisy):
                    print("nj: "+str(nj)+" nb: "+str(nb)+" kin: "+str(kin))
                    print str(EventCount.GetBinContent(1))+" +- "+str(EventCount.GetBinError(1))

                plot.Fill(Bin, EventCount.GetBinContent(1))
                plot.SetBinError(Bin, EventCount.GetBinError(1))
                Bin+=1

    #################################################################################
    # end loop and fill bins with observed number of events
    #################################################################################

    runtime = time.time()-start_time

    if(doNoisy):
        print "runtime "+str(runtime)

    return plot

def getFudgeHist(val=None, error=None, njRange=None, njSplit=None, nbRange=None, kinRange=None, njFudge=None, nbFudge=None, kinFudge=None):
 
    if(val==None):
        val = 1.
    if(error==None):
        error= 1.
    if(njSplit==None):
        njSplit=True
    if(njSplit==True):
        if(njRange==None):
            njRange = range(1,9)
    if(njRange==None):
        njRange = range(1,6)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    elif(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)
    if(njFudge==None):
        njFudge=njRange
    elif(type(njFudge) is int):
        njFudge=range(njFudge,njFudge+1)
    if(nbFudge==None):
        nbFudge=nbRange
    elif(type(nbFudge) is int):
        nbFudge = range(nbFudge,nbFudge+1)
    if(kinFudge==None):
        kinFudge=kinRange
    elif(type(kinFudge) is int):
        kinFudge = range(kinFudge,kinFudge+1)
    #############################################################################
    # define unique root TH1F 
    #############################################################################
    removedBins = getRemovedBins(njRange, nbRange, kinRange, njSplit)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]    

    nbins = len(kinRange)*len(njRange)*len(nbRange)-subtractBins

    plotlabel = "fudgeplot_nj"+str(len(njRange))+"_nb"+str(len(nbRange))+"_kin"+str(len(kinRange))
    
    hIter = 1
    while(type(ROOT.gROOT.FindObject(plotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1
        
    plotlabel+="_"+str(hIter)
  
    hFudge = ROOT.TH1F(plotlabel,plotlabel,nbins,0.5,nbins+0.5)

    #################################################################################
    # end define unique root TH1F 
    #################################################################################

    Bin=1
    uncutBin=0
    for nj in njRange:
        for nb in nbRange:
            for kin in kinRange:
                uncutBin+=1
                if(uncutBin in avoidBins):
                    continue

                if((nj in njFudge) and (nb in nbFudge) and (kin in kinFudge)):
                    
                    hFudge.SetBinContent(Bin,val)
                    hFudge.SetBinError(Bin,error)

                else:

                    hFudge.SetBinContent(Bin,0)
                    hFudge.SetBinError(Bin,0)

                Bin+=1

    return hFudge

def get0bPrediction(samples, njRange=None, njSplit=None, nbRange=None, kinRange=None, doLumi=None, removeZkfactor=None, doNoisy=None, doProof=None, dphiCut=None, treeLoc=None, applyEffs=None, treeName=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, applySF=None,applyPuWeight=None):

    if(njSplit==None):
        njSplit=True
    if(njSplit==True):
        if(njRange==None):
            njRange = range(1,9)
    if(njRange==None):
        njRange = range(1,6)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    elif(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)
    if(type(samples) is not list):
        samples = [samples]
    #############################################################################
    # define unique root TH1F 
    #############################################################################

    removedBins = getRemovedBins(njRange, nbRange, kinRange, njSplit)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]
    
    removedBins0b = getRemovedBins(njRange, 0, kinRange, njSplit)
    avoidBins0b = removedBins0b[1]
    
    nbins = len(kinRange)*len(njRange)*len(nbRange)-subtractBins

    plotlabel = "extrapplot_"+str(samples[0])+"_nj"+str(len(njRange))+"_nb"+str(len(nbRange))+"_kin"+str(len(kinRange))
    
    hIter = 1
    while(type(ROOT.gROOT.FindObject(plotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1
        
    plotlabel+="_"+str(hIter)
  
    hPred = ROOT.TH1F(plotlabel,plotlabel,nbins,0.5,nbins+0.5)

    #################################################################################
    # end define unique root TH1F 
    #################################################################################

    hlist = []
    for sample in samples:
        hlist.append(getHist(sample, njRange, 0, kinRange, doLumi, removeZkfactor, doNoisy, doProof, dphiCut, treeLoc, applyEffs, treeName, applyMassCut, applyPtCut, applyHTCut, applyMHTCut, applyNJetsCut, applySF, njSplit=njSplit,applyPuWeight=applyPuWeight))

    h1 = hlist[0].Clone()
    for i in range(1,len(samples)):
        h1.Add(hlist[i])
    
    Bin=1
    uncutBin=0
    for nj in njRange:
        for nb in nbRange:
            for kin in kinRange:
                uncutBin+=1
                if(uncutBin in avoidBins):
                    continue

                kin = max(kin,1)
                nj = max(nj,1)

                bin0b = kin + (nj-1)*len(kinRange)
                subtractBins0b = sum(i<=bin0b for i in avoidBins0b)
                bin0b-=subtractBins0b

                pred0b = h1.GetBinContent(bin0b)
                epred0b = h1.GetBinError(bin0b)

                hPred.SetBinContent(Bin,pred0b)
                hPred.SetBinError(Bin,epred0b)

                Bin += 1


    return hPred


def getExtrapolation(samples, njFactorRange=None, doFactorization=None, njRange=None, njSplit=None, nbRange=None, kinRange=None, doLumi=None, removeZkfactor=None, doNoisy=None, doProof=None, dphiCut=None, treeLoc=None, applyEffs=None, treeName=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, applySF=None, njExtrapRange=None, applyPuWeight=None):

    if(doFactorization==None):
        doFactorization=False
    if(njSplit==None):
        njSplit=True
    if(doFactorization==True):
        if(njFactorRange==None):
            if(njSplit==True):
                njFactorRange=range(6,8)
            else:
                njFactorRange=range(4,5)
        if(njExtrapRange==None):
            if(njSplit==True):
                njExtrapRange=range(8,9)
            else:
                njExtrapRange=range(5,6)
    if(njSplit==True):
        if(njRange==None):
            njRange = range(1,9)
    if(njRange==None):
        njRange = range(1,6)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    elif(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)
        
    if(type(samples) is not list):
        samples = [samples]
    #############################################################################
    # define unique root TH1F 
    #############################################################################

    removedBins = getRemovedBins(njRange, nbRange, kinRange, njSplit)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]

    nbins = len(kinRange)*len(njRange)*len(nbRange)-subtractBins

    plotlabel = "extrapplot_"+str(samples[0])+"_nj"+str(len(njRange))+"_nb"+str(len(nbRange))+"_kin"+str(len(kinRange))
    
    hIter = 1
    while(type(ROOT.gROOT.FindObject(plotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1
        
    plotlabel+="_"+str(hIter)
  
    hExtrap = ROOT.TH1F(plotlabel,plotlabel,nbins,0.5,nbins+0.5)

    #################################################################################
    # end define unique root TH1F 
    #################################################################################

    hlist = []
    for sample in samples:
        hlist.append(getHist(sample, njRange, nbRange, -1, doLumi, removeZkfactor, doNoisy, doProof, dphiCut, treeLoc, applyEffs, treeName, applyMassCut, applyPtCut, applyHTCut, applyMHTCut, applyNJetsCut, applySF, njSplit=njSplit, applyPuWeight=applyPuWeight))

    h1 = hlist[0].Clone()
    for i in range(1,len(samples)):
        h1.Add(hlist[i])

    Bin = 1
    uncutBin = 0
    binSubtract = 0
    binNb = [0]
    bin0b = [0]
    for nj in njRange:
        for nb in nbRange:
            if(not(nj==1 and nb==3)):
                binNb[0]+=1
                if(nb==0):
                    bin0b[0] = binNb[0]
            else:
                binSubtract+=1
            for kin in kinRange:
                uncutBin+=1
                if(nj==1 and nb==3):
                    continue
                if(uncutBin in avoidBins):
                    continue
                kin = max(kin,1)
                nj = max(nj,1)

                if(doFactorization==True):
                    if((njSplit and nj in njExtrapRange) or (not njSplit and nj in njFactorRange or nj in njExtrapRange)):
                        njExtrap = njFactorRange
                        binNb = []
                        bin0b = []
                        for njE in njExtrap:
                            binNb.append(1 + nb + (njE-1)*len(nbRange)-binSubtract)
                            bin0b.append(1 + (njE-1)*len(nbRange)-binSubtract)
                    else:
                        bin0b = []
                        binNb = []
                        #Bin = kin + nb*len(kinRange) + (nj-1)*len(nbRange)*len(kinRange)
                        binNb.append(1 + nb + (nj-1)*len(nbRange)-binSubtract)
                        bin0b.append(1 + (nj-1)*len(nbRange)-binSubtract)
                
                b0 = 0
                eb0 = 0.
                for Ob in bin0b:
                    b0+=h1.GetBinContent(Ob)
                    eb0 += h1.GetBinError(Ob)
                if(b0==0):
                    b0 = 1.
                    eb0 = 1.
                else:
                    eb0 = eb0/b0

                b = 0.
                eb = 0.
                for Nb in binNb:
                    b += h1.GetBinContent(Nb)
                    eb = h1.GetBinError(Nb)
                if(b>0.):
                    eb = eb/b
                else:
                    eb = 0
                
                hExtrap.SetBinContent(Bin,b/b0)
                eExtrap = hExtrap.GetBinContent(Bin)*ROOT.TMath.Sqrt(eb**2+eb0**2)

                if(round(b/b0,4)==1):
                    eExtrap = 0.

                hExtrap.SetBinError(Bin,eExtrap)

                Bin +=1
                
    return hExtrap

def getBinomialCorrection(prob=None, njFactorRange=None, njRange=None, njSplit=None, nbRange=None, kinRange=None,njData=None):

    if(prob==None):
        prob = 0.0600
    if(njSplit==None):
        njSplit = True
    if(njFactorRange==None):
        if(njSplit==True):
            njFactorRange=range(6,8)
        else:
            njFactorRange=range(4,5)
    if(njSplit==True):
        njRange = range(1,9)
    if(njRange==None):
        njRange = range(1,6)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    elif(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)
    if(njData==None):
        if(njSplit==True):
            njData=range(1,8)
        else:
            njData=range(1,5)
    #############################################################################
    # define unique root TH1F 
    #############################################################################

    removedBins = getRemovedBins(njRange, nbRange, kinRange, njSplit)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]

    nbins = len(kinRange)*len(njRange)*len(nbRange)-subtractBins
    
    plotlabel = "binomplot_nj"+str(len(njRange))+"_nb"+str(len(nbRange))+"_kin"+str(len(kinRange))
    
    hIter = 1
    while(type(ROOT.gROOT.FindObject(plotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1
        
    plotlabel+="_"+str(hIter)
  
    hCorr = ROOT.TH1F(plotlabel,plotlabel,nbins,0.5,nbins+0.5)

    #################################################################################
    # end define unique root TH1F 
    #################################################################################


    h_nj = getDist('zinv','NJets',distRange=[2.5,12.5])
    NJets4o3 = h_nj.GetBinContent(2)/ h_nj.GetBinContent(1)
    NJets5o3 = h_nj.GetBinContent(3)/ h_nj.GetBinContent(1)
    NJets6o3 = h_nj.GetBinContent(4)/ h_nj.GetBinContent(1)
    NJets6o5 = h_nj.GetBinContent(4)/ h_nj.GetBinContent(3)
    NJets8o7 = h_nj.GetBinContent(6)/ h_nj.GetBinContent(5)
    NJets10o9 = h_nj.GetBinContent(8)/ h_nj.GetBinContent(7)
    NJets11o9 = h_nj.GetBinContent(9)/ h_nj.GetBinContent(7)
    NJets12o9 = h_nj.GetBinContent(10)/ h_nj.GetBinContent(7)

    uncutBin = 0
    Bin = 1
    for nj in njRange:
        for nb in nbRange:
            for kin in kinRange:
                uncutBin+=1
                if(uncutBin in avoidBins):
                    continue

                kin = max(kin,1)
                nj = max(nj,1)

                J = nj+2

                if(njSplit==True):
                    #corrTerm = ROOT.TMath.Binomial(J,nb)*(1-prob)**(J-3)/(ROOT.TMath.Binomial(3,nb))
                    #normTerm = (1+NJets4o3*(1-prob)+NJets5o3*(1-prob)**2+NJets6o3*(1-prob)**3)/(ROOT.TMath.Binomial(3,nb)+NJets4o3*ROOT.TMath.Binomial(4,nb)*(1-prob)+NJets5o3*ROOT.TMath.Binomial(5,nb)*(1-prob)**2+NJets6o3*ROOT.TMath.Binomial(6,nb)*(1-prob)**3)
                    #normTerm = (1+NJets6o5*(1-prob))/(ROOT.TMath.Binomial(5,nb)+NJets6o5*ROOT.TMath.Binomial(6,nb)*(1-prob))
                    normTerm = (1+NJets8o7*(1-prob))/(ROOT.TMath.Binomial(7,nb)+NJets8o7*ROOT.TMath.Binomial(8,nb)*(1-prob))
                    if(nj==8):
                        corrTerm = (ROOT.TMath.Binomial(9,nb)+NJets10o9*ROOT.TMath.Binomial(10,nb)*(1-prob)+NJets11o9*ROOT.TMath.Binomial(11,nb)*(1-prob)**2+NJets12o9*ROOT.TMath.Binomial(12,nb)*(1-prob)**3)*normTerm/(1+NJets10o9*(1-prob)+NJets11o9*(1-prob)**2+NJets12o9*(1-prob)**3)
                    else:
                        corrTerm = ROOT.TMath.Binomial(J,nb)*normTerm
                    if(nb==0):
                        corrTerm = 1.
                else:
                    normTerm = (1+NJets8o7*(1-prob))/(ROOT.TMath.Binomial(7,nb)+NJets8o7*ROOT.TMath.Binomial(8,nb)*(1-prob))

                    # if(nj==3):
                    #     corrTerm = (ROOT.TMath.Binomial(7,nb)+NJets8o7*ROOT.TMath.Binomial(8,nb)*(1-prob))*normTerm/(1+NJets8o7*(1-prob))
                    if(nj==5):
                        corrTerm = (ROOT.TMath.Binomial(9,nb)+NJets10o9*ROOT.TMath.Binomial(10,nb)*(1-prob)+NJets11o9*ROOT.TMath.Binomial(11,nb)*(1-prob)**2+NJets12o9*ROOT.TMath.Binomial(12,nb)*(1-prob)**3)*normTerm/(1+NJets10o9*(1-prob)+NJets11o9*(1-prob)**2+NJets12o9*(1-prob)**3)
                    else:
                        corrTerm = 1.0

                    if(nb==0):
                        corrTerm = 1.0

                if(nj in njData):
                    corrTerm=1.
                hCorr.SetBinContent(Bin,corrTerm)
                hCorr.SetBinError(Bin,0.) # get this from toy study if needed
                Bin+=1

    return hCorr

def getMonteCarloCorrection(sample, njFactorRange=None, njRange=None, njSplit=None, nbRange=None, kinRange=None, doLumi=None, removeZkfactor=None, doNoisy=None, doProof=None, dphiCut=None, treeLoc=None, applyEffs=None, treeName=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, applySF=None, applyPuWeight=None):

    if(njSplit==None):
        njSplit=True
    if(njSplit==True):
        if(njRange==None):
            njRange = range(1,9)
    if(njFactorRange==None):
        if(njSplit==True):
            njFactorRange=range(6,8)
        else:
            njFactorRange=range(4,5)
    if(njRange==None):
        njRange = range(1,6)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    elif(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)

    h1 = getExtrapolation(sample, njFactorRange, False, njRange, njSplit, nbRange, kinRange, doLumi, removeZkfactor, doNoisy, doProof, dphiCut, treeLoc, applyEffs, treeName, applyMassCut, applyPtCut, applyHTCut, applyMHTCut, applyNJetsCut, applySF, applyPuWeight=applyPuWeight)
    #h1 =  getExtrapolation(sample, False)
    h2 = getExtrapolation(sample, njFactorRange, True, njRange, njSplit, nbRange, kinRange, doLumi, removeZkfactor, doNoisy, doProof, dphiCut, treeLoc, applyEffs, treeName, applyMassCut, applyPtCut, applyHTCut, applyMHTCut, applyNJetsCut, applySF, applyPuWeight=applyPuWeight)
    #h2 = getExtrapolation(sample, True)

    h1.Divide(h2)

    removedBins = getRemovedBins(njRange, nbRange, kinRange, njSplit)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]

    uncutBin=0
    Bin = 1
    for nj in njRange:
        for nb in nbRange:
            for kin in kinRange:
                uncutBin+=1
                if(uncutBin in avoidBins):
                    continue
                if(nj<=max(njFactorRange)):
                    h1.SetBinError(Bin,0.)
                if(nb==0):
                    h1.SetBinError(Bin,0.)
                Bin+=1

    return h1

def mergeHist(hist,njRange=None,nbRange=None,kinRange=None):

    if(njRange==None):
        njRange = range(1,6)
    if(nbRange==None):
        nbRange = range(4)
    if(kinRange==None):
        kinRange = range(1,11)

    plotlabel = hist.GetTitle()+'_merge'
    hIter = 1
    while(type(ROOT.gROOT.FindObject(plotlabel+"_"+str(hIter)))==ROOT.TH1F):
        hIter+=1
        
    plotlabel+="_"+str(hIter)

    removedBins = getRemovedBins(njRange, nbRange, kinRange)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]

    removedBinsSplit =  getRemovedBins(nbRange=nbRange, kinRange=kinRange, njSplit=True)
    avoidBinsSplit = removedBinsSplit[1]

    nbins = len(kinRange)*len(njRange)*len(nbRange)-subtractBins

    h1 = ROOT.TH1F(plotlabel,plotlabel,nbins,0.5,nbins+0.5)
    
    splitDict = {
        1: [1],
        2: [2,3],
        3: [4,5],
        4: [6,7],
        5: [8],
    }

    uncutBin = 0
    Bin = 1
    for nj in njRange:
        for nb in nbRange:
            for kin in kinRange:
                uncutBin+=1
                if(uncutBin in avoidBins):
                    continue

                nEvents = 0.
                eEvents = 0.
                for j in splitDict[nj]:
                    splitBin = kin + nb*len(kinRange) + (j-1)*len(nbRange)*len(kinRange)
                    subtractBinsSplit = sum(i<splitBin for i in avoidBinsSplit)
                    splitBin-=subtractBinsSplit
                    nEvents+=hist.GetBinContent(splitBin)
                    eEvents+=hist.GetBinError(splitBin)**2
                eEvents = eEvents**0.5
                h1.SetBinContent(Bin,nEvents)
                h1.SetBinError(Bin,eEvents)
                Bin+=1
    return h1


def setError(hist, error=None, njRange=None, nbRange=None, kinRange=None, njSet=None, nbSet=None, kinSet=None):

    h1 = hist.Clone()

    if(error==None):
        error = 0.
    if(njRange==None and nbRange==None and kinRange==None and njSet==None and nbSet==None and kinSet==None):
        for i in range(1,hist.GetNbinsX()+1):
            h1.SetBinError(i, error)
    else:
        if(njRange==None):
            njRange = range(1,6)
        elif(type(njRange) is int):
            njRange = range(njRange,njRange+1)
        if(nbRange==None):
            nbRange = range(4)
        elif(type(nbRange) is int):
            nbRange = range(nbRange,nbRange+1)
        if(kinRange==None):
            kinRange = range(1,11)
        elif(type(kinRange) is int):
            kinRange = range(kinRange,kinRange+1)
        if(njSet==None):
            njSet=njRange
        elif(type(njSet) is int):
            njSet=range(njSet,njSet+1)
        if(nbSet==None):
            nbSet=nbRange
        elif(type(nbSet) is int):
            nbSet=range(nbSet,nbSet+1)
        if(kinSet==None):
            kinSet=kinRange
        elif(type(kinSet) is int):
            kinSet=range(kinSet,kinSet+1)

        for nj in njRange:
            for nb in nbRange:
                for kin in kinRange:
                    Bin = kin + nb*len(kinRange) + (nj-1)*len(nbRange)*len(kinRange)
                    if( (nj in njSet) and (nb in nbSet) and (kin in kinSet)):
                        h1.SetBinError(Bin,error)
    return h1

def addFractionalError(hist, njSplit=None, njRange=None, nbRange=None, kinRange=None, njError=None, nbError=None, kinError=None):
       
    if(njSplit==None):
        njSplit=False
    if(njSplit==True):
        if(njRange==None):
            njRange = range(1,8)
    if(njRange==None):
        njRange = range(1,5)
    elif(type(njRange) is int):
        njRange = range(njRange,njRange+1)
    if(nbRange==None):
        nbRange = range(4)
    elif(type(nbRange) is int):
        nbRange = range(nbRange,nbRange+1)
    if(kinRange==None):
        kinRange = range(1,11)
    elif(type(kinRange) is int):
        kinRange = range(kinRange,kinRange+1)
    if(njError==None):
        njError = numpy.zeros(len(njRange))
    if(nbError==None):
        nbError = numpy.zeros(len(nbRange))
    if(kinError==None):
        kinError = numpy.zeros(len(kinRange))

    h1 = hist.Clone()

    for nj in njRange:
        for nb in nbRange:
            for kin in kinRange:
                Bin = kin + nb*len(kinRange) + (nj-1)*len(nbRange)*len(kinRange)
                fracError = 0.
                if(h1.GetBinContent(Bin)>0.):
                    fracError = h1.GetBinError(Bin)/h1.GetBinContent(Bin)
                error = h1.GetBinContent(Bin)*(fracError**2+njError[nj-1]**2+nbError[nb]**2+kinError[kin-1]**2)**0.5
                h1.SetBinError(Bin,error)
  
    
    return h1

# from Jim Hirschauer
# modified by TM to return fractional error
def getClopperPearsonError(rNi, rNo, fracError=None):

    if(fracError==None):
        fracError=True

    # tail = (1-CL)/2.  If you want 95% CL, the tail = (1.-0.95)/2
    # 95%
    #tail = 0.025
    # 68%
    tail = 0.16
    rvrat = rNi/float(rNo)

    qHi = ROOT.Math.fdistribution_quantile_c(tail   ,2.*(rNi+1.),2.*(rNo)   )
    qLo = ROOT.Math.fdistribution_quantile_c(1.-tail,2.*(rNi)   ,2.*(rNo+1.))

    Rhi = 1./(rNo/(rNi+1.)/qHi)
    if rNi>0.:
        Rlo = 1./((rNo+1.)/rNi/qLo)
    else :
        Rlo = 0.

    #eGaus  = rvrat*math.sqrt(1./rNi + 1./rNo)

    eHi = Rhi-rvrat
    eLo = rvrat-Rlo
    if(fracError):
        eHi *= 1./rvrat
        eLo *= 1./rvrat

    error = [eLo, eHi]

    return error

def getRatioHist(numHist, denomHist, doFlip=None):
    if(doFlip==None):
        doFlip=False
        
    h = ROOT.TH1F()

    if(doFlip):
        h = denomHist.Clone()
        h.Divide(numHist)
    else:
        h = numHist.Clone()
        h.Divide(denomHist)

    return h

def getDiffHist(obsHist, expHist, doFlip=None, divide=None):
    if(doFlip==None):
        doFlip=False
        
    if(divide==None):
        divide=True
    h = ROOT.TH1F()

    if(doFlip):
        h = expHist.Clone()
        h.Add(obsHist,-1)
        if(divide):
            h.Divide(obsHist)
    else:
        h = obsHist.Clone()
        h.Add(expHist,-1)
        if(divide):
            h.Divide(expHist)
    return h

    #h_dr = getRatioHist(zll_r,pho_r)
    #g_dr = ROOT.TGraphAsymmErrors(h_dr)

def getPullHist(postHist, preHist, doFlip=None):
    if(doFlip==None):
        doFlip=False
        
    h = ROOT.TH1F()
    sig = ROOT.TH1F('sig','sig',preHist.GetNbinsX(),preHist.GetBinLowEdge(1),preHist.GetBinLowEdge(preHist.GetNbinsX()+1))
    for i in range(1,preHist.GetNbinsX()+1):
        sig.SetBinContent(i,preHist.GetBinError(i))

    if(doFlip):
        h = preHist.Clone()
        h.Add(postHist,-1)
        h.Divide(sig)
    else:
        h = postHist.Clone()
        h.Add(preHist,-1)
        h.Divide(sig)

    for i in range(1,preHist.GetNbinsX()+1):
        h.SetBinError(i,0.)

    return h

    #h_dr = getRatioHist(zll_r,pho_r)
    #g_dr = ROOT.TGraphAsymmErrors(h_dr)

def getRatioGraph(numHist,denomHist,ratioHist=None,doFlip=None,cg=None,doClopperPearsonError=None):

    if(doFlip==None):
        doFlip=False
    if(cg==None):
        cg = []
    if(doClopperPearsonError==None):
        doClopperPearsonError=True

    h_r = getRatioHist(numHist,denomHist)

    if(ratioHist is not None):
        h_r = ratioHist

    g_r = ROOT.TGraphAsymmErrors(h_r)


    if(cg==[]):
        # center of gravity not provided
        # using bin centers
        for i in range(1,numHist.GetNbinsX()+1):
            cg.append(numHist.GetBinCenter(i))

    for i in range(1,numHist.GetNbinsX()+1):

        if(doClopperPearsonError):
            cpError = getClopperPearsonError(numHist.GetBinContent(i),denomHist.GetBinContent(i))
            eLo = float()
            eHi = float()

            eLo = h_r.GetBinContent(i)*ROOT.TMath.Sqrt(cpError[0]**2+(h_r.GetBinError(i)/h_r.GetBinContent(i))**2)
            eHi = h_r.GetBinContent(i)*ROOT.TMath.Sqrt(cpError[1]**2+(h_r.GetBinError(i)/h_r.GetBinContent(i))**2)
            # eLo = h_r.GetBinContent(i)*cpError[0]**2
            # eHi = h_r.GetBinContent(i)*cpError[1]**2
        else:
            eLo = h_r.GetBinError(i)/2.
            eHi = eLo

        cv = h_r.GetBinContent(i)

        if(doFlip):
            eLo = (eLo/h_r.GetBinContent(i)**2)
            eHi = (eHi/h_r.GetBinContent(i)**2)
            cv = 1/cv

        g_r.SetPointError(i-1,0,0,eLo,eHi)
        g_r.SetPoint(i-1,cg[i-1],cv)

    return g_r

def getPlotAndRatio(numHists, denomHists, bottomPlots=None, doStack=None, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, ratioTitle=None, ratioMin=None, ratioMax=None, doLogy=None, doFlip=None, doDiff=None, doPull=None, makeLeg=None, legList=None, legCoords=None, textCoords=None, canvasSize=None, canvasName=None, numColors=None, denomColor=None, numMarkers=None, denomMarker=None, markerSize=None, lineWidth=None, numDrawStyles=None, denomDrawStyle=None, drawErrorBand=None, stackColors=None, axisTitleSize=None, drawVerticalLines=None, drawHorizontalLine=None, statBox=None, drawText=None, text=None, setMax=None, setMin=None, doClosureStyle=None,errorBandColor=None,errorBandFillStyle=None,legHeader=None,nDivRatio=None,doNumFill=None, hLineVal=None, hLineColors=None,nDivX=None,ratioGridx=None,ratioGridy=None,topGridx=None,topGridy=None):

    ROOT.gROOT.Reset()

    setTdrStyle()

    denomHist = ROOT.TH1F

    #####################################################################
    # set defaults and initialize 
    #####################################################################
    if(doClosureStyle==None):
        doClosureStyle=False
    if type(numHists) is ROOT.TH1F:
        numHists = [numHists]
    if type(denomHists) is ROOT.TH1F:
        denomHist = denomHists
        doStack=False
    else:
        denomHist = denomHists[0].Clone()
        for i in range(1,len(denomHists)):
            denomHist.Add(denomHists[i])
        if(doStack==None):
            doStack=True
    if(Title==None):
        Title=""
    if(xTitle==None):
        xTitle = "Bin"
    if(yTitle==None):
        yTitle = "Events"
    if(doCMSlumi==None):
        doCMSlumi=True
    if(iPos==None):
        iPos=0
    if(doLogy==None):
        doLogy=True
    if(doDiff==None):
        doDiff = False
    if(doDiff==True and ratioTitle==None):
        ratioTitle = "#frac{obs-exp}{exp}"
    if(ratioTitle==None):
        ratioTitle = "Ratio"
    if(doDiff==True and ratioMin==None):
        ratioMin=-1.01
    if(doDiff==True and ratioMax==None):
        ratioMax=1.01
    if(ratioMin==None):
        ratioMin=0.01
    if(ratioMax==None):
        ratioMax=1.9
    if(doFlip==None):
        doFlip = False
    if(makeLeg==None):
        makeLeg=True
    if(makeLeg==False and (legCoords!=None or legList!=None)):
        print "****WARNING****"
        print "****Legend coordinates and/or legend list were given but makeLeg was set to false"
        print "****Legend will not be drawn*************"
    if(legList==None):
        legList=[]
    if(legCoords==None):
        legCoords = [0.6,0.6,0.89,0.89]
    if(textCoords==None):
        textCoords = [0.20,0.70,.40,.83]
    if(canvasSize==None):
        canvasSize = [900,600]
    if(canvasName==None):
        canvasName = "canvas"
    if(numColors==None):
        numColors = [1,4,2,8,12,7,5,13,46,41,6]
    if(denomColor==None):
        denomColor=4
    if(numMarkers==None):
        numMarkers = [20,20,20,20,20,20,20,20,20]
    if(denomMarker==None):
        denomMarker=20
    if(markerSize==None):
        if(denomHist.GetNbinsX()>=30):
            markerSize=1.1
        else:
            markerSize=1.2
    if(lineWidth==None):
        lineWidth=2
    if(denomDrawStyle==None):
        denomDrawStyle="HIST"
    if(denomDrawStyle=="HIST"):
        denomMarker=1
    if(numDrawStyles==None):
        numDrawStyles=["e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0"]
    if(drawErrorBand==None):
        drawErrorBand=True
    if(stackColors==None):
        stackColors = [3, 46, 12, 7, 4, 5, 40, 5, 42, 8, 38]
    if(axisTitleSize==None):
        axisTitleSize=0.1
    if(drawVerticalLines==None):
        drawVerticalLines=True
    if(statBox==None):
        statBox = 0
    if(drawText==None):
        if(doStack==True and len(numHists)==1):
            drawText=True
        else:
            drawText=False
    if(text==None):
        if(doStack==True and len(numHists)==1):
            dataOverMC = numHists[0].Integral()/denomHist.Integral()
            text = "Data/MC = "+str(round(dataOverMC,2))
        else:
            text = "arXiv:1602.06581"
    if(extraText==None):
        extraText = "Preliminary"
    if(drawHorizontalLine==None):
        if(doStack==True and len(numHists)==1):
            drawHorizontalLine=True
        else:
            drawHorizontalLine=False
    if(hLineVal==None):
        if(doStack==True and len(numHists)==1):
            dataOverMC = numHists[0].Integral()/denomHist.Integral()
            hLineVal=[1.,dataOverMC]
        else:
            hLineVal=[1.]
        if(doDiff==True):
            if(doStack==True and len(numHists)==1):
                dataOverMC = (numHists[0].Integral()-denomHist.Integral())/denomHist.Integral()
                hLineVal=[0.,dataOverMC]
            else:
                hLineVal=[0.]
    if(type(hLineVal) is not list):
        hLineVal=[hLineVal]
    if(doPull==None):
        doPull=False
    if(errorBandColor==None):
        errorBandColor=denomColor
    if(errorBandFillStyle==None):
        errorBandFillStyle=3002
    if(nDivRatio==None):
        nDivRatio=5
    if(nDivX==None):
        nDivX=510
    if(ratioGridx==None):
        ratioGridx=False
    if(ratioGridy==None):
        ratioGridy=False
    if(topGridx==None):
        topGridx=False
    if(topGridy==None):
        topGridy=False
    if(doClosureStyle==True):
        errorBandFillStyle=3144
        errorBandColor = ROOT.kRed-10
        ratioTitle = "#frac{Direct}{Prediction} "
        extraText = "Simulation"
        xTitle = "Search region bin number"
        #legCoords = [.65, .95, .54, .79]
        if(legHeader==None):
            legHeader = "Z #rightarrow #nu#bar{#nu} background"
        ratioMax=2.15
        ratioMin=0.001
        nDivRatio=505
        markerSize=1.1
        legList = ['Direct from simulation','Treat simulation like data']
    if(legHeader==None):
        legHeader=""
    if(doNumFill==None):
        doNumFill=False
    if(hLineColors==None):
        hLineColors = [4,2,1,1]
    if(type(hLineColors) is not list):
        hLineColors=[hLineColors]
    #####################################################################
    # end set defaults
    #####################################################################

    #####################################################################
    # get ratio/diff plots
    #####################################################################

    numHistsRatio = []

    if(bottomPlots==None):
        if(doClosureStyle==True):
            denomHistNoError = denomHist.Clone()
            for i in range(1,denomHistNoError.GetNbinsX()+1):
                denomHistNoError.SetBinError(i,0.)
            for i in range(len(numHists)):
                numHistsRatio.append(getRatioHist(numHists[i],denomHistNoError,doFlip))
            numHistsRatio.append(getRatioHist(denomHist,denomHistNoError,doFlip))
        else:
            for i in range(len(numHists)):
                if(doDiff):
                    numHistsRatio.append(getDiffHist(numHists[i],denomHist,doFlip))
                elif(doPull):
                    numHistsRatio.append(getPullHist(numHists[i],denomHist,doFlip))
                else: #default
                    numHistsRatio.append(getRatioHist(numHists[i],denomHist,doFlip))
    else:
        if(type(bottomPlots) is not list):
            bottomPlots = [bottomPlots]
        for p in bottomPlots:
            numHistsRatio.append(p)
    
    #####################################################################
    # end get ratio/diff plots
    #####################################################################
    
    #####################################################################
    # make unique root canvas object 
    #####################################################################

    cIter = 1
    while(type(ROOT.gROOT.FindObject(canvasName+str(cIter)))==ROOT.TCanvas):
        cIter+=1
    canvasName+=str(cIter)
  
    canvas = ROOT.TCanvas(canvasName, canvasName, 0,0, canvasSize[0],canvasSize[1])
    canvas.Range(0,0,1,1)
    canvas.cd()
 
    #####################################################################
    # end make unique root canvas object 
    #####################################################################

    n=0
    seg=0

    if(denomHist.GetSize()==74):
        seg = 6
        n = 12
    elif(denomHist.GetSize()==20):
        seg = 6
        n = 3
    elif(denomHist.GetSize()==30):
        seg = 4
        n = 7
    elif(denomHist.GetSize()==18):
        seg = 4
        n = 4
    elif(denomHist.GetSize()==162):
        seg = 10
        n = 16
        
    #####################################################################
    # bottom plot
    #####################################################################

    pad_1 = ROOT.TPad("pad_1","pad_1",0.00,0.01,0.99,0.32)
    pad_1.Draw()
    pad_1.cd()
    pad_1.SetTopMargin(0.01)
    pad_1.SetBottomMargin(0.4)
    pad_1.SetRightMargin(0.1)
    pad_1.SetGridx(ratioGridx)
    pad_1.SetGridy(ratioGridy)

    for i in range(len(numHistsRatio)):
        numHistsRatio[i].SetMinimum(ratioMin)
        numHistsRatio[i].SetMaximum(ratioMax)
        numHistsRatio[i].SetTitle(";"+xTitle+";"+ratioTitle)
        numHistsRatio[i].GetXaxis().SetTitleSize((2*axisTitleSize))
        numHistsRatio[i].GetXaxis().SetTitleOffset(0.043/(0.5*axisTitleSize))
        numHistsRatio[i].GetXaxis().SetLabelSize(1.5*axisTitleSize)
        numHistsRatio[i].GetYaxis().SetLabelSize(1.5*axisTitleSize)
        numHistsRatio[i].GetYaxis().SetTitleSize(2*axisTitleSize)
        numHistsRatio[i].GetYaxis().SetTitleOffset(0.033/(1.0*axisTitleSize))
        numHistsRatio[i].GetYaxis().SetNdivisions(nDivRatio)
        numHistsRatio[i].GetXaxis().SetNdivisions(nDivX)
        numHistsRatio[i].SetLineWidth(lineWidth)
        numHistsRatio[i].SetLineColor(numColors[i])
        numHistsRatio[i].SetMarkerStyle(numMarkers[i])
        numHistsRatio[i].SetMarkerSize(markerSize)
        numHistsRatio[i].SetMarkerColor(numColors[i])
        if(doNumFill):
            if('p' not in numDrawStyles[i]):
               numHistsRatio[i].SetFillColor(numColors[i])
               numHistsRatio[i].SetFillStyle(3004)
               numHistsRatio[i].SetMarkerSize(0.001)
               numHistsRatio[i].SetMarkerStyle(1)
               numHistsRatio[i].SetLineWidth(0)
        if(doClosureStyle==True):
            lineAtOne = numHistsRatio[len(numHistsRatio)-1].Clone()
            for j in range(1,lineAtOne.GetNbinsX()+1):
                lineAtOne.SetBinContent(j,1)
                lineAtOne.SetBinError(j,0.)
            numHistsRatio[len(numHistsRatio)-1].SetMarkerSize(0.001)
            numHistsRatio[len(numHistsRatio)-1].SetMarkerStyle(1)
            numHistsRatio[len(numHistsRatio)-1].SetFillColor(errorBandColor)
            numHistsRatio[len(numHistsRatio)-1].SetFillStyle(errorBandFillStyle)
            numHistsRatio[len(numHistsRatio)-1].SetLineColor(errorBandColor)
            numHistsRatio[len(numHistsRatio)-1].SetMarkerColor(errorBandColor)
            if(i==0):
                numHistsRatio[len(numHistsRatio)-1].Draw("e2")
                lineAtOne.Draw("hist same")
            ROOT.SetOwnership(lineAtOne, 0)
            if(i is not (len(numHistsRatio)-1)):
               numHistsRatio[i].SetStats(statBox)
               numHistsRatio[i].Draw("same"+numDrawStyles[i])
        elif(type(numHistsRatio[i]) is not ROOT.TGraphAsymmErrors):
            numHistsRatio[i].SetStats(statBox)
            #numHistsRatio[i].Draw("same"+numDrawStyles[i])
            numHistsRatio[i].Draw(numDrawStyles[i]+'same')
        else:
            numHistsRatio[i].GetXaxis().SetLimits(denomHist.GetBinLowEdge(1),denomHist.GetBinLowEdge(denomHist.GetNbinsX()+1))
            numHistsRatio[i].Draw('ap same')
            

    # draw lines 

    horizontalLine = []
    hLineWidths = [1,2,1,1]
    hLineStyles = [1,2,1,1]
    for i in range(len(hLineVal)):
        horizontalLine.append(ROOT.TLine(denomHist.GetBinLowEdge(1),hLineVal[i],denomHist.GetBinLowEdge(denomHist.GetNbinsX()+1),hLineVal[i]))
        horizontalLine[i].SetLineWidth(hLineWidths[i])
        horizontalLine[i].SetLineStyle(hLineStyles[i])
        horizontalLine[i].SetLineColor(hLineColors[i])
        if(drawHorizontalLine==True):
            horizontalLine[i].Draw()
        ROOT.SetOwnership(horizontalLine[i],0)
        
    verticalLineB = []
    if(drawVerticalLines==True):
        for i in range(1,n):
            verticalLineB.append(ROOT.TLine(i*seg+0.5,ratioMin,i*seg+0.5,ratioMax))
            verticalLineB[i-1].SetLineWidth(3)
            verticalLineB[i-1].SetLineStyle(2)
            verticalLineB[i-1].Draw()
    
    pad_1.Draw()
    pad_1.Update()
    canvas.Update()
    #####################################################################
    # end bottom plot
    #####################################################################

    #####################################################################
    # top plot
    #####################################################################
    canvas.cd()
    pad_2 = ROOT.TPad("pad_2","pad_2",0.00,0.33,0.99,0.99)
    pad_2.Draw()
    pad_2.cd()
    pad_2.SetTopMargin(0.1)
    pad_2.SetBottomMargin(0.03)
    pad_2.SetRightMargin(0.1)
    pad_2.SetGridx(topGridx)
    pad_2.SetGridy(topGridy)

    denomHist.SetStats(statBox)
    denomHist.SetTitleSize(axisTitleSize)
    denomHist.GetXaxis().SetTitleSize(0.)
    denomHist.GetYaxis().SetTitleSize(axisTitleSize)
    denomHist.GetYaxis().SetTitleOffset(0.033/(0.5*axisTitleSize))
    denomHist.SetTitle(Title+";;"+yTitle)
    denomHist.SetLineWidth(lineWidth)
    denomHist.SetLineColor(denomColor)
    denomHist.SetMarkerStyle(denomMarker)
    denomHist.SetMarkerSize(markerSize)
    denomHist.SetMarkerColor(denomColor)
    denomHist.GetYaxis().SetLabelSize(0.8*axisTitleSize)
    denomHist.GetXaxis().SetLabelSize(0.0)
    
    if(not(doClosureStyle)):
        denomHist.Draw(denomDrawStyle)
        
    hs = ROOT.THStack("hs","hs")
    if(doStack):
        for i in range(len(denomHists)):
            denomHists[i].SetLineColor(1)
            denomHists[i].SetFillColor(stackColors[len(stackColors)-len(denomHists)+i])
            hs.Add(denomHists[i])
        hs.Draw("HIST")
        hs.GetHistogram().SetTitleSize(axisTitleSize)
        hs.GetHistogram().GetXaxis().SetTitleSize(0.)
        hs.GetHistogram().GetXaxis().SetLabelSize(0.)
        hs.GetHistogram().GetYaxis().SetTitleSize(axisTitleSize)
        hs.GetHistogram().GetYaxis().SetTitleOffset(0.033/(0.5*axisTitleSize))
        hs.GetHistogram().SetTitle(Title+";;"+yTitle)
        hs.GetHistogram().GetYaxis().SetLabelSize(0.8*axisTitleSize)

    errorBandHist = denomHist.Clone()
    errorBandHist.SetLineColor(denomColor)
    errorBandHist.SetFillColor(errorBandColor)
    #errorBandHist.SetFillStyle(3018)
    errorBandHist.SetFillStyle(errorBandFillStyle)

    if(drawErrorBand==True):
        errorBandHist.Draw("e2 same")

    if(doClosureStyle):
        denomHist.Draw(denomDrawStyle+"same")

    Min = denomHist.GetBinContent(denomHist.FindLastBinAbove(0))
    Max = denomHist.GetBinContent(denomHist.FindFirstBinAbove(0))
    for i in range(len(numHists)):
        numHists[i].SetStats(statBox)
        numHists[i].SetLineWidth(lineWidth)
        numHists[i].SetLineColor(numColors[i])
        numHists[i].SetMarkerStyle(numMarkers[i])
        numHists[i].SetMarkerSize(markerSize)
        if(doNumFill):
            if('p' not in numDrawStyles[i]):
               numHists[i].SetFillColor(numColors[i])
               numHists[i].SetFillStyle(3004)
               numHists[i].SetMarkerSize(0.001)
               numHists[i].SetMarkerStyle(1)
               numHists[i].SetLineWidth(0)
        numHists[i].SetMarkerColor(numColors[i])
        numHists[i].Draw("same"+numDrawStyles[i])
        if(Min>numHists[i].GetBinContent(numHists[i].FindLastBinAbove(0))):
            Min=numHists[i].GetBinContent(numHists[i].FindLastBinAbove(0))
        if(Max<numHists[i].GetMaximum()):
            Max=numHists[i].GetMaximum()

    if(doLogy):
        Min*=0.5
        Max*=5
        if(doStack==True): # goofy ROOT THStack max min issue
            hs.SetMinimum(Min*(1+0.5*ROOT.TMath.Log10(Max/Min)))
            hs.SetMaximum(Max/(1+0.2*ROOT.TMath.Log10(Max/Min)))
    else:
        Min=0
        Max*=1.2
        if(doStack==True): # goofy ROOT THStack max min issue
            hs.SetMinimum(Min)
            hs.SetMaximum(Max/(1+ROOT.gStyle.GetHistTopMargin()))

    
    denomHist.SetMinimum(Min)
    denomHist.SetMaximum(Max)
    errorBandHist.SetMinimum(Min)
    errorBandHist.SetMaximum(Max)

    if(setMax is not None):
        Max = setMax
        denomHist.SetMaximum(setMax)
        #errorBandHist.SetMaximum(Max)
    if(setMin is not None):
        Min = setMin
        denomHist.SetMinimum(setMin)
        #errorBandHist.SetMinimum(setMin)
    #####################################################################
    # build legend
    #####################################################################

    leg = pad_2.BuildLegend(legCoords[0],legCoords[1],legCoords[2],legCoords[3], "")  

    if(len(legList)>0):
        leg.Delete()
        leg = ROOT.TLegend(legCoords[0],legCoords[1],legCoords[2],legCoords[3], "","brNDC")
        for i in range(len(legList)):
            if(i>=len(numHists)):
                if(doStack==True):
                    leg.AddEntry(denomHists[i-len(denomHists)], legList[i],'f')
                elif(drawErrorBand==True):
                    leg.AddEntry(errorBandHist, legList[i], 'fl')
                elif(denomDrawStyle=='HIST'):
                    leg.AddEntry(denomHist, legList[i], 'l')
                else:
                    leg.AddEntry(denomHist, legList[i], 'ep')
            else:
                if(numDrawStyles[i]=='HIST' or numDrawStyles[i]=='eHIST'):
                    leg.AddEntry(numHists[i], legList[i],'l')
                else:
                    leg.AddEntry(numHists[i], legList[i],numDrawStyles[i])
    leg.SetBorderSize(1)
    leg.SetFillStyle(1001) 
    leg.SetFillColor(0) 
    leg.SetLineWidth(0)
    leg.SetHeader(legHeader)

    pad_2.Draw()
    pad_2.SetLogy(doLogy)

    #####################################################################
    # end build legend
    #####################################################################

    # draw lines
    verticalLineT = []
    for i in range(1,n):
        verticalLineT.append(ROOT.TLine(i*seg+0.5,Min,i*seg+0.5,Max))
        verticalLineT[i-1].SetLineWidth(3)
        verticalLineT[i-1].SetLineStyle(2)
        if(drawVerticalLines==True):
            verticalLineT[i-1].Draw()


    #####################################################################
    # end top plot
    #####################################################################


    pt = ROOT.TPaveText(textCoords[0],textCoords[1],textCoords[2],textCoords[3],"brNDC")
    #pt = ROOT.TPaveText(0.65,0.80,.85,.85,"brNDC")
    if(type(text) is not list):
        text = [text]
    for txt in text:
        pt.AddText(txt)

    pt.SetFillColor(0)
    if(drawText==True):
        pt.Draw("NB")

    if(makeLeg):
        leg.Draw("NB")

    # deal with pyROOT scope issues
    ROOT.SetOwnership(errorBandHist, 0)
    ROOT.SetOwnership(hs, 0)
    ROOT.SetOwnership(leg, 0)
    ROOT.SetOwnership(pt, 0)
    ROOT.SetOwnership(pad_1,0)
    ROOT.SetOwnership(pad_2,0)

    for line in verticalLineB:
        ROOT.SetOwnership(line, 0)
    for line in verticalLineT:
        ROOT.SetOwnership(line, 0)
    for hist in numHistsRatio:
        ROOT.SetOwnership(hist, 0)

    # apply cms lumi info to pad

    if(doCMSlumi):
        cmsLumi(pad_2, iPeriod, iPos, extraText)

    return canvas

def getDevSyst(hists, binSplit=None, doFlip=None):

    if(type(hists) is not list):
        hists = [hists]
    if(binSplit==None):
        binSplit=[range(1,hists[0].GetNbinsX()+1)]
    if(type(binSplit) is not list):
        binSplit = [[binSplit]]
    if(type(binSplit[0]) is not list):
        binSplit = [binSplit]
    if(doFlip==None):
        doFlip=False

    if(len(hists)>2):
        print "**** cannot accept more than 2 hists **** using the first 2 hists"
    if(len(hists)>1):
        h1 = getRatioHist(hists[0],hists[1],doFlip)
    else:
        h1 = hists[0].Clone()

    systList = []
    canvasList = []
    fitList = []
    hList = []
    eHistList = []
    ptList = []

    canvasName = 'canvas_'
    histName = 'h_syst_'

    for i in range(len(binSplit)):

        cIter = 1
        while(type(ROOT.gROOT.FindObject(canvasName+str(cIter)))==ROOT.TCanvas):
            cIter+=1
            
        canvasName+=str(cIter)
        histName+=str(cIter)
        canvasList.append(ROOT.TCanvas(canvasName, canvasName, 0, 0, 600, 600))

        hList.append(ROOT.TH1F(histName+str(i),histName+str(i),len(binSplit[i]),binSplit[i][0]-0.5,binSplit[i][len(binSplit[i])-1]+0.5))

        for j in range(1,len(binSplit[i])+1):
            hList[i].SetBinContent(j,h1.GetBinContent(binSplit[i][j-1]))
            hList[i].SetBinError(j,h1.GetBinError(binSplit[i][j-1]))

        fitList.append(ROOT.TF1("fit_"+str(cIter), 'pol0(0)', hList[i].GetBinLowEdge(1), hList[i].GetBinLowEdge(hList[i].GetNbinsX()+1)+0.5))

        hList[i].Fit(fitList[i])
                    
        rms = 0.
        for j in range(1,len(binSplit[i])+1):
            rms += (hList[i].GetBinContent(j)-fitList[i].GetParameter(0))**2
        
        rms = (rms/fitList[i].GetNDF())**0.5

        if((fitList[i].GetChisquare()/fitList[i].GetNDF())>1.):
            systList.append(rms*(1.-1./(fitList[i].GetChisquare()/fitList[i].GetNDF())**2)**0.5)
        else:
            systList.append(rms)

        eHistList.append(hList[i].Clone())
        for j in range(1,len(binSplit[i])+1):
            eHistList[i].SetBinContent(j,fitList[i].GetParameter(0))
            eHistList[i].SetBinError(j,systList[i])
            
        eHistList[i].SetFillStyle(3004)
        eHistList[i].SetFillColor(ROOT.kRed)

        eHistList[i].SetMarkerSize(0.0001)

        eHistList[i].SetMinimum(0.0)
        eHistList[i].SetMaximum(1.5)
        
        eHistList[i].Draw('e2')
        hList[i].Draw('ep same')

        eHistList[i].GetXaxis().SetTitle('RA2b Bin')
        eHistList[i].GetYaxis().SetTitle('Sim/Pred')

        ptList.append(ROOT.TPaveText(0.2,0.2,0.5,0.5,"brNDC"))
        if((fitList[i].GetChisquare()/fitList[i].GetNDF())>1.):
            ptList[i].AddText('#sqrt{1-#frac{1}{(#chi^{2}/ndf)^{2}}}*(x_{rms}) = '+str(round(systList[i],3)))
        else:
            ptList[i].AddText('x_{rms} = '+str(round(systList[i],3)))
        ptList[i].SetFillColor(0)
        ptList[i].Draw("NB")

        ROOT.SetOwnership(hList[i], 0)
        ROOT.SetOwnership(canvasList[i], 0)
        ROOT.SetOwnership(fitList[i], 0)
        ROOT.SetOwnership(eHistList[i], 0)
        ROOT.SetOwnership(ptList[i], 0)


    return systList
        


def getDoubleRatioGraph(dist, binning=None, extraCuts=None, nJetBin=None, bJetBin=None, kinBin=None, dphiCut=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, doLumi=None, removeZkfactor=None, doNoisy=None, doProof=None, treeLoc=None, applyEffs=None, treeName=None, applySF=None, doCG=None, do20=None, doClopperPearsonError=None, applyPuWeight=None, doDiMu=None, doDiEl=None):

    #####################################################################
    # set defaults and initialize 
    #####################################################################
    if(doClopperPearsonError==None):
        doClopperPearsonError=True
    if(do20==None):
        do20=False
    if(do20):
        if(doLumi==None):
            doLumi=24.50
    if(dist=='HT'):
        if(binning==None):
            #binning=[300,400,500.,650.,800.,1000.,1200.,1600.]
            #binning=[300,350,425,550,800.,1600.]

            #binning=[300.,400.,500,750.,1000.,1600.]
            binning=[300,400,450,550,650,750,1000,1200,1600]
    if(dist=='MHT'):
        if(binning==None):
            #binning=[300.,400.,500.,600.,750.,900.]
            binning=[300.,350.,400.,450.,600.,750,900.]
    if(dist=='NJets'):
        if(applyNJetsCut==None):
            applyNJetsCut=False
        if(binning==None):
            binning=[1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]
            #binning=[1.5,2.5,4.5,6.5,9.5]
            #binning=[2.5,3.5,4.5,5.5,6.5,7.5]
    if(dist=='BTags'):
        bJetBin=-1
        if(binning==None):
            binning=[-0.5,0.5,1.5,2.5,3.5]
    if(doLumi==None):
        doLumi = 24.5
    if(removeZkfactor==None):
       removeZkfactor=True
    if(doNoisy==None):
        doNoisy = False
    if(doProof==None):
        doProof = False
    if(extraCuts==None):
        extraCuts=""
    if(nJetBin==None):
        nJetBin=-1
    if(bJetBin==None):
        bJetBin=0
    if(kinBin==None):
        kinBin=-1
    if(applyMassCut==None):
        applyMassCut=True
    if(applyPtCut==None):
        applyPtCut=True
    if(applyHTCut==None):
        applyHTCut=True
    if(applyMHTCut==None):
        applyMHTCut=True
    if(applyNJetsCut==None):
        applyNJetsCut=False
    if(dphiCut==None):
        dphiCut = 'sig'
    if(treeLoc==None):
        treeLoc = "/home/ww/work/data/lpcTrees/Skims/Run2ProductionV10"
    if(applyEffs==None):
        applyEffs=True
    if(treeName==None):
        treeName = "tree"
    if(applySF==None):
        applySF=False
    if(doCG==None):
        doCG=True
    if(applyPuWeight==None):
        applyPuWeight=True
    sampleSuffix = ''
    if(dphiCut=='ldp'):
        sampleSuffix='LDP'
    if(doDiMu==None):
        doDiMu=True
    if(doDiEl==None):
        doDiEl=True

    ########################################################
    # get bin center of gravity from photon distribution
    ########################################################
    cg = []
    if(doCG==True):
        for i in range(len(binning)-1):
            h_cg = getDist('photon'+sampleSuffix, dist, distRange=[binning[i],binning[i+1]], nBins=100, doVarBinning=False,
                           extraCuts=extraCuts, nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                           applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                           doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName)
            cg.append(h_cg.GetMean())

    ########################################################
    # end get bin center of gravity from photon distribution
    ########################################################


    ########################################################
    # get zll,dyll,photon, and gjets histograms
    ########################################################
    
    overflowmax = (binning[len(binning)-1]+binning[len(binning)-2])/2.
    dist = "min("+dist+","+str(overflowmax)+")"

    pho_da = getDist('photon'+sampleSuffix, dist, doVarBinning=True, binning=binning, extraCuts=extraCuts, 
                     nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                     applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                     doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName)

    doCombo=False
    if(doCombo==True):
        extraCuts = '@GenJets.size()>=5'
    pho_mc = getDist('gjets'+sampleSuffix, dist, doVarBinning=True, binning=binning, extraCuts=extraCuts, 
                     nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                     applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                     doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName, applySF=applySF,
                     applyPuWeight=applyPuWeight)
    
    if(doCombo==True):
        extraCuts = '@GenJets.size()<5'
        pho_mc_old = getDist('gjetsold'+sampleSuffix, dist, doVarBinning=True, binning=binning, extraCuts=extraCuts, 
                         nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                         applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                         doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName, applySF=applySF,
                         applyPuWeight=applyPuWeight)
        pho_mc.Add(pho_mc_old)
        extraCuts = None

    zmm_da = getDist('zmm'+sampleSuffix, dist, doVarBinning=True, binning=binning, extraCuts=extraCuts, 
                     nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                     applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                     doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName)

    zmm_mc = getDist('dymm'+sampleSuffix, dist, doVarBinning=True, binning=binning, extraCuts=extraCuts, 
                     nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                     applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                     doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName, applySF=applySF,
                     applyPuWeight=applyPuWeight)

    zee_da = getDist('zee'+sampleSuffix, dist, doVarBinning=True, binning=binning, extraCuts=extraCuts, 
                     nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                     applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                     doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName)

    zee_mc = getDist('dyee'+sampleSuffix, dist, doVarBinning=True, binning=binning, extraCuts=extraCuts, 
                     nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                     applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                     doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName, applySF=applySF,
                     applyPuWeight=applyPuWeight)

    # if(do20):
    #     pho_da.Sumw2(False)
    #     zmm_da.Sumw2(False)
    #     zee_da.Sumw2(False)

    if(doDiMu==True):
        zll_da = zmm_da.Clone()
        zll_mc = zmm_mc.Clone()
        if(doDiEl==True):
            zll_da.Add(zee_da)
            zll_mc.Add(zee_mc)
    else:
        zll_da = zee_da.Clone()
        zll_mc = zee_mc.Clone()
        

    zll_r = getRatioHist(zll_da,zll_mc)
    pho_r = getRatioHist(pho_da,pho_mc)
    h_dr = getRatioHist(zll_r,pho_r)

    g_dr = getRatioGraph(zll_da,pho_da,ratioHist=h_dr,cg=cg,doClopperPearsonError=doClopperPearsonError)

    g_dr.GetXaxis().SetLimits(binning[0],binning[len(binning)-1])
    
    return g_dr

# from Jim Hirschauer
def getNewPars(parsList, covMatSize, sqrtCovMat):
    # Get new 2nd order poly
    # parsList = [0.4015, 0.0001412, -0.00000004877]
    size = len(parsList)
    parsArray = array('d',parsList)
    parsVec = ROOT.TVectorD(covMatSize, parsArray)
    unitGaussRand = ROOT.TVectorD(covMatSize)
    
    for loop in range(covMatSize):
        unitGaussRand[loop] = random.gauss(0., 1.)
        
    for outer in range(covMatSize):
        for inner in range(covMatSize):
            parsVec[outer] += sqrtCovMat(outer, inner) * unitGaussRand(inner)

    return parsVec

# from Jim Hirschauer
# modified by TM to accept tgraph instead of hist
def getDoubleRatioFit(hist, func, returnDiff=False):
    hist.Fit(func,"RNQ")

    loLim = ROOT.Double(-0.5)
    hiLim = ROOT.Double(0.)
    func.GetRange(loLim, hiLim)
    curveDef = ROOT.TString()
    curveDef = func.GetExpFormula()
    # Put fit results in parsList
    fixParsList = []
    freeParsList = []
    fixParsPairList = []
    freeParsPairList = []
    size = func.GetNpar()
    covMatSize = 0
    for ipar in range(size):
        # Free params will have finite error
        if func.GetParError(ipar) > 0:
            covMatSize+=1
            freeParsList.append(func.GetParameter(ipar))
            freeParsPairList.append((ipar,func.GetParameter(ipar)))
        else:
            fixParsList.append(func.GetParameter(ipar))
            fixParsPairList.append((ipar,func.GetParameter(ipar)))
            
    # Get fit covariance matrix
    fitter = ROOT.TVirtualFitter.GetFitter()
    covMat = ROOT.TMatrixD(covMatSize,covMatSize,fitter.GetCovarianceMatrix())
    #covMat.Print()

    # Get square root of cov matrix
    valMat = ROOT.TMatrixD(covMatSize,covMatSize)
    sqrtCovMat = ROOT.TMatrixD(covMatSize,covMatSize)
    tmpMat = ROOT.TMatrixD(covMatSize,covMatSize)
    eigenvalues = ROOT.TVectorD(covMatSize)
    eigenvectors = covMat.EigenVectors(eigenvalues)
    ev_invert = eigenvectors.Clone()
    ev_invert.Invert()
    for loop in range(covMatSize):
        valMat[loop][loop] = math.sqrt(eigenvalues(loop))
    tmpMat.Mult(valMat,ev_invert)
    sqrtCovMat.Mult(eigenvectors,tmpMat)

    # Now generate ncurves random curves by throwing params based on the
    # fit cov matrix

    #    print "Starting loop over curves"
    curves = {}
    hists = {}

    # First get binning of hist that I originally fit, dijet ratio uses
    # variable binning
    # nbins = hist.GetNbinsX()
    nbins = hist.GetN()
    binCenters = []
    binLowEdges = []
    ptsForGrph = []

    #########################################
    # Edit mbins for granularity of output
    #########################################
    # JFH 09-Oct-15
    # loEnd = hist.GetBinLowEdge(1)
    #loEnd = hist.GetX().__getitem__(0)*0.001
    loEnd = -0.5
    # hiEnd = hist.GetBinLowEdge(nbins)+hist.GetBinWidth(nbins)
    hiEnd = hist.GetX().__getitem__(nbins-1)*2
    mbins = 100
    bwidth = (hiEnd-loEnd)/float(mbins)
    for ibin in range(mbins+1):
        ptsForGrph.append(loEnd + ibin*bwidth)
    #############################

    #    for bin in range(nbins+2): # JFH 09-Oct-15    
    for bin in range(mbins+1): # JFH 09-Oct-15    
        hists[bin] = ROOT.TH1F("bin"+str(bin), "bin"+str(bin), 1000,-5.0,5.0)
        if bin == 0:
            # binLowEdge = hist.GetBinLowEdge(bin+1)
            binLowEdge = hist.GetX().__getitem__(bin)
            #ptsForGrph.append(binLowEdge)  # JFH 09-Oct-15
        else:
            # binCenter = hist.GetBinCenter(bin)
            binCenter = hist.GetX().__getitem__(bin)
            if binCenter > loLim  and binCenter < hiLim:
                binCenters.append(binCenter)
                #ptsForGrph.append(binCenter)  # JFH 09-Oct-15
                binLowEdges.append(binLowEdge)
    
    #    nbins = len(binCenters)
    nbins = len(ptsForGrph)

    # Do actual throwing and fill curves dictionary
    ncurves = 10000
    for curve in range(ncurves):
        # Get new params using getNewPars()
        parsVec = getNewPars(freeParsList, covMatSize, sqrtCovMat)
        #        curves[curve] = ROOT.TF1("curve"+str(curve),curveDef,loLim, hiLim)
        curves[curve] = func.Clone()
        for el in fixParsPairList:
            curves[curve].SetParameter(el[0],el[1])
        npar = 0
        for el in freeParsPairList:
            curves[curve].SetParameter(el[0],parsVec[npar])
            npar += 1

    # To get high and low systematic error curves, make histograms of
    # curve values in each bin.  Fill histograms, one for each bin, of
    # curve values:

    curve_vals = {}
    for bin in range(mbins+1): # JFH 09-Oct-15
    # for bin in range(nbins):  # JFH 09-Oct-15
        curve_vals[bin] = []
        for curve in range(ncurves):
            #            funcVal = curves[curve].Eval(binCenters[bin])
            funcVal = curves[curve].Eval(ptsForGrph[bin])
            hists[bin].Fill(funcVal)
            curve_vals[bin].append(funcVal)
        
    # For high and low systematic error curves, fill lists of rms values from
    # bin-by-bin histograms
    list_mid = []
    list_hi = []
    list_lo = []
    list_rms = []

    gr = ROOT.TGraphErrors(nbins)
    # Using RMS
    useRMS = False
    usePercentile = True
    #useRMS = True
    #usePercentile = False
    #    returnDiff = False
    returnRMS = False
    if useRMS:
        for bin in range(nbins):
            funcVal = func.Eval(ptsForGrph[bin])
            rms = hists[bin].GetRMS()
            list_mid.append(funcVal)
            if returnDiff: 
                list_hi.append(rms)
                list_lo.append(rms)
            else:
                list_hi.append(funcVal+rms)
                list_lo.append(funcVal-rms)

            list_rms.append(rms)
    if usePercentile:
        hiQ = int(ncurves*0.84) # 84 = 50+68/2
        loQ = int(ncurves*0.26) # 26 = 50-68/2
        if (ncurves*0.84)%1 != 0. or (ncurves*0.26)%1 != 0.:
            print "Percentile calculation is not using integers, just thought you should know."
            print "Use 100 or 1000 curves to fix this."
        for bin in range(nbins):
            funcVal = func.Eval(ptsForGrph[bin]) # JFH 09-Oct-15
            list_mid.append(funcVal)             # JFH 09-Oct-15
            curve_vals[bin].sort()
            # Return difference
            if returnDiff:
                list_hi.append(curve_vals[bin][hiQ]-funcVal)
                list_lo.append(funcVal-curve_vals[bin][loQ])
            # Return absolute points
            else:
                list_hi.append(curve_vals[bin][hiQ])
                list_lo.append(curve_vals[bin][loQ])
        

    # Make arrays for ROOT
    binCenterArray = array ('d', ptsForGrph)
    ylist = []
    xerr = []
    for ibin in range(nbins):
        ylist.append(0.)
        xerr.append(0.)
    array_y  = array ('d', ylist)
    array_ex = array ('d', xerr)
        
    array_mid = array ('d', list_mid)
    array_hi  = array ('d', list_hi )
    array_lo  = array ('d', list_lo )
    if returnRMS:
        array_rms = array ('d', list_rms)

    # Convert arrays to TGraphs
    #    print "len(binCenterArray), binCenterArray, array_mid = ", len(binCenterArray), binCenterArray, array_mid
    graph_mid = ROOT.TGraph(len(binCenterArray), binCenterArray, array_mid)
    graph_hi  = ROOT.TGraph(len(binCenterArray), binCenterArray, array_hi)
    graph_lo  = ROOT.TGraph(len(binCenterArray), binCenterArray, array_lo)
    if returnRMS:
        graph_rms = ROOT.TGraphErrors(len(binCenterArray), binCenterArray, array_y, array_ex, array_rms)

    if returnRMS:
        return graph_rms
    else:
        return (graph_mid, graph_lo, graph_hi)

def getDoubleRatioPlot(dr_graphs, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, fitFunc=None, returnDiff=None, drawText=None, text=None, textCoords=None, addDeviationInQuad=None):
    """takes as input a TGraph (or list of TGraphs) returns canvas with fit and error
    default is linear fit"""
    
    setTdrStyle()

    ###############################################################################        
    # set default parameters
    ###############################################################################        

    if(Title==None):
        Title=""
    if(yTitle==None):
        yTitle = "R^{data}_{Z_{ll}/#gamma}/R^{sim}_{Z_{ll}/#gamma}"
    if(doCMSlumi==None):
        doCMSlumi=True
    if(iPos==None):
        iPos=11
    if(extraText==None):
        extraText="Preliminary"
    if(fitFunc==None):
        fitFunc = "pol1(0)"
    if(returnDiff==None):
        returnDiff = False
    if((type(dr_graphs) is ROOT.TGraphAsymmErrors) or (type(dr_graphs) is ROOT.TGraph)):
        dr_graphs = [dr_graphs]
    if(drawText==False):
        drawText = False
    if(text==None):
        text = "arXiv:1602.06581"
    if(textCoords==None):
        textCoords = [0.65,0.85,.85,.90]
    if(addDeviationInQuad==None):
        addDeviationInQuad=True
    ###############################################################################        
    # end set default parameters
    ###############################################################################        

    it = 0
    c = []
    line = []

    errorDict = {}

    for graph in dr_graphs:

        ROOT.gStyle.SetOptFit(0000)
        nbins = graph.GetXaxis().GetNbins()
        xlow = graph.GetXaxis().GetBinLowEdge(1)
        if(xlow==0.0):
            xlow = -0.5
        xhigh = graph.GetXaxis().GetBinUpEdge(nbins)

        #meanDR = graph.GetMean(2)

        func_mid  = ROOT.TF1("func_mid", fitFunc, xlow, xhigh)

        line.append(ROOT.TF1("ave", 'pol0(0)', xlow, xhigh))
        graph.Fit(line[it],"QN")

        #line.append(ROOT.TLine(xlow,meanDR,xhigh,meanDR))

        fitGraph = getDoubleRatioFit(graph, func_mid, returnDiff)

        cIter = it
        while(type(ROOT.gROOT.FindObject("c"+str(cIter+1)))==ROOT.TCanvas):
            cIter+=1
        c.append(ROOT.TCanvas("c"+str(cIter+1),"c"+str(cIter+1),640*it+40,500,600,600))

        fitGraph[0].Draw('al')
        fitGraph[1].Draw('l same')
        fitGraph[2].Draw('l same')
        line[it].Draw('l same')
        graph.Draw('p same')

        graph.SetMarkerStyle(20)
        graph.SetLineColor(1)

        fitGraph[0].SetMaximum(2.)
        fitGraph[0].SetMinimum(0.)
        fitGraph[0].GetXaxis().SetLimits(xlow,xhigh)
        
        fitGraph[1].SetLineStyle(2)
        fitGraph[2].SetLineStyle(2)

        fitGraph[0].SetLineWidth(2)
        fitGraph[1].SetLineWidth(2)
        fitGraph[2].SetLineWidth(2)

        fitGraph[0].SetLineColor(4)
        fitGraph[1].SetLineColor(4)
        fitGraph[2].SetLineColor(4)

        line[it].SetLineWidth(2)
        line[it].SetLineStyle(2)
        line[it].SetLineColor(2)

        error1D = []
        errorDev = []
        if(xTitle==None and 'NJets' in graph.GetTitle()):
            xTitle = "N_{jet}"

            bin1 = 2.            
            h_bin2 = getDist("photon","NJets",distRange=[2.5,4.5],doVarBinning=False)
            bin2 = h_bin2.GetMean()
            h_bin3 = getDist("photon","NJets",distRange=[4.5,6.5],doVarBinning=False)
            bin3 = h_bin3.GetMean()
            h_bin4 = getDist("photon","NJets",distRange=[6.5,8.5],doVarBinning=False)
            bin4 = h_bin4.GetMean()
            h_bin5 = getDist("photon","NJets",distRange=[8.5,12.5],doVarBinning=False)
            bin5 = h_bin5.GetMean()
            error1D.append(((fitGraph[2].Eval(bin1)-fitGraph[0].Eval(bin1))/line[it].Eval(bin1),(fitGraph[0].Eval(bin1)-fitGraph[1].Eval(bin1))/line[it].Eval(bin1)))
            error1D.append(((fitGraph[2].Eval(bin2)-fitGraph[0].Eval(bin2))/line[it].Eval(bin2),(fitGraph[0].Eval(bin2)-fitGraph[1].Eval(bin2))/line[it].Eval(bin2)))
            error1D.append(((fitGraph[2].Eval(bin3)-fitGraph[0].Eval(bin3))/line[it].Eval(bin3),(fitGraph[0].Eval(bin3)-fitGraph[1].Eval(bin3))/line[it].Eval(bin3)))
            error1D.append(((fitGraph[2].Eval(bin4)-fitGraph[0].Eval(bin4))/line[it].Eval(bin4),(fitGraph[0].Eval(bin4)-fitGraph[1].Eval(bin4))/line[it].Eval(bin4)))
            error1D.append(((fitGraph[2].Eval(bin5)-fitGraph[0].Eval(bin5))/line[it].Eval(bin5),(fitGraph[0].Eval(bin5)-fitGraph[1].Eval(bin5))/line[it].Eval(bin1)))
            errorDev.append(abs(fitGraph[0].Eval(bin1)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin2)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin3)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin4)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin5)-line[it].Eval(bin1))/line[it].Eval(bin1))

            if(addDeviationInQuad):
                for i in range(len(error1D)):
                    error1D[i] = ((errorDev[i]**2+error1D[i][0]**2)**0.5,(errorDev[i]**2+error1D[i][1]**2)**0.5)                    

            errorDict['NJets'] = error1D

            drVal = line[it].Eval(bin1)
            drErr = line[it].GetParError(0)/drVal

        if(xTitle==None and 'BTags' in graph.GetTitle()):
            xTitle = "N_{b-Jets}"
            fitGraph[0].SetMaximum(2.1)

            
            bin1 = 0.0
            bin2 = 1.0
            bin3 = 2.0
            h_bin4 = getDist("photon","BTags",distRange=[2.5,6.5],doVarBinning=False)
            bin4 = h_bin4.GetMean()
            
            errorDev.append(abs(fitGraph[0].Eval(bin1)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin2)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin3)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin4)-line[it].Eval(bin1))/line[it].Eval(bin1))


            if(addDeviationInQuad):
                for i in range(len(error1D)):
                    error1D[i] = ((errorDev[i]**2+error1D[i][0]**2)**0.5,(errorDev[i]**2+error1D[i][1]**2)**0.5)
                    
            errorDict['BTags'] = error1D
            drVal = line[it].Eval(bin1)
            drErr = line[it].GetParError(0)/drVal

        if(xTitle==None and 'MHT' in graph.GetTitle()):
            xTitle = "H_{T}^{miss} [GeV]"
            
            h_bin1 = getDist("photon","MHT",distRange=[300,350],doVarBinning=False)
            bin1 = h_bin1.GetMean()
            h_bin2 = getDist("photon","MHT",distRange=[350,500],doVarBinning=False)
            bin2 = h_bin2.GetMean()
            h_bin3 = getDist("photon","MHT",distRange=[500,750],doVarBinning=False)
            bin3 = h_bin3.GetMean()
            h_bin4 = getDist("photon","MHT",distRange=[750,1500],doVarBinning=False)
            bin4 = h_bin4.GetMean()
            h_bin5 = getDist("photon","MHT",distRange=[250,300],doVarBinning=False,applyMHTCut=False)
            bin5 = h_bin4.GetMean()

            error1D.append(((fitGraph[2].Eval(bin1)-fitGraph[0].Eval(bin1))/line[it].Eval(bin1),(fitGraph[0].Eval(bin1)-fitGraph[1].Eval(bin1))/line[it].Eval(bin1)))
            error1D.append(((fitGraph[2].Eval(bin2)-fitGraph[0].Eval(bin2))/line[it].Eval(bin2),(fitGraph[0].Eval(bin2)-fitGraph[1].Eval(bin2))/line[it].Eval(bin2)))
            error1D.append(((fitGraph[2].Eval(bin3)-fitGraph[0].Eval(bin3))/line[it].Eval(bin3),(fitGraph[0].Eval(bin3)-fitGraph[1].Eval(bin3))/line[it].Eval(bin3)))
            error1D.append(((fitGraph[2].Eval(bin4)-fitGraph[0].Eval(bin4))/line[it].Eval(bin4),(fitGraph[0].Eval(bin4)-fitGraph[1].Eval(bin4))/line[it].Eval(bin4)))
            error1D.append(((fitGraph[2].Eval(bin5)-fitGraph[0].Eval(bin5))/line[it].Eval(bin5),(fitGraph[0].Eval(bin5)-fitGraph[1].Eval(bin5))/line[it].Eval(bin5)))
            
            errorDev.append(abs(fitGraph[0].Eval(bin1)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin2)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin3)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin4)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin5)-line[it].Eval(bin1))/line[it].Eval(bin1))

            if(addDeviationInQuad):
                for i in range(len(error1D)):
                    error1D[i] = ((errorDev[i]**2+error1D[i][0]**2)**0.5,(errorDev[i]**2+error1D[i][1]**2)**0.5)
                    
            errorDict['MHT'] = error1D
            drVal = line[it].Eval(bin1)
            drErr = line[it].GetParError(0)/drVal

        if(xTitle==None and 'HT' in graph.GetTitle()):
            xTitle = "H_{T} [GeV]"
            
            h_bin1 = getDist("photon","HT",distRange=[300,500],doVarBinning=False)
            bin1 = h_bin1.GetMean()
            h_bin2 = getDist("photon","HT",distRange=[500,1000],doVarBinning=False)
            bin2 = h_bin2.GetMean()
            h_bin3 = getDist("photon","HT",distRange=[1000,2500],doVarBinning=False)
            bin3 = h_bin3.GetMean()
            h_bin4 = getDist("photon","HT",distRange=[350,500],doVarBinning=False)
            bin4 = h_bin4.GetMean()
            h_bin5 = getDist("photon","HT",distRange=[750,1500],doVarBinning=False)
            bin5 = h_bin5.GetMean()
            h_bin6 = getDist("photon","HT",distRange=[1500,2500],doVarBinning=False)
            bin6 = h_bin6.GetMean()
            
            error1D.append(((fitGraph[2].Eval(bin1)-fitGraph[0].Eval(bin1))/line[it].Eval(bin1),(fitGraph[0].Eval(bin1)-fitGraph[1].Eval(bin1))/line[it].Eval(bin1)))
            error1D.append(((fitGraph[2].Eval(bin2)-fitGraph[0].Eval(bin2))/line[it].Eval(bin2),(fitGraph[0].Eval(bin2)-fitGraph[1].Eval(bin2))/line[it].Eval(bin2)))
            error1D.append(((fitGraph[2].Eval(bin3)-fitGraph[0].Eval(bin3))/line[it].Eval(bin3),(fitGraph[0].Eval(bin3)-fitGraph[1].Eval(bin3))/line[it].Eval(bin3)))
            error1D.append(((fitGraph[2].Eval(bin4)-fitGraph[0].Eval(bin4))/line[it].Eval(bin4),(fitGraph[0].Eval(bin4)-fitGraph[1].Eval(bin4))/line[it].Eval(bin4)))
            error1D.append(((fitGraph[2].Eval(bin5)-fitGraph[0].Eval(bin5))/line[it].Eval(bin5),(fitGraph[0].Eval(bin5)-fitGraph[1].Eval(bin5))/line[it].Eval(bin5)))
            error1D.append(((fitGraph[2].Eval(bin6)-fitGraph[0].Eval(bin6))/line[it].Eval(bin6),(fitGraph[0].Eval(bin6)-fitGraph[1].Eval(bin6))/line[it].Eval(bin6)))

            errorDev.append(abs(fitGraph[0].Eval(bin1)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin2)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin3)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin4)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin5)-line[it].Eval(bin1))/line[it].Eval(bin1))
            errorDev.append(abs(fitGraph[0].Eval(bin6)-line[it].Eval(bin1))/line[it].Eval(bin1))

            if(addDeviationInQuad):
                for i in range(len(error1D)):
                    error1D[i] = ((errorDev[i]**2+error1D[i][0]**2)**0.5,(errorDev[i]**2+error1D[i][1]**2)**0.5)
                    
            errorDict['HT'] = error1D
            drVal = line[it].Eval(bin1)
            drErr = line[it].GetParError(0)/drVal

        if(xTitle==None):
            xTitle = graph.GetTitle()

        fitGraph[0].GetXaxis().SetTitle(xTitle)
        xTitle=None
        fitGraph[0].GetYaxis().SetTitle(yTitle)

        fitGraph[0].SetTitle("")

        c[it].SetRightMargin(0.05)

        cmsLumi(c[it], iPeriod, iPos, extraText)

        pt = ROOT.TPaveText(textCoords[0],textCoords[1],textCoords[2],textCoords[3],"brNDC")
        pt.AddText(text)
        pt.SetFillColor(0)
        if(drawText==True):
            pt.Draw("NB")

        ROOT.SetOwnership(c[it],0)
        ROOT.SetOwnership(pt,0)

        for g in fitGraph:
            ROOT.SetOwnership(g,0)

        it+=1

    for l in line:
        ROOT.SetOwnership(l, 0)

    

    return [(drVal,drErr),errorDict]


def getZmassFitPlot(fitFunc=None, dataSet=None, mcSet=None, plotMC=None, doDiMu=None, doDiEl=None, doDiLep=None, getShapeFromLoose=None, nBins=None, distRange=None, nJetBin=None, bJetBin=None, kinBin=None, doVarBinning=None, binning=None, extraCuts=None, dphiCut=None, doLumi=None, do20=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, keepCanvas=None, drawText=None, text=None, textCoords=None, drawText2=None, text2=None, textCoords2=None):
    """getZmassFitPlot() returns diMu and diEl mass plots (canvas's) incuding 
    fits to dataset bin with bJetBin, nJetBin, or kinBin, default is -1 (inclusive)
    set doDiLep=True if you want to merge mumu ee mass distributions
    predefined Z fit functions: (fitFunc parameter) 
    'CB'-crystal ball
    'V'-Voigtian i.e. breit-wigner convolved with Gaussian
    'DG'-Double Gaussian
    'TG'-Triple Gaussian"""

    setTdrStyle()
    
    if(fitFunc==None):
        fitFunc='V'
    if(keepCanvas==None):
        keepCanvas=True
    if(do20==None):
        do20=False
    sampleSuffix=''
    if(do20):
        sampleSuffix='20'
        if(doLumi==None):
            doLumi=24.50
        if(extraText==None):
            extraText="Simulation"
    if(text==None):
        text = "arXiv:1602.06581"
    elif(drawText==None):
        drawText = True
    if(drawText==None):
        drawText=False
    if(textCoords==None):
        textCoords = [0.20,0.70,.40,.75]
    if(text2==None):
        text2 = "N_{b}>=0"
    elif(drawText2==None):
        drawText2 = True
    if(drawText2==None):
        drawText2=False
    if(textCoords2==None):
        textCoords2 = [0.20,0.77,.40,.82]
    if(iPeriod==None):
        iPeriod=6
    if(doCMSlumi==None):
        doCMSlumi=True
    if(iPos==None):
        iPos=11
    if(extraText==None):
        extraText="Preliminary"
    if(getShapeFromLoose==None):
        getShapeFromLoose=True
    if(doDiMu==None):
        doDiMu=True
    if(doDiEl==None):
        doDiEl=True
    if(doDiLep==None):
        doDiLep=False
    if(doVarBinning==None):
        doVarBinning=False
    if(doLumi==None):
        doLumi=24.5
    if(binning==None):
        binning=[60.,70.,80.,90.,100.,110.,120.]
    if(extraCuts==None):
        extraCuts=""
    if(dphiCut==None):
        dphiCut='sig'
    if(distRange==None):
        distRange = [60.,120.]
    if(nBins==None):
        nBins = 30
    if(nJetBin==None):
        nJetBin = -1
    if(bJetBin==None):
        bJetBin = -1
    if(kinBin==None):
        kinBin = -1
    if(dataSet==None):
        d1 = {}
        d2 = {}
        d3 = {}
        dataSet=[]
        nEvents=[]
        if(doDiMu==True):
            h_zmmL = getDist("zmm"+sampleSuffix,"ZCandidates.M()",applyNJetsCut=False,  nBins=nBins, distRange=distRange)
            d1['loose'] = h_zmmL
        if(doDiEl==True):
            h_zeeL = getDist("zee"+sampleSuffix,"ZCandidates.M()",applyNJetsCut=False,  nBins=nBins, distRange=distRange)
            d2['loose'] = h_zeeL
        if(doDiLep==True):
            h_zllL = getDist("zll"+sampleSuffix,"ZCandidates.M()",applyNJetsCut=False,  nBins=nBins, distRange=distRange)
            d3['loose'] = h_zllL
        if(doDiMu==True):
            h_zmmS = getDist("zmm"+sampleSuffix,"ZCandidates.M()", distRange=distRange, nBins=nBins, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts)
            d1['sig'] = h_zmmS
            h_zmmE = getDist("zmm"+sampleSuffix, "ZCandidates.M()", distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts,applyMassCut=True)
            nEvents.append(h_zmmE)
            dataSet.append(d1)
        if(doDiEl==True):
            h_zeeS = getDist("zee"+sampleSuffix,"ZCandidates.M()", distRange=distRange, nBins=nBins, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts)
            d2['sig'] = h_zeeS
            h_zeeE = getDist("zee"+sampleSuffix, "ZCandidates.M()", distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts,applyMassCut=True)
            nEvents.append(h_zeeE)
            dataSet.append(d2)
        if(doDiLep==True):
            h_llS = getDist("zll"+sampleSuffix,"ZCandidates.M()", distRange=distRange, nBins=nBins, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts)
            d3['sig'] = h_llS
            h_zllE = getDist("zll"+sampleSuffix, "ZCandidates.M()", distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts,applyMassCut=True)
            nEvents.append(h_zllE)
            dataSet.append(d3)
    elif(mcSet==None):
        if(plotMC==None):
            plotMC=False ## user gives dataSet, can't predict corresponding MC set

    if(type(dataSet) is not list):
        dataSet = [dataSet]
    for i in range(len(dataSet)):
        if(type(dataSet[i]) is not dict):
            key = "k"+str(i)
            dataSet[i] = {key: dataSet[i]}
    if(plotMC==None):
        plotMC=True
    mcList = []
    if(plotMC):
        if(mcSet==None):
            mcSet = []
            if(doDiMu):
                h_ttmm = getDist("ttmm","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_ttzmm = getDist("ttzmm","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_dibosmm = getDist("dibosonmm","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_dymm = getDist("dymm","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                hsMM = ROOT.THStack("hsMM","dimuon mass [GeV]")
                h_dymm.SetFillColor(2)
                h_ttzmm.SetFillColor(5)
                h_dibosmm.SetFillColor(8)
                h_ttmm.SetFillColor(38)
                hsMM.Add(h_ttmm)
                hsMM.Add(h_ttzmm)
                hsMM.Add(h_dibosmm)
                hsMM.Add(h_dymm)
                mcSet.append(hsMM)
                mcList.append(h_dymm)
                mcList.append(h_ttzmm)
                mcList.append(h_dibosmm)
                mcList.append(h_ttmm)

            if(doDiEl):
                h_ttee = getDist("ttee","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_ttzee = getDist("ttzee","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_dibosee = getDist("dibosonee","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_dyee = getDist("dyee","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                hsEE = ROOT.THStack("hsEE","dielectron mass [GeV]")
                h_dyee.SetFillColor(2)
                h_ttzee.SetFillColor(5)
                h_dibosee.SetFillColor(8)
                h_ttee.SetFillColor(38)
                hsEE.Add(h_ttee)
                hsEE.Add(h_ttzee)
                hsEE.Add(h_dibosee)
                hsEE.Add(h_dyee)
                mcSet.append(hsEE)
                mcList.append(h_dyee)
                mcList.append(h_ttzee)
                mcList.append(h_dibosee)
                mcList.append(h_ttee)

            if(doDiLep):
                h_ttll = getDist("ttll","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_ttzll = getDist("ttzll","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_dibosll = getDist("dibosonll","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                h_dyll = getDist("dyll","ZCandidates.M()", nBins=nBins, distRange=distRange, nJetBin=nJetBin, bJetBin=bJetBin, kinBin=kinBin, extraCuts=extraCuts, doLumi=doLumi)
                hsLL = ROOT.THStack("hsLL","dilepton mass [GeV]")            
                h_dyll.SetFillColor(2)
                h_ttzll.SetFillColor(5)
                h_dibosll.SetFillColor(8)
                h_ttll.SetFillColor(38)
                hsLL.Add(h_ttll)
                hsLL.Add(h_ttzll)
                hsLL.Add(h_dibosll)
                hsLL.Add(h_dyll)
                mcSet.append(hsLL)
                mcList.append(h_dyll)
                mcList.append(h_ttzll)
                mcList.append(h_dibosll)
                mcList.append(h_ttll)

    canvas = []
    purityList = []

    j = -1 
    for data in dataSet:
        j+=1

        keyList = data.keys()

        it = 0
        frame = {}

        eventsPerBin = int((distRange[1]-distRange[0])/nBins)

        xTitle = 'rooVarX'
        
        yTitle = 'Events / '+str(eventsPerBin)+' GeV'

        if(plotMC):
            xTitle = mcSet[j].GetTitle()

        rooVarX = ROOT.RooRealVar("rooVarX",xTitle,distRange[0],distRange[1])
        rooArgL = ROOT.RooArgList(rooVarX)
        
        # Cheby polynomial for bg
        c0 = ROOT.RooRealVar("c0","c0",-0.599,-1,1)
        c1 = ROOT.RooRealVar("c1","c1",-0.099,-.1,.1)
        c2 = ROOT.RooRealVar("c2","c2",0.099,-.2,.2)
        bfit = ROOT.RooChebychev("background","Chebychev PDF",rooVarX,ROOT.RooArgList(c0,c1,c2))
                
        # Crystal Ball
        cbm = ROOT.RooRealVar("cbm","m in crystal ball",93.,88.,96.)
        cbs = ROOT.RooRealVar("cbs","s in crystal ball",3.38,2.0,4.0)
        cba = ROOT.RooRealVar("cba","a in crystal ball",3.38,0.1,5)
        cbn = ROOT.RooRealVar("cbn","n in crystal ball",10.0,1.0,100.)
        zfitCB = ROOT.RooCBShape("CBShape","CB PDF",rooVarX,cbm,cbs,cba,cbn)

        # BW convolved with Gaussian
        vm = ROOT.RooRealVar("vm","m in voigtian",91.19,90.,92.)
        vg = ROOT.RooRealVar("vg","g in voigtian",2.64,1.0,3.0)
        vs = ROOT.RooRealVar("vs","s in voigtian",1.45,1.0,2.0)
        zfitV = ROOT.RooVoigtian("voigtian","Voigtian PDF",rooVarX,vm,vg,vs)

        # double Gaussian
        meanCore = ROOT.RooRealVar("meanCore","meanCore",90.9,80,100)
        widthCore = ROOT.RooRealVar("widthCore","widthCore",2.0,1.0,5.0)
        meanTail = ROOT.RooRealVar("meanTail","meanTail",86.5,80,100)
        widthTail= ROOT.RooRealVar("widthTail","widthTail",5.0,2.0,15.0)
        fracZ = ROOT.RooRealVar("fracZ","fracZ",0.65,0.40,1.0)
        zfitCore = ROOT.RooGaussian("zfitCore","zfitCore",rooVarX,meanCore,widthCore)
        zfitTail = ROOT.RooGaussian("zfitTail","zfitTail",rooVarX,meanTail,widthTail)
        zfitDG = ROOT.RooAddPdf("zfitDG","zfitDG",ROOT.RooArgList(zfitCore,zfitTail),ROOT.RooArgList(fracZ))

        # triple Gaussian
        meanOut = ROOT.RooRealVar("meanOut","meanOut",90.1,80,100)
        widthOut = ROOT.RooRealVar("widthOut","widthOut",15.0,3.0,40.0)
        zfitOut = ROOT.RooGaussian("zfitOut","zfitOut",rooVarX,meanOut,widthOut)
        fracZout= ROOT.RooRealVar("fracZout","fracZout",0.05,0.00,0.30)
        zfitTG = ROOT.RooAddPdf("zfitTG","zfitTG",ROOT.RooArgList(zfitCore,zfitTail,zfitOut),ROOT.RooArgList(fracZ,fracZout))

        fitFuncDict = {
            'CB': zfitCB,
            'V': zfitV,
            'DG': zfitDG,
            'TG': zfitTG,
        }

        correctionFactor = {
            'CB': 0.9246,
            'V': 0.9780,
            'DG': 0.9478,
            'TG': 0.9505,
        }

        for key in data:

            Data = ROOT.RooDataHist('d','d',rooArgL,data[key])
        

            if(key == keyList[0]):
                c0.setConstant(False)
                c1.setConstant(False)
                c2.setConstant(False)
                cbm.setConstant(False)
                cbs.setConstant(False)
                cba.setConstant(False)
                cbn.setConstant(False)
                meanCore.setConstant(False)
                widthCore.setConstant(False)
                meanTail.setConstant(False)
                widthTail.setConstant(False)
                meanOut.setConstant(False)
                widthOut.setConstant(False)
                fracZ.setConstant(False)
                fracZout.setConstant(False)
                vg.setConstant(False)
                vs.setConstant(False)
            elif(getShapeFromLoose):
                c0.setConstant(True)
                c1.setConstant(True)
                c2.setConstant(True)
                cbm.setConstant(True)
                cbs.setConstant(True)
                cba.setConstant(True)
                cbn.setConstant(True)
                meanCore.setConstant(True)
                widthCore.setConstant(True)
                meanTail.setConstant(True)
                widthTail.setConstant(True)
                meanOut.setConstant(True)
                widthOut.setConstant(True)
                fracZ.setConstant(True)
                fracZout.setConstant(True)
                vm.setConstant(True)
                vg.setConstant(True)
                vs.setConstant(True)


            nsig = ROOT.RooRealVar("nsig","nsig",10000,0,10000000)
            nbkg = ROOT.RooRealVar("nbkg","nbkg",1000,0,10000000)
            tot = ROOT.RooRealVar("tot","tot",10000,0,10000000)

            frac = ROOT.RooRealVar("frac","frac",0.9,0.0,1.0)
            fsig = ROOT.RooRealVar("fsig","fsig",0.7,0.2,1.0)

            model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(fitFuncDict[fitFunc],bfit),ROOT.RooArgList(nsig,nbkg))
            model2 = ROOT.RooAddPdf("model2","model2",ROOT.RooArgList(fitFuncDict[fitFunc],bfit),ROOT.RooArgList(frac))

            if(it==0):
                model.fitTo(Data)
            else:
                model.fitTo(Data)
                model2.fitTo(Data)
                
                rooVarX.setRange("all",60,120)
                rooVarX.setRange("selection",76.19,106.19)

                selFraction = frac.getVal()
                selFractionE = frac.getError()
                allSel = fitFuncDict[fitFunc].createIntegral(ROOT.RooArgSet(rooArgL),"all").getVal()
                selectionSel = fitFuncDict[fitFunc].createIntegral(ROOT.RooArgSet(rooArgL),"selection").getVal()
                allBGSel = bfit.createIntegral(ROOT.RooArgSet(rooArgL),"all").getVal()
                selectionBGSel = bfit.createIntegral(ROOT.RooArgSet(rooArgL),"selection").getVal()

                nsigLVLB = nsig.getVal()
                nbkgLVLB = nbkg.getVal()

                integralSigSel = ((selectionSel/allSel)*nsigLVLB)
                integralBgSel = ((selectionBGSel/allBGSel)*nbkgLVLB)

                selentries = data[key].GetEntries()
                selSelentries = nEvents[j].GetEntries()
                selFractionENaive = ROOT.TMath.Sqrt( (selFraction*(1-selFraction))/selentries )
                selFractionSel = integralSigSel/(integralSigSel+integralBgSel)
                selFractionESelNaive = ROOT.TMath.Sqrt( (selFractionSel*(1-selFractionSel))/(float(selSelentries) ))
                selFractionESel = selFractionESelNaive*(selFractionE/selFractionENaive)

                purity = (1/correctionFactor[fitFunc])*integralSigSel/(integralSigSel+integralBgSel)
                Epurity= (1/correctionFactor[fitFunc])*selFractionESel 
                
                purityList.append((purity,Epurity))

            ###############################
            # unique canvas name
            ##############################
            canvasName="canvas"

            cIter = 1
            while(type(ROOT.gROOT.FindObject(canvasName+str(cIter)))==ROOT.TCanvas):
                cIter+=1
            canvasName+=str(cIter)
                
            canvas.append(ROOT.TCanvas(canvasName, canvasName, 20+it*600, j*600, 500,500))

            canvas[it+2*j].SetRightMargin(0.05)

            if(plotMC):
                if(key == keyList[len(keyList)-1]):
                    mcSet[j].Draw("HIST")
                    mcSet[j].GetXaxis().SetTitle(xTitle)
                    mcSet[j].GetYaxis().SetTitle(yTitle)
                    if(data[key].GetMaximum()>mcSet[j].GetMaximum()):
                        mcSet[j].SetMaximum(1.25*data[key].GetMaximum())
                    ROOT.SetOwnership(mcSet[j], 0)
                    

            frame[key] = rooVarX.frame()
            nbins = data[key].GetNbinsX()
            interval = int((data[key].GetBinLowEdge(nbins+1)-data[key].GetBinLowEdge(1))/nbins)
            #frame[key].GetYaxis.SetTitle("Events / "+str(interval)+" GeV")
            Data.plotOn(frame[key])
            #model.plotOn(frame[key],ROOT.RooFit.Name("fit"),ROOT.RooFit.LineColor(4))
   
            dataStr = 'Data'
            if(do20):
                dataStr = '\"Data\"'
            if(it==0):
                model.plotOn(frame[key],ROOT.RooFit.Name("fit"))
                model.plotOn(frame[key],ROOT.RooFit.Components("background"),ROOT.RooFit.LineColor(4),ROOT.RooFit.LineStyle(2))
                # leg0 = ROOT.TLegend(0.23,0.79,0.65,0.88, "","brNDC")
                # leg0.SetBorderSize(2)
                # leg0.SetFillStyle(1001)
                # leg0.SetFillColor(ROOT.kWhite) 
   
                # leg0.AddEntry(data[key], dataStr+" (loosened selection)", "p")
                # leg0.AddEntry(frame[key].findObject("fit"), "Fit to "+dataStr, "l")
   
                # leg0.Draw("NB")
                # ROOT.SetOwnership(leg0,0)

            else:
                model2.plotOn(frame[key],ROOT.RooFit.Name("fit"))
                model2.plotOn(frame[key],ROOT.RooFit.Components("background"),ROOT.RooFit.LineColor(4),ROOT.RooFit.LineStyle(2))

                leg1 = ROOT.TLegend(0.65,0.60,0.85,0.88, "","brNDC")
                leg1.SetBorderSize(2)
                leg1.SetFillStyle(1001)
                leg1.SetFillColor(ROOT.kWhite) 
                leg1.AddEntry(mcList[0],"Drell-Yan MC","f")
                leg1.AddEntry(mcList[1],"t#bar{t}Z MC","f")
                leg1.AddEntry(mcList[2],"VV MC","f")
                leg1.AddEntry(mcList[3],"t#bar{t} MC","f")
                leg1.AddEntry(data[key], dataStr, "p")
                leg1.AddEntry(frame[key].findObject("fit"), "Fit to "+dataStr, "l")

                leg1.Draw("NB")
                ROOT.SetOwnership(leg1,0)
            frame[key].Draw("same")

            pt = ROOT.TPaveText(textCoords[0],textCoords[1],textCoords[2],textCoords[3],"brNDC")
            pt.AddText(text)
            pt.SetFillColor(0)
            if(drawText==True):
                pt.Draw("NB")
                ROOT.SetOwnership(pt,0)
            pt2 = ROOT.TPaveText(textCoords2[0],textCoords2[1],textCoords2[2],textCoords2[3],"brNDC")
            pt2.AddText(text2)
            pt2.SetFillColor(0)
            if(drawText2==True):
                pt2.Draw("NB")
                ROOT.SetOwnership(pt2,0)

            cmsLumi(canvas[it+2*j], iPeriod, iPos, extraText)
 
            it+=1

    for c in canvas:
        ROOT.SetOwnership(c, int(not keepCanvas))

    return purityList
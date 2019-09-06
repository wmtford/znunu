"""Module containing tools for Z->nunu estimation studies.
Designed to run on LPC produced skims
T. Mulholland: troy.mulholland@cern.ch
This is a subset of RA2b.py; Functions that project ntuple onto
histograms have been stripped out.  -- W. T. Ford"""


import time
import re
import math
import random
import numpy
from array import array
import ROOT

def cmsLumi(pad,  iPeriod=None,  iPosX=None, extraText=None, cmsText=None):

    """ modifed by TM 
    modified form of CMS_lumi
    Initiated by: Gautier Hamel de Monchenault (Saclay)
    Translated in Python by: Joshua Hardenbrook (Princeton)
    Updated by:   Dinko Ferencek (Rutgers) """

    if(cmsText==None):
        cmsText     = "CMS"
    cmsTextFont   = 61  

    writeExtraText = True
    if(extraText==None):
        extraText   = "Preliminary"
    if(iPeriod==None):
        iPeriod = 5
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
    lumi_13TeV_V8 = "4.0 fb^{-1}"
    lumi_13TeV_2016 = "35.9 fb^{-1}"
    lumi_13TeV_2017 = "41.5 fb^{-1}"
    lumi_13TeV_2018 = "59.6 fb^{-1}"
    lumi_13TeV_Run2 = "137 fb^{-1}"
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
        lumiText += lumi_13TeV_2017
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
    elif ( iPeriod==8 ):
        lumiText += lumi_13TeV_2018
        lumiText += " (13 TeV)"
    elif ( iPeriod==9 ):
        lumiText += lumi_13TeV_Run2
        lumiText += " (13 TeV)"
    elif ( iPeriod==12 ):
        lumiText += "8 TeV"
    elif ( iPeriod==0 ):
        lumiText += lumi_sqrtS
            
    #print lumiText

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    # latex.SetTextColor(ROOT.kBlack)    
    latex.SetTextColor(1)    
    
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
    # tdrStyle.SetCanvasColor(ROOT.kWhite)
    tdrStyle.SetCanvasColor(0)
    tdrStyle.SetCanvasDefH(600) #Height of canvas
    tdrStyle.SetCanvasDefW(600) #Width of canvas
    tdrStyle.SetCanvasDefX(0)   #POsition on screen
    tdrStyle.SetCanvasDefY(0)
    
    tdrStyle.SetPadBorderMode(0)
    # tdrStyle.SetPadColor(ROOT.kWhite)
    tdrStyle.SetPadColor(0)
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
    # tdrStyle.SetStatColor(ROOT.kWhite)
    tdrStyle.SetStatColor(0)
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

        removedBins = getRemovedBins()
        subtractBins = removedBins[0]
        avoidBins = removedBins[1]

        Bin = 1
        uncutBin=0
        for nj in njRange:
            for nb in nbRange:
                for kin in kinRange:
                    uncutBin+=1
                    if(uncutBin in avoidBins):
                        continue
                    if( (nj in njSet) and (nb in nbSet) and (kin in kinSet)):
                        h1.SetBinError(Bin,error)
                    Bin+=1
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
        
    h = ROOT.TH1D()

    if(doFlip):
        h = denomHist.Clone()
        h.Divide(numHist)
    else:
        h = numHist.Clone()
        h.Divide(denomHist)

    return h

def getDiffHist(obsHist, expHist, doFlip=None, divide=None, setZero=None):
    if(doFlip==None):
        doFlip=False
        
    if(divide==None):
        divide=True
    if('TGraph' in str(type(expHist))):
        h = ROOT.TGraphAsymmErrors()
    else:
        h = ROOT.TH1D()

    if(doFlip):
        h = expHist.Clone()
        h.Add(obsHist,-1)
        if(divide):
            h.Divide(obsHist)
    else:
        if('TGraph' in str(type(expHist))):
            for i in range(1,obsHist.GetNbinsX()+1):
                binWidth = obsHist.GetBinWidth(1)
                binStart = (obsHist.GetBinLowEdge(1)+obsHist.GetBinLowEdge(2))/2.
                expVal = expHist.Eval(binWidth*(i-1)+binStart)
                if expVal==0:
                    expVal=1.
                val = obsHist.GetBinContent(i)/expVal-1
                if(obsHist.GetBinContent(i)>0):
                    errDn = (val+1)*( (obsHist.GetBinError(i)/obsHist.GetBinContent(i))**2+(expHist.GetErrorYhigh(i-1)/expVal)**2 )**0.5
                    errUp = (val+1)*( (obsHist.GetBinError(i)/obsHist.GetBinContent(i))**2+(expHist.GetErrorYlow(i-1)/expVal)**2 )**0.5
                if(setZero==True and val==-1):
                    continue
                h.SetPointEYhigh(i-1, errUp)
                h.SetPointEYlow(i-1,  errDn)
                h.SetPoint(i-1,binWidth*(i-1)+binStart,val)
               
        else:
            h = obsHist.Clone()
            h2 = setError(expHist)
            h.Add(h2,-1)
            if(divide):
                h.Divide(expHist)

    return h

    #h_dr = getRatioHist(zll_r,pho_r)
    #g_dr = ROOT.TGraphAsymmErrors(h_dr)

def getPullHist(dataHist, predHist, doFlip=None, setZero=None):
    if(doFlip==None):
        doFlip=False
    if(setZero==None):
        setZero=False
        
    # h = ROOT.TH1D()
    # sig = ROOT.TH1D('sig','sig',preHist.GetNbinsX(),preHist.GetBinLowEdge(1),preHist.GetBinLowEdge(preHist.GetNbinsX()+1))
    # for i in range(1,preHist.GetNbinsX()+1):
    #     sig.SetBinContent(i,preHist.GetBinError(i))

    # if(doFlip):
    #     h = preHist.Clone()
    #     h.Add(postHist,-1)
    #     h.Divide(sig)
    # else:
    #     h = postHist.Clone()
    #     h.Add(preHist,-1)
    #     h.Divide(sig)

    # for i in range(1,preHist.GetNbinsX()+1):
    #     h.SetBinError(i,0.)

    h = dataHist.Clone()
    for i in range(1, h.GetNbinsX()+1):
        if 'TGraph' in str(type(predHist)):
            binWidth = h.GetBinWidth(1)
            binStart = (h.GetBinLowEdge(1)+h.GetBinLowEdge(2))/2.
            predVal = predHist.Eval(binWidth*(i-1)+binStart)
            num = dataHist.GetBinContent(i)-predVal
            if(num>0):
                predErr = predHist.GetErrorYhigh(i-1)
            else:
                predErr = predHist.GetErrorYlow(i-1)
        else:
            predVal = predHist.GetBinContent(i)
            predErr = predHist.GetBinError(i)
            num = dataHist.GetBinContent(i)-predVal

        denom = (predVal+predErr**2)**0.5
        
        if setZero:
            if(dataHist.GetBinContent(i)==0):
                num = 0

        h.SetBinContent(i, num/denom)
        h.SetBinError(i, 0)

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

def addGraphs(grList, binning):
    
    graph = ROOT.TGraphAsymmErrors()

    for i in range(grList[0].GetN()):
        val = 0.
        ehi = 0.
        elo = 0.
        for gr in grList:
            val+=gr.Eval(binning[i])
            ehi+=gr.GetErrorYhigh(i)**2
            elo+=gr.GetErrorYlow(i)**2
        graph.SetPoint(i, binning[i], val)
        graph.SetPointEYhigh(i, ehi**0.5)
        graph.SetPointEYlow(i,  elo**0.5)
        graph.SetPointEXhigh(i, binning[0])
        graph.SetPointEXlow(i, binning[0])

    return graph

def getPlotAndRatio(numHists, denomHists=None, bottomPlots=None, doStack=None, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, ratioTitle=None, ratioMin=None, ratioMax=None, doLogy=None, doFlip=None, doDiff=None, doPull=None, makeLeg=None, legList=None, legCoords=None, textCoords=None, canvasSize=None, canvasName=None, numColors=None, denomColor=None, numMarkers=None, denomMarker=None, markerSize=None, lineWidth=None, numDrawStyles=None, denomDrawStyle=None, drawErrorBand=None, stackColors=None, axisTitleSize=None, drawVerticalLines=None, drawHorizontalLine=None, statBox=None, drawText=None, text=None, setMax=None, setMin=None, doClosureStyle=None,errorBandColor=None,errorBandFillStyle=None,legHeader=None,nDivRatio=None,doNumFill=None, hLineVal=None, hLineColors=None,nDivX=None,ratioGridx=None,ratioGridy=None,topGridx=None,topGridy=None,doRatio=None,numFillStyles=None,numFillColors=None, drawSigText=None, errorBandHist=None,nDivY=None,numLineWidth=None,order=None):

    ROOT.gROOT.Reset()

    setTdrStyle()

    denomHist = ROOT.TH1D

    #####################################################################
    # set defaults and initialize 
    #####################################################################
    if(denomHists==None):
        doRatio=False
    if(doRatio==None):
        doRatio=True
    if(doClosureStyle==None):
        doClosureStyle=False
    if (type(numHists) is ROOT.TH1D or type(denomHists) is ROOT.TH1F):
        numHists = [numHists]
    if (type(denomHists) is ROOT.TH1D or type(denomHists) is ROOT.TH1F):
        denomHist = denomHists
        doStack=False
    elif type(denomHists) is list:
        if (type(denomHists[0]) is ROOT.TH1D or type(denomHists[0]) is ROOT.TH1F):
            denomHist = denomHists[0].Clone()
            for i in range(1,len(denomHists)):
                denomHist.Add(denomHists[i])
        elif 'TGraph' in str(type(denomHists[0])):
            denomHist = numHists[0].Clone()
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
        if(doRatio==False):
            doLogy=False
        else:
            doLogy=True
    if(doDiff==None):
        doDiff = False
    if(doDiff==True and ratioTitle==None):
        ratioTitle = "#frac{Obs.-Exp.}{Exp.}"
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
        textCoords = [0.35,0.75,.55,.88]  # wtf
        # textCoords = [0.20,0.70,.40,.83]
    if(canvasSize==None):
        canvasSize = [900,600]
    if(canvasName==None):
        canvasName = "canvas"
    if(numColors==None):
        if(doRatio==False):
            numColors = [4,8,2,1,12,7,5,13,46,41,6]
        else:
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
    if(numLineWidth==None):
        numLineWidth=1
    if(denomDrawStyle==None):
        denomDrawStyle="HIST"
    if(denomDrawStyle=="HIST"):
        denomMarker=1
    if(numDrawStyles==None):
        if(doRatio==False):
            numDrawStyles=["HIST","HIST","HIST","HIST","HIST","HIST","HIST","HIST","HIST"]
        else:
            numDrawStyles=["e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0","e0p X0"]
    if(numFillStyles==None):
        numFillStyles=[3001,3001,3001,3001,3001,3001]
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
        errorBandFillStyle=3005  # Troy mod (was 3002)
    if(nDivRatio==None):
        nDivRatio=5
    if(nDivX==None):
        nDivX=510
    if(nDivY==None):
        nDivY=510
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
        # errorBandColor = ROOT.kRed-10
        errorBandColor = 632-10
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
    if(numFillColors==None):
        numFillColors = numColors
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

    if(doRatio):
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
            numHistsRatio[i].SetLineWidth(numLineWidth)
            numHistsRatio[i].SetLineColor(numColors[i])
            numHistsRatio[i].SetMarkerStyle(numMarkers[i])
            numHistsRatio[i].SetMarkerSize(markerSize)
            numHistsRatio[i].SetMarkerColor(numColors[i])
            if(doNumFill):
                if('p' not in numDrawStyles[i]):
                    numHistsRatio[i].SetFillColor(numFillColors[i])
                    numHistsRatio[i].SetFillStyle(numFillStyles[i])
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
                numHistsRatio[i].Draw('ae0p same')
                # numHistsRatio[i].Draw('ap same')
            

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
    if(doRatio):
        canvas.cd()
        pad_2 = ROOT.TPad("pad_2","pad_2",0.00,0.33,0.99,0.99)
    else:
        pad_2 = ROOT.TPad("pad_2","pad_2",0.00,0.01,0.99,0.99)
    pad_2.Draw()
    pad_2.cd()
    pad_2.SetTopMargin(0.1)
    pad_2.SetBottomMargin(0.03)
    pad_2.SetRightMargin(0.1)
    pad_2.SetGridx(topGridx)
    pad_2.SetGridy(topGridy)

    denomHist.SetStats(statBox)
    denomHist.SetTitleSize(axisTitleSize)
    if(doRatio):
        denomHist.GetXaxis().SetTitleSize(0.)
    else:
        denomHist.GetXaxis().SetTitleSize((2*axisTitleSize))
        denomHist.GetXaxis().SetTitleOffset(0.043/(0.5*axisTitleSize))
    denomHist.GetYaxis().SetTitleSize(axisTitleSize)
    denomHist.GetYaxis().SetTitleOffset(0.033/(0.5*axisTitleSize))
    denomHist.SetTitle(Title+";;"+yTitle)
    denomHist.SetLineWidth(lineWidth)
    denomHist.SetLineColor(denomColor)
    denomHist.SetMarkerStyle(denomMarker)
    denomHist.GetYaxis().SetNdivisions(nDivY)
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
        if(doRatio):
            hs.GetHistogram().GetXaxis().SetTitleSize(0.)
            hs.GetHistogram().GetXaxis().SetLabelSize(0.)
        else:
            hs.GetXaxis().SetTitleSize((2*axisTitleSize))
            hs.GetXaxis().SetTitleOffset(0.043/(0.5*axisTitleSize))
        hs.GetHistogram().GetYaxis().SetTitleSize(axisTitleSize)
        hs.GetHistogram().GetYaxis().SetTitleOffset(0.033/(0.5*axisTitleSize))
        hs.GetHistogram().SetTitle(Title+";;"+yTitle)
        hs.GetHistogram().GetYaxis().SetLabelSize(0.8*axisTitleSize)

    if(errorBandHist == None):
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
        numHists[i].SetLineWidth(numLineWidth)
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
            # hs.SetMaximum(Max/(1+ROOT.gStyle.GetHistTopMargin()))  # wtf causes a segv

    
    denomHist.SetMinimum(Min)
    denomHist.SetMaximum(Max)
    errorBandHist.SetMinimum(Min)
    errorBandHist.SetMaximum(Max)

    if(setMax is not None):
        Max = setMax
        denomHist.SetMaximum(setMax)
        #errorBandHist.SetMaximum(Max)
        if(doStack==True): # goofy ROOT THStack max min issue
            if(not doLogy):
                hs.SetMaximum(Max/(1+ROOT.gStyle.GetHistTopMargin()))
            else:
                hs.SetMaximum(Max/(1+0.2*ROOT.TMath.Log10(Max/Min)))
    if(setMin is not None):
        Min = setMin
        denomHist.SetMinimum(setMin)
        if(doStack==True): # goofy ROOT THStack max min issue
            if(not doLogy):
                hs.SetMinimum(Min)
            else:
                hs.SetMinimum(Min*(1+0.5*ROOT.TMath.Log10(Max/Min)))

        #errorBandHist.SetMinimum(setMin)
    #####################################################################
    # build legend
    #####################################################################

    leg = pad_2.BuildLegend(legCoords[0],legCoords[1],legCoords[2],legCoords[3], "")  

    if(len(legList)>0):
        if order==None:
            order = range(len(legList))
        leg.Delete()
        leg = ROOT.TLegend(legCoords[0],legCoords[1],legCoords[2],legCoords[3], "","brNDC")
        leg.SetTextFont(42)  # wtf
        for i in order:
            if('skip' in legList[i] or 'Skip' in legList[i] or 'SKIP' in legList[i]):
                continue
            if(i>=len(numHists)):
                if(doStack==True):
                #     leg.AddEntry(denomHists[i-len(numHists)], legList[i],'f')  # wtf
                    inv = len(legList) - i
                    leg.AddEntry(denomHists[inv-len(numHists)], legList[inv],'f')  # wtf
                elif(drawErrorBand==True):
                    leg.AddEntry(errorBandHist, legList[i], 'fl')
                elif(denomDrawStyle=='HIST' or denomDrawStyle=='hist'):
                    leg.AddEntry(denomHist, legList[i], 'l')
                else:
                    leg.AddEntry(denomHist, legList[i], 'ep')
            else:
                if(numDrawStyles[i]=='HIST' or numDrawStyles[i]=='eHIST' or numDrawStyles[i]=='hist'):
                    leg.AddEntry(numHists[i], legList[i],'l')
                else:
                    leg.AddEntry(numHists[i], legList[i],'ep')
    leg.SetBorderSize(1)
    leg.SetFillStyle(1001) 
    leg.SetFillColor(0) 
    leg.SetLineWidth(0)
    leg.SetHeader(legHeader)
    # leg.SetTextSize(0.07)  # wtf

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
    if(doRatio):
        ROOT.SetOwnership(pad_1,0)
    ROOT.SetOwnership(pad_2,0)

    if(doRatio):
        for line in verticalLineB:
            ROOT.SetOwnership(line, 0)
    for line in verticalLineT:
        ROOT.SetOwnership(line, 0)
    for hist in numHistsRatio:
        ROOT.SetOwnership(hist, 0)

    # apply cms lumi info to pad

    if(doCMSlumi):
        cmsLumi(pad_2, iPeriod, iPos, extraText)

    return (canvas,pad_1,pad_2)

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

        hList.append(ROOT.TH1D(histName+str(i),histName+str(i),len(binSplit[i]),binSplit[i][0]-0.5,binSplit[i][len(binSplit[i])-1]+0.5))

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
        # eHistList[i].SetFillColor(ROOT.kRed)
        eHistList[i].SetFillColor(632)

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
        


def getDoubleRatioGraph(dist, histos=None, binning=None, extraCuts=None, nJetBin=None, bJetBin=None, kinBin=None, dphiCut=None, applyMassCut=None, applyPtCut=None, applyHTCut=None, applyMHTCut=None, applyNJetsCut=None, doLumi=None, removeZkfactor=None, doNoisy=None, doProof=None, treeLoc=None, applyEffs=None, treeName=None, applySF=None, doCG=None, do20=None, doClopperPearsonError=None, applyPuWeight=None, doDiMu=None, doDiEl=None):

    #####################################################################
    # set defaults and initialize 
    #####################################################################
    if(doClopperPearsonError==None):
        doClopperPearsonError=True
    if(do20==None):
        do20=False
    if(do20):
        if(doLumi==None):
            doLumi=35.9
    if (histos is None):
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
                binning=[1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]
                #binning=[1.5,2.5,4.5,6.5,9.5]
                #binning=[2.5,3.5,4.5,5.5,6.5,7.5]
        if(dist=='BTags'):
            bJetBin=-1
            if(binning==None):
                binning=[-0.5,0.5,1.5,2.5,3.5]
    else:
        binning = []
    if(doLumi==None):
        doLumi = 35.9
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
        treeLoc = "/nfs/data38/cms/wtford/lpcTrees/Skims/Run2ProductionV12"
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
    if(dphiCut=='none'):
        sampleSuffix='IDP'
    if(doDiMu==None):
        doDiMu=True
    if(doDiEl==None):
        doDiEl=True

    ########################################################
    # get bin center of gravity from photon distribution
    ########################################################
    cg = []
    if(doCG==True):
        if (histos is None):
            for i in range(len(binning)-1):
                h_cg = getDist('photon'+sampleSuffix, dist, distRange=[binning[i],binning[i+1]], nBins=100, doVarBinning=False,
                               extraCuts=extraCuts, nJetBin=nJetBin, bJetBin=bJetBin, dphiCut=dphiCut, applyMassCut=applyMassCut,
                               applyPtCut=applyPtCut, applyHTCut=applyHTCut, applyMHTCut=applyMHTCut, applyNJetsCut=applyNJetsCut, 
                               doLumi=doLumi, applyEffs=applyEffs, treeLoc=treeLoc, treeName=treeName)
                cg.append(h_cg.GetMean())
        else:
            for i in range(1, histos['pho_cg'].GetNbinsX()+1):
                cg.append(histos['pho_cg'].GetBinContent(i))
                binning.append(histos['pho_cg'].GetBinLowEdge(i))
            cgAxis = histos['pho_cg'].GetXaxis()
            binning.append(cgAxis.GetBinUpEdge(cgAxis.GetNbins()))
            

    ########################################################
    # end get bin center of gravity from photon distribution
    ########################################################


    ########################################################
    # get zll,dyll,photon, and gjets histograms
    ########################################################
    
    if (histos is None):
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

    else:
        pho_da = histos['pho_da']
        pho_mc = histos['pho_mc']
        zmm_da = histos['zmm_da']
        zmm_mc = histos['zmm_mc']
        zee_da = histos['zee_da']
        zee_mc = histos['zee_mc']

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
    print str(hist.GetName())
    for ipar in range(size):
        print "Function parameter "+str(ipar)+":  "+str(func.GetParameter(ipar))+" +/- "+str(func.GetParError(ipar))
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
    # mbins = 100
    mbins = nbins -1  # wtf
    bwidth = (hiEnd-loEnd)/float(mbins)
    for ibin in range(mbins+1):
        ptsForGrph.append(loEnd + ibin*bwidth)
    #############################

    #    for bin in range(nbins+2): # JFH 09-Oct-15    
    for bin in range(mbins+1): # JFH 09-Oct-15    
        hists[bin] = ROOT.TH1D("bin"+str(bin), "bin"+str(bin), 1000,-5.0,5.0)
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
    graph_mid.SetName('fitLine')
    graph_hi  = ROOT.TGraph(len(binCenterArray), binCenterArray, array_hi)
    graph_hi.SetName('fitLineUp')
    graph_lo  = ROOT.TGraph(len(binCenterArray), binCenterArray, array_lo)
    graph_lo.SetName('fitLineLow')
    if returnRMS:
        graph_rms = ROOT.TGraphErrors(len(binCenterArray), binCenterArray, array_y, array_ex, array_rms)

    if returnRMS:
        return graph_rms
    else:
        return (graph_mid, graph_lo, graph_hi)

def getDoubleRatioPlot(dr_graphs, binMeans=None, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, fitFunc=None, returnDiff=None, drawText=None, text=None, textCoords=None, addDeviationInQuad=None, Ymin=None, Ymax=None, plotRefLine=None):
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
    if(Ymin is None):
        Ymin = 0.
    if(Ymax is None):
        Ymax = 2.
    if(plotRefLine is None):
        plotRefLine = False
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
        # if(xlow==0.0):
        #     xlow = -0.5
        xhigh = graph.GetXaxis().GetBinUpEdge(nbins)

        #meanDR = graph.GetMean(2)

        func_mid  = ROOT.TF1("func_mid", fitFunc, xlow, xhigh)

        line.append(ROOT.TF1("ave", 'pol0(0)', xlow, xhigh))
        graph.Fit(line[it],"QN")

        #line.append(ROOT.TLine(xlow,meanDR,xhigh,meanDR))
        refLine = None

        fitGraph = getDoubleRatioFit(graph, func_mid, returnDiff)

        cIter = it
        while(type(ROOT.gROOT.FindObject("c"+str(cIter+1)))==ROOT.TCanvas):
            cIter+=1
        c.append(ROOT.TCanvas("c"+str(cIter+1),"c"+str(cIter+1),640*it+40,500,600,600))

        fitGraph[0].Draw('al')
        fitGraph[1].Draw('l same')
        fitGraph[2].Draw('l same')
        line[it].Draw('l same')
        if (plotRefLine and 'HT' in graph.GetTitle() and not 'MHT' in graph.GetTitle()):
            refLine = ROOT.TF1("refLine", 'pol1(0)', xlow, xhigh)
            refLine.SetParameters(0.8566, 0.0001950)
            refLine.Draw('l same')
            refLine.SetLineWidth(2)
            refLine.SetLineStyle(4)
            refLine.SetLineColor(8)
        graph.Draw('e0p same')

        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1)  # wtf
        graph.SetLineColor(1)

        fitGraph[0].SetMaximum(Ymax)
        fitGraph[0].SetMinimum(Ymin)
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

            if (binMeans is None):
                bin1 = 2.
                h_bin2 = getDist("photon","NJets",distRange=[2.5,4.5],doVarBinning=False)
                bin2 = h_bin2.GetMean()
                h_bin3 = getDist("photon","NJets",distRange=[4.5,6.5],doVarBinning=False)
                bin3 = h_bin3.GetMean()
                h_bin4 = getDist("photon","NJets",distRange=[6.5,8.5],doVarBinning=False)
                bin4 = h_bin4.GetMean()
                h_bin5 = getDist("photon","NJets",distRange=[8.5,12.5],doVarBinning=False)
                bin5 = h_bin5.GetMean()
            else:
                bin1 = binMeans['NJets'][0]
                bin2 = binMeans['NJets'][1]
                bin3 = binMeans['NJets'][2]
                bin4 = binMeans['NJets'][3]
                bin5 = binMeans['NJets'][4]

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
            
            if (binMeans is None):
                h_bin1 = getDist("photon","MHT",distRange=[300,350],doVarBinning=False)
                bin1 = h_bin1.GetMean()
                h_bin2 = getDist("photon","MHT",distRange=[350,500],doVarBinning=False)
                bin2 = h_bin2.GetMean()
                h_bin3 = getDist("photon","MHT",distRange=[500,750],doVarBinning=False)
                bin3 = h_bin3.GetMean()
                h_bin4 = getDist("photon","MHT",distRange=[750,1500],doVarBinning=False)
                bin4 = h_bin4.GetMean()
                h_bin5 = getDist("photon","MHT",distRange=[250,300],doVarBinning=False,applyMHTCut=False)
                bin5 = h_bin5.GetMean()
            else:
                bin1 = binMeans['MHT'][1]
                bin2 = binMeans['MHT'][2]
                bin3 = binMeans['MHT'][3]
                bin4 = binMeans['MHT'][4]
                bin5 = binMeans['MHT'][0]

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

            if (binMeans is None):
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
            else:
                bin1 = binMeans['HT'][0]
                bin2 = binMeans['HT'][1]
                bin3 = binMeans['HT'][2]
                bin4 = binMeans['HT'][3+1]
                bin5 = binMeans['HT'][3+3]
                bin6 = binMeans['HT'][3+4]
            
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

    if (type(refLine) == ROOT.TF1):
        ROOT.SetOwnership(refLine, 0)

    return [(drVal,drErr),errorDict]


def getZmassFitPlot(fitFunc=None, dataSet=None, mcSet=None, plotMC=None, doDiMu=None, doDiEl=None, doDiLep=None, getShapeFromLoose=None, nBins=None, distRange=None, nJetBin=None, bJetBin=None, kinBin=None, doVarBinning=None, binning=None, extraCuts=None, dphiCut=None, doLumi=None, do20=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, keepCanvas=None, drawText=None, text=None, textCoords=None, drawText2=None, text2=None, textCoords2=None, savePlot=None):
    """getZmassFitPlot() returns diMu and diEl mass plots (canvas's) incuding 
    fits to dataset bin with bJetBin, nJetBin, or kinBin, default is -1 (inclusive)
    set doDiLep=True if you want to merge mumu ee mass distributions
    predefined Z fit functions: (fitFunc parameter) 
    'CB'-crystal ball
    'V'-Voigtian i.e. breit-wigner convolved with Gaussian
    'DG'-Double Gaussian
    'TG'-Triple Gaussian"""

    setTdrStyle()
    
    if(savePlot==None):
        savePlot=False
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
            doLumi=35.9
        if(extraText==None):
            extraText="Simulation"
    if(text==None):
        text = "arXiv:1704.07781"
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
        iPeriod=5
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
        doLumi=35.9
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
    else:
        # nEvents=[]
        for dataSetDict in dataSet:
            dataSetDict["sig"].SetMarkerStyle(20)
            dataSetDict["sig"].SetMarkerSize(1.1)
            hS = dataSetDict["sig"]
            # Here we kluge the Z-mass-cut histogram; assumes current binning.
            # hE = hS.Clone()
            # for bin in range(9):  # underflow to 76.0
            #     hE.SetBinContent(bin, 0)
            # for bin in range(24, hE.GetXaxis().GetNbins()+2):  # 106.0 to overflow
            #     hE.SetBinContent(bin, 0)
            # nEvents.append(hE)
        if(mcSet==None):
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

        else:
            # The MC stacked histograms are provided
            # Assume stack order is tt, diboson, ttz, dy
            for stackhist in mcSet:
                if ((doDiMu and ("hsMM" in stackhist.GetName())) or 
                    (doDiEl and ("hsEE" in stackhist.GetName())) or 
                    (doDiLep and ("hsLL" in stackhist.GetName()))):
                    hlist = stackhist.GetHists()
                    if (len(hlist) >= 4):
                        hlist[3].SetFillColor(2)
                        mcList.append(hlist[3])
                    if (len(hlist) >= 2):
                        hlist[1].SetFillColor(5)
                        mcList.append(hlist[1])
                    if (len(hlist) >= 3):
                        hlist[2].SetFillColor(8)
                        mcList.append(hlist[2])
                    if (len(hlist) >= 1):
                        hlist[0].SetFillColor(38)
                        mcList.append(hlist[0])

    canvas = []
    purityList = []

    jDict = {0: 'dimu', 1: 'diel'}

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
        vg = ROOT.RooRealVar("vg","g in voigtian",2.64,1.0,5.0)  # wtf
        vs = ROOT.RooRealVar("vs","s in voigtian",1.45,1.0,3.0)  # wtf
        # vg = ROOT.RooRealVar("vg","g in voigtian",2.64,1.0,3.0)
        # vs = ROOT.RooRealVar("vs","s in voigtian",1.45,1.0,2.0)
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
            # 'V': 0.9780,
            'V': 0.9892,  # 2017v16
            # 'V': 1.0,  # 2017v16
            'DG': 0.9478,
            'TG': 0.9505,
        }
        # correctionFactor = {
        #     'CB': 1.,
        #     'V': 1.,
        #     'DG': 1.,
        #     'TG': 1.,
        # }

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
            # fsig = ROOT.RooRealVar("fsig","fsig",0.7,0.2,1.0)  # wtf

            model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(fitFuncDict[fitFunc],bfit),ROOT.RooArgList(nsig,nbkg))
            model2 = ROOT.RooAddPdf("model2","model2",ROOT.RooArgList(fitFuncDict[fitFunc],bfit),ROOT.RooArgList(frac))

            if(it==0):
                model.fitTo(Data, ROOT.RooFit.SumW2Error(0))  # for fitting weighted MC; set 1 for wt-corrected errors
            else:
                # model.fitTo(Data)
                model2.fitTo(Data, ROOT.RooFit.SumW2Error(0))
                
                rooVarX.setRange("all",60,120)
                rooVarX.setRange("selection",76.19,106.19)

                selFraction = frac.getVal()
                selFractionE = frac.getError()
                allSel = fitFuncDict[fitFunc].createIntegral(ROOT.RooArgSet(rooArgL),"all").getVal()
                selectionSel = fitFuncDict[fitFunc].createIntegral(ROOT.RooArgSet(rooArgL),"selection").getVal()
                allBGSel = bfit.createIntegral(ROOT.RooArgSet(rooArgL),"all").getVal()
                selectionBGSel = bfit.createIntegral(ROOT.RooArgSet(rooArgL),"selection").getVal()
                # print "PDF integral:  sig, bg = "+str(selectionSel/allSel)+", "+str(selectionBGSel/allBGSel)

                # nsigLVLB = nsig.getVal()
                # nbkgLVLB = nbkg.getVal()

                # integralSigSel = ((selectionSel/allSel)*nsigLVLB)
                # integralBgSel = ((selectionBGSel/allBGSel)*nbkgLVLB)

                # # print "Nbins, (Entries, Sum) in sig, E = "+str(nEvents[j].GetXaxis().GetNbins())+", ("+str(data[key].GetEntries())+", "+str(data[key].GetSumOfWeights())+"), ("+str(nEvents[j].GetEntries())+", "+str(nEvents[j].GetSumOfWeights())+")"
                # # selentries = data[key].GetEntries()
                # # selSelentries = nEvents[j].GetEntries()
                # selentries = data[key].GetSumOfWeights()  # wtf
                # selSelentries = nEvents[j].GetSumOfWeights()  # wtf
                # selFractionENaive = ROOT.TMath.Sqrt( (selFraction*(1-selFraction))/selentries )
                # selFractionSel = integralSigSel/(integralSigSel+integralBgSel)
                # selFractionESelNaive = ROOT.TMath.Sqrt( (selFractionSel*(1-selFractionSel))/(float(selSelentries) ))
                # selFractionESel = selFractionESelNaive*(selFractionE/selFractionENaive)

                # purity_orig = (1/correctionFactor[fitFunc])*integralSigSel/(integralSigSel+integralBgSel)
                # Epurity_orig = (1/correctionFactor[fitFunc])*selFractionESel

                # Alternate purity calculation (wtf)
                Rbs = (selectionBGSel/allBGSel) / (selectionSel/allSel)
                purity = (1/correctionFactor[fitFunc]) * selFraction / (Rbs + (1-Rbs)*selFraction)
                Epurity = (1/correctionFactor[fitFunc]) * Rbs*selFractionE / (Rbs + (1-Rbs)*selFraction)**2
                # print "purity orig, new = "+str(purity_orig)+"+/-"+str(Epurity_orig)+", "+str(purity)+"+/-"+str(Epurity)
                # Take 1-correctionFactor as a systematic:  # wtf
                Epurity = math.sqrt(Epurity**2 + (1-correctionFactor[fitFunc])**2)
                
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
            canvas[it+2*j].SetLeftMargin(0.18)

            if(plotMC):
                if(key == keyList[len(keyList)-1]):
                    mcSet[j].Draw("HIST")
                    mcSet[j].GetXaxis().SetTitle(xTitle)
                    mcSet[j].GetYaxis().SetTitle(yTitle)
                    mcSet[j].GetYaxis().SetTitleOffset(1.5)
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

                if(plotMC):
                    leg1 = ROOT.TLegend(0.65,0.60,0.85,0.88, "","brNDC")
                    leg1.SetBorderSize(2)
                    leg1.SetFillStyle(1001)
                    # leg1.SetFillColor(ROOT.kWhite) 
                    leg1.SetFillColor(0) 
                    if (len(mcList) / len(mcSet) >= 1):
                        leg1.AddEntry(mcList[0],"Drell-Yan Sim","f")
                    if (len(mcList) / len(mcSet) >= 2):
                        leg1.AddEntry(mcList[1],"t#bar{t}Z Sim","f")
                    if (len(mcList) / len(mcSet) >= 3):
                        leg1.AddEntry(mcList[2],"VV Sim","f")
                    if (len(mcList) / len(mcSet) >= 4):
                        leg1.AddEntry(mcList[3],"t#bar{t} Sim","f")
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
 
            if(savePlot and it==1):
                canvas[it+2*j].SaveAs('../plots/massPlots/'+jDict[j]+'_nj'+str(nJetBin)+'_nb'+str(bJetBin)+'.pdf')

            it+=1

    for c in canvas:
        ROOT.SetOwnership(c, int(not keepCanvas))

    return purityList

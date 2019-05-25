#!/usr/bin/python
"""
Make two-panel figure with histograms
"""
from array import array
import ROOT
import RA2b

def textCMS(margin):
  CMStxt = ROOT.TPaveText(margin, 0.98, 0.23, 1.0, "NDC")
  CMStxt.SetBorderSize(0)
  CMStxt.SetFillColor(0)
  CMStxt.SetTextFont(42)
  CMStxt.SetTextAlign(13)
  CMStxt.SetTextSize(0.11)
  CMStxt.SetMargin(0.)
  CMStxt.AddText("#bf{CMS}")
  return CMStxt

def textLumi(margin):
  Lumitxt = ROOT.TPaveText(0.80, 0.98, 1-margin, 1.0, "NDC")
  Lumitxt.SetBorderSize(0)
  Lumitxt.SetFillColor(0)
  Lumitxt.SetTextFont(42)
  Lumitxt.SetTextAlign(33)
  Lumitxt.SetTextSize(0.09)
  Lumitxt.SetMargin(0.)
  Lumitxt.AddText("137 fb^{-1} (13 TeV)")
  return Lumitxt

def addDivisionsUp(ymin, ymax, axisTitleSize):
  # Putting lines and labels explaining search region definitions
  ymax2 = ymax/10
  ls = 0.5
  # Njet labels
  ttext_njet = ROOT.TLatex();
  ttext_njet.SetTextFont(42);
  ttext_njet.SetTextSize(1.3*axisTitleSize);
  ttext_njet.SetTextAlign(22);
  ttext_njet.DrawLatex(2.1 , ymax/1.5 , "2 #leq N_{#scale[0.2]{ }jet} #leq 3");
  ttext_njet.DrawLatex(5.5 , ymax/1.5 , "4 #leq N_{#scale[0.2]{ }jet} #leq 5");
  ttext_njet.DrawLatex(9.5 , ymax/1.5 , "6 #leq N_{#scale[0.2]{ }jet} #leq 7");
  ttext_njet.DrawLatex(13.5 , ymax/1.5 , "8 #leq N_{#scale[0.2]{ }jet} #leq 9");
  ttext_njet.DrawLatex(17.5 , ymax/1.5 , "N_{#scale[0.2]{ }jet} #geq 10");

  # Njet separation lines
  tl_njet = ROOT.TLine();
  tl_njet.SetLineStyle(2);
  tl_njet.DrawLine(3.+ls,ymin/10,3.+ls,1.2*ymax);
  tl_njet.DrawLine(7.+ls,ymin/10,7.+ls,1.2*ymax); 
  tl_njet.DrawLine(11.+ls,ymin/10,11.+ls,1.2*ymax);
  tl_njet.DrawLine(15.+ls,ymin/10,15.+ls,1.2*ymax);

  # Nb labels
  ttext_nb = ROOT.TLatex();
  ttext_nb.SetTextFont(42);
  ttext_nb.SetTextSize(1.3*axisTitleSize);
  ttext_nb.SetTextAlign(22);
    
  ttext_nb.DrawLatex(2., ymax/5.5 , "N_{#scale[0.2]{ }b-jet}");
  ttext_nb.DrawLatex(1., ymax/20. , "0");
  ttext_nb.DrawLatex(2., ymax/20. , "1");
  ttext_nb.DrawLatex(3., ymax/20. , "#geq 2");
  ttext_nb.DrawLatex(4., ymax/20. , "0");
  ttext_nb.DrawLatex(5., ymax/20. , "1");
  ttext_nb.DrawLatex(6., ymax/20. , "2");
  ttext_nb.DrawLatex(7., ymax/20. , "#geq 3");

  # Nb separation lines
  tl_nb = ROOT.TLine();
  tl_nb.SetLineStyle(3);
  tl_nb.DrawLine(1.+ls,ymin/10,1.+ls,ymax2); 
  tl_nb.DrawLine(2.+ls,ymin/10,2.+ls,ymax2); 

  tl_nb.DrawLine(4.+ls,ymin/10,4.+ls,ymax2); 
  tl_nb.DrawLine(5.+ls,ymin/10,5.+ls,ymax2); 
  tl_nb.DrawLine(6.+ls,ymin/10,6.+ls,ymax2); 

  tl_nb.DrawLine(8.+ls,ymin/10,8.+ls,ymax2); 
  tl_nb.DrawLine(9.+ls,ymin/10,9.+ls,ymax2); 
  tl_nb.DrawLine(10.+ls,ymin/10,10.+ls,ymax2); 

  tl_nb.DrawLine(11.+ls,ymin/10,11.+ls,ymax2); 
  tl_nb.DrawLine(12.+ls,ymin/10,12.+ls,ymax2); 
  tl_nb.DrawLine(13.+ls,ymin/10,13.+ls,ymax2); 

  tl_nb.DrawLine(14.+ls,ymin/10,14.+ls,ymax2); 
  tl_nb.DrawLine(15.+ls,ymin/10,15.+ls,ymax2); 
  tl_nb.DrawLine(16.+ls,ymin/10,16.+ls,ymax2); 

  tl_nb.DrawLine(17.+ls,ymin/10,17.+ls,ymax2); 
  tl_nb.DrawLine(18.+ls,ymin/10,18.+ls,ymax2); 

def addDivisionsLow(ymin, ymax, axisTitleSize):
  # Putting lines and labels explaining search region definitions
  ls = 0.5

  # Njet separation lines
  tl_njet = ROOT.TLine();
  tl_njet.SetLineStyle(2);
  tl_njet.DrawLine(3.+ls,ymin,3.+ls,10*ymax);
  tl_njet.DrawLine(7.+ls,ymin,7.+ls,10*ymax); 
  tl_njet.DrawLine(11.+ls,ymin,11.+ls,10*ymax);
  tl_njet.DrawLine(15.+ls,ymin,15.+ls,10*ymax);

  # Nb separation lines
  tl_nb = ROOT.TLine();
  tl_nb.SetLineStyle(3);
  tl_nb.DrawLine(1.+ls,ymin,1.+ls,ymax); 
  tl_nb.DrawLine(2.+ls,ymin,2.+ls,ymax); 

  tl_nb.DrawLine(4.+ls,ymin,4.+ls,ymax); 
  tl_nb.DrawLine(5.+ls,ymin,5.+ls,ymax); 
  tl_nb.DrawLine(6.+ls,ymin,6.+ls,ymax); 

  tl_nb.DrawLine(8.+ls,ymin,8.+ls,ymax); 
  tl_nb.DrawLine(9.+ls,ymin,9.+ls,ymax); 
  tl_nb.DrawLine(10.+ls,ymin,10.+ls,ymax); 

  tl_nb.DrawLine(11.+ls,ymin,11.+ls,ymax); 
  tl_nb.DrawLine(12.+ls,ymin,12.+ls,ymax); 
  tl_nb.DrawLine(13.+ls,ymin,13.+ls,ymax); 

  tl_nb.DrawLine(14.+ls,ymin,14.+ls,ymax); 
  tl_nb.DrawLine(15.+ls,ymin,15.+ls,ymax); 
  tl_nb.DrawLine(16.+ls,ymin,16.+ls,ymax); 

  tl_nb.DrawLine(17.+ls,ymin,17.+ls,ymax); 
  tl_nb.DrawLine(18.+ls,ymin,18.+ls,ymax); 

def panelHists():
  ROOT.gROOT.Reset()
  ROOT.gROOT.SetBatch(1)
  
  doMumu = False
  doEe = False
  doLl = True
  doPhoton = False
  nScale = [1, 1, 1, 1]
  dScale = [1, 1, 1, 1]
  MZmmMax = 0
  MZeeMax = 0
  setMin = None
  setMax = None
  ratioMin = [None, None, None, None]
  ratioMax = [None, None, None, None]
  extraText = ""
  
  legHeader = None
  legList = []
  # legCoordsDefault = [0.6,0.6,0.89,0.89]
  legCoordsDefault = [0.68,0.6,0.89,0.89]
  legCoords = None
  drawText = False
  textCoordsDefault = [0.35,0.75,.55,.88]
  textCoords = None
  
  def histNamesBuild(namesRootList, suffix):
    histNames = []
    for nameRoot in namesRootList:
      histNames.append(nameRoot+"_"+suffix)
    return histNames
  
  hists = []
  filehistsM = {}
  filehistsE = {}
  filehistsL = {}
  filehistsP = {}
  
  #  ========================================================================================
  
  varNames = [['hCCjb']]
  sampleSuffixN = ['zll']
  sampleSuffixD = [['ttll', 'VVll', 'ttzll', 'dyll']]
  legendsD = [['  t #bar{t}', '  Diboson', '  t #bar{t} Z', '  Z + jets']]
  iPeriod = 9
  NfileZll = ROOT.TFile('../outputs/histsDY_Run2v17.root')
  # NfilePhoton = ROOT.TFile('../outputs/histsPhoton_Run2v17.root')
  DfileZll = ROOT.TFile('../outputs/histsDYMC_Run2v17.root')
  # DfilePhoton = ROOT.TFile('../outputs/histsGjets_Run2v17_DRr2wt.root')
  extrapFile = ROOT.TFile('../outputs/hExtrap_Run2v17.root')
  legHeader = 'Z #rightarrow l #bar{l}'
  legendsN = ['  Data']
  
  
  
  # =============================================================================
  reaction = -1
  if (doMumu):
    reaction += 1
    legend = []
    histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
    legend.append(legendsN[reaction])
    histnamesD = []
    for hnd in sampleSuffixD[reaction]:
      histnamesD.append(histNamesBuild(varNames[reaction], hnd))
    for ld in legendsD[reaction]:
      legend.append(ld)
    filehistsM['N'] = (NfileZll, nScale[reaction], histnamesN)
    filehistsM['D'] = (DfileZll, dScale[reaction], histnamesD)
    hists.append(filehistsM)
    legList.append(legend)
  if (doEe):
    reaction += 1
    legend = []
    histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
    legend.append(legendsN[reaction])
    histnamesD = []
    for hnd in sampleSuffixD[reaction]:
      histnamesD.append(histNamesBuild(varNames[reaction], hnd))
    for ld in legendsD[reaction]:
      legend.append(ld)
    filehistsE['N'] = (NfileZll, nScale[reaction], histnamesN)
    filehistsE['D'] = (DfileZll, dScale[reaction], histnamesD)
    hists.append(filehistsE)
    legList.append(legend)
  if (doLl):
    reaction += 1
    legend = []
    histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
    legend.append(legendsN[reaction])
    histnamesD = []
    for hnd in sampleSuffixD[reaction]:
      histnamesD.append(histNamesBuild(varNames[reaction], hnd))
    for ld in legendsD[reaction]:
      legend.append(ld)
    filehistsL['N'] = (NfileZll, nScale[reaction], histnamesN)
    filehistsL['D'] = (DfileZll, dScale[reaction], histnamesD)
    hists.append(filehistsL)
    legList.append(legend)
  if (doPhoton):
    reaction += 1
    legend = []
    histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
    legend.append(legendsN[reaction])
    histnamesD = []
    for hnd in sampleSuffixD[reaction]:
      histnamesD.append(histNamesBuild(varNames[reaction], hnd))
    for ld in legendsD[reaction]:
      legend.append(ld)
    filehistsP['N'] = (NfilePhoton, nScale[reaction], histnamesN)
    filehistsP['D'] = (DfilePhoton, dScale[reaction], histnamesD)
    hists.append(filehistsP)
    legList.append(legend)
  
  #  ========================================================================================
  
  reaction = -1
  for samples in hists:
    reaction += 1
    Nfile = samples['N'][0]
    Dfile = samples['D'][0]
    for i in range(len(samples['N'][2])):
      nhName = samples['N'][2][i]
      print "nhName = "+str(nhName)+" in file "+str(Nfile)
      hNumer = Nfile.Get(nhName)
      hNumer.SetName(str(nhName)+"N")
      hNumer.Scale(samples['N'][1])
      hDenList = []
      for c in range(len(samples['D'][2])):  # Denominator may have multiple components
        dhName = samples['D'][2][c][i]  # scale is stored in samples['D'][1]
        print "dhName = "+str(dhName)+" in file "+str(Dfile)
        hDen = Dfile.Get(dhName)
        hDen.SetName(str(dhName)+"D")
        hDen.Scale(samples['D'][1])
        hDenList.append(hDen)
      if ("hZmass" in nhName):
        doLogy = False
        if ("mm" in nhName and MZmmMax != 0):
          hNumer.SetMaximum(MZmmMax)
        elif ("ee" in nhName and MZeeMax != 0):
          hNumer.SetMaximum(MZeeMax)
      else:
        doLogy = True
      if ("Eta" in nhName and not "Photon" in nhName):
        setMin = 1
        # legCoords = legCoordsDefault
        # textCoords = textCoordsDefault
        # legCoords = [0.6,0.6,0.89,0.89]
        # drawText = True
        # textCoords = [0.35,0.75,.55,.88]
      else:
        setMin = None
        # drawText = True
        # textCoords = textCoordsDefault
        # legCoords = legCoordsDefault
      # print "legCoords = "+str(legCoords[0])+", "+str(legCoords[1])+", "+str(legCoords[2])+", "+str(legCoords[3])
      doClosureStyle = False
      # if ("hCC" in nhName):
      #   doClosureStyle = True
      canvTuple = RA2b.getPlotAndRatio(
        numHists=hNumer, denomHists=hDenList, doRatio=True,
        doLogy=doLogy, doCMSlumi=True, iPeriod=iPeriod, drawHorizontalLine=False,
        xTitle=hNumer.GetXaxis().GetTitle(), yTitle=hNumer.GetYaxis().GetTitle(),
        ratioMin=ratioMin[reaction], ratioMax=ratioMax[reaction], setMin=setMin, setMax=setMax,
        legHeader = legHeader, legList = legList[reaction],
        # legCoords = legCoords,
        drawText = drawText, textCoords = textCoords, extraText = extraText,
        doClosureStyle = doClosureStyle, markerSize=0.8, canvasSize = [900,600]
        )
  
      hExtrap = extrapFile.Get("hExtrap")
      # Extract objects from canvTuple
      padHi = canvTuple[2]
      for obj in padHi.GetListOfPrimitives():
        objName = obj.GetName()
        print "Object in pad:  "+objName
        # Object in pad:  TFrame
        # Object in pad:  hs  THStack
        # Object in pad:  hCCjb_ttllD  TH1
        # Object in pad:  hCCjb_zllN  TH1
        # Object in pad:  TPave  TPave, Tlegend
        # Object in pad:
        # Object in pad:
        # Object in pad:
        if (obj.InheritsFrom("TH1")):
          if (objName == "hCCjb_zllN"):
            hData = obj.Clone("hData")
          # elif (objName == "hCCjb_ttllD"):
          #   denHist = obj.Clone("denHist")
        if (obj.InheritsFrom("THStack")):
          hMC = obj.Clone("hMC")
        if (obj.InheritsFrom("TLegend")):
          legend = obj.Clone("legend")
        if (obj.InheritsFrom("TFrame")):
          frame2 = obj.Clone("frame2")
          frame2.SetX1(0)
          frame2.SetX2(20)
      
      # Set up canvas dimensions
      W = 900
      H = 720
      newCanv = ROOT.TCanvas("newCanv", "newCanv", 0,0, W, H)
      newCanv.Range(0,0,1,1)
      newCanv.cd()

      pml = 0.13  # pad horizontal margins
      pmr = 0.01

      # Lower panel
      p1b = 0.01  # pad1 vertical dimensions
      p1t = 0.443
      p1mt = 0.005  # pad1 vertical margins
      p1mb = 0.25
      p1min = 0.001
      p1max = 2
      axisTitleSize = 0.055
      markerSize = 1.4
      pad1 = ROOT.TPad("pad1", "pad1", 0.00, p1b, 0.99, p1t)
      pad1.Draw()
      pad1.cd()
      pad1.SetTopMargin(p1mt)
      pad1.SetBottomMargin(p1mb)
      pad1.SetLeftMargin(pml)
      pad1.SetRightMargin(pmr)
      pad1.SetLogy(1)
      hExtrap.Draw("e0p X0")
      hExtrap.GetXaxis().SetTitleSize(2*axisTitleSize)
      hExtrap.GetXaxis().SetLabelSize(1.5*axisTitleSize)
      hExtrap.GetYaxis().SetLabelSize(1.5*axisTitleSize)
      hExtrap.GetXaxis().SetTitle('N_{#scale[0.5]{ }jet}, N_{#scale[0.5]{ }b-jet} bin index')

      hExtrap.GetYaxis().SetTitle('F^{#scale[0.5]{ }data}_{j,b}')
      hExtrap.GetYaxis().SetTitleSize(2*axisTitleSize)
      hExtrap.GetYaxis().SetTitleOffset(0.025/(1.0*axisTitleSize))

      # Attempt to get \mathCal{F}
      # hExtrap.GetYaxis().SetTitle('')
      # Ytitle = ROOT.TMathText()
      # Ytitle.SetTextSize(2*axisTitleSize)
      # Ytitle.SetTextAlign(32)
      # # Ytitle.SetTextAngle(90)
      # Ytitle.DrawMathText(0.5, 0.6, "F")
      # # Ytitle.DrawMathText(0.5, 0.6, "\\mathscr{F}")

      hExtrap.SetMinimum(p1min)
      hExtrap.SetMaximum(p1max)
      hExtrap.SetLineWidth(1)
      hExtrap.SetLineColor(1)
      hExtrap.SetMarkerStyle(20)
      hExtrap.SetMarkerSize(markerSize)
      hExtrap.SetMarkerColor(1)
      addDivisionsLow(p1min, p1max, axisTitleSize)
      pad1.Draw()
      pad1.Update()
      newCanv.Update()
      newCanv.cd()

      # Upper panel
      p2b = 0.44  # pad2 vertical dimensions
      p2t = 0.99
      p2mt = 0.11  # pad2 vertical margins
      p2mb = 0.005
      p2min = .5
      p2max = 5.0e6
      titleScale = (p1t-p1b)/(p2t-p2b)
      print "title scale factor = "+str(titleScale)
      axisTitleSize *= titleScale
      pad2 = ROOT.TPad("pad2", "pad2", 0.00, p2b, 0.99, p2t)
      pad2.Draw()
      pad2.cd()
      pad2.SetTopMargin(p2mt)
      pad2.SetBottomMargin(p2mb)
      pad2.SetLeftMargin(pml)
      pad2.SetRightMargin(pmr)
      pad2.SetLogy(1)
      pad2.SetTickx(0)
      hMC.SetMinimum(p2min)
      hMC.SetMaximum(p2max)
      hMC.GetYaxis().SetTitleSize((2*axisTitleSize))
      hMC.GetYaxis().SetLabelSize(1.5*axisTitleSize)
      hMC.GetYaxis().SetTitleOffset(0.03/(1.0*axisTitleSize))
      hMC.GetYaxis().SetTitle('Events')
      hMC.GetXaxis().SetTicks('')
      print "Ordinate title offset = "+str(hMC.GetYaxis().GetTitleOffset())
      hMC.Draw()
      hData.SetMarkerSize(markerSize)
      hData.Draw("e0p X0 same")
      CMStxt = textCMS(pml)
      CMStxt.Draw()
      Lumitxt = textLumi(pmr)
      Lumitxt.Draw()
      addDivisionsUp(p2min, p2max, axisTitleSize)
      legend.SetX1NDC(0.78)
      legend.SetY1NDC(0.32)
      legend.SetX2NDC(0.93)
      legend.SetY2NDC(0.75)
      legend.SetTextSize(1.3*axisTitleSize)
      legend.SetLineWidth(1)
      legend.SetFillStyle(1001)
      legend.SetFillColor(0)
      legend.SetLineColor(1)
      legend.SetBorderSize(1)
      # legend.SetMargin(0.3)  # this 0.3 seems to be the default
      legend.Draw()

      newCanv.cd()

      newCanv.SaveAs(str(nhName)+"-Fextrap.pdf")

      exit
    
  # =================================================================================
  
  
  # def getPlotAndRatio(numHists, denomHists=None, bottomPlots=None, doStack=None, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, ratioTitle=None, ratioMin=None, ratioMax=None, doLogy=None, doFlip=None, doDiff=None, doPull=None, makeLeg=None, legList=None, legCoords=None, textCoords=None, canvasSize=None, canvasName=None, numColors=None, denomColor=None, numMarkers=None, denomMarker=None, markerSize=None, lineWidth=None, numDrawStyles=None, denomDrawStyle=None, drawErrorBand=None, stackColors=None, axisTitleSize=None, drawVerticalLines=None, drawHorizontalLine=None, statBox=None, drawText=None, text=None, setMax=None, setMin=None, doClosureStyle=None,errorBandColor=None,errorBandFillStyle=None,legHeader=None,nDivRatio=None,doNumFill=None, hLineVal=None, hLineColors=None,nDivX=None,ratioGridx=None,ratioGridy=None,topGridx=None,topGridy=None,doRatio=None,numFillStyles=None,numFillColors=None)


def main():
  ROOT.gROOT.Reset()
  ROOT.gROOT.SetBatch(1)
  panelHists()
  
if __name__ == "__main__":
  main()

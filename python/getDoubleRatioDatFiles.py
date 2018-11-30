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
    doSample = ['sig','hdp','ldp']
#################################################

haveHistograms = True
runBlock = "2016"

DYkfactor = 1.23

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
        if (runBlock is "2016"):
            dataZllFile = ROOT.TFile('../outputs/histsDY_2016v15.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_2016v15.root')
            # dataZllFile = ROOT.TFile('../outputs/histsDY_2016v12.root')
            # dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_2016v12.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMC_2016v12.root')
            mcPhotonFile = ROOT.TFile('../outputs/histsGjets_2016v12.root')
        histoNJets['pho_da'] = dataPhotonFile.Get("hNJets_DR_photon")
        histoNJets['pho_cg'] = histoNJets['pho_da'].Clone()
        for i in range(1, histoNJets['pho_cg'].GetNbinsX()+1):
            histoNJets['pho_cg'].SetBinContent(i, histoNJets['pho_cg'].GetBinCenter(i))
        histoNJets['pho_mc'] = mcPhotonFile.Get("hNJets_DR_gjets")
        histoNJets['zmm_da'] = dataZllFile.Get("hNJets_DR_zmm")
        histoNJets['zmm_mc'] = mcZllFile.Get("hNJets_DR_dymm")
        histoNJets['zmm_mc'].Scale(1/DYkfactor)
        histoNJets['zee_da'] = dataZllFile.Get("hNJets_DR_zee")
        histoNJets['zee_mc'] = mcZllFile.Get("hNJets_DR_dyee")
        histoNJets['zee_mc'].Scale(1/DYkfactor)
        histoMHT['pho_da'] = dataPhotonFile.Get("hMHT_DR_photon")
        histoMHT['pho_cg'] = dataPhotonFile.Get("hMHT_DR_xWt_photon")
        histoMHT['pho_cg'].Divide(histoMHT['pho_da'])
        histoMHT['pho_mc'] = mcPhotonFile.Get("hMHT_DR_gjets")
        histoMHT['zmm_da'] = dataZllFile.Get("hMHT_DR_zmm")
        histoMHT['zmm_mc'] = mcZllFile.Get("hMHT_DR_dymm")
        histoMHT['zmm_mc'].Scale(1/DYkfactor)
        histoMHT['zee_da'] = dataZllFile.Get("hMHT_DR_zee")
        histoMHT['zee_mc'] = mcZllFile.Get("hMHT_DR_dyee")
        histoMHT['zee_mc'].Scale(1/DYkfactor)
        histoHT['pho_da'] = dataPhotonFile.Get("hHT_DR_photon")
        histoHT['pho_cg'] = dataPhotonFile.Get("hHT_DR_xWt_photon")
        histoHT['pho_cg'].Divide(histoHT['pho_da'])
        histoHT['pho_mc'] = mcPhotonFile.Get("hHT_DR_gjets")
        histoHT['zmm_da'] = dataZllFile.Get("hHT_DR_zmm")
        histoHT['zmm_mc'] = mcZllFile.Get("hHT_DR_dymm")
        histoHT['zmm_mc'].Scale(1/DYkfactor)
        histoHT['zee_da'] = dataZllFile.Get("hHT_DR_zee")
        histoHT['zee_mc'] = mcZllFile.Get("hHT_DR_dyee")
        histoHT['zee_mc'].Scale(1/DYkfactor)

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
        dr_out = RA2b.getDoubleRatioPlot([nj_dr,mht_dr,ht_dr], binMeans = binMeans)

    # raw_input("Press the <ENTER> key to continue...")
    if (sample == "sig"):
        c1 = ROOT.gROOT.FindObject("c1")
        if (type(c1)==ROOT.TCanvas):
            fn = "DR_NJets_"+sample+".pdf"
            c1.SaveAs(fn)
            c2 = ROOT.gROOT.FindObject("c2")
        if (type(c2)==ROOT.TCanvas):
            fn = "DR_MHT_"+sample+".pdf"
            c2.SaveAs(fn)
            c3 = ROOT.gROOT.FindObject("c3")
        if (type(c3)==ROOT.TCanvas):
            fn = "DR_HT_"+sample+".pdf"
            c3.SaveAs(fn)

    ## define kinematic range ##
    kinRange = [] 
    ## if qcd binning add 11-13 to beginning of kinRange
    if(sample=='ldp' or sample=='hdp'):
        kinRange+=range(11,14)
    kinRange+=range(1,11)

    ## define the dictionaries to transcribe the 
    ## double ratio bins to the analysis bins    
    mhtDict = {
            1: 0,
            2: 0,
            3: 0,
            4: 1,
            5: 1,
            6: 1,
            7: 2,
            8: 2,
            9: 3,
            10: 3,
            11: 4,
            12: 4,
            13: 4,
            -1: 0,
    }
    htDict = {
            1: 0,
            2: 1,
            3: 2,
            4: 0,
            5: 1,
            6: 2,
            7: 3,
            8: 1,
            9: 2,
            10: 1,
            11: 2,
            12: 4,
            13: 5,
            -1: 0,
    }
    
    dr = dr_out[0][0]
    edr = dr_out[0][1]
    
    effFile = ROOT.TFile("../plots/histograms/effHists.root")
        
    h_pur_m = effFile.Get("h_pur_m")
    h_pur_e = effFile.Get("h_pur_e")
    h_trig_m = effFile.Get("h_trig_m1")    
    h_trig_e = effFile.Get("h_trig_e1")    
    
    ## hard code the btag SF and error for now
    btagSF = 1.
    btagSFerror = 0.005
    
    ## hard coded lepton SFs
    lepSF = 0.95
    eLepSF = 0.05
    
    ## trigger
    trig = (h_trig_m.GetBinContent(1)+h_trig_e.GetBinContent(1))/2.
    trigError = (h_trig_m.GetBinError(1)+h_trig_e.GetBinError(1))/2.
    
    ## purity just take 3-4 jet 0b
    pur = (h_pur_m.GetBinContent(4)+h_pur_e.GetBinContent(4))/2.
    purError = (h_pur_m.GetBinError(4)+h_pur_e.GetBinError(4))/2.
        
    ## had to add this funny business to 
    ## remove the bins that don't have
    ## any events
    removedBins = RA2b.getRemovedBins(nbRange=0,kinRange=kinRange)
    subtractBins = removedBins[0]
    avoidBins = removedBins[1]
    uncutBin = 0
    Bin = 1
    
    print "DR | DR cv error | DR shape up | DR shape down | 0b purity | purity error | trig eff | trig eff error | lepton SF | lepton SF error | btag SF | btag SF error |"
    for nj in range(1,6):
        for kin in kinRange:
            uncutBin+=1
            if(uncutBin in avoidBins):
                continue
            mht = mhtDict[kin]
            ht = htDict[kin]
            EdrUp = max(dr_out[1]['NJets'][min(max(nj,1),5)-1][0], dr_out[1]['HT'][ht][0], dr_out[1]['MHT'][mht][0])
            EdrUp = (max(EdrUp**2-edr**2,0))**0.5
            EdrDown = max(dr_out[1]['NJets'][min(max(nj,1),5)-1][1], dr_out[1]['HT'][ht][1], dr_out[1]['MHT'][mht][1])
            EdrDown = (max(EdrDown**2-edr**2,0))**0.5
            print str(round(dr,4))+" | "+str(round(edr,4))+" | "+str(round(EdrUp,4))+" | "+str(round(EdrDown,4))+" | "+str(round(pur,4))+" | "+str(round(purError/pur,4))+" | "+str(round(trig,4))+" | "+str(round(trigError/trig,4))+" | "+str(round(lepSF,4))+" | "+str(round(eLepSF,4))+" | "+str(round(btagSF,4))+" | "+str(round(btagSFerror,4)) + " | "
            Bin+=1

    print "\a"

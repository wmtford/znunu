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
runYear = "Run2ldpnominal"
postDRfit = True
plotRefLine = False

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
            dataZllFile = ROOT.TFile('../outputs/histsDY_Run2v161617.root')
            dataPhotonFile = ROOT.TFile('../outputs/histsPhoton_Run2v161617.root')
            mcZllFile = ROOT.TFile('../outputs/histsDYMC_Run2v161617.root')
            if (postDRfit):
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_Run2v161617_DRr2wt.root')
            else:
                mcPhotonFile = ROOT.TFile('../outputs/histsGjets_Run2v161617.root')
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
                                         iPeriod = iPeriod, plotRefLine = plotRefLine)

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
            # c3.SaveAs("DR_HT_"+sample+".root")

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
            4: 3,
            5: 1,
            6: 2,
            7: 1,
            8: 2,
            9: 4,
            10: 5,
            11: 0,
            12: 1,
            13: 2,
            -1: 0,
    }

    dr = dr_out[0][0]
    edr = dr_out[0][1]
    
    ## hard code the btag SF and error for now
    btagSF = 1.
    btagSFerror = 0.005
    
    ## lepton SFs
    # lepSF = 0.95
    # eLepSF = 0.05
    # mcZllFile = ROOT.TFile('../src/histsDYMC2018.root')  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # dataZllFile = ROOT.TFile('../src/histsZll2018AB.root')  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    h_SF_m = mcZllFile.Get("hSFwt_DR_dymm")
    h_SFerr_m = mcZllFile.Get("hSFsys_DR_dymm")
    h_SF_e = mcZllFile.Get("hSFwt_DR_dyee")
    h_SFerr_e = mcZllFile.Get("hSFsys_DR_dyee")
    mcy_m = h_SF_m.GetSumOfWeights()
    mcy_e = h_SF_e.GetSumOfWeights()
    lepSF = (mcy_m*h_SF_m.GetMean() + mcy_e*h_SF_e.GetMean()) / (mcy_m + mcy_e)
    eLepSF = (mcy_m*h_SFerr_m.GetMean() + mcy_e*h_SFerr_e.GetMean()) / (mcy_m + mcy_e)
    
    ## purity as function of Njets bin
    #  (Actually, we don't use the central values, and the errors should be
    #  zeroed out here and accounted for in the extraplation file.)
    #  This assumes Run2 era binning
    #  First get mumu, ee yields for each NJets bin
    hNJetsCC_zmm_da = dataZllFile.Get("hNJets_DRCC_zmm")
    day_m = []
    day_m.append(hNJetsCC_zmm_da.GetBinContent(1) + hNJetsCC_zmm_da.GetBinContent(2))
    day_m.append(hNJetsCC_zmm_da.GetBinContent(3) + hNJetsCC_zmm_da.GetBinContent(4))
    day_m.append(hNJetsCC_zmm_da.GetBinContent(5) + hNJetsCC_zmm_da.GetBinContent(6))
    day_m.append(hNJetsCC_zmm_da.GetBinContent(7) + hNJetsCC_zmm_da.GetBinContent(8))
    day_m.append(hNJetsCC_zmm_da.GetBinContent(9))
    hNJetsCC_zee_da = dataZllFile.Get("hNJets_DRCC_zee")
    day_e = []
    day_e.append(hNJetsCC_zee_da.GetBinContent(1) + hNJetsCC_zee_da.GetBinContent(2))
    day_e.append(hNJetsCC_zee_da.GetBinContent(3) + hNJetsCC_zee_da.GetBinContent(4))
    day_e.append(hNJetsCC_zee_da.GetBinContent(5) + hNJetsCC_zee_da.GetBinContent(6))
    day_e.append(hNJetsCC_zee_da.GetBinContent(7) + hNJetsCC_zee_da.GetBinContent(8))
    day_e.append(hNJetsCC_zee_da.GetBinContent(9))
    effFile = ROOT.TFile("../plots/histograms/effHists.root")
    h_pur_m = effFile.Get("h_pur_m")
    h_pur_e = effFile.Get("h_pur_e")
    pur_m = []
    pur_m.append(h_pur_m.GetBinContent(1))
    pur_m.append(h_pur_m.GetBinContent(4))
    pur_m.append(h_pur_m.GetBinContent(8))
    purError_m = []
    purError_m.append(h_pur_m.GetBinError(1))
    purError_m.append(h_pur_m.GetBinError(4))
    purError_m.append(h_pur_m.GetBinError(8))
    pur_e = []
    pur_e.append(h_pur_e.GetBinContent(1))
    pur_e.append(h_pur_e.GetBinContent(4))
    pur_e.append(h_pur_e.GetBinContent(8))
    purError_e = []
    purError_e.append(h_pur_e.GetBinError(1))
    purError_e.append(h_pur_e.GetBinError(4))
    purError_e.append(h_pur_e.GetBinError(8))

    ## trigger
    # h_trig_m = effFile.Get("h_trig_m")    
    # h_trig_e = effFile.Get("h_trig_e")    
    # trig = (h_trig_m.GetBinContent(1)+h_trig_e.GetBinContent(1))/2.
    # trigError = (h_trig_m.GetBinError(1)+h_trig_e.GetBinError(1))/2.
    ## Now absorb trigger efficiency into SF
    trig = 1
    trigError = 0
    
    ## had to add this funny business to 
    ## remove the bins that don't have
    ## any events
    kinSkip = [1,4,11]
    njSkip = [4,5]
    Bin = 1
    
    print "DR | DR cv error | DR shape up | DR shape down | 0b purity | purity error | trig eff | trig eff error | lepton SF | lepton SF error | btag SF | btag SF error |"
    for nj in range(1,6):
        indpur = min(3, nj)-1
        for kin in kinRange:
            if (nj in njSkip and kin in kinSkip):
                continue
            mht = mhtDict[kin]
            ht = htDict[kin]
            EdrUp = max(dr_out[1]['NJets'][min(max(nj,1),5)-1][0], dr_out[1]['HT'][ht][0], dr_out[1]['MHT'][mht][0])
            EdrUp = (max(EdrUp**2-edr**2,0))**0.5
            EdrDown = max(dr_out[1]['NJets'][min(max(nj,1),5)-1][1], dr_out[1]['HT'][ht][1], dr_out[1]['MHT'][mht][1])
            EdrDown = (max(EdrDown**2-edr**2,0))**0.5
            daytot = day_m[nj-1] + day_e[nj-1]
            if (daytot > 0):
                purav = (day_m[nj-1]*pur_m[indpur] + day_e[nj-1]*pur_e[indpur]) / daytot
                purErrav = (day_m[nj-1]*purError_m[indpur] + day_e[nj-1]*purError_e[indpur]) / (daytot*purav)
            else:
                purav = (pur_m[indpur] + pur_e[indpur]) / 2
                purErrav = (purError_m[indpur] + purError_e[indpur]) / (2*purav)
            print ("%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |%7.4f |" %
                   (dr, edr, EdrUp, EdrDown, purav, 0*purErrav, trig, trigError/trig, lepSF, eLepSF, btagSF, btagSFerror)
                   )
            Bin+=1

    print "\a"

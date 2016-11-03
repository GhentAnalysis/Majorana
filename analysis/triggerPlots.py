#! /usr/bin/env python

from ROOT import TFile, TTree, TH2D, TCanvas, gROOT
import ROOT
from collections import defaultdict
import os

ROOT.gROOT.SetBatch(True)


def makeCumul(hist):
  cumulHist = hist.Clone()
  for i in range(hist.GetNbinsX()+1):
    for j in range(hist.GetNbinsX()+1):
      cumulHist.SetBinContent(hist.GetBin(i, j), 0)
      for ii in range(i, hist.GetNbinsX()+2):
        for jj in range(j, hist.GetNbisY()+2):
           cumulHist.SetBinContent(hist.GetBin(i, j), cumulHist.GetBinContent(i, j) + hist.GetBinContent(ii, jj))
  return cumulHist


def passTriggers(tree, triggerList):
  for trigger in triggerList:
    if not hasattr(tree, trigger): continue
    if getattr(tree, trigger + '_prescale') != 1:
      print 'Warning: prescaled trigger: ' + trigger + ' (' + hasattr(tree, trigger + '_prescale') + ')'
    if getattr(tree, trigger):
      return True
  return False


def run(inputTree, triggers, dir):
  try: os.mkdir(dir)
  except: pass
  hists = {}
  ptThresholds = ['10to15','15to20','20to30','30to40','40to50']
  for ptThreshold in ptThresholds:
    hists[ptThreshold] = {}
    for triggerApplied in [True, False]:
      hists[ptThreshold][triggerApplied] = {}
      for channel in ['eee', 'eemu', 'emumu', 'mumumu']: 
        name = channel + "_pt" + str(ptThreshold) + ("" if triggerApplied else "_noTrig")
        hists[ptThreshold][triggerApplied][channel] = TH2D(name, name, 9, 5, 50, 9, 5, 50)

  # loop over events
  for i in range(inputTree.GetEntries()):
    inputTree.GetEntry(i)

    if inputTree._nLeptons != 3: continue
    for i in range(3):
      if inputTree._flavors[i] == 0:
        if inputTree._miniisolation[i] > 0.2: continue
        if not inputTree._passedMVA90[i]: continue
      else:
        if not inputTree._istight[i]: continue
     
    # flavor of muon is 1, flavor of electron is 0 
    nMuon = sum(inputTree._flavors[i] for i in range(3))
    if   nMuon == 3: channel = "mumumu"
    elif nMuon == 2: channel = "eemu"
    elif nMuon == 1: channel = "emumu"
    elif nMuon == 0: channel = "eee"

    pt = [i for i in inputTree._lPt]
    pt.sort(reverse=True)
    for ptThreshold in ptThresholds:
      if isinstance(ptThreshold, int): 
        if pt[0] < ptThreshold: continue
      else:
        ptMin = int(ptThreshold.split('to')[0])
        ptMax = int(ptThreshold.split('to')[1])
        if pt[0] < ptMin or pt[0] > ptMax: continue

      hists[ptThreshold][False][channel].Fill(pt[2], pt[1])
      if passTriggers(inputTree, triggers[nMuon]): 
        hists[ptThreshold][True][channel].Fill(pt[2], pt[1])

  ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/Majorana/analysis/tdrstyle.C")
  ROOT.setTDRStyle()
  ROOT.gROOT.SetStyle('tdrStyle')
  ROOT.gStyle.SetPaintTextFormat("3.2f")
  ROOT.gStyle.SetPadRightMargin(0.12)
  ROOT.gStyle.SetPadTopMargin(0.15)
  ROOT.gStyle.SetTitleX(0.005)
  ROOT.gStyle.SetTitleY(0.985)
  ROOT.gStyle.SetOptTitle(1)
  ROOT.gROOT.ForceStyle()
 

  for ptThreshold in ptThresholds:
    for channel in ['eee', 'eemu', 'emumu', 'mumumu']:
      for useCumul in [False]: # do not look at cumulative distributions, because they make it sample-dependent
        if useCumul:
	  hists[ptThreshold][False][channel] = makeCumul(hists[ptThreshold][False][channel])
	  hists[ptThreshold][True][channel]  = makeCumul(hists[ptThreshold][True][channel])
	eff = hists[ptThreshold][True][channel].Clone()
	eff.Divide(hists[ptThreshold][False][channel])
   
	if hasattr("ROOT","c1"): del ROOT.c1 
	c1 = TCanvas(eff.GetName() + ('_cumul' if useCumul else ''), eff.GetName())
	c1.cd()
	eff.GetZaxis().SetRangeUser(0, 1);
	eff.GetXaxis().SetTitle("trailing lepton p_{T} [Gev]")
	eff.GetYaxis().SetTitle("subleading lepton p_{T} [Gev]")
        if isinstance(ptThreshold, int): ptRange = "p_{T} > " + str(ptThreshold) + " Gev"
        else:                            ptRange = ptThreshold.split('to')[0] + ' GeV < p_{T} < ' + ptThreshold.split('to')[1] + ' GeV'
	eff.SetTitle(channel + " channel, leading lepton " + ptRange + (' (cumulative distributions)' if useCumul else ''))
	eff.Draw("COLZ TEXT")
	c1.RedrawAxis()
	c1.Print(os.path.join(dir, eff.GetName() + ('_cumul' if useCumul else '') + '.pdf'))
	c1.Print(os.path.join(dir, eff.GetName() + ('_cumul' if useCumul else '') + '.png'))

  return True

import glob
chain = ROOT.TChain('trileptonProducer/trileptonTree')
# listOfFiles = glob.glob('/pnfs/iihe/cms/store/user/tomc/majorana/WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8/crab_test2/161030_121836/0000/trilepton*.root') # of course again pnfs problems
listOfFiles = glob.glob('/user/tomc/public/majorana/WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8/local_test2/trilepton*.root')
for i in listOfFiles:
  chain.Add(i)

# 3l triggers only
triggers_3l = {3: ['HLT_TripleMu_12_10_5'],
               2: ['HLT_DiMu9_Ele9_CaloIdL_TrackIdL'],
               1: ['HLT_Mu8_DiEle12_CaloIdL_TrackIdL'],
               0: ['HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL']}
run(chain, triggers_3l, '3l')


# 2l triggers
triggers_2l = {2: ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ','HLT_Mu30_TkMu11','HLT_Mu40_TkMu11'],
               1: ['HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL'],
               0: ['HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL','HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ','HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL']}
triggers_3l2l = {3: triggers_3l[3] + triggers_2l[2],
                 2: triggers_3l[2] + triggers_2l[2] + triggers_2l[1],
                 1: triggers_3l[1] + triggers_2l[1] + triggers_2l[0],
                 0: triggers_3l[0] + triggers_2l[0]}
run(chain, triggers_3l2l, '3l2l')


# 1l triggers (removed HLT_IsoMu24_eta2p1, HLT_IsoTkMu24_eta2p1, HLT_Ele30_WPTight_Gsf, HLT_Ele30_eta2p1_WPTight_Gsf)
triggers_1l = {1: ['HLT_Mu50','HLT_TkMu50','HLT_IsoMu27','HLT_IsoTkMu27','HLT_IsoMu24','HLT_IsoMu24','HLT_IsoMu22_eta2p1','HLT_IsoTkMu22_eta2p1','HLT_IsoTkMu22','HLT_Mu45_eta2p1'],
               0: ['HLT_Ele32_WPTight_Gsf','HLT_Ele32_eta2p1_WPTight_Gsf','HLT_Ele27_WPTight_Gsf','HLT_Ele27_eta2p1_WPTight_Gsf','HLT_Ele27_eta2p1_WPTight_Gsf','HLT_Ele25_eta2p1_WPTight_Gsf','HLT_Ele25_WPTight_Gsf']}
triggers_3l2l1l = {3: triggers_3l2l[3] + triggers_1l[1],
                   2: triggers_3l2l[2] + triggers_1l[1] + triggers_1l[0],
                   1: triggers_3l2l[1] + triggers_1l[1] + triggers_1l[0],
                   0: triggers_3l2l[0] + triggers_1l[0]}
run(chain, triggers_3l2l1l, '3l2l1l')

triggers_only1l = {3: triggers_1l[1],
                   2: triggers_1l[1] + triggers_1l[0],
                   1: triggers_1l[1] + triggers_1l[0],
                   0: triggers_1l[0]}
run(chain, triggers_only1l, '1l')

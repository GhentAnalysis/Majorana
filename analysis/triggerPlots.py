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


def run(inputTree):
  hists = {}
  ptThresholds = [20,25,30,35,40,50]
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
     
    # flavor of muon is 1, flavor of electron is 0 
    nMuon = sum(inputTree._flavors[i] for i in range(3))
    if   nMuon == 3: channel = "mumumu"
    elif nMuon == 2: channel = "eemu"
    elif nMuon == 1: channel = "emumu"
    elif nMuon == 0: channel = "eee"

    if   nMuon == 3: trigger = inputTree.HLT_TripleMu_12_10_5
    elif nMuon == 2: trigger = inputTree.HLT_DiMu9_Ele9_CaloIdL_TrackIdL
    elif nMuon == 1: trigger = inputTree.HLT_Mu8_DiEle12_CaloIdL_TrackIdL
    elif nMuon == 0: trigger = inputTree.HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL

    pt = [i for i in inputTree._lPt]
    pt.sort(reverse=True)
    for ptThreshold in ptThresholds:
      if pt[0] < ptThreshold: continue
      hists[ptThreshold][False][channel].Fill(pt[2], pt[1])
      if trigger: hists[ptThreshold][True][channel].Fill(pt[2], pt[1])

  ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/Majorana/analysis/tdrstyle.C")
  ROOT.setTDRStyle()
  ROOT.gROOT.SetStyle('tdrStyle')
  ROOT.gStyle.SetPaintTextFormat("3.2f")
  ROOT.gStyle.SetPadRightMargin(0.12)
  ROOT.gStyle.SetPadTopMargin(0.15)
  ROOT.gStyle.SetTitleX(0.005)
  ROOT.gStyle.SetTitleY(0.985)
  ROOT.gStyle.SetOptTitle(1);
  
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
	eff.SetTitle(channel + " channel, leading lepton p_{T} > " + str(ptThreshold) + " GeV" + (' (cumulative distributions)' if useCumul else ''))
	eff.Draw("COLZ TEXT")
	c1.RedrawAxis()
	c1.Print(eff.GetName() + ('_cumul' if useCumul else '') + '.pdf')

  return True

import glob
chain = ROOT.TChain('trileptonProducer/trileptonTree')
listOfFiles = glob.glob('/pnfs/iihe/cms/store/user/tomc/majorana/WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8/crab_20161027_134409/161027_114444/0000/trilepton*.root')
for i in listOfFiles:
  chain.Add(i)

run(chain)

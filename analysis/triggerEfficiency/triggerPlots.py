#! /usr/bin/env python

import ROOT
from collections import defaultdict
import os
import Majorana.tools.samples      as samples
import Majorana.tools.style        as style
import Majorana.analysis.selectors as selectors

samples.init(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Majorana', 'analysis', 'samples.txt'))

triggers_3l     = {3: ['HLT_TripleMu_12_10_5'],
                   2: ['HLT_DiMu9_Ele9_CaloIdL_TrackIdL'],
                   1: ['HLT_Mu8_DiEle12_CaloIdL_TrackIdL'],
                   0: ['HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL']}

triggers_2l     = {2: ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL','HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL','HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL','HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ','HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ','HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ'],
                   1: ['HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ','HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ','HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ','HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ'],
                   0: ['HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ']}

triggers_1l     = {1: ['HLT_IsoMu24','HLT_IsoTkMu24'],
                   0: ['HLT_Ele27_WPTight_Gsf']}

triggers_3l2l   = {3: triggers_3l[3] + triggers_2l[2],
                   2: triggers_3l[2] + triggers_2l[2] + triggers_2l[1],
                   1: triggers_3l[1] + triggers_2l[1] + triggers_2l[0],
                   0: triggers_3l[0] + triggers_2l[0]}

triggers_3l2l1l = {3: triggers_3l2l[3] + triggers_1l[1],
                   2: triggers_3l2l[2] + triggers_1l[1] + triggers_1l[0],
                   1: triggers_3l2l[1] + triggers_1l[1] + triggers_1l[0],
                   0: triggers_3l2l[0] + triggers_1l[0]}


def makeCumul(hist):
  cumulHist = hist.Clone()
  for i in range(hist.GetNbinsX()+1):
    for j in range(hist.GetNbinsX()+1):
      cumulHist.SetBinContent(hist.GetBin(i, j), 0)
      for ii in range(i, hist.GetNbinsX()+2):
        for jj in range(j, hist.GetNbisY()+2):
           cumulHist.SetBinContent(hist.GetBin(i, j), cumulHist.GetBinContent(i, j) + hist.GetBinContent(ii, jj))
  return cumulHist


listOfWarnings = []
def passTriggers(tree, triggerList):
  for trigger in triggerList:
    if not hasattr(tree, trigger):                  listOfWarnings.append('Warning: trigger ' + trigger + ' not in tree!')
    elif getattr(tree, trigger + '_prescale') != 1: listOfWarnings.append('Warning: prescaled trigger: ' + trigger + ' (' + str(getattr(tree, trigger + '_prescale')) + ') in run ' + str(tree._runNb))
    elif getattr(tree, trigger):                    return True
  return False


def run(inputTree, triggers, dir, nLepton):
  try: os.makedirs(dir)
  except: pass
  hists = {}

  if nLepton ==1: style.setDefault()
  else:           style.setDefault2D()

  if nLepton == 3:   channels = ['eee', 'eemu', 'emumu', 'mumumu']
  elif nLepton == 2: channels = ['ee', 'emu', 'mumu']
  elif nLepton == 1: channels = ['e', 'mu']

  ptThresholds = ['10to15','15to20','20to30','30to40','40to50'] if nLepton == 3 else [None]
  for ptThreshold in ptThresholds:
    hists[ptThreshold] = {}
    for triggerApplied in [True, False]:
      hists[ptThreshold][triggerApplied] = {}
      for channel in channels:
	if nLepton==3: name = channel + "_pt" + str(ptThreshold) + ("" if triggerApplied else "_noTrig")
	else:          name = channel + ("" if triggerApplied else "_noTrig")
	if nLepton==1: hists[ptThreshold][triggerApplied][channel] = ROOT.TH1D(name, name, 26, 5, 70)
        else:          hists[ptThreshold][triggerApplied][channel] = ROOT.TH2D(name, name, 13, 5, 70, 13, 5, 70)

  # loop over events
  print inputTree.GetEntries()
  for i in range(inputTree.GetEntries()):
    inputTree.GetEntry(i)
    if i%100000==0: print i, '/', inputTree.GetEntries()

    if inputTree._nLeptons != nLepton: continue
    if not all(selectors.leptonSelector(inputTree, i) for i in range(nLepton)): continue
     
    # flavor of muon is 1, flavor of electron is 0 
    nMuon = sum(inputTree._flavors[i] for i in range(nLepton))
    if nLepton == 3:
      if   nMuon == 3: channel = "mumumu"
      elif nMuon == 2: channel = "emumu"
      elif nMuon == 1: channel = "eemu"
      elif nMuon == 0: channel = "eee"
    elif nLepton == 2:
      if   nMuon == 2: channel = "mumu"
      elif nMuon == 1: channel = "emu"
      elif nMuon == 0: channel = "ee"
    elif nLepton == 1:
      if   nMuon == 1: channel = "mu"
      elif nMuon == 0: channel = "e"

    pt = [i for i in inputTree._lPt]
    pt.sort(reverse=True)

    if nLepton == 3:
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
    elif nLepton == 2:
	hists[None][False][channel].Fill(pt[1], pt[0])
	if passTriggers(inputTree, triggers[nMuon]): 
	  hists[None][True][channel].Fill(pt[1], pt[0])
    elif nLepton == 1:
	hists[None][False][channel].Fill(pt[0])
	if passTriggers(inputTree, triggers[nMuon]):
	  hists[None][True][channel].Fill(pt[0])

  for ptThreshold in ptThresholds:
    for channel in channels:
      for useCumul in [False]: # do not look at cumulative distributions, because they make it sample-dependent
        if useCumul:
	  hists[ptThreshold][False][channel] = makeCumul(hists[ptThreshold][False][channel])
	  hists[ptThreshold][True][channel]  = makeCumul(hists[ptThreshold][True][channel])
	hists[ptThreshold][False][channel].Sumw2()
	hists[ptThreshold][True][channel].Sumw2()
	eff = hists[ptThreshold][True][channel].Clone()
	eff.Divide(hists[ptThreshold][False][channel])
   
	if hasattr("ROOT","c1"): del ROOT.c1 
	c1 = ROOT.TCanvas(eff.GetName() + ('_cumul' if useCumul else ''), eff.GetName())
	c1.cd()

        if nLepton == 1:
	  eff.GetXaxis().SetTitle("lepton p_{T} [Gev]")
        else:
	  eff.GetZaxis().SetRangeUser(0, 1)
	  eff.GetXaxis().SetTitle("trailing lepton p_{T} [Gev]")
	  eff.GetYaxis().SetTitle(("sub" if nLepton==3 else "") + "leading lepton p_{T} [Gev]")

        if ptThreshold:
	  if isinstance(ptThreshold, int): ptRange = "p_{T} > " + str(ptThreshold) + " Gev"
	  else:                            ptRange = ptThreshold.split('to')[0] + ' GeV < p_{T} < ' + ptThreshold.split('to')[1] + ' GeV'
	  eff.SetTitle(channel + " channel, leading lepton " + ptRange + (' (cumulative distributions)' if useCumul else ''))
        else:
	  eff.SetTitle(channel + " channel")

	if nLepton==1: eff.Draw("E")
        else:          eff.Draw("COLZ TEXT")

	c1.RedrawAxis()
	c1.Print(os.path.join(dir, eff.GetName() + ('_cumul' if useCumul else '') + '.pdf'))
	c1.Print(os.path.join(dir, eff.GetName() + ('_cumul' if useCumul else '') + '.png'))

  return True


for pd in ['WZ','JetHT','MET']:
  chain = samples.getTree(pd, productionLabel='triggerEfficiency_v3', shortDebug=False)

  run(chain, triggers_3l, pd + '/3l', 3)
  run(chain, triggers_2l, pd + '/2l', 2)
  run(chain, triggers_1l, pd + '/1l', 1)
 #run(chain, triggers_3l2l1l, '3l2l1l', 3)

for w in sorted(set(listOfWarnings)): print w

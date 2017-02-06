#! /usr/bin/env python

import ROOT
from collections import defaultdict
import os
import Majorana.tools.samples      as samples
import Majorana.tools.style        as style
import Majorana.analysis.selectors as selectors


triggerSamples      = ['WZ','WJets','MET','JetHT']
triggerCombinations = ['1l','2l','3l','2l1l','3l2l1l','1l1l']

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--sample',         action='store',      default=None, nargs='?', choices=triggerSamples)
argParser.add_argument('--type',           action='store',      default=None, nargs='?', choices=triggerCombinations)
argParser.add_argument('--isChild',        action='store_true', default=False)
argParser.add_argument('--dryRun',         action='store_true', default=False)
argParser.add_argument('--runLocal',       action='store_true', default=False)
args = argParser.parse_args()



samples.init(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Majorana', 'analysis', 'samples.txt'))

triggers_3l     = {3: ['HLT_TripleMu_12_10_5'],
                   2: ['HLT_DiMu9_Ele9_CaloIdL_TrackIdL'],
                   1: ['HLT_Mu8_DiEle12_CaloIdL_TrackIdL'],
                   0: ['HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL']}

triggers_2l     = {2: ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ','HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ','HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ','HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL','HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL','HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL'],
                   1: ['HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ','HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ','HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ','HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL','HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL'],
                   0: ['HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ']}

triggers_1l     = {1: ['HLT_IsoMu24','HLT_IsoTkMu24'],
                   0: ['HLT_Ele27_WPTight_Gsf']}

triggers_1l1l   = {2: triggers_1l[1],
                   1: triggers_1l[1] + triggers_1l[0],
                   0: triggers_1l[0]}

triggers_2l1l   = {2: triggers_2l[2] + triggers_1l[1],
                   1: triggers_2l[1] + triggers_1l[1] + triggers_1l[0],
                   0: triggers_2l[0] + triggers_1l[0]}

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

  if nLepton == 3:   channels = ['eee', 'eemu', 'emue', 'muee', 'emumu', 'muemu', 'mumue', 'mumumu']
  elif nLepton == 2: channels = ['ee', 'emu', 'mue', 'mumu']
  elif nLepton == 1: channels = ['e', 'mu']

  ptThresholds = ['10to15','15to20','20to30','30to40','40to50','50to70','70to100'] if nLepton == 3 else [None]
  for ptThreshold in ptThresholds:
    hists[ptThreshold] = {}
    for triggerApplied in [True, False]:
      hists[ptThreshold][triggerApplied] = {}
      for channel in channels:
	if nLepton==3: name = channel + "_pt" + str(ptThreshold) + ("" if triggerApplied else "_noTrig")
	else:          name = channel + ("" if triggerApplied else "_noTrig")
	if nLepton==1: hists[ptThreshold][triggerApplied][channel] = ROOT.TH1D(name, name, 38, 5, 100)
  else:          hists[ptThreshold][triggerApplied][channel] = ROOT.TH2D(name, name, 19, 5, 100, 19, 5, 100)

  def getPtThreshold(pt):
    for ptThreshold in ptThresholds:
      ptMin = int(ptThreshold.split('to')[0])
      ptMax = int(ptThreshold.split('to')[1])
      if pt > ptMin or pt < ptMax: return ptThreshold

  # loop over events
  print inputTree.GetEntries()
  for i in range(inputTree.GetEntries()):
    inputTree.GetEntry(i)
    if i%100000==0: print i, '/', inputTree.GetEntries()

    if inputTree._nLeptons != nLepton: continue
    if not all(selectors.leptonSelector(inputTree, i) for i in range(nLepton)): continue
 
    # flavor of muon is 1, flavor of electron is 0 
    ptAndFlavour = [(inputTree._lPt[i], inputTree._flavors[i]) for i in range(nLepton)]
    def getSortKey(item): return item[0]
    ptAndFlavour.sort(reverse=True, key=getSortKey)
    pt      = [i[0] for i in ptAndFlavour]
    flavour = [i[1] for i in ptAndFlavour]

    nMuon   = sum(i for i in flavour)
    channel = "".join(["mu" if i==1 else "e" for i in flavour])

    if nLepton == 3: ptThreshold = getPtThreshold(pt[0])
    else:            ptThreshold = None

    def fill(*args):
      hists[ptThreshold][False][channel].Fill(*args)
      if passTriggers(inputTree, triggers[nMuon]):
	      hists[ptThreshold][True][channel].Fill(*args)

    if   nLepton == 3: fill(pt[2], pt[1])
    elif nLepton == 2: fill(pt[1], pt[0])
    elif nLepton == 1: fill(pt[0])

  for ptThreshold in ptThresholds:
    for channel in channels:
      hists[ptThreshold][False][channel].Sumw2()
      hists[ptThreshold][True][channel].Sumw2()
      eff = hists[ptThreshold][True][channel].Clone()
      eff.Divide(hists[ptThreshold][False][channel])
 
      c = ROOT.TCanvas(eff.GetName(), eff.GetName())
      c.cd()

      if nLepton == 1:
        eff.GetXaxis().SetTitle("lepton p_{T} [Gev]")
        eff.GetYaxis().SetRangeUser(0, 1)
      else:
        eff.GetZaxis().SetRangeUser(0, 1)
        eff.GetXaxis().SetTitle("trailing lepton p_{T} [Gev]")
        eff.GetYaxis().SetTitle(("sub" if nLepton==3 else "") + "leading lepton p_{T} [Gev]")

      if ptThreshold: eff.SetTitle(channel + " channel, leading lepton " + ptThreshold.split('to')[0] + ' GeV < p_{T} < ' + ptThreshold.split('to')[1] + ' GeV')
      else:           eff.SetTitle(channel + " channel")

      if nLepton==1: eff.Draw("E")
      else:          eff.Draw("COLZ TEXT")

      c.RedrawAxis()
      c.Print(os.path.join(dir, eff.GetName() + '.pdf'))
      c.Print(os.path.join(dir, eff.GetName() + '.png'))
      file = ROOT.TFile(os.path.join("triggerEfficiencies.root"),"UPDATE")
      eff.Write(dir + '_' + eff.GetName(), ROOT.TObject.kOverwrite)
      file.Close()

  return True


def launch(command, logfile):
  if args.dryRun:     print command
  elif args.runLocal: os.system(command + " --isChild &> " + logfile)
  else:               os.system("qsub -v command=\"" + command + " --isChild\" -q localgrid@cream02 -o " + logfile + " -e " + logfile + " -l walltime=10:00:00 runTriggerPlotsOnCream02.sh")

if not args.sample:
  try:    os.makedirs('log')
  except: pass
  for pd in triggerSamples:
    for type in triggerCombinations:
      launch('./triggerPlots.py --sample=' + pd + ' --type=' + type + ' --isChild', 'log/' + pd + '_' + type + '.log')
elif args.isChild:
  pd = args.sample
  chain = samples.getTree(pd, treeType='singleLep', productionLabel='triggerEfficiency_v6', shortDebug=False)

  extraLabel = ''
#  extraLabel = 'test_multiIsolation_'
  if args.type == '1l':       run(chain, triggers_1l, extraLabel + pd + '/1l', 1)
  elif args.type == '3l':     run(chain, triggers_3l, extraLabel + pd + '/3l', 3)
  elif args.type == '2l':     run(chain, triggers_2l, extraLabel + pd + '/2l', 2)
  elif args.type == '1l1l':   run(chain, triggers_1l1l, extraLabel + pd + '/1l1l', 2)
  elif args.type == '2l1l':   run(chain, triggers_2l1l, extraLabel + pd + '/2l1l', 2)
  elif args.type == '3l2l1l': run(chain, triggers_3l2l1l, extraLabel + pd +'/3l2l1l', 3)

with open('triggerWarnings.txt','a') as f:
  for w in sorted(set(listOfWarnings)): f.write(w + '\n')

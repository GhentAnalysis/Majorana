#! /usr/bin/env python
import ROOT,os

import Majorana.tools.samples      as samples
import Majorana.tools.style        as style
import Majorana.analysis.selectors as selectors
import Majorana.tools.plot         as plot

style.setDefault2D()

file = ROOT.TFile(os.path.join("triggerEfficiencies.root"))
file.cd()

channels = ['eee', 'eemu', 'emue', 'muee', 'emumu', 'muemu', 'mumue', 'mumumu']

def getEfficiency(sample, l, type):
  keys = [k for k in file.GetListOfKeys() if k.GetName().count(sample + '/' + type + '_' + l)]
  if len(keys) == 1: return keys[0].ReadObj()
  else:              
    print sample, l, type, " missing"
    return None

for lep in range(3): # probevar
  graphs = {c: ROOT.TMultiGraph() for c in channels}
  for c in channels: graphs[c].forLegend = []
  marker = 21

  for map in ['MET','JetHT']:
    marker += 1
    maps = {0: getEfficiency(map, 'e_etaBinned',  '1l_eta'), 1: getEfficiency(map, 'mu_etaBinned', '1l_eta')}

    signals = ["Majorana_M40", "Majorana_M100"]

    samples.init(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Majorana', 'analysis', 'samples.txt'))

    color = 0
    for signal in signals:
      color += 1
      name        = map + " map applied to " + signal

      effHists = {}
      for channel in channels: 
        effHists[channel] = {}
        for err in [None,'up','down']:
          effHists[channel][err] = {}
          for effWeight in [True, False]:
            effHists[channel][err][effWeight] = ROOT.TH1D(name+channel+(err if err else '')+('pass' if effWeight else 'total'), name, 19, 5, 100)

      chain = samples.getTree(signal, productionLabel='triggerEfficiency_v7', shortDebug=False)
      for i in range(chain.GetEntries()):
        chain.GetEntry(i)
        if chain._nLeptons < 3: continue
        if sum(selectors.leptonSelector(chain, i) for i in range(chain._nLeptons)) < 3: continue

        ptAndFlavour = [(chain._lPt[i], chain._flavors[i], chain._lEta[i]) for i in range(chain._nLeptons)]
        def getSortKey(item): return item[0]
        ptAndFlavour.sort(reverse=True, key=getSortKey)

        pt      = [j[0] for j in ptAndFlavour]
        flavour = [j[1] for j in ptAndFlavour]
        eta     = [j[2] for j in ptAndFlavour]

        channel = "".join(["mu" if i==1 else "e" for i in flavour[:3]]) if channels[0] else None

        for err in [None, "up", "down"]:
          eff = 1.
          for k in range(3):
            legEff                      = maps[flavour[k]].GetEfficiency(maps[flavour[k]].FindFixBin(abs(eta[k]),pt[k]))
            if err == "up":     legEff += maps[flavour[k]].GetEfficiencyErrorUp(maps[flavour[k]].FindFixBin(abs(eta[k]),pt[k]))
            elif err == "down": legEff -= maps[flavour[k]].GetEfficiencyErrorLow(maps[flavour[k]].FindFixBin(abs(eta[k]),pt[k]))
            eff *= (1-legEff)
          eff = 1-eff

          effHists[channel][err][False].Fill(pt[lep], 1.)
          effHists[channel][err][True].Fill(pt[lep], eff)

      for channel in channels:
        effErrors = {err: ROOT.TGraphAsymmErrors(effHists[channel][err][True], effHists[channel][err][False]) for err in [None, "up", "down"]}
        for bin in range(effErrors[None].GetN()):
          effErrors[None].SetPointEYhigh(bin, effErrors["up"].GetY()[bin]-effErrors[None].GetY()[bin])
          effErrors[None].SetPointEYlow(bin,  effErrors[None].GetY()[bin]-effErrors["down"].GetY()[bin])

        effErrors[None].SetMarkerStyle(marker)
        effErrors[None].SetMarkerColor(color)
        effErrors[None].SetLineColor(color)
        graphs[channel].forLegend.append((effErrors[None], name))
        graphs[channel].Add(effErrors[None],"P")


  for (channel, graph) in graphs.iteritems():
    c = ROOT.TCanvas(channel, channel)
    graph.SetTitle(channel + ';p_{T}(l_{' + str(lep+1) + '}) [GeV];3l efficiency')
    graph.Draw('a')

    legendCoordinates = (0.40,0.15,0.92,0.4)
    legend = ROOT.TLegend(*legendCoordinates)
    legend.SetFillStyle(0)
    legend.SetShadowColor(ROOT.kWhite)
    legend.SetBorderSize(0)
    for (g, n) in graph.forLegend:
      legend.AddEntry(g, n)
      legend.Draw()

    dir = 'plots/3l_eff_vs_lep' + str(lep+1) + '/'
    try:    os.makedirs(dir)
    except: pass
    c.Print(dir + channel + '.pdf')
    c.Print(dir + channel + '.png')

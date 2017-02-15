#! /usr/bin/env python
import ROOT,os
import Majorana.tools.style as style
import Majorana.tools.plot  as plot

style.setDefault()

file = ROOT.TFile(os.path.join("triggerEfficiencies.root"))
file.cd()
 
samples = ['WZ','WJets','MET','JetHT','JetHT-met30']
for l in ['e','mu']:
  graphs = ROOT.TMultiGraph()
  graphs.forLegend = []
  color = 0
  for sample in samples:
    keys = [k for k in file.GetListOfKeys() if (k.GetName().count(sample + '/1l_' + l)) and not k.GetName().count('etaBinned')]
    if len(keys) == 1:
      print 'Adding ' + sample
      color += 1
      graph = keys[0].ReadObj().CreateGraph()
      graph.SetMarkerStyle(20)
      graph.SetMarkerColor(color)
      graph.SetLineColor(color)
      graphs.forLegend.append((graph, sample))
      graphs.Add(graph, 'p')

  c = ROOT.TCanvas(l, l)
  graphs.SetTitle(';p_{T}(l) [GeV];1l efficiency')
  graphs.Draw('a')

  legendCoordinates = (0.40,0.15,0.92,0.4)
  legend = ROOT.TLegend(*legendCoordinates)
  legend.SetFillStyle(0)
  legend.SetShadowColor(ROOT.kWhite)
  legend.SetBorderSize(0)
  for (g, n) in graphs.forLegend:
    legend.AddEntry(g, n)
    legend.Draw()

  dir = 'plots/combine/'
  try:    os.makedirs(dir)
  except: pass
  c.Print(dir + l + '.pdf')
  c.Print(dir + l + '.png')

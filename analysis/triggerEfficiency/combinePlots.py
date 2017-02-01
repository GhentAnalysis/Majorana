#! /usr/bin/env python
import ROOT,os
import Majorana.tools.style as style
import Majorana.tools.plot  as plot

style.setDefault()

file = ROOT.TFile(os.path.join("triggerEfficiencies.root"))
file.cd()
 
samples = ['WZ','MET','JetHT']
for l in ['e','mu']:
  c1 = ROOT.TCanvas(l, l)
  histos = []
  i = 0
  for sample in samples:
    keys = [k for k in file.GetListOfKeys() if (k.GetName() == sample + '/1l_' + l)]
    if len(keys) == 1:
      print 'Adding ' + sample
      i += 1
      histos.append(keys[0].ReadObj())
      histos[-1].texName = sample
      histos[-1].style   = style.errorStyle(i)
      histos[-1].xTitle  = 'p_T(l)'
      histos[-1].yTitle  = 'eff'
  histos[0].SetName('combine/eff_1l_' + l)
  plot.draw1D(histos)

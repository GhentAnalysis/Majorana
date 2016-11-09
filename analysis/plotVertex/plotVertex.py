#! /usr/bin/env python

import ROOT
from collections import defaultdict
import os
import Majorana.tools.sample as sample
import Majorana.tools.style  as style

style.setDefault()


def getHist(inputTree):
  hist_dxy = ROOT.TH1D("dxy","dxy", 300, -1.5, 1.5)

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

      hist_dxy.Fill(inputTree._ipPV[i], inputTree._weight)
  return hist_dxy

hist_displaced = getHist(sample.getTree('Majorana_M10_ctau10mm'))
hist_prompt    = getHist(sample.getTree('Majorana_M10'))

for hist in [hist_displaced, hist_prompt]:
  hist.Scale(1/hist.Integral())

if hasattr("ROOT","c1"): del ROOT.c1 
c1 = ROOT.TCanvas("c1","c1")
c1.cd()
hist_prompt.GetXaxis().SetTitle('d_{xy} [cm]')
hist_prompt.GetYaxis().SetTitle('1/N dN/dx')
hist_prompt.SetLineColor(ROOT.kRed)
hist_prompt.Draw("HIST")
hist_displaced.Draw("HIST SAME")
c1.RedrawAxis()
c1.Print('dxy.pdf')
c1.Print('dxy.png')


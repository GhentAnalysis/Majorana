#! /usr/bin/env python

import ROOT
from collections import defaultdict
import os
import Majorana.tools.sample as sample
import Majorana.tools.style  as style

style.setDefault()


def run(inputTree):
  hist_dxy = ROOT.TH1D("dxy","dxy", 50, -10, 10)

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

      hist_dxy.Fill(inputTree._ipPV[i])

  if hasattr("ROOT","c1"): del ROOT.c1 
  c1 = ROOT.TCanvas("c1","c1")
  c1.cd()
  hist_dxy.GetXaxis().SetTitle('d_{xy}')
  hist_dxy.GetYaxis().SetTitle('events')
  hist_dxy.Draw("HIST")
  c1.RedrawAxis()
  c1.Print('dxy.pdf')
  c1.Print('dxy.png')

run(sample.getTree('Majorana_M10_ctau10mm'))

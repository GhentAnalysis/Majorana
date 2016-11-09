#! /usr/bin/env python

import ROOT
from collections import defaultdict
import os
import Majorana.tools.sample as sample
import Majorana.tools.style  as style
import Majorana.tools.plot   as plot

style.setDefault()


def getHist(inputTree):
  hist_dxy = ROOT.TH1D("dxy","dxy", 300, -1.5, 1.5)
  hist_dxy.xTitle = 'd_{xy} [cm]'
  hist_dxy.yTitle = '1/N dN/dx'

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
        if inputTree._miniisolation[i] > 0.2: continue

      hist_dxy.Fill(inputTree._ipPV[i], inputTree._weight)
  return hist_dxy

hist_displaced = getHist(sample.getTree('Majorana_M10_ctau10mm'))
hist_prompt    = getHist(sample.getTree('Majorana_M10'))
hist_displaced.style = style.lineStyle(ROOT.kRed)
hist_displaced.texName = 'displaced'
hist_prompt.style = style.lineStyle(ROOT.kBlack)
hist_prompt.texName = 'prompt'

for hist in [hist_displaced, hist_prompt]:
  hist.Scale(1/hist.Integral())

plot.draw1D([hist_displaced, hist_prompt])

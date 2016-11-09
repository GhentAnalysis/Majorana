import ROOT

def setDefault():
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/Majorana/tools/tdrstyle.C")
  ROOT.setTDRStyle()
  ROOT.gROOT.SetStyle('tdrStyle')
  ROOT.gStyle.SetPaintTextFormat("3.2f")
  ROOT.gStyle.SetPadRightMargin(0.12)
  ROOT.gStyle.SetPadTopMargin(0.15)
  ROOT.gStyle.SetTitleX(0.005)
  ROOT.gStyle.SetTitleY(0.985)
  ROOT.gStyle.SetOptTitle(1)
  ROOT.gROOT.ForceStyle()
 



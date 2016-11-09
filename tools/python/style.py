import ROOT

def setDefault():
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/Majorana/tools/tdrstyle.C")
  ROOT.setTDRStyle()
  ROOT.gROOT.SetStyle('tdrStyle')
  ROOT.gStyle.SetPaintTextFormat("3.2f")
  ROOT.gROOT.ForceStyle()
 
def setDefault2D():
  setDefault()
  ROOT.gStyle.SetPadRightMargin(0.12)
  ROOT.gStyle.SetPadTopMargin(0.15)
  ROOT.gStyle.SetTitleX(0.005)
  ROOT.gStyle.SetTitleY(0.985)
  ROOT.gStyle.SetOptTitle(1)
  ROOT.gROOT.ForceStyle()

def errorStyle(color, markerStyle = 20, markerSize = 1):
  def func(histo):
    histo.SetLineColor(color)
    histo.SetMarkerSize(markerSize)
    histo.SetMarkerStyle(20)
    histo.SetMarkerColor(color)
    histo.SetLineWidth(1)
    histo.drawOption = "e1"
    return 
  return func

def lineStyle(color, width = None, dotted=False, dashed=False):
  def func(histo):
    histo.SetLineColor(color)
    histo.SetMarkerSize(0)
    histo.SetMarkerStyle(0)
    histo.SetFillColor(0)
    if dotted: histo.SetLineStyle(3)
    if dashed: histo.SetLineStyle(7)
    if width: histo.SetLineWidth(width)
    histo.drawOption = "hist"
    return 
  return func

def fillStyle(color, lineColor = None):
  def func(histo):
    histo.SetLineColor(lineColor if lineColor is not None else color)
    histo.SetMarkerSize(0)
    histo.SetMarkerStyle(0)
    histo.SetFillColor(color)
    histo.drawOption = "hist"
    return 
  return func

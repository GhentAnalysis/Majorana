import ROOT



def draw1D(histos, legend='auto'):
  if hasattr("ROOT","c1"): del ROOT.c1 
  c1 = ROOT.TCanvas("c1","c1")
  c1.cd()
  for h in histos:
    h.SetTitle('')
    if hasattr(h, 'style'):  h.style(h)
    if hasattr(h, 'xTitle'): h.GetXaxis().SetTitle(h.xTitle)
    if hasattr(h, 'yTitle'): h.GetYaxis().SetTitle(h.yTitle)

  histos[0].SetMaximum(max(h.GetMaximum() for h in histos)*1.1)
  histos[0].SetMinimum(min(h.GetMinimum() for h in histos))
  histos[0].Draw(histos[0].drawOption)
  for h in histos[1:]: h.Draw(h.drawOption + " SAME")
  c1.RedrawAxis()

  if legend is not None:
    legendCoordinates = (0.50,0.93-0.05*len(histos),0.92,0.93) if legend == 'auto' else legend
    legend = ROOT.TLegend(*legendCoordinates)
    legend.SetFillStyle(0)
    legend.SetShadowColor(ROOT.kWhite)
    legend.SetBorderSize(0)
    for h in histos:
      if hasattr(h, 'texName'):
        legend.AddEntry(h, h.texName)
    legend.Draw()

  c1.Print(histos[0].GetName() + '.pdf')
  c1.Print(histos[0].GetName() + '.png')

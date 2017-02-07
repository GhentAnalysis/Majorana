import ROOT,os


def draw1D(histos, legend='auto'):
  c = ROOT.TCanvas(histos[0].GetName(), histos[0].GetName())
  c.cd()
  for h in histos:
    h.SetTitle('')
    if hasattr(h, 'style'):  h.style(h)
    if hasattr(h, 'xTitle'): h.GetXaxis().SetTitle(h.xTitle)
    if hasattr(h, 'yTitle'): h.GetYaxis().SetTitle(h.yTitle)

  histos[0].SetMaximum(max(h.GetMaximum() for h in histos)*1.1)
  histos[0].SetMinimum(min(h.GetMinimum() for h in histos))
  histos[0].Draw(histos[0].drawOption)
  for h in histos[1:]: h.Draw(h.drawOption + " SAME")
  c.RedrawAxis()

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
  
  try:    os.makedirs(os.path.dirname(histos[0].GetName()))
  except: pass
  c.Print(histos[0].GetName() + '.pdf')
  c.Print(histos[0].GetName() + '.png')

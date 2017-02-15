import ROOT

def matchedGenParticle(tree, index):
  gen = ROOT.TLorentzVector()
  lep = ROOT.TLorentzVector()
  gen.SetPtEtaPhiE(tree._lPtmc[index], tree._lEtamc[index], tree._lPhimc[index], tree._lEmc[index])
  lep.SetPtEtaPhiE(tree._lPt[index],   tree._lEta[index],   tree._lPhi[index],   tree._lE[index])
  if gen.DeltaR(lep) < 0.1: return gen
  else:                     return None


# Current 'tight' selections
def electronSelector(tree, index, ipCuts=True):
  if tree._isolation[index] > 0.1:          return False
  if not tree._trigEmulator[index]:         return False
  if not tree._isotrigEmulator[index]:      return False
  if not tree._passedMVA_SUSY[index*3]:     return False   # yes, times three, this is a result of the horrible SUSY-style coding we currently have in our trilepton.cc because someone says speed is more important than quality
  if ipCuts and tree._3dIPsig[index] > 4:   return False
  if ipCuts and tree._ipPV[index] > 0.05:   return False
  if ipCuts and tree._ipZPV[index] > 0.1:   return False
  return True

def muonSelector(tree, index, ipCuts=True):
  if tree._isolation[index] > 0.1:          return False
  if not tree._ismedium[index]:             return False
  if ipCuts and tree._3dIPsig[index] > 4:   return False
  if ipCuts and tree._ipPV[index] > 0.05:   return False
  if ipCuts and tree._ipZPV[index] > 0.1:   return False
  return True

def leptonSelector(tree, index, ipCuts=True):
  if   tree._flavors[index] == 0: return electronSelector(tree, index, ipCuts)
  elif tree._flavors[index] == 1: return muonSelector(tree, index, ipCuts)
  else:                           return False

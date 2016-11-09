#
# Keep track of used tuples here
#
import glob, os
import ROOT

defaultProductionLabel = 'test3'
paths = {'Majorana_M10_ctau10mm' : '/user/tomc/public/majorana/MajoranaNeutrino_trilepton_M-10_5f_NLO_ctau10mm',
         'WZ'                    : '/user/tomc/public/majorana/WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8'}

def getTree(sample, treeType='trilepton', productionLabel=defaultProductionLabel):
  chain           = ROOT.TChain('trileptonProducer/' + treeType + 'Tree')
  productionLabel = ('crab' if 'pnfs' in paths[sample] else 'local') + '_' + productionLabel
  listOfFiles     = glob.glob(os.path.join(paths[sample], productionLabel, treeType + '*.root'))
  for i in listOfFiles:
    chain.Add(i)
  return chain

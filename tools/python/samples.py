import glob, os
import ROOT

defaultProductionLabel = 'test3'
sampleInfos = None
lumiWeight = 1000

def init(file):
  global sampleInfos									# Access the global variable
  sampleInfos = [line.split('#')[0].strip() for line in open(file)] 			# Strip comments and \n charachters
  sampleInfos = [line.split() for line in sampleInfos if line]				# Get (name, path, x-sec) tuples
  for info in sampleInfos:
#    if 'pnfs' in info[1]: info[1] = 'dcap://maite.iihe.ac.be' + info[1]
    if   len(info) == 2: info += (None,)						# Set no xsec for data
    elif len(info) == 3: info[2] = eval(info[2])					# Transform multiplication string into xsec
    else:                exit(1)
  sampleInfos = {name : (path, xsec) for name, path, xsec in sampleInfos}		# Transform list into dictionary


def getTree(name, treeType='trilepton', productionLabel=defaultProductionLabel, shortDebug=False):
  if not sampleInfos or name not in sampleInfos.keys(): 
    print 'No dataset know with name ', name 
    exit(1)

  chain           = ROOT.TChain('trileptonProducer/' + treeType + 'Tree')
  print os.path.join(sampleInfos[name][0], '*' + productionLabel, '*', '*', treeType + '*.root')
  listOfFiles     = glob.glob(os.path.join(sampleInfos[name][0], '*' + productionLabel, treeType + '*.root'))
  listOfFiles    += glob.glob(os.path.join(sampleInfos[name][0], '*' + productionLabel, '*', '*', treeType + '*.root'))
  for i in (listOfFiles[:20] if shortDebug else listOfFiles):
    print "Adding " + i
    chain.Add(i)
  chain.sampleWeight = sampleInfos[name][1]*lumiWeight if sampleInfos[name][1] else 1.
  return chain

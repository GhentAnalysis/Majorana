#!/usr/bin/env python
import os, glob, sys

productionLabel = 'last_FR'                                                                               # Label to keep track of the tuple versio
outDir          = '/user/' + os.environ['USER'] + '/public/FR_runD_good'                                  # Output directory in case of local submission
datasets        = [dataset.strip() for dataset in open(sys.argv[1])]                                      # Get list of datasets from file given as first argument
datasets        = [dataset.split()[0] for dataset in datasets if dataset and not dataset.startswith('#')] # Clean empty and comment lines
groupFiles      = 5                                                                                       # Group files together when running locally


for dataset in datasets:
  if dataset.startswith('FAKERATE:'):
    outputName      = 'fakeRate'
    treeForFakeRate = True
    singleLep       = False
    dataset         = dataset.split(':')[-1]
  elif dataset.startswith('SINGLELEP:'):
    outputName      = 'trilepton'
    treeForFakeRate = False
    singleLep       = True
    dataset         = dataset.split(':')[-1]
  else:
    outputName      = 'trilepton'
    treeForFakeRate = False
    singleLep       = False

  if 'pnfs' in dataset or 'user' in dataset:
    if 'pnfs' in dataset: datasetName = dataset.split('/MINIAOD')[0].split('/')[-1]
    else:                 datasetName = dataset.split('/')[-1]
    print dataset

    i = 0
    j = 0
    inputFiles = []
    for file in glob.glob(dataset + ('/*/*.root' if 'pnfs' in dataset else '/*.root')):
      j          += 1
      inputFiles += [('dcap://maite.iihe.ac.be' if 'pnfs' in dataset else 'file://') + file]
      if j%groupFiles!=0: continue

      dir        = os.getcwd()
      wallTime   = '05:00:00'
      inputFile  = ','.join(inputFiles)
      outputFile = os.path.join(outDir, datasetName, 'local_' + productionLabel, outputName + '_' + str(i) + '.root')
      logFile    = os.path.join(outDir, datasetName, 'local_' + productionLabel, 'log', outputName + '_' + str(i) + '.log')

      for filename in [outputFile, logFile]:
        try:    os.makedirs(os.path.dirname(filename))
        except: pass
      
      print 'Submitting ' + inputFile + ' to cream02:'
      args  = 'dir=' + dir + ',inputFile=\"' + inputFile + '\",outputFile=' + outputFile + ',events=-1'
      args += ',isData='          + ('True' if 'Run2016' in dataset else 'False')
      args += ',treeForFakeRate=' + ('True' if treeForFakeRate      else 'False')
      args += ',singleLep='       + ('True' if singleLep            else 'False')
      os.system('qsub -v ' + args + ' -q localgrid@cream02 -o ' + logFile + ' -e ' + logFile + ' -l walltime=' + wallTime + ' runOnCream02.sh')
      i += 1
      inputFiles = []

  else: # use crab
    print 'Submitting ' + dataset + ' using crab:'
    os.environ['CRAB_PRODUCTIONLABEL'] = productionLabel
    os.environ['CRAB_DATASET']         = dataset
    os.environ['CRAB_TREEFORFAKERATE'] = 'True' if treeForFakeRate else 'False'
    os.environ['CRAB_SINGLELEP']       = 'True' if singleLep else 'False'
    os.system('crab submit -c crab.py')


#!/usr/bin/env python
import os, glob, sys

productionLabel = 'fakeRate_v3'		                                					# Label to keep track of the tuple versio
outDir          = '/user/' + os.environ['USER'] + '/public/majorana'						# Output directory in case of local submission
datasets        = [dataset.strip() for dataset in open(sys.argv[1])]						# Get list of datasets from file given as first argument
datasets        = [dataset.split()[0] for dataset in datasets if dataset and not dataset.startswith('#')]	# Clean empty and comment lines

for dataset in datasets:
  if dataset.startswith('FAKERATE:'):
    outputName      = 'fakeRate'
    treeForFakeRate = True
    dataset         = dataset.split(':')[-1]
  else:
    outputName      = 'trilepton'
    treeForFakeRate = False

  if 'pnfs' in dataset or 'user' in dataset:
    if 'pnfs' in dataset: datasetName = dataset.split('/MINIAOD')[0].split('/')[-1]
    else:                 datasetName = dataset.split('/')[-1]

    i = 0
    for file in glob.glob(dataset + ('/*/*.root' if 'pnfs' in dataset else '/*.root')):
      dir        = os.getcwd()
      wallTime   = '05:00:00'
      inputFile  = ('dcap://maite.iihe.ac.be' if 'pnfs' in dataset else 'file://') + file
      outputFile = os.path.join(outDir, datasetName, 'local_' + productionLabel, outputName + '_' + str(i) + '.root')
      logFile    = os.path.join(outDir, datasetName, 'local_' + productionLabel, 'log', outputName + '_' + str(i) + '.log')

      for filename in [outputFile, logFile]:
        try:    os.makedirs(os.path.dirname(filename))
        except: pass
      
      print 'Submitting ' + inputFile + ' to cream02:'
      args  = 'dir=' + dir + ',inputFile=' + inputFile + ',outputFile=' + outputFile + ',events=-1'
      args += ',isData='          + ('True' if 'Run2016' in dataset else 'False')
      args += ',treeForFakeRate=' + ('True' if treeForFakeRate      else 'False')
      os.system('qsub -v ' + args + ' -q localgrid@cream02 -o ' + logFile + ' -e ' + logFile + ' -l walltime=' + wallTime + ' runOnCream02.sh')
      i += 1

  else: # use crab
    print 'Submitting ' + dataset + ' using crab:'
    os.environ['CRAB_PRODUCTIONLABEL'] = productionLabel
    os.environ['CRAB_DATASET']         = dataset
    os.environ['CRAB_TREEFORFAKERATE'] = 'True' if treeForFakeRate else 'False'
    os.system('crab submit -c crab.py')


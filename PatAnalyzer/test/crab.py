from WMCore.Configuration import Configuration
import os

# Get parameters from submit script
productionLabel = os.environ['CRAB_PRODUCTIONLABEL']
dataset         = os.environ['CRAB_DATASET']
treeForFakeRate = os.environ['CRAB_TREEFORFAKERATE']
singleLep       = os.environ['CRAB_SINGLELEP']


# Crab configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = True
config.General.requestName  = dataset.split('/')[2] + '_' + productionLabel
config.General.workArea     = os.path.join('crab', productionLabel, dataset.split('/')[1] + dataset.split('/')[2])

config.section_('JobType')
config.JobType.psetName                = 'trilepton.py'
config.JobType.pyCfgParams             = (['isData=True'] if not 'SIM' in dataset else []) + ['treeForFakeRate='+treeForFakeRate, 'singleLep='+singleLep, 'events=-1']
config.JobType.pluginName              = 'analysis'
config.JobType.outputFiles             = ['fakeRate.root' if treeForFakeRate=='True' else ('singleLep.root' if singleLep=='True' else 'trilepton.root')]
config.JobType.sendExternalFolder      = True
config.JobType.allowUndistributedCMSSW = True 

config.section_('Data')
config.Data.inputDataset               = dataset
config.Data.unitsPerJob                = 100
config.Data.splitting                  = 'LumiBased'
config.Data.outLFNDirBase              = '/store/user/' + os.environ['USER'] + '/majorana/'
config.Data.publication                = False

config.section_('Site')
config.Site.storageSite                = 'T2_BE_IIHE'

config.section_('User')
config.User.voGroup                    = 'becms'

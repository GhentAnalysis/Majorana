from WMCore.Configuration import Configuration
import os

# Get parameters from submit script
productionLabel = os.environ['CRAB_PRODUCTIONLABEL']
dataset         = os.environ['CRAB_DATASET']
treeForFakeRate = os.environ['CRAB_TREEFORFAKERATE']


# Crab configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = True
config.General.requestName  = dataset.split('/')[2] + '_' + productionLabel
config.General.workArea     = os.path.join('crab', productionLabel, dataset.split('/')[1] + dataset.split('/')[2])

config.section_('JobType')
config.JobType.psetName                = 'trilepton.py'
config.JobType.pyCfgParams             = (['isData=True'] if not 'SIM' in dataset else []) + (['treeForFakeRate=True'] if treeForFakeRate else []) + ['events=-1']
config.JobType.pluginName              = 'analysis'
config.JobType.outputFiles             = ['trilepton.root' if treeForFakeRate=='False' else 'fakeRate.root']
config.JobType.allowUndistributedCMSSW = True 

config.section_('Data')
config.Data.inputDataset               = dataset
config.Data.unitsPerJob                = 1000
config.Data.splitting                  = 'LumiBased'
config.Data.outLFNDirBase              = '/store/user/' + os.environ['USER'] + '/majorana/'
config.Data.publication                = False

config.section_('Site')
config.Site.storageSite                = 'T2_BE_IIHE'

config.section_('User')
config.User.voGroup                    = 'becms'

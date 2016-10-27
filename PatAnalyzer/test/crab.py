from WMCore.Configuration import Configuration
import shutil, os

config = Configuration()
config.section_('General')
config.General.transferLogs = True

config.section_('JobType')
config.JobType.psetName    = 'trilepton.py'
#config.JobType.pyCfgParams  = ['isData=True']
config.JobType.pluginName  = 'analysis'
config.JobType.outputFiles = ['trilepton.root']
config.JobType.allowUndistributedCMSSW = True 

config.section_('Data')
config.Data.inputDataset  = '/WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
#config.Data.inputDataset  = '/DoubleMuon/Run2016B-PromptReco-v2/MINIAOD'
config.Data.unitsPerJob   = 30
config.Data.splitting     = 'LumiBased'
config.Data.outLFNDirBase = '/store/user/' + os.environ['USER'] + '/majorana/'
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_BE_IIHE'

config.section_('User')
config.User.voGroup = 'becms'

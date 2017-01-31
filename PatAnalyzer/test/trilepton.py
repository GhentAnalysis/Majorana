import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts

isData          = True
treeForFakeRate = False
singleLep       = False
#inputFile       = "file:///user/mvit/public/Majorana/MajoranaNeutrino_trilepton_M-10_5f_NLO/Majorana_trilepton_RunIISpring16MiniAODv2_96.root"
inputFile       = "/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v3/000/284/037/00000/54D1BAE3-BD9F-E611-9E45-02163E012706.root"
nEvents         = 10
outputFile      = None

def getVal(arg):
    return arg.split('=')[-1]

## loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if   "isData"          in sys.argv[i]: isData          = (getVal(sys.argv[i]) == "True")
    elif "treeForFakeRate" in sys.argv[i]: treeForFakeRate = (getVal(sys.argv[i]) == "True")
    elif "singleLep"       in sys.argv[i]: singleLep       = (getVal(sys.argv[i]) == "True")
    elif "output"          in sys.argv[i]: outputFile      = getVal(sys.argv[i])
    elif "input"           in sys.argv[i]: inputFile       = getVal(sys.argv[i])
    elif "events"          in sys.argv[i]: nEvents         = int(getVal(sys.argv[i]))

process = cms.Process("trilepton")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO' # Options: INFO, WARNING, ERROR
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source    = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nEvents)) # for debugging

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if isData: process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4'
#else:      process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
else:	   process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')


# Set up electron ID (VID framework)
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


if not outputFile:
  if treeForFakeRate: outputFile = 'fakeRate.root'
  if singleLep:       outputFile = 'singleLep.root'
  else:               outputFile = 'trilepton.root'
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))


process.trileptonProducer = cms.EDAnalyzer("trilepton",
                                          #METFilter                              = cms.InputTag("TriggerResults::PAT"),
                                           qualityCuts                            = PFTauQualityCuts,
					   isData                                 = cms.untracked.bool(isData),
					   SampleName                             = cms.untracked.string("default"),
					   electronMvaIdMap                       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
					   mvaValuesMap_HZZ                       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
	 	                           electronMvaId90Map 			  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
					   electronMvaId80Map 			  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
					   electronCutBasedIdTightMap             = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
					   electronCutBasedIdMediumMap            = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
					   genPartsLabel                          = cms.InputTag("prunedGenParticles"),
					   pdfvariablesLabel                      = cms.InputTag("generator"),
					   BeamSpotLabel                          = cms.InputTag("offlineBeamSpot"),
					   slimmedAddPileupInfoLabel              = cms.InputTag("slimmedAddPileupInfo"),
					   goodOfflinePrimaryVerticesLabel        = cms.InputTag("goodOfflinePrimaryVertices"),
					   packedPFCandidatesLabel                = cms.InputTag("packedPFCandidates"),
					   MuonLabel                              = cms.InputTag("slimmedMuons"),
					   ElectronLabel                          = cms.InputTag("slimmedElectrons"),
					   JetLabel                               = cms.InputTag("slimmedJets"),
					   METLabel                               = cms.InputTag("slimmedMETs"),
					   reducedEgammaLabel                     = cms.InputTag("reducedEgamma:reducedConversions"),
					   TauLabel                               = cms.InputTag("slimmedTaus"),
					   fixedGridRhoFastjetCentralNeutralLabel = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
					   fixedGridRhoFastjetAllLabel            = cms.InputTag("fixedGridRhoFastjetAll"),
					   prescales                              = cms.InputTag("patTrigger"),
                                           triggerResultsHLT                      = cms.InputTag("TriggerResults","","HLT"), # reHLT samples have their triggers stored in HLT2
                                           triggerResultsRECO                     = cms.InputTag("TriggerResults::RECO"),
					   exernalLHEPLabel                       = cms.InputTag("externalLHEProducer"),
                                           treeForFakeRate                        = cms.untracked.bool(treeForFakeRate),
                                           singleLep                              = cms.untracked.bool(singleLep),
					   )

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
                                                  src          = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                                  filterParams = cms.PSet(minNdof = cms.double(4.),
                                                                          maxZ    = cms.double(24.),
                                                                          maxRho  = cms.double(2.)),
                                                  filter       = cms.bool( False )
)

if isData:
   import FWCore.PythonUtilities.LumiList as LumiList
   process.source.lumisToProcess = LumiList.LumiList(filename = 'JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices
		    *process.egmGsfElectronIDSequence
		    *process.trileptonProducer)

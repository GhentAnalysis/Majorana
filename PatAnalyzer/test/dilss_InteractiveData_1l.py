import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts


skim     ="NOSkim"
isMC     = bool(False)
inputFile=""
#outputFile="results/FakeMuons_Interactive.root"
outputFile=""
def getVal(arg):
    i=0
    while i < len(arg) and arg[i] != "=": i+=1
    return arg[i+1:]

## loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "isMC" in sys.argv[i] :
        isMC=False if getVal(sys.argv[i]) == "False" else True
    elif "isMSUGRA" in sys.argv[i]:
        isMSUGRA=True
    elif "isSMS" in sys.argv[i]:
        isSMS=True
    elif "skim" in sys.argv[i]:
        skim=getVal(sys.argv[i])
    elif "output" in sys.argv[i]:
        outputFile=getVal(sys.argv[i])
    elif "input" in sys.argv[i] :
        inputFile=getVal(sys.argv[i])


if skim=="" : print "WARNING: No Skim Conditions have been provided \n"

#process = cms.Process("FakeElectrons")
process = cms.Process("pippo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'WARNING' # Options: INFO, WARNING, ERROR
process.MessageLogger.cerr.FwkReport.reportEvery = 10

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


if isMC:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    process.GlobalTag.globaltag="MCRUN2_74_V9" # tag fo 53X 2012A/B 13Jul2012 rereco
else:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    process.GlobalTag.globaltag="74X_dataRun2_reMiniAOD_v1" # tag fo 53X 2012A/B 13Jul2012 rereco

# RECO AND GEN SETUP
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')


process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             fileNames = cms.untracked.vstring(),
                             )

# input files
for i in inputFile.split(","):
    print "Adding: ", i
    process.source.fileNames.append(i)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
print "TWO"

#FakeElectrons
process.FakeElectrons = cms.EDAnalyzer("dilss13_1l", #
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults::HLT"),
                                       METLabel = cms.InputTag("slimmedMETs"),

                                       mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       SampleLabel = cms.untracked.string("ElectronsData"), # defines a piece of code to run; helps to avoid code recompilation
                                       SampleName = cms.untracked.string("DoubleLeptons"), # defines a piece of code to run; helps to avoid code recompilation
                                       eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                       eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                       eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                       eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")


                                       )

print "FOUR"


process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
                                                  , filterParams =cms.PSet(
                                                                           minNdof = cms.double( 4. )
                                                                           , maxZ    = cms.double( 24. )
                                                                           , maxRho  = cms.double( 2. ) )
                                                  , filter       = cms.bool( False ) # use only as producer
                                                  , src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
                                                 )

#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                    minimumNDOF = cms.uint32(4) ,
#                                    maxAbsZ = cms.double(24),
#                                    maxd0 = cms.double(2)
#                                  )

#process.load('RecoMET.METFilters.eeBadScFilter_cfi')

##___________________________HCAL_Noise_Filter________________________________||
#process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
#process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
#process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
#        inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
#        reverseDecision = cms.bool(False)
#        )

#process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
#        inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
#        reverseDecision = cms.bool(False)
#        )

from RecoMET.METProducers.CSCHaloData_cfi import *
from RecoMET.METProducers.EcalHaloData_cfi import *
from RecoMET.METProducers.HcalHaloData_cfi import *
from RecoMET.METProducers.GlobalHaloData_cfi import *
from RecoMET.METProducers.BeamHaloSummary_cfi import *

process.load('RecoMET.METFilters.CSCTightHalo2015Filter_cfi')
process.CSCTightHalo2015Filter.taggingMode = cms.bool(True)

process.BeamHaloId = cms.Sequence(CSCHaloData*EcalHaloData*HcalHaloData*GlobalHaloData*BeamHaloSummary)


if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices # before previously, changed 4Feb 2016
    #primaryVertexFilter
    #*process.eeBadScFilter
    #*process.HBHENoiseFilterResultProducer #produces HBHE baseline bools
    #*process.ApplyBaselineHBHENoiseFilter  #reject events based
    #*process.BeamHaloId
    #*process.CSCTightHalo2015Filter
    *process.egmGsfElectronIDSequence
    *process.FakeElectrons
    )
else:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '/localgrid/ikhvastu/CMSSW_7_4_14/src/Majorana/PatAnalyzer/test/JSON/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.json').getVLuminosityBlockRange()
    process.p = cms.Path(
    process.goodOfflinePrimaryVertices
    #primaryVertexFilter
    #*process.eeBadScFilter
    #*process.HBHENoiseFilterResultProducer #produces HBHE baseline bools
    #*process.ApplyBaselineHBHENoiseFilter  #reject events based
    *process.BeamHaloId
    *process.CSCTightHalo2015Filter
    *process.egmGsfElectronIDSequence
    *process.FakeElectrons
    )

print "SIX"



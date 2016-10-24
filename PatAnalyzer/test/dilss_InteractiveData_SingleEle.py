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
    process.GlobalTag.globaltag="74X_dataRun2_reMiniAOD_v0" # tag fo 53X 2012A/B 13Jul2012 rereco

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.Reconstruction_cff')
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

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

#if isMC:

# input files
for i in inputFile.split(","):
    print "Adding: ", i
    process.source.fileNames.append(i)
#else:
#    process.source.fileNames = cms.untracked.vstring(
#            'dcap://maite.iihe.ac.be:/pnfs/iihe/cms/ph/sc4/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v3/000/256/629/00000/8EA4C10E-F35E-E511-ABF9-02163E014108.root')
            #'dcap://maite.iihe.ac.be:/pnfs/iihe/cms/ph/sc4/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/D8C225F4-5D6B-E411-A778-0025901D42C0.root',
            #'/cms/data/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/484D51C6-2673-E411-8AB0-001E67398412.root'


#else:
   #import PhysicsTools.PythonAnalysis.LumiList as LumiList
   #process.source.lumisToProcess = LumiList.LumiList(filename = '/localgrid/depoyraz/CMSSW_7_4_12/src/SUSYAnalyzer/PatAnalyzer/test/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

# input files
   #for j in inputFile.split(","):
   #    print "Adding: ", j
   #    process.source.fileNames.append(j)
   #process.source.fileNames = cms.untracked.vstring(
   # '/scratch/lfs/lesya/FakeSync_jets_test_data_Electrons.root'
   # )
#    runOnData(process)


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
process.FakeElectrons = cms.EDAnalyzer("dilss13", #
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults::HLT"),
                                       METLabel = cms.InputTag("slimmedMETs"),

                                       #ElMediumMVALabel = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                       #ElTightMVALabel = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                       #ElMVALabel = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                       mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                       #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                       #ElMVACategoryLabel = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       #qualityCuts = PFTauQualityCuts,
                                       SampleLabel = cms.untracked.string("ElectronsData"), # defines a piece of code to run; helps to avoid code recompilation
                                       SampleName = cms.untracked.string("SingleEle"), # defines a piece of code to run; helps to avoid code recompilation
                                       eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                       eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                       eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                       eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")
                                       #eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                       #electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                       #genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                       #
                                       # ID decisions (common to all formats)
                                       #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                       # ValueMaps with MVA results

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

#process.GlobalTag.toGet = cms.VPSet(
#cms.PSet(record = cms.string('EcalLinearCorrectionsRcd'),
#tag = cms.string('EcalLinearCorrections_mc'),
#connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
#)
#)



if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
    *process.egmGsfElectronIDSequence
    *process.FakeElectrons
        )
else:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '/localgrid/ikhvastu/CMSSW_7_4_14/src/SUSYAnalyzer/PatAnalyzer/test/JSON/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.json').getVLuminosityBlockRange()
    process.p = cms.Path(process.goodOfflinePrimaryVertices * process.egmGsfElectronIDSequence * process.FakeElectrons)

print "SIX"



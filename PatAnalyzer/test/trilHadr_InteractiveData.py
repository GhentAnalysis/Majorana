import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts


skim     ="NOSkim"
isMC     = bool(False)
isMSUGRA = bool(False)
isSMS    = bool(False)
useCrab  = bool(False)
doSusyTopProjection = bool(False)
inputFile=""
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

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if isMC:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    #process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
else:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    process.GlobalTag.globaltag="80X_dataRun2_Prompt_ICHEP16JEC_v0" # tag fo 53X 2012A/B 13Jul2012 rereco

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             fileNames = cms.untracked.vstring(),
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#if isMC:

# input files
for i in inputFile.split(","):
    print "Adding: ", i
    process.source.fileNames.append(i)


process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )

#FakeElectrons
process.FakeElectrons = cms.EDAnalyzer("trilHadr", #dilss13
                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
                                       #HLTResultsLabel = cms.InputTag("TriggerResults"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults::HLT"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       qualityCuts = PFTauQualityCuts,
                                       SampleLabel = cms.untracked.string("ElectronsData"), # defines a piece of code to run; helps to avoid code recompilation
                                       SampleName = cms.untracked.string("default"),

                                       # Objects specific to AOD format
                                       #
                                       #electrons    = cms.InputTag("gedGsfElectrons"),
                                       #genParticles = cms.InputTag("genParticles"),
                                       #
                                       # Objects specific to MiniAOD format
                                       #
                                       #electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                       #genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                       #
                                       # ID decisions (common to all formats)
                                       #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                       #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                       # ValueMaps with MVA results
                                       mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                       genPartsLabel    = cms.InputTag("prunedGenParticles"),
                                       pdfvariablesLabel = cms.InputTag("generator"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       slimmedAddPileupInfoLabel = cms.InputTag("slimmedAddPileupInfo"),
                                       goodOfflinePrimaryVerticesLabel = cms.InputTag("goodOfflinePrimaryVertices"),
                                       packedPFCandidatesLabel = cms.InputTag("packedPFCandidates"),
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       METLabel = cms.InputTag("slimmedMETs"),
                                       reducedEgammaLabel = cms.InputTag("reducedEgamma:reducedConversions"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       fixedGridRhoFastjetCentralNeutralLabel = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                                       fixedGridRhoFastjetAllLabel = cms.InputTag("fixedGridRhoFastjetAll"),
                                       prescales = cms.InputTag("patTrigger"),
                                       bits = cms.InputTag("TriggerResults","","HLT"),
                                       #mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),

                                       #eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                       #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                       #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                       #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")

                                       )

#process.tauMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                     src = cms.InputTag("slimmedTaus"),
#                                     matched = cms.InputTag("prunedGenParticles"),
#                                     distMin = cms.double(0.15),
#                                     matchPDGId = cms.vint32()
#                                     )

#process.tauMCMatch = cms.EDProducer("MCMatcher",
#    				src     = cms.InputTag("slimmedTaus"),
#    				matched = cms.InputTag("prunedGenParticles"),
#                   mcPdgId     = cms.vint32(),
#                   checkCharge = cms.bool(False),
#                   mcStatus = cms.vint32(),
#                   maxDeltaR = cms.double(0.15),
#                   maxDPtRel = cms.double(3.0),
#                   resolveAmbiguities = cms.bool(True),
#                   resolveByMatchQuality = cms.bool(False),
#)


#
# Example for a configuration of the MC match
# for taus (cuts are NOT tuned)
# (using old values from TQAF, january 2008)
#
process.tauMatch = cms.EDProducer("MCMatcher",
                          src         = cms.InputTag("slimmedTaus"),          # RECO objects to match
                          matched     = cms.InputTag("prunedGenParticles"),              # mc-truth particle collection
                          mcPdgId     = cms.vint32(15),                            # one or more PDG ID (15 = tau); absolute values (see below)
                          checkCharge = cms.bool(True),                            # True = require RECO and MC objects to have the same charge
                          mcStatus    = cms.vint32(2),                             # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                          # NOTE that Taus can only be status 3 or 2, never 1!
                          maxDeltaR   = cms.double(999.9),                         # Minimum deltaR for the match.     By default any deltaR is allowed (why??)
                          maxDPtRel   = cms.double(999.9),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                          resolveAmbiguities    = cms.bool(True),                  # Forbid two RECO objects to match to the same GEN object
                          resolveByMatchQuality = cms.bool(False),                 # False = just match input in order; True = pick lowest deltaR pair first
                          )

process.tauGenJetMatch = cms.EDProducer("GenJetMatcher",             # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                src         = cms.InputTag("slimmedTaus"),          # RECO jets (any View<Jet> is ok)
                                matched     = cms.InputTag("tauGenJets"),  # GEN jets  (must be GenJetCollection) "tauGenJetsSelectorAllHadrons"
                                mcPdgId     = cms.vint32(),                              # n/a
                                mcStatus    = cms.vint32(),                              # n/a
                                checkCharge = cms.bool(False),                           # n/a
                                maxDeltaR   = cms.double(0.1),                           # Minimum deltaR for the match
                                maxDPtRel   = cms.double(3.0),                           # Minimum deltaPt/Pt for the match
                                resolveAmbiguities    = cms.bool(True),                  # Forbid two RECO objects to match to the same GEN object
                                resolveByMatchQuality = cms.bool(False),                 # False = just match input in order; True = pick lowest deltaR pair first
                                )
process.tauGenJets = cms.EDProducer(
                            "TauGenJetProducer",
                            prunedGenParticles =  cms.InputTag('prunedGenParticles'),
                            includeNeutrinos = cms.bool( False ),
                            verbose = cms.untracked.bool( False )
                            )

process.patMCTruth_Tau =  cms.Sequence ( process.tauMatch+
                                process.tauGenJets*
                                process.tauGenJetMatch )


process.electronMatch = cms.EDProducer("MCMatcher",
                                  src         = cms.InputTag("slimmedElectrons"),          # RECO objects to match
                                  matched     = cms.InputTag("prunedGenParticles"),              # mc-truth particle collection
                                  mcPdgId     = cms.vint32(),                            # one or more PDG ID (15 = tau); absolute values (see below)
                                  checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
                                  mcStatus    = cms.vint32(),                             # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                  # NOTE that Taus can only be status 3 or 2, never 1!
                                  maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed (why??)
                                  maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                  resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                  resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                  )

process.muonMatch = cms.EDProducer("MCMatcher",
                                   src         = cms.InputTag("slimmedMuons"),              # RECO objects to match
                                   matched     = cms.InputTag("prunedGenParticles"),          # mc-truth particle collection
                                   mcPdgId     = cms.vint32(),                          # one or more PDG ID (15 = tau); absolute values (see below)
                                   checkCharge = cms.bool(False),                       # True = require RECO and MC objects to have the same charge
                                   mcStatus    = cms.vint32(),                          # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   # NOTE that Taus can only be status 3 or 2, never 1!
                                   maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed
                                   maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                   resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                   )

## reco-generator(parton) matching for jets
process.jetPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, -1, -2, -3, -4, -5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.jetPartonMatch2 = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, -1, -2, -3, -4, -5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(2),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(3.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.jetPartonMatchJets = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("ak5GenJets"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(3.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.electronPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                             src = cms.InputTag("slimmedElectrons"),      # RECO objects to match
                                             matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                             mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                             mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                             checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                             maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                             maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                             resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                             resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                             )

process.muonPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                         src = cms.InputTag("slimmedMuons"),      # RECO objects to match
                                         matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                         mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                         mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                         checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                         maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                         maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                         resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                         resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                         )

process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
                                                  , filterParams =cms.PSet(
                                                                           minNdof = cms.double( 4. )
                                                                           , maxZ    = cms.double( 24. )
                                                                           , maxRho  = cms.double( 2. ) )
                                                  , filter       = cms.bool( False ) # use only as producer
                                                  , src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
)

if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
#   *process.tauMatch
#   *process.muonMatch
#   *process.electronMatch
#   *process.jetPartonMatch
#   *process.jetPartonMatch2
#   *process.electronPartonMatch
#   *process.muonPartonMatch
    *process.egmGsfElectronIDSequence
    *process.FakeElectrons
        )
else:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '/user/ikhvastu/CMSSW_8_0_5/src/Majorana/PatAnalyzer/test/JSON/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON_NoL1T.json').getVLuminosityBlockRange()
    process.p = cms.Path(
    process.goodOfflinePrimaryVertices
    #primaryVertexFilter
    #*process.eeBadScFilter
    #*process.HBHENoiseFilterResultProducer #produces HBHE baseline bools
    #*process.ApplyBaselineHBHENoiseFilter  #reject events based
    #*process.BeamHaloId
    #*process.CSCTightHalo2015Filter
    *process.egmGsfElectronIDSequence
    *process.FakeElectrons
)



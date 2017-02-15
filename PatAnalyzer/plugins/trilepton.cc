#include "trilepton.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/CandAlgos/interface/CandMatcher.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthPairSelector.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "PhysicsTools/CandUtils/interface/CandMatcherNew.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

//btagging
#include "Majorana/PatAnalyzer/interface/BTagCalibrationStandalone.h"

#include <iostream>

trilepton::trilepton(const edm::ParameterSet & iConfig) :
  genparticleToken(                      consumes<reco::GenParticleCollection>(   iConfig.getParameter<edm::InputTag>("genPartsLabel"))),
  electronMvaIdMapToken(                 consumes<edm::ValueMap<float>>(          iConfig.getParameter<edm::InputTag>("electronMvaIdMap"))),
  mvaValuesMapToken_HZZ_(                consumes<edm::ValueMap<float>>(          iConfig.getParameter<edm::InputTag>("mvaValuesMap_HZZ"))),
  electronMvaIdMap80Token(               consumes<edm::ValueMap<bool>>(           iConfig.getParameter<edm::InputTag>("electronMvaId80Map"))),
  electronMvaIdMap90Token(               consumes<edm::ValueMap<bool>>(           iConfig.getParameter<edm::InputTag>("electronMvaId90Map"))),
  electronCutBasedIdMapTightToken(       consumes<edm::ValueMap<bool>>(           iConfig.getParameter<edm::InputTag>("electronCutBasedIdTightMap"))),
  electronCutBasedIdMapMediumToken(      consumes<edm::ValueMap<bool>>(           iConfig.getParameter<edm::InputTag>("electronCutBasedIdMediumMap"))),
  pdfvariablesToken(                     consumes<GenEventInfoProduct>(           iConfig.getParameter<edm::InputTag>("pdfvariablesLabel"))),
  IT_beamspot(                           consumes<reco::BeamSpot>(                iConfig.getParameter<edm::InputTag>("BeamSpotLabel"))),
  PileUpToken(                           consumes<vector<PileupSummaryInfo>>(     iConfig.getParameter<edm::InputTag>("slimmedAddPileupInfoLabel"))),
  goodOfflinePrimaryVerticesToken(       consumes<std::vector<Vertex>>(           iConfig.getParameter<edm::InputTag>("goodOfflinePrimaryVerticesLabel"))),
  packedPFCandidatesToken(               consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidatesLabel"))),
  IT_muon(                               consumes<std::vector<pat::Muon>>(        iConfig.getParameter<edm::InputTag>("MuonLabel"))),
  IT_electron(                           consumes<edm::View<pat::Electron>>(      iConfig.getParameter<edm::InputTag>("ElectronLabel"))),
  IT_jet(                                consumes<std::vector<pat::Jet>>(         iConfig.getParameter<edm::InputTag>("JetLabel"))),
  IT_pfmet(                              consumes<std::vector<pat::MET>>(         iConfig.getParameter<edm::InputTag>("METLabel"))),
  reducedEgammaToken(                    consumes<std::vector<reco::Conversion>>( iConfig.getParameter<edm::InputTag>("reducedEgammaLabel"))),
  IT_tau(                                consumes<std::vector<pat::Tau>>(         iConfig.getParameter<edm::InputTag>("TauLabel"))),
  fixedGridRhoFastjetCentralNeutralToken(consumes<double>(                        iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetCentralNeutralLabel"))),
  fixedGridRhoFastjetAllToken(           consumes<double>(                        iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetAllLabel"))),
  triggerResultsHLTToken(                consumes<edm::TriggerResults>(           iConfig.getParameter<edm::InputTag>("triggerResultsHLT"))),
  triggerResultsRECOToken(               consumes<edm::TriggerResults>(           iConfig.getParameter<edm::InputTag>("triggerResultsRECO"))),
  triggerPrescalesToken(                 consumes<pat::PackedTriggerPrescales>(   iConfig.getParameter<edm::InputTag>("prescales"))),
  IT_externalLHEProducer(                consumes<LHEEventProduct>(               iConfig.getParameter<edm::InputTag>("exernalLHEPLabel"))),
  isData(                                                                         iConfig.getUntrackedParameter<bool>("isData")),
  treeForFakeRate(                                                                iConfig.getUntrackedParameter<bool>("treeForFakeRate")),
  singleLep(                                                                      iConfig.getUntrackedParameter<bool>("singleLep")),
  SampleName(                                                                     iConfig.getUntrackedParameter<std::string>("SampleName")),
  _relIsoCutE(999.),//0.15
  _relIsoCutMu(999.),//0.15
  _relIsoCutEloose(999.), //0.5
  _relIsoCutMuloose(999.), //0.5
  _chargeConsistency(true),
  _minPt0(5),
  _minPt1(10.),
  _tightD0Mu(0.05),
  _tightD0E(0.05),
  _looseD0Mu(0.05), // 0.05 for sync   // Not sure if this was some temporary change?
  _looseD0E(0.05), // 0.05 for sync
  //_looseD0Mu(0.2),
  //_looseD0E(9999999.),
  //_jetPtCut(40.),
  _jetPtCut(20.),
  _jetEtaCut(2.4),
  _tauPt(20),
  _tauEta(2.3),
  _regression(false)
{
    // eventually move this to config level
    // just a bunch of lepton triggers, might need a closer look for cleanup or additions
    triggersToSave = {"HLT_TripleMu_12_10_5", "HLT_DiMu9_Ele9_CaloIdL_TrackIdL", "HLT_Mu8_DiEle12_CaloIdL_TrackIdL", "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",  // 3l
                      "HLT_TripleMu_5_3_3_DZ_Mass3p8",                                                                                                         // 3l (with mass requirement)
                      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",               // mumu
                      "HLT_Mu30_TkMu11", "HLT_Mu40_TkMu11",
                      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",    "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL",                  // mumu prescaled
                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",                                                                                     // mue
                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",                                     // mue prescaled
                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                      "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",                                                                                        // mue partially off
                      "HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf",        // ee
                      "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
                      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",        // ee prescaled
                      "HLT_Mu50", "HLT_TkMu50", "HLT_IsoMu27", "HLT_IsoTkMu27",                                                                                // mu
                      "HLT_IsoMu24", "HLT_IsoTkMu24", "HLT_IsoMu24_eta2p1", "HLT_IsoTkMu24_eta2p1",
                      "HLT_IsoMu22_eta2p1", "HLT_IsoTkMu22_eta2p1",
                      "HLT_IsoTkMu22", "HLT_Mu45_eta2p1",                                                                                                      // mu partially off
                      "HLT_Ele32_WPTight_Gsf", "HLT_Ele32_eta2p1_WPTight_Gsf", "HLT_Ele30_WPTight_Gsf", "HLT_Ele30_eta2p1_WPTight_Gsf",                        // e
                      "HLT_Ele27_WPTight_Gsf", "HLT_Ele27_eta2p1_WPTight_Gsf", "HLT_Ele27_eta2p1_WPLoose_Gsf", "HLT_Ele25_eta2p1_WPTight_Gsf",
                      "HLT_Ele25_WPTight_Gsf",                                                                                                                 // e partially off
                      "HLT_Ele8_CaloIdM_TrackIdM_PFJet30","HLT_Ele12_CaloIdM_TrackIdM_PFJet30",                                                                // e + PF jet (for fake rate measurment)
                      "HLT_Ele8_CaloIdM_TrackIdM_IsoVL_PFJet30","HLT_Ele12_CaloIdM_TrackIdM_IsoVL_PFJet30",
                      "HLT_Mu8_TrkIsoVVL","HLT_Mu17_TrkIsoVVL","HLT_Mu24_TrkIsoVVL","HLT_Mu34_TrkIsoVVL",                                                      // FR
                      "HLT_Ele18_CaloIdM_TrackIdM_PFJet30",
                      "HLT_Ele23_CaloIdM_TrackIdM_PFJet30","HLT_Ele33_CaloIdM_TrackIdM_PFJet30",
                      "HLT_Mu3_PFJet40",
                      "HLT_Mu8","HLT_Mu17","HLT_Mu24","HLT_Mu34",
                      "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30",
                      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30",
                      "HLT_Ele12_CaloIdL_TrkIdL_IsoVL","HLT_Ele17_CaloIdL_TrackIdL_IsoVL",
                      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL"};
    filtersToSave  = {"Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_eeBadScFilter", "Flag_globalTightHalo2016Filter"};
}


void trilepton::beginJob()
{
    // Read in effective areas from text files, to be used for lepton isolation calculations
    tools::readEffAreas(edm::FileInPath("Majorana/PatAnalyzer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt").fullPath(), 11);
    tools::readEffAreas(edm::FileInPath("Majorana/PatAnalyzer/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt").fullPath(), 13);
    tools::initMultiIsoConstants(); // that's the problem when using namespaces instead of classes, you need to have separate init functions

    if(singleLep)            outputTree = fs->make<TTree>("singleLepTree", "singleLepTree");
    else if(treeForFakeRate) outputTree = fs->make<TTree>("fakeRateTree",  "fakeRateTree");
    else                     outputTree = fs->make<TTree>("trileptonTree","trileptonTree");
    Nvtx = fs->make<TH1F>("N_{vtx}", "Number of vertices;N_{vtx};events / 1"  ,    40, 0., 40.);
    
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);

   
    outputTree->Branch("_eventNb",                    &_eventNb,                    "_eventNb/l");
    outputTree->Branch("_runNb",                      &_runNb,                      "_runNb/l");
    outputTree->Branch("_lumiBlock",                  &_lumiBlock,                  "_lumiBlock/l");

    outputTree->Branch("_nLeptons",                   &_nLeptons,                   "_nLeptons/I");
    outputTree->Branch("_lPt",                        &_lPt,                        "_lPt[_nLeptons]/D");
    outputTree->Branch("_lEta",                       &_lEta,                       "_lEta[_nLeptons]/D");
    outputTree->Branch("_lPhi",                       &_lPhi,                       "_lPhi[_nLeptons]/D");
    outputTree->Branch("_lE",                         &_lE,                         "_lE[_nLeptons]/D");
    outputTree->Branch("_lPtmc",                      &_lPtmc,                      "_lPtmc[_nLeptons]/D");
    outputTree->Branch("_lEtamc",                     &_lEtamc,                     "_lEtamc[_nLeptons]/D");
    outputTree->Branch("_lPhimc",                     &_lPhimc,                     "_lPhimc[_nLeptons]/D");
    outputTree->Branch("_lEmc",                       &_lEmc,                       "_lEmc[_nLeptons]/D");
    outputTree->Branch("_lpdgmc",                     &_lpdgmc,                     "_lpdgmc[_nLeptons]/D");
    outputTree->Branch("_lchargemc",                  &_lchargemc,                  "_lchargemc[_nLeptons]/D");

    outputTree->Branch("_nuPtmc",                     &_nuPtmc,                     "_nuPtmc[_nLeptons]/D");
    outputTree->Branch("_nuEtamc",                    &_nuEtamc,                    "_nuEtamc[_nLeptons]/D");
    outputTree->Branch("_nuPhimc",                    &_nuPhimc,                    "_nuPhimc[_nLeptons]/D");
    outputTree->Branch("_nuEmc",                      &_nuEmc,                      "_nuEmc[_nLeptons]/D");

    outputTree->Branch("_mtmc",                       &_mtmc,                       "_mtmc[_nLeptons]/D");

    outputTree->Branch("_nZboson",                    &_nZboson,                    "_nZboson/I");

    outputTree->Branch("_nEle",                       &_nEle,                       "_nEle/I");
    outputTree->Branch("_nMu",                        &_nMu,                        "_nMu/I");
    outputTree->Branch("_nTau",                       &_nTau,                       "_nTau/I");

    outputTree->Branch("_flavors",                    &_flavors,                    "_flavors[_nLeptons]/I");
    outputTree->Branch("_charges",                    &_charges,                    "_charges[_nLeptons]/D");
    outputTree->Branch("_isolation",                  &_isolation,                  "_isolation[_nLeptons]/D");
    outputTree->Branch("_isolation_absolute",         &_isolation_absolute,         "_isolation_absolute[_nLeptons]/D");
    outputTree->Branch("_miniisolation",              &_miniisolation,              "_miniisolation[_nLeptons]/D");
    outputTree->Branch("_miniisolationCharged",       &_miniisolationCharged,       "_miniisolationCharged[_nLeptons]/D");
    outputTree->Branch("_ptrel",                      &_ptrel,                      "_ptrel[_nLeptons]/D");
    outputTree->Branch("_ptratio",                    &_ptratio,                    "_ptratio[_nLeptons]/D");
    outputTree->Branch("_muonSegmentComp",            &_muonSegmentComp,            "_muonSegmentComp[_nLeptons]/D");

    outputTree->Branch("_mll",                        &_mll,                        "_mll[3]/D");

    outputTree->Branch("_origin",                     &_origin,                     "_origin[_nLeptons]/I");
    outputTree->Branch("_originReduced",              &_originReduced,              "_originReduced[_nLeptons]/I");
    outputTree->Branch("_originPhot",                &_originPhot,                "_originPhot[_nLeptons]/I");
    outputTree->Branch("_originDetailed",                &_originDetailed,                "_originDetailed[_nLeptons]/I");

    outputTree->Branch("_isPromptFinalState",         &_isPromptFinalState,         "_isPromptFinalState[_nLeptons]/O");
    outputTree->Branch("_fromHardProcessFinalState",  &_fromHardProcessFinalState,  "_fromHardProcessFinalState[_nLeptons]/O");

    outputTree->Branch("_PVchi2",                     &_PVchi2,                     "_PVchi2/D");
    outputTree->Branch("_PVerr",                      &_PVerr,                      "_PVerr[3]/D");
    outputTree->Branch("_ipPV",                       &_ipPV,                       "_ipPV[_nLeptons]/D");
    outputTree->Branch("_ipPVerr",                    &_ipPVerr,                    "_ipPVerr[_nLeptons]/D");
    outputTree->Branch("_ipZPV",                      &_ipZPV,                      "_ipZPV[_nLeptons]/D");
    outputTree->Branch("_ipZPVerr",                   &_ipZPVerr,                   "_ipZPVerr[_nLeptons]/D");

    outputTree->Branch("_ipPVmc",                     &_ipPVmc,                     "_ipPVmc[_nLeptons]/D");

    outputTree->Branch("_3dIP",                       &_3dIP,                       "_3dIP[_nLeptons]/D");
    outputTree->Branch("_3dIPerr",                    &_3dIPerr,                    "_3dIPerr[_nLeptons]/D");
    outputTree->Branch("_3dIPsig",                    &_3dIPsig,                    "_3dIPsig[_nLeptons]/D");


    outputTree->Branch("_mt",                         &_mt,                         "_mt[_nLeptons]/D");
    outputTree->Branch("_isloose",                    &_isloose,                    "_isloose[_nLeptons]/O");
    outputTree->Branch("_ismedium",                   &_ismedium,                   "_ismedium[_nLeptons]/O");
    outputTree->Branch("_istight",                    &_istight,                    "_istight[_nLeptons]/O");

    outputTree->Branch("_trigEmulator",               &_trigEmulator,               "_trigEmulator[_nLeptons]/O");
    outputTree->Branch("_isotrigEmulator",            &_isotrigEmulator,            "_isotrigEmulator[_nLeptons]/O");

    outputTree->Branch("_closeJetPtAll",              &_closeJetPtAll,              "_closeJetPtAll[_nLeptons]/D");
    outputTree->Branch("_closeJetAngAll",             &_closeJetAngAll,             "_closeJetAngAll[_nLeptons]/D");

    outputTree->Branch("_trackSelectionMultiplicity", &_trackSelectionMultiplicity, "_trackSelectionMultiplicity[_nLeptons]/I");
    outputTree->Branch("_closeJetCSVAll",             &_closeJetCSVAll,             "_closeJetCSVAll[_nLeptons]/D");

    outputTree->Branch("_chargeConst",                &_chargeConst,                "_chargeConst[_nLeptons]/O");
    outputTree->Branch("_hitsNumber",                 &_hitsNumber,                 "_hitsNumber[_nLeptons]/I");
    outputTree->Branch("_vtxFitConversion",           &_vtxFitConversion,           "_vtxFitConversion[_nLeptons]/O");

    outputTree->Branch("_mvaValue",                   &_mvaValue,                   "_mvaValue[_nLeptons]/D");
    outputTree->Branch("_mvaValue_HZZ",               &_mvaValue_HZZ,               "_mvaValue_HZZ[_nLeptons]/D");

    outputTree->Branch("_passedCutBasedIdTight",      &_passedCutBasedIdTight,      "_passedCutBasedIdTight[_nLeptons]/O");
    outputTree->Branch("_passedCutBasedIdMedium",     &_passedCutBasedIdMedium,     "_passedCutBasedIdMedium[_nLeptons]/O");
    outputTree->Branch("_passedMVA80",                &_passedMVA80,                "_passedMVA80[_nLeptons]/O");
    outputTree->Branch("_passedMVA90",                &_passedMVA90,                "_passedMVA90[_nLeptons]/O");
    outputTree->Branch("_passedMVA_SUSY",             &_passedMVA_SUSY,             "_passedMVA_SUSY[_nLeptons][3]/O"); // I really don't like this SUSY-style of saving things in the tree, but lets keep it for now to not break the current workflow
 
    outputTree->Branch("_n_PV",                       &_n_PV,                       "_n_PV/I");
    outputTree->Branch("_n_MCTruth_PV",               &_n_MCTruth_PV,               "_n_MCTruth_PV/D");
    outputTree->Branch("_met",                        &_met,                        "_met/D");
    outputTree->Branch("_met_phi",                    &_met_phi,                    "_met_phi/D");
    outputTree->Branch("HT",                          &HT,                          "HT/D");
 
    outputTree->Branch("_genmet",                     &_genmet,                     "_genmet/D");
    outputTree->Branch("_genmet_phi",                 &_genmet_phi,                 "_genmet_phi/D");
    outputTree->Branch("_genqpt",                     &_genqpt,                     "_genqpt/D");
    outputTree->Branch("_genqpt20",                   &_genqpt20,                   "_genqpt20/D");

    outputTree->Branch("_mompt",                      &_mompt,                      "_mompt[_nLeptons]/D");
    outputTree->Branch("_momphi",                     &_momphi,                     "_momphi[_nLeptons]/D");
    outputTree->Branch("_mometa",                     &_mometa,                     "_mometa[_nLeptons]/D");
    outputTree->Branch("_mompdg",                     &_mompdg,                     "_mompdg[_nLeptons]/I");

    outputTree->Branch("_n_bJets",                    &_n_bJets,                    "_n_bJets/I");
    outputTree->Branch("_n_Jets",                     &_n_Jets,                     "_n_Jets/I");
    outputTree->Branch("_bTagged",                    &_bTagged,                    "_bTagged[_n_Jets]/O");
    outputTree->Branch("_jetEta",                     &_jetEta,                     "_jetEta[_n_Jets]/D");
    outputTree->Branch("_jetPhi",                     &_jetPhi,                     "_jetPhi[_n_Jets]/D");
    outputTree->Branch("_jetPt",                      &_jetPt,                      "_jetPt[_n_Jets]/D");
    outputTree->Branch("_jetFlavour",                 &_jetFlavour,                 "_jetFlavour[_n_Jets]/I");
    outputTree->Branch("_jetE",                       &_jetE,                       "_jetE[_n_Jets]/D");
    outputTree->Branch("_csv",                        &_csv,                        "_csv[_n_Jets]/D");
    outputTree->Branch("_clean",                      &_clean,                      "_clean[_n_Jets]/I");

    // JEC
    outputTree->Branch("_jecUnc",                     &_jecUnc,                     "_jecUnc[_n_Jets]/D");
    outputTree->Branch("_jetPtUp",                    &_jetPtUp,                    "_jetPtUp[_n_Jets]/D");
    outputTree->Branch("_jetPtDown",                  &_jetPtDown,                  "_jetPtDown[_n_Jets]/D");

    outputTree->Branch("_matchedjetPt",               &_matchedjetPt,               "_matchedjetPt[_n_Jets]/D");
    outputTree->Branch("_matchedjetEta",              &_matchedjetEta,              "_matchedjetEta[_n_Jets]/D");
    outputTree->Branch("_matchedjetPhi",              &_matchedjetPhi,              "_matchedjetPhi[_n_Jets]/D");
    outputTree->Branch("_matchedjetE",                &_matchedjetE,                "_matchedjetE[_n_Jets]/D");
    outputTree->Branch("_matchedjetM",                &_matchedjetM,                "_matchedjetM[_n_Jets]/D");
    outputTree->Branch("_matchGjet",                  &_matchGjet,                  "_matchGjet[_n_Jets]/O");

    // b-tag SF
    outputTree->Branch("_btagSF",                     &_btagSF,                     "_btagSF[19][_n_Jets]/D");


    // theoretical
    outputTree->Branch("_weight",                     &_weight,                     "_weight/D");
    outputTree->Branch("_LHEweight",                  &_LHEweight,                  "_LHEweight[111]/D");
    outputTree->Branch("_LHEweightID",                &_LHEweightID,                "_LHEweightID[111]/D");

    outputTree->Branch("_nMajorana",                  &_nMajorana,                  "_nMajorana/I");
    outputTree->Branch("_findMatched",                &_findMatched,                "_findMatched/I");
        

  

    // trigger
    for(TString triggerName : triggersToSave){
      outputTree->Branch(triggerName, &triggerFlags[triggerName], triggerName + "/O");
      outputTree->Branch(triggerName + "_prescale", &triggerPrescales[triggerName], triggerName + "_prescale/I");
    }

    // MET filters
    for(TString filterName : filtersToSave){
      outputTree->Branch(filterName, &triggerFlags[filterName], filterName + "/O");
    }


    //Gen Level outputTrees////////////////////////////////////////////////////
    if(not isData){
      outputTree->Branch("_gen_nL",       &_gen_nL,       "_gen_nL/I"); 
      outputTree->Branch("_gen_lPt",      &_gen_lPt,      "_gen_lPt[_gen_nL]/D");
      outputTree->Branch("_gen_lE",       &_gen_lE,       "_gen_lE[_gen_nL]/D");
      outputTree->Branch("_gen_lEta",     &_gen_lEta,     "_gen_lEta[_gen_nL]/D");
      outputTree->Branch("_gen_lPhi",     &_gen_lPhi,     "_gen_lPhi[_gen_nL]/D");
      outputTree->Branch("_gen_lmompdg",  &_gen_lmompdg,  "_gen_lmompdg[_gen_nL]/I");
      outputTree->Branch("_gen_charges",  &_gen_charges,  "_gen_charges[_gen_nL]/D");
      outputTree->Branch("_gen_flavors",  &_gen_flavors,  "_gen_flavors[_gen_nL]/I");
      //Added neutrinos  ////////////////////////////////////////////////////
      outputTree->Branch("_gen_nNu",      &_gen_nNu,      "_gen_nNu/I");
      outputTree->Branch("_gen_nuPt",     &_gen_nuPt,     "_gen_nuPt[_gen_nNu]/D");
      outputTree->Branch("_gen_nuE",      &_gen_nuE,      "_gen_nuE[_gen_nNu]/D");
      outputTree->Branch("_gen_nuEta",    &_gen_nuEta,    "_gen_nuEta[_gen_nNu]/D");
      outputTree->Branch("_gen_nuPhi",    &_gen_nuPhi,    "_gen_nuPhi[_gen_nNu]/D");
      outputTree->Branch("_gen_numompdg", &_gen_numompdg, "_gen_numompdg[_gen_nNu]/I");
      //End of neutrino Edit////////////////////////////////////////////////////////
      //Added neutrinos Majorana  ////////////////////////////////////////////////////
      outputTree->Branch("_gen_nMajo",    &_gen_nMajo,    "_gen_nMajo/I");
      outputTree->Branch("_gen_majoPt",   &_gen_majoPt,   "_gen_majoPt[_gen_nMajo]/D");
      outputTree->Branch("_gen_majoE",    &_gen_majoE,    "_gen_majoE[_gen_nMajo]/D");
      outputTree->Branch("_gen_majoEta",  &_gen_majoEta,  "_gen_majoEta[_gen_nMajo]/D");
      outputTree->Branch("_gen_majoPhi",  &_gen_majoPhi,  "_gen_majoPhi[_gen_nMajo]/D");
      //Added W  ////////////////////////////////////////////////////
      outputTree->Branch("_gen_nW",       &_gen_nW,       "_gen_nW/I");
      outputTree->Branch("_gen_wPt",      &_gen_wPt,      "_gen_wPt[_gen_nW]/D");
      outputTree->Branch("_gen_wE",       &_gen_wE,       "_gen_wE[_gen_nW]/D");
      outputTree->Branch("_gen_wEta",     &_gen_wEta,     "_gen_wEta[_gen_nW]/D");
      outputTree->Branch("_gen_wPhi",     &_gen_wPhi,     "_gen_wPhi[_gen_nW]/D");
      outputTree->Branch("_gen_wmompdg",  &_gen_wmompdg,  "_gen_wmompdg[_gen_nW]/I");

    }

    
    GPM = GenParticleManager();

    fMetCorrector = new OnTheFlyCorrections("Majorana/PatAnalyzer/data/", "Summer16_23Sep2016", isData);
    if (isData) _corrLevel = "L2L3Residual";
    else        _corrLevel = "L3Absolute";

    jecUnc = new JetCorrectionUncertainty(edm::FileInPath("Majorana/PatAnalyzer/data/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt").fullPath());
    
    for(TString wp : {"VT","T","M","L","VL"}){
      leptonWorkingPoints["multiIsolation" + wp] = new std::vector<bool>();
      outputTree->Branch("multiIsolation" + wp, "vector<bool>", leptonWorkingPoints["multiIsolation" + wp]);
    }

    for(TString wp : {"VT","T","M","L","VL","RelIso04","MiniIso04"}){
      leptonConeCorrectedPt[wp] = new std::vector<double>();
      outputTree->Branch("leptonConeCorrectedPt" + wp, "vector<double>", leptonConeCorrectedPt[wp]);
    }
        
    _nEventsTotal = 0;
    _nEventsFiltered = 0;
    _nEventsTotalCounted = 0;
    
    firstEvent_ = true;
}

void trilepton::endJob() {
    for(TString wp : {"VT","T","M","L","VL"})                        delete leptonWorkingPoints["multiIsolation" + wp];
    for(TString wp : {"VT","T","M","L","VL","RelIso04","MiniIso04"}) delete leptonConeCorrectedPt[wp];

    std::cout<<_nEventsTotal<<std::endl;
    std::cout<<_nEventsFiltered<<std::endl;
    
    delete fMetCorrector;
    
    
}


void trilepton::getTriggerResults(const edm::Event& iEvent, bool isHLT, edm::EDGetTokenT<edm::TriggerResults> token, std::vector<TString> toSave){
  edm::Handle<edm::TriggerResults> trigResults;       iEvent.getByToken(token,                 trigResults);
  edm::Handle<pat::PackedTriggerPrescales> prescales; iEvent.getByToken(triggerPrescalesToken, prescales);

  if(trigResults.failedToGet()){
    std::cout << "WARNING: no trigger results for " << (isHLT? "HLT" : "RECO") << "!" << std::endl;
    return;
  }

  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*trigResults);

  // Get full trigger list, and remember the indices of triggers we need to save
  if(firstEvent_){
    for(TString triggerName : toSave) triggerIndices[triggerName] = -1;
    std::cout << "Available triggers:" << std::endl;
    for (unsigned int i = 0; i < trigResults->size(); ++i){
      std::cout << "  " << triggerNames.triggerName(i) << std::endl;
      for(TString triggerName : toSave){
        if(TString(triggerNames.triggerName(i)).Contains(triggerName + "_v")){
          triggerIndices[triggerName] = i;
          std::cout << "     --> index: " << i << std::endl; 
        }
      }
    }
    std::cout << std::endl;
  }

  // Save our triggers/flags
  for(TString triggerName : toSave){
    if(triggerIndices[triggerName] == -1){
      triggerFlags[triggerName] = false;
      triggerPrescales[triggerName] = -1;
      continue;
    }
    triggerFlags[triggerName] = trigResults->accept(triggerIndices[triggerName]);
    if(isHLT){
      triggerPrescales[triggerName] = prescales->getPrescaleForIndex(triggerIndices[triggerName]);
    }
  }
}



void trilepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{
    for(TString wp : {"VT","T","M","L","VL"})                        leptonWorkingPoints["multiIsolation" + wp]->clear();
    for(TString wp : {"VT","T","M","L","VL","RelIso04","MiniIso04"}) leptonConeCorrectedPt[wp]->clear();

    _nZboson = 0;

    if(not isData) {
        //******************************************************************************************************************
        // Gen level particles                  ****************************************************************************
        //******************************************************************************************************************
        //iEvent.getByLabel("packedGenParticles", TheGenParticles);
        iEvent.getByToken(genparticleToken, TheGenParticles);
        std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW, vGenMajorana;
        if( TheGenParticles.isValid() ) {
            GPM.SetCollection(TheGenParticles);
            GPM.Classify();
            vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),1);
            vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),1);
            vGenNPMuons = GPM.filterByStatus(GPM.getNonPromptMuons(),1);
            vGenNPElectrons = GPM.filterByStatus(GPM.getNonPromptElectrons(),1);
          //  std::cout<<"********************************************************************"<<std::endl;
          //              std::cout<<"********************************************************************"<<std::endl;

            TLorentzVector Gen0;
            Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
            _genqpt = 0;
            _genqpt20 = 0;
            //cout << "Particle ids in the event" << endl;
            for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ )
            {
                int id = std::abs(p->pdgId());
                
                if(p->status() == 1){
                    const GenParticle *mom = GPM.getMother(&*p);
                    //cout << id << " " << p->pt() << " " << std::abs(mom->pdgId()) << endl;
                }
                
                if(id == 9900012) _nMajorana++;
                if(id == 23)      _nZboson++;
                if ( (id == 12 || id == 14 || id == 16 ) && (p->status() == 1) ) {
                    TLorentzVector Gen;
                    Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );
                    Gen0 += Gen;
                }
                if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 21 || id == 22 ) && (p->status() == 23)){
                    _genqpt += p->pt();
                    if(p->pt() > 20)
                        _genqpt20 += p->pt();
                }
            }
            cout << endl;
            if (Gen0.E()!=0) {
              _genmet = Gen0.Pt();
              _genmet_phi = Gen0.Phi();
            } else {
              _genmet = 0;
              _genmet_phi = 0;
            }

            for(int i = 0; i <20 ; ++i) _gen_lmompdg[i] = 0;
            int nL = 0,  nNu = 0, nMajo=0, nw=0;;

            for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; ++p ){
              //if(nL > 10 || nPh > 10 || nNu > 10) break;
              if(p->status() == 0) continue;  //status 0 contains no useful information and should be skipped
              unsigned id = std::abs(p->pdgId()); 
              const reco::GenStatusFlags StatusFlags;
              //MCTruthHelper::fillGenStatusFlags(*p,StatusFlags);
              //if( id == 11 || id ==13 || id == 15){
               //if(p->status() == 1 || (id == 15 && p->status() == 2)){        //Tau's always decay! 
              if(id == 11 || id == 13){
                if(p->status() == 1){
                  if (nL > 29) continue;
                  _gen_lPt[nL] = p->pt();
                  _gen_lEta[nL] = p->eta();
                  _gen_lPhi[nL] = p->phi();
                  _gen_lE[nL] = p->energy();
                  const GenParticle *mom = GPM.getMother(&*p);        //Use getmother or getmotherparton here?
                  while( fabs(mom->pdgId()) == id && GPM.getMother(mom) != nullptr) mom = GPM.getMother(mom);
                  if(mom != nullptr) _gen_lmompdg[nL] = mom->pdgId();
                  else               _gen_lmompdg[nL] = 0;  //Particles for which mompdg is 0 have no mother.
                  if(id == 11)       _gen_flavors[nL] = 0;
                  else if(id ==13)   _gen_flavors[nL] = 1;
                //else if(id ==15) {_gen_flavors[nL] = 2;}
                  _gen_charges[nL] = p->charge();

                  _gen_lPt[nL] = p->pt();
                  _gen_lEta[nL] = p->eta();
                  _gen_lPhi[nL] = p->phi();
                  _gen_lE[nL] = p->energy();
                  ++nL;
                }
              } else if(id == 12 || id == 14 || id == 16){
                if(p->status() == 1){
                  if (nNu > 29) continue;
                  _gen_nuPt[nNu] = p->pt();
                  _gen_nuE[nNu] = p->energy();
                  _gen_nuEta[nNu] = p->eta();
                  _gen_nuPhi[nNu] = p->phi();
                  
                  const GenParticle *mom = GPM.getMother(&*p);        //Use getmother or getmotherparton here?
                  while( fabs(mom->pdgId()) == id && GPM.getMother(mom) != nullptr) mom = GPM.getMother(mom);
                  if(mom != nullptr) _gen_numompdg[nNu] = mom->pdgId();
                  else               _gen_numompdg[nNu] = 0;  //Particles for which mompdg is 0 have no mother.
                  
                  ++nNu;
                }
              } else if(  id == 24 ){
                if(p->status() == 22 ||p->status() == 52){
                  if (nw > 29) continue;
                  _gen_wPt[nw] = p->pt();
                  _gen_wE[nw] = p->energy();
                  _gen_wEta[nw] = p->eta();
                  _gen_wPhi[nw] = p->phi();
                  if ( p->pt() == 0 &&  p->phi() == 0 && p->eta()== 0 ) continue;
                  
                  const GenParticle *mom = GPM.getMother(&*p);         //Use getmother or getmotherparton here?
                  while( fabs(mom->pdgId()) == id && GPM.getMother(mom) != nullptr) mom = GPM.getMother(mom);
                  if(mom != nullptr) _gen_wmompdg[nw] = mom->pdgId();
                  else               _gen_wmompdg[nw] = 0;  //Particles for which mompdg is 0 have no mother.
                  ++nw;
                }
              } else if(id == 9900012){
                _gen_majoPt[nMajo] = p->pt();
                _gen_majoE[nMajo] = p->energy();
                _gen_majoEta[nMajo] = p->eta();
                _gen_majoPhi[nMajo] = p->phi();
              
                ++nMajo;
              } 
            }
            _gen_nL = nL;
            _gen_nNu = nNu;
            _gen_nMajo = nMajo;
            _gen_nW= nw;
        }
    }

    getTriggerResults(iEvent, true,  triggerResultsHLTToken,  triggersToSave);
    getTriggerResults(iEvent, false, triggerResultsRECOToken, filtersToSave);
    firstEvent_ = false;

    realdata_ = iEvent.isRealData();

    _weight = 1.;
    // This means we are not getting the weights from NLO samples????? NEED CHECK:
  /*
    if(not isData) {
        edm::Handle<GenEventInfoProduct> pdfvariables;
        iEvent.getByToken(pdfvariablesToken, pdfvariables);
        _weight = pdfvariables->weight();
    } 
  */
    _runNb     = iEvent.id().run();
    _eventNb   = iEvent.id().event();
    _lumiBlock = iEvent.luminosityBlock();
   
    //============ Total number of events is the sum of the events ============
    //============ in each of these luminosity blocks ============
    _nEventsTotalCounted++;
    _hCounter->Fill(0., _weight);
    
   
    
    //============ Beamspot ============
    edm::Handle< reco::BeamSpot > theBeamSpot;
    iEvent.getByToken( IT_beamspot, theBeamSpot );
    BeamSpot::Point  BS= theBeamSpot->position();;
    //==================================
    
    //============ MC Truth Primary vertices ============
    //edm::InputTag IT_MCtruthVtx = edm::InputTag("addPileupInfo");
    //edm::Handle<std::vector<PileupSummaryInfo> > theMCTruthVertices;
    //iEvent.getByLabel("addPileupInfo", theMCTruthVertices) ;
    //if( ! theVertices.isValid() ) ERR(IT_MCTruthVtx ) ;
    //int nMCTruthvertex = theMCTruthVertices->size();

    if (not isData) {
        edm::Handle<vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByToken(PileUpToken, PupInfo);
        //std::cout << "Get Bunch Crossing: ";
        for(vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            //std::cout << PVI->getBunchCrossing() << " ";
            if (PVI->getBunchCrossing() == 0) {
                double ntruth = PVI->getTrueNumInteractions();
                _n_MCTruth_PV = ntruth;
            }
        }
        //std::cout << std::endl;
    }

    //============ Primary vertices ============
    //edm::InputTag IT_goodVtx = edm::InputTag("offlineSlimmedPrimaryVertices");
    edm::Handle<std::vector<Vertex> > theVertices;
    iEvent.getByToken( goodOfflinePrimaryVerticesToken, theVertices) ;
    _n_PV = theVertices->size();
    
    Nvtx->Fill(_n_PV);
    if(! _n_PV ){
        cout << "[WARNING]: No candidate primary vertices passed the quality cuts, so skipping event" << endl;
        return ;
    }

    //std::cout << "Truth and reco vertices: " << _n_MCTruth_PV << " " << _n_PV << std::endl;
    
    Vertex::Point PV = theVertices->begin()->position();
    //std::cout<<PV.x()<<" "<<PV.y()<<" "<<PV.z()<<std::endl;
    const Vertex* PVtx = &((*theVertices)[0]);
    _PVchi2 = PVtx->chi2();
    _PVerr[0] = PVtx->xError();
    _PVerr[1] = PVtx->yError();
    _PVerr[2] = PVtx->zError();
    //==================================
    
    edm::Handle<pat::PackedCandidateCollection> pfcands;       iEvent.getByToken(packedPFCandidatesToken, pfcands);
    edm::Handle<std::vector<pat::Muon>> muons;                 iEvent.getByToken(IT_muon, muons);
    edm::Handle<edm::View<pat::Electron>> electrons;           iEvent.getByToken(IT_electron, electrons);
    edm::Handle<pat::TauCollection> PatTaus;                   iEvent.getByToken(IT_tau, PatTaus);
    edm::Handle<edm::ValueMap<float>> electronMvaIdMap;        iEvent.getByToken(electronMvaIdMapToken, electronMvaIdMap);
    edm::Handle<edm::ValueMap<float>> electronMvaIdHZZ;        iEvent.getByToken(mvaValuesMapToken_HZZ_, electronMvaIdHZZ);
    edm::Handle<edm::ValueMap<bool>> electronMvaIdMap90;       iEvent.getByToken(electronMvaIdMap90Token, electronMvaIdMap90);
    edm::Handle<edm::ValueMap<bool>> electronMvaIdMap80;       iEvent.getByToken(electronMvaIdMap80Token, electronMvaIdMap80);
    edm::Handle<edm::ValueMap<bool>> electronCutBasedIdMapT;   iEvent.getByToken(electronCutBasedIdMapTightToken, electronCutBasedIdMapT);
    edm::Handle<edm::ValueMap<bool>> electronCutBasedIdMapM;   iEvent.getByToken(electronCutBasedIdMapMediumToken, electronCutBasedIdMapM);
    edm::Handle<std::vector<reco::Conversion>> theConversions; iEvent.getByToken(reducedEgammaToken, theConversions);
    edm::Handle<std::vector<pat::Jet>> thePatJets;             iEvent.getByToken(IT_jet , thePatJets );
    edm::Handle<vector<pat::MET>> ThePFMET;                    iEvent.getByToken(IT_pfmet, ThePFMET);
    edm::Handle<double> rhoJets;                               iEvent.getByToken(fixedGridRhoFastjetCentralNeutralToken , rhoJets);
    edm::Handle<double> rhoJECJets;                            iEvent.getByToken(fixedGridRhoFastjetAllToken , rhoJECJets);//kt6PFJets

    myRhoJets    = *rhoJets;
    myRhoJECJets = *rhoJECJets;

    const vector<pat::MET> *pfmetcol = ThePFMET.product();
    const pat::MET *pfmet;
    pfmet = &(pfmetcol->front());
    // Old MET implementation, 2 FEB 2016
//    _met = pfmet->pt();
//    _met_phi = pfmet->phi();

    // New MET implementation
    //==================================

    //in miniAOD v1
    //double rawmetX = pfmet->shiftedP2_74x(pat::MET::METUncertainty(12),pat::MET::Raw).px;
    //double rawmetY = pfmet->shiftedP2_74x(pat::MET::METUncertainty(12),pat::MET::Raw).py;

    //in miniAOD v2
    double rawmetX = pfmet->uncorPx();
    double rawmetY = pfmet->uncorPy();

    //std::cout << "Raw met: " << rawmetX << " " << rawmetY << std::endl;

    
    //double rawmetX = shiftedPt( pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
    //double rawmetY = shiftedPhi( pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);

    
    double corrMetx = rawmetX;
    double corrMety = rawmetY;

    // Ok, this whole part down here is a bit strange and I am not sure if we can trust, but it only affects met which we do not use yet 
    for(auto jet = (*thePatJets).begin(); jet != (*thePatJets).end(); jet++){
        TLorentzVector pJet;
        pJet.SetPtEtaPhiE(((&*jet)->correctedP4("Uncorrected")).Pt(), ((&*jet)->correctedP4("Uncorrected")).Eta(),((&*jet)->correctedP4("Uncorrected")).Phi(), ((&*jet)->correctedP4("Uncorrected")).E());
        const std::vector<reco::CandidatePtr> & cands = (&*jet)->daughterPtrVector();
        for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();cand != cands.end(); ++cand ) {
            const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
            const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
            if ( mu != 0  && ( mu->isGlobalMuon() || mu->isStandAloneMuon() ) ) {
                  //reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
                  TLorentzVector pMu;
                  pMu.SetPtEtaPhiE((*cand)->pt(),(*cand)->eta(),(*cand)->phi(),(*cand)->energy());
                  //rawJetP4 -= muonP4;
                  //  cout << "mu found new method: "<< muonP4.Pt()<< ", "<<muonP4.Eta()<< ", "<<muonP4.Phi()<<endl; 
                  pJet-= pMu;
             }
        }

       
        std::pair <double, double> corr = fMetCorrector->getCorrections(
                                                                      ((&*jet)->correctedP4("Uncorrected")).Pt(),
                                                                      ((&*jet)->correctedP4("Uncorrected")).Eta(),
                                                                      pJet.Pt(),
                                                                      pJet.Phi(),
                                                                      (&*jet)->neutralEmEnergyFraction()+(&*jet)->chargedEmEnergyFraction(),
                                                                      myRhoJECJets,
                                                                      (&*jet)->jetArea(), _runNb);
        corrMetx += corr.first ;
        corrMety += corr.second;
    }

    _met     = sqrt(corrMetx*corrMetx + corrMety*corrMety);
    _met_phi = atan2(corrMety, corrMetx);


    enum decay {
        W_L,  // 0
        W_T_L, // 1
        W_B_L, // 2
        W_B_D_L, //3
        W_B_D_T_L, // 4
        W_B_T_L, // 5
        W_D_L, // 6
        W_D_T_L, //7
        B_L, // 8
        B_D_L, //9
        B_D_T_L, //10
        B_T_L,  // 11
        D_L, //12
        D_T_L, //13
        B_Baryon, // 14
        C_Baryon, //15
        pi_0, //16
        photon_, //17
        F_L, //18
        N_U_L_L, // 19
        W_N_W // 20
    };
    

    //Taus
    std::vector<const pat::Tau* > sTau;
    for (const pat::Tau &tau : *PatTaus) {
        if(tau.pt() < _tauPt) continue;
        if(fabs(tau.eta()) > _tauEta) continue;
        sTau.push_back(&tau);
    }
    
    // This are actually not really "selected" jets, it are basically all jets in the central region, only to be used for lepton variables
    SelectedJetsAll.clear();
    for(auto jet = thePatJets->begin(); jet != thePatJets->end(); ++jet){
      if(jet->pt() < 5) continue;
      if(jet->eta() > 3.0) continue;
      SelectedJetsAll.push_back(&*jet);
    }
    _n_JetsAll = SelectedJetsAll.size();



    /*
     * Leptons
     */
    int leptonCounter = 0;
//  bool _islooseCut[leptonCounter];
    for(auto muon = muons->begin(); muon != muons->end(); ++muon){
      if(muon->pt() < _minPt0)        continue;
      if(std::abs(muon->eta()) > 2.4) continue;
//    if(!muon->isLooseMuon())        continue;  // Store only loose muons, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
      bool goodGlb      = muon->isGlobalMuon() and muon->globalTrack()->normalizedChi2() < 3 and muon->combinedQuality().chi2LocalPosition < 12 and muon->combinedQuality().trkKink < 20;
      bool isMedium     = muon->isLooseMuon() and muon->innerTrack()->validFraction() > 0.8 and muon->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451); // temporary ICHEP recommendation    
      //if (!isMedium) continue;
      if(muon->innerTrack().isNull()) continue;

      if (leptonCounter == 10) continue;  // using arrays, so do not store more than 10 muons, they are sorted in pt anyway
      double relIso = tools::pfRelIso(&*muon,myRhoJets);

      // might be preferable to change to vectors instead of arrays because we do not know the length
      _flavors[leptonCounter]     = 1;
      _charges[leptonCounter]     = muon->charge();
      _isolation[leptonCounter]   = relIso;
      _isolation_absolute[leptonCounter] = tools::pfAbsIso(&*muon);

      _miniisolation[leptonCounter]        = tools::getMiniIsolation(pfcands, &*muon, 0.05, 0.2, 10., myRhoJets, false);
      _miniisolationCharged[leptonCounter] = tools::getMiniIsolation(pfcands, &*muon, 0.05, 0.2, 10., myRhoJets, false);

      _ipPV[leptonCounter]     = muon->innerTrack()->dxy(PV);
      _ipPVerr[leptonCounter]  = muon->innerTrack()->dxyError();
      _ipZPV[leptonCounter]    = muon->innerTrack()->dz(PV);
      _ipZPVerr[leptonCounter] = muon->innerTrack()->dzError();
      _3dIP[leptonCounter]     = muon->dB(pat::Muon::PV3D);
      _3dIPerr[leptonCounter]  = muon->edB(pat::Muon::PV3D);
      _3dIPsig[leptonCounter]  = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);

      // Important: we are not using POG definitions anymore, NEVER call this a "loose muon" in a presentation
      _isloose[leptonCounter] = std::abs(muon->eta()) < 2.4 && muon->pt() > 5 && fabs(_ipPV[leptonCounter]) < 0.05 && fabs(_ipZPV[leptonCounter]) < 0.1 && muon->isPFMuon() && (muon->isTrackerMuon() || muon->isGlobalMuon() ) && _isolation[leptonCounter] < 0.6;

      if(std::abs(muon->innerTrack()->dxy(PV)) > 0.05) continue;
      if(std::abs(muon->innerTrack()->dz(PV)) > 0.1  ) continue;

      if(relIso > 0.6) continue; // loose selection for FR


     
      //bool goodGlb      = muon->isGlobalMuon() and muon->globalTrack()->normalizedChi2() < 3 and muon->combinedQuality().chi2LocalPosition < 12 and muon->combinedQuality().trkKink < 20;
      //bool isMedium     = muon->isLooseMuon() and muon->innerTrack()->validFraction() > 0.49 and muon->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451); // temporary ICHEP recommendation
//    bool isMedium     = muon->isMediumMuon();  // default definition with  muon->innerTrack()->validFraction() > 0.8 instead
 
      _muonSegmentComp[leptonCounter] = muon->segmentCompatibility(); // do we really use this variable in our trees? should be enough to select loose or medium
 
      //_isloose[leptonCounter]  = muon->isLooseMuon();
      _ismedium[leptonCounter] = isMedium;
      _istight[leptonCounter]  = muon->isTightMuon(*PVtx);

      //if(!_isloose[leptonCounter]) continue;
      TLorentzVector vect;
      vect.SetPtEtaPhiE(muon->pt(),muon->eta(),muon->phi(), muon->energy());
      _mt[leptonCounter] = tools::MT_calc(vect, _met, _met_phi);

      // So for some reason, here we are adding again the same stuff to out tree
      _lPt[leptonCounter]  = muon->pt();
      _lEta[leptonCounter] = muon->eta(); 
      _lPhi[leptonCounter] = muon->phi();
      _lE[leptonCounter]   = muon->energy();

      fillCloseJetVars(leptonCounter, PV);
      for(TString wp : {"VT","T","M","L","VL"}){
        leptonWorkingPoints["multiIsolation" + wp]->push_back(tools::passMultiIsolation(wp, _miniisolation[leptonCounter], _ptratio[leptonCounter], _ptrel[leptonCounter]));
        leptonConeCorrectedPt[wp]->push_back(tools::leptonConeCorrectedPt(_lPt[leptonCounter], wp, _miniisolation[leptonCounter], _ptratio[leptonCounter], _ptrel[leptonCounter]));
      }
      leptonConeCorrectedPt["MiniIso04"]->push_back(tools::leptonConeCorrectedPt(_lPt[leptonCounter], "MiniIso04", _miniisolation[leptonCounter], 0, 0));
      leptonConeCorrectedPt["RelIso04"]->push_back(tools::leptonConeCorrectedPt(_lPt[leptonCounter], "RelIso04", _isolation[leptonCounter], 0, 0));

      if(not isData){  
        const GenParticle* mc = GPM.matchedMC(&*muon);
        if(mc!=0){
          fillMCVars(mc, leptonCounter);
          _ipPVmc[leptonCounter] = std::abs(muon->innerTrack()->dxy(PVmc));
        } else {
          _origin[leptonCounter] = -1;
          _originReduced[leptonCounter] = -1;
          
          _mompt[leptonCounter]  = 0;
          _momphi[leptonCounter] = 0;
          _mometa[leptonCounter] = 0;
          _mompdg[leptonCounter] = 0;
        }
      }
      leptonCounter++;
    }
    _nMu = leptonCounter;

    for(auto electron = electrons->begin(); electron != electrons->end(); ++electron){
      if (leptonCounter == 10) continue; // This will go terribly wrong when there are about 8-10 muons or more (even though this is an improbable situation)

      _passedCutBasedIdTight[leptonCounter] = false; // unitialized for muons? probably not a big issue though
      _passedMVA90[leptonCounter]           = false;
      _passedMVA80[leptonCounter]           = false;
      _findMatched[leptonCounter]           = -1;

      if(electron->pt() < _minPt1)          continue;
      if(std::abs(electron->eta()) > 2.5)   continue;
      if(!electron->gsfTrack().isNonnull()) continue;
      //if( TMath::Abs(gsfTrack->dxy(PV)) > 0.05  )  continue;
      //if( TMath::Abs(gsfTrack->dz(PV)) > 0.1  ) continue;

      _lCleanCut[leptonCounter] = true;
        for(int m = 0; m < _nMu; ++m){
          if(!_isloose[m]) continue;
          TLorentzVector mu, el;
          mu.SetPtEtaPhiE(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
          el.SetPtEtaPhiE(electron->pt(), electron->eta(), electron->phi(), electron->energy());
          if(el.DeltaR(mu) < 0.05){
            _lCleanCut[leptonCounter] = false;
            break;
          }
        }    
      if(!_lCleanCut[leptonCounter]) continue;

      //if(!tools::isLooseCutBasedElectronWithoutIsolation(&*electron)) continue; //only store those passing the loose cut based id, without isolation requirement
      double relIso = tools::pfRelIso(&*electron, myRhoJets);
      if(relIso > 0.6) continue;
      // There will be a new electron MVA soon
      edm::RefToBase<pat::Electron> electronRef(edm::Ref<edm::View<pat::Electron>>(electrons, (electron - electrons->begin())));
      _mvaValue[leptonCounter]               = (*electronMvaIdMap)[electronRef];
      _passedCutBasedIdTight[leptonCounter]  = (*electronCutBasedIdMapT)[electronRef];
      _passedCutBasedIdMedium[leptonCounter] = (*electronCutBasedIdMapM)[electronRef];
      _passedMVA80[leptonCounter]            = (*electronMvaIdMap80)[electronRef];
      _passedMVA90[leptonCounter]            = (*electronMvaIdMap90)[electronRef];

      _passedMVA_SUSY[leptonCounter][0] = tools::passed_loose_MVA_FR_slidingCut( &*electron, _mvaValue[leptonCounter], _mvaValue_HZZ[leptonCounter]);
      _passedMVA_SUSY[leptonCounter][1] = tools::passed_medium_MVA_FR_slidingCut(&*electron, _mvaValue[leptonCounter]);
      _passedMVA_SUSY[leptonCounter][2] = tools::passed_tight_MVA_FR_slidingCut( &*electron, _mvaValue[leptonCounter]);

      //bool crossCheckTight = tools::isTightCutBasedElectronWithoutIsolation(&*electron, false) and tools::pfRelIso(&*electron, myRhoJECJets) < (electron->isEB() ? 0.0588 : 0.0571);

      _flavors[leptonCounter]            = 0;
      _charges[leptonCounter]            = electron->charge();
      _isolation[leptonCounter]          = relIso;
      _isolation_absolute[leptonCounter] = tools::pfAbsIso(&*electron, myRhoJECJets);

      _ipPV[leptonCounter]  = electron->gsfTrack()->dxy(PV);
      _ipZPV[leptonCounter] = electron->gsfTrack()->dz(PV);

      if(std::abs(_ipPV[leptonCounter]) > 0.05) continue;
      if(std::abs(_ipZPV[leptonCounter]) > 0.1) continue;
      
      _miniisolation[leptonCounter]        = tools::getMiniIsolation(pfcands, &*electron, 0.05, 0.2, 10., myRhoJets, false);
      _miniisolationCharged[leptonCounter] = tools::getMiniIsolation(pfcands, &*electron, 0.05, 0.2, 10., myRhoJets, false);

      _chargeConst[leptonCounter]      = electron->isGsfCtfScPixChargeConsistent();
      _vtxFitConversion[leptonCounter] = ConversionTools::hasMatchedConversion(reco::GsfElectron (*&*electron), theConversions, BS);
      _hitsNumber[leptonCounter]       = electron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

      _trigEmulator[leptonCounter]    = tools::triggerEmulatorReturned(&*electron);
      _isotrigEmulator[leptonCounter] = tools::isoTriggerEmulator(&*electron);

      _3dIP[leptonCounter]    = electron->dB(pat::Electron::PV3D);
      _3dIPerr[leptonCounter] = electron->edB(pat::Electron::PV3D);
      _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);

      if( _vtxFitConversion[leptonCounter] )  continue;
      if (_hitsNumber[leptonCounter] > 0) continue;   
     
      //if(_isolation[leptonCounter] > 1) continue;
      //if(_3dIPsig[leptonCounter] > 4) continue;

      TLorentzVector vect_e;
      vect_e.SetPtEtaPhiE(electron->pt(),electron->eta(),electron->phi(), electron->energy());
      _mt[leptonCounter] = tools::MT_calc(vect_e, _met, _met_phi);
      
      _lPt[leptonCounter]  = electron->pt();
      _lEta[leptonCounter] = electron->eta();
      _lPhi[leptonCounter] = electron->phi();
      _lE[leptonCounter]   = electron->energy();

     
      fillCloseJetVars(leptonCounter, PV);

      for(TString wp : {"VT","T","M","L","VL"}){
        leptonWorkingPoints["multiIsolation" + wp]->push_back(tools::passMultiIsolation(wp, _miniisolation[leptonCounter], _ptratio[leptonCounter], _ptrel[leptonCounter]));
        leptonConeCorrectedPt[wp]->push_back(tools::leptonConeCorrectedPt(_lPt[leptonCounter], wp, _miniisolation[leptonCounter], _ptratio[leptonCounter], _ptrel[leptonCounter]));
      }
      leptonConeCorrectedPt["MiniIso04"]->push_back(tools::leptonConeCorrectedPt(_lPt[leptonCounter], "MiniIso04", _miniisolation[leptonCounter], 0, 0));
      leptonConeCorrectedPt["RelIso04"]->push_back(tools::leptonConeCorrectedPt(_lPt[leptonCounter], "RelIso04", _isolation[leptonCounter], 0, 0));

      if (not isData) {
        _findMatched[leptonCounter]=-1;
        const GenParticle* mc = GPM.matchedMC(&*electron,11);
        _originPhot[leptonCounter] = -2;
        if ( mc!=0 ) {
          fillMCVars(mc, leptonCounter);
            if (mc->pdgId() != 22) //not a photon => should be a fake, and GenParticleManager does the job
              _origin[leptonCounter] = GPM.originReduced(_originDetailed[leptonCounter]);
            std::cout<<" != 22"<<std::endl;
            else {
              std::cout<<"------------>     == 22"<<std::endl;
              _originPhot[leptonCounter] = photonOrigin(mc); //is from conversion, so I store info on where a photon come from (fragmentation, FSR etc)
              GPM.printInheritance(&(*mc));

            } 
          
            //Vertex::Point PVmc = mcMom->vertex();
          _ipPVmc[leptonCounter] = std::abs(electron->gsfTrack()->dxy(PVmc));
          _findMatched[leptonCounter]=1;
        } else {
          _originReduced[leptonCounter] = -1;
          _originDetailed[leptonCounter] = -1;
          _origin[leptonCounter] = 4;
         _mompt[leptonCounter] = 0;
         _momphi[leptonCounter] = 0;
          _mometa[leptonCounter] = 0;
          _mompdg[leptonCounter] = 0;
          _findMatched[leptonCounter]=0;
        }
      }
      
      

      
      
      leptonCounter++;
    }
    _nEle = leptonCounter-_nMu;
    _nTau = 0;

    _nLeptons = leptonCounter;
    _sb = false;
    _n_Jets = 0;
    _n_bJets = 0;
    HT = 0;

    std::vector<const pat::Jet*> SelectedJets = tools::JetSelector(*thePatJets, _jetPtCut, _jetEtaCut);
    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){

        double uncPt = (SelectedJets[i]->correctedP4("Uncorrected")).Pt();
        double uncEta = (SelectedJets[i]->correctedP4("Uncorrected")).Eta();
         
        //double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),_corrLevel);
        double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJECJets, SelectedJets[i]->jetArea(),_corrLevel, _runNb);
        //std::cout << "Before correction pt of jet: " << uncPt << " " << uncPt*corr << std::endl;
                                 
        if (uncPt*corr < 25) continue;
        //if (uncPt*corr < 20) continue;
                                         
        _jetEta[_n_Jets] = SelectedJets[i]->eta();
        _jetPhi[_n_Jets] = SelectedJets[i]->phi();
        _jetPt[_n_Jets] = uncPt*corr;
        _jetE[_n_Jets] = (SelectedJets[i]->correctedP4("Uncorrected")).E()*corr;

        _jetFlavour[_n_Jets] = SelectedJets[i]->hadronFlavour();
        
        TLorentzVector jt; jt.SetPtEtaPhiE(_jetPt[_n_Jets],_jetEta[_n_Jets],_jetPhi[_n_Jets], _jetE[_n_Jets]);
        //std::cout<<_jetPt[_n_Jets]<<" "<<_jetEta[_n_Jets]<<" "<<_jetPhi[_n_Jets]<<std::endl;
        //TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[_n_Jets],_jetEta[_n_Jets],_jetPhi[_n_Jets],0);

        bool clean = true;
        for (int k=0; k!=_nLeptons; ++k) {
          if (_isloose[k] && _lPt[k] > 10){
            TLorentzVector vect_k;
            vect_k.SetPtEtaPhiE(_lPt[k],_lEta[k],_lPhi[k], _lE[k]);

            double dR1 = (vect_k.DeltaR( jt ));
            //std::cout << "jet cleaning: " << dR1 << " " << ((TLorentzVector *)_leptonP4->At(k))->Pt()  << std::endl;
            clean = clean && (dR1 > 0.4) ;
          }
        }
        _clean[_n_Jets]=0;
        if (clean) _clean[_n_Jets] =1;
        else       _clean[_n_Jets] =-1;
                
        for(int j=0; j != _nLeptons; ++j){
          TLorentzVector vect_i;
          vect_i.SetPtEtaPhiE(_lPt[j],_lEta[j],_lPhi[j], _lE[j]);
          _jetDeltaR[_n_Jets][j] = (vect_i.DeltaR( jt ) );
        }

        _csv[_n_Jets] = SelectedJets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        if(!realdata_){
          jecUnc->setJetEta(SelectedJets[i]->eta());
          jecUnc->setJetPt(SelectedJets[i]->pt());
          double unc = jecUnc->getUncertainty(true);
          _jecUnc[_n_Jets] = unc;
          //cout << " JEC unc " << unc << endl;
          _jetPtUp[_n_Jets] = (1+unc)*SelectedJets[i]->pt();
          _jetPtDown[_n_Jets] = (1-unc)*SelectedJets[i]->pt();
        }

        if(_csv[_n_Jets] > 0.460) {
          _bTagged[_n_Jets] = true;
          _n_bJets++;
        } else _bTagged[_n_Jets] = false;
        
        if (_jetPt[_n_Jets] > 25) HT+= _jetPt[_n_Jets];
        _n_Jets++;
    }

    // Check if this event has triggered at least one of our save triggers
    bool triggerFound = false;
    for(auto trigger : triggersToSave){
      if(triggerFlags[trigger]) triggerFound = true;
    }

    // Here we make the decision on what to save
    if(treeForFakeRate){
      if(!triggerFound)        return;
      //if(_nLeptons != 1)       return;
      //if(_n_Jets < 1)          return;        // For fake rate tree: exactly 1 loose lepton + at least 1 jet
      //if(_jetPt[0] < 30)       return;        // with deltaR(j, l) > 1 (back-to-back)
      //if(_jetDeltaR[0][0] < 1) return;
      if(_nLeptons < 3) return; // fake ttbar

    } else if(singleLep){                     // Important: do not require trigger for singleLep trees which we use to measure trigger efficiencies
      if(_nLeptons < 1) return;
    } else {
      if(!triggerFound) return;
      if(_nLeptons < 3) return;
    }
    outputTree->Fill();
}

void trilepton::fillMCVars(const GenParticle* mc, const int leptonCounter) {
    
    _lPtmc[leptonCounter] = mc->pt();
    _lEmc[leptonCounter] = mc->energy();
    _lPhimc[leptonCounter] = mc->phi();
    _lEtamc[leptonCounter] = mc->eta();
    _lpdgmc[leptonCounter] = mc->pdgId();
    _lchargemc[leptonCounter] = mc->charge();
    //GPM.printInheritance(mc);
    //std::cout << "pdg of tau truth: " << mc->pdgId() << std::endl;
    
  
    _originDetailed[leptonCounter] = GPM.origin(mc);
    _origin[leptonCounter] = GPM.originReduced(_originDetailed[leptonCounter]);
    _originReduced[leptonCounter] = GPM.originReduced(_originDetailed[leptonCounter]);
    _isPromptFinalState[leptonCounter] = mc->isPromptFinalState();
    _fromHardProcessFinalState[leptonCounter] = mc->fromHardProcessFinalState();
  
  if (_isPromptFinalState[leptonCounter] || _fromHardProcessFinalState[leptonCounter]) {
        _origin[leptonCounter] = 0;
    }

    

    //std::cout << "Origin of tau: " << _originReduced[leptonCounter] << std::endl;

    unsigned int n = mc->numberOfDaughters();
    _flag_lepton[leptonCounter] = false;
    for(unsigned int j = 0; j < n; ++j) {
        const Candidate * d = mc->daughter( j );
        int dauId = d->pdgId();
        if(dauId == 11 || dauId == 13)
            _flag_lepton[leptonCounter] = true;
    }
    //std::cout << "Tau leptonic: " << _flag_lepton[leptonCounter] << std::endl;

    const GenParticle* mcMom = GPM.getMotherParton(mc);
    if (mcMom!=NULL) {
        
        //std::cout<<"Mother: "<<std::endl;
        //std::cout<<mcMom->pdgId()<<" "<<mcMom->pt()<<std::endl;
        //std::cout<<mc->pt()<<std::endl;
        
        _mompt[leptonCounter] = mcMom->pt();
        _momphi[leptonCounter] = mcMom->phi();
        _mometa[leptonCounter] = mcMom->eta();
        _mompdg[leptonCounter] = mcMom->pdgId();
        
        PVmc = mcMom->vertex();
        
        //std::cout<<"d0: "<<_ipPVmc<<" "<<_ipPV<<std::endl;
        
        TLorentzVector Gen0;
        Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
        
        
        _isolationMC[leptonCounter][0] = 0;
        _isolationMC[leptonCounter][1] = 0;
        std::vector<GenParticleRef> dauts;
        for (unsigned int dau=0; dau!=mcMom->numberOfDaughters(); ++dau) {
            GenParticleRef daut = mcMom->daughterRef(dau);
            dauts.push_back(daut);
        }
        unsigned int counterD = 0;
        //if (_isolation == 0)
        //    std::cout<<"==== "<<chargedHadronIso*iM->pt()<<" "<<neutralHadronIso*iM->pt()<<" "<<photonIso*iM->pt()<<" "<<beta*iM->pt()<<std::endl;
        while (counterD < dauts.size()) {
            if (dauts.at(counterD)->status() == 1) {
                _isolationMC[leptonCounter][0]+=dauts.at(counterD)->pt();
                if ((fabs(dauts.at(counterD)->pdgId())!= 12) && (fabs(dauts.at(counterD)->pdgId())!= 14) && (fabs(dauts.at(counterD)->pdgId())!= 16)) {
                    _isolationMC[leptonCounter][1]+=dauts.at(counterD)->pt();
                } else {
                    TLorentzVector Gen;
                    Gen.SetPtEtaPhiE( dauts.at(counterD)->pt(), dauts.at(counterD)->eta(), dauts.at(counterD)->phi(), dauts.at(counterD)->energy() );
                    Gen0 += Gen;
                }
            } else {
                for (unsigned int dau=0; dau!=dauts.at(counterD)->numberOfDaughters(); ++dau) {
                    GenParticleRef daut = dauts.at(counterD)->daughterRef(dau);
                    dauts.push_back(daut);
                }
            }
            counterD++;
        }
        
        //std::cout<<"PTs: "<<mc->pt()<<" "<<iM->pt()<<std::endl;
        _isolationMC[leptonCounter][0]/=mc->pt();
        _isolationMC[leptonCounter][0]-=1;
        _isolationMC[leptonCounter][1]/=mc->pt();
        _isolationMC[leptonCounter][1]-=1;
        
        _nuPtmc[leptonCounter] = Gen0.Pt();
        _nuPhimc[leptonCounter] = Gen0.Phi();
        _nuEtamc[leptonCounter] = Gen0.Eta();
        _nuEmc[leptonCounter] = Gen0.E();
        
        TLorentzVector lmc;
        lmc.SetPtEtaPhiE(_lPtmc[leptonCounter],_lEtamc[leptonCounter],_lPhimc[leptonCounter],_lEmc[leptonCounter]);
        _mtmc[leptonCounter] = tools::MT_calc(lmc, _nuPtmc[leptonCounter], _nuPhimc[leptonCounter]);
    } else {
        _mompt[leptonCounter] = 0;
        _momphi[leptonCounter] = 0;
        _mometa[leptonCounter] = 0;
        _mompdg[leptonCounter] = 0;
    }
    
}


int trilepton::photonOrigin(const GenParticle* photon) {
    //implementation of a flag for the photon truth matching:
    //https://hypernews.cern.ch/HyperNews/CMS/get/susy-interpretations/192.html
    //-1: not a photon
    //0: direct prompt photons (prompt and delta R > 0.4)
    //1: fragmentation photons (prompt and delta R < 0.4)
    //2: non-prompt photons

    
    if (photon->pdgId() != 22) return -1;
    if (!photon->isPromptFinalState()) return 2;
    if (photon->pt() < 10) return 1;

    
    TLorentzVector photonP4; photonP4.SetPtEtaPhiE(photon->pt(),photon->eta(),photon->phi(),photon->energy());
    TLorentzVector partonP4;
    bool smallDR = false;
    for( unsigned p=0; p<TheGenParticles->size(); ++p ){
        if ((*TheGenParticles)[p].status() != 23) continue;
        if (!((*TheGenParticles)[p].pdgId() == 21 || (fabs((*TheGenParticles)[p].pdgId()) > 0 && fabs((*TheGenParticles)[p].pdgId()) <7))) continue;
        partonP4.SetPtEtaPhiE((*TheGenParticles)[p].pt(),(*TheGenParticles)[p].eta(),(*TheGenParticles)[p].phi(),(*TheGenParticles)[p].energy());
        double deltaR = photonP4.DeltaR(partonP4);
        smallDR |= (deltaR < 0.05);
    }
    if (smallDR) return 1;
    else return 0;
}











void trilepton::fillCloseJetVars(const int leptonCounter, Vertex::Point PV) {

    _closeJetPtAll[leptonCounter] = 0;
    _closeJetAngAll[leptonCounter] = 10000;
    _ptRelAll[leptonCounter] = 0;
    _ptrel[leptonCounter] = 0;
    _ptratio[leptonCounter] = 1.;

    _closeIndex[leptonCounter] = 0;
    TLorentzVector pJet;
    //std::cout << "Selected Jets Number:\t"<< SelectedJetsAll.size() << std::endl;
    for(unsigned int k = 0 ; k < SelectedJetsAll.size() ;k++ ){
      double uncorrPt = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Pt();
      TLorentzVector vect_j;
      vect_j.SetPtEtaPhiE(_lPt[leptonCounter],_lEta[leptonCounter],_lPhi[leptonCounter], _lE[leptonCounter]);

      double corr = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJECJets, SelectedJetsAll[k]->jetArea(),"L1FastJet", _runNb);
      double corr2 = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJECJets, SelectedJetsAll[k]->jetArea(),_corrLevel, _runNb);
        
      pJet.SetPtEtaPhiE( corr*uncorrPt, SelectedJetsAll[k]->eta(), SelectedJetsAll[k]->phi(), corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E());
      pJet-=vect_j;
       
      pJet*=corr2/corr;
      pJet+=vect_j;

              
      double ang = vect_j.DeltaR( pJet );
      if (ang < _closeJetAngAll[leptonCounter]) {
        _closeJetAngAll[leptonCounter] = vect_j.DeltaR( pJet );
        _closeJetPtAll[leptonCounter] = pJet.Pt();
                                
      //_closeJetEtaAll[leptonCounter] = pJet.Eta();
      //_closeJetPhiAll[leptonCounter] = pJet.Phi();
      //_closeJetEAll[leptonCounter] = pJet.E();
      //_closeJetMAll[leptonCounter] = pJet.M();
      //_closeJetNconstAll[leptonCounter] = SelectedJetsAll[k]->numberOfDaughters();
        _closeJetCSVAll[leptonCounter] = SelectedJetsAll[k]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    
        _ptrel[leptonCounter] = vect_j.Perp(pJet.Vect() - (vect_j.Vect()));

        _ptratio[leptonCounter] = (vect_j.Pt())/(_closeJetPtAll[leptonCounter]);

        _closeIndex[leptonCounter] = k;

        int trackSelectionMult = 0;
        for(int unsigned it=0; it!=SelectedJetsAll[_closeIndex[leptonCounter]]->numberOfDaughters(); ++it) {
          const pat::PackedCandidate * icand = dynamic_cast<const pat::PackedCandidate *> (SelectedJetsAll[_closeIndex[leptonCounter]]->daughter(it));
          const reco::Track& trk = icand->pseudoTrack();
          Double_t deta = trk.eta() - pJet.Eta();
          Double_t dphi = TVector2::Phi_mpi_pi(trk.phi() - pJet.Phi());
          bool trackSelection =  (trk.pt() > 1) && trk.charge() != 0 && (TMath::Sqrt( deta*deta+dphi*dphi ) < 0.4) && (icand->fromPV() > 1) && (trk.hitPattern().numberOfValidHits()>=8) && (trk.hitPattern().numberOfValidPixelHits()>=2) && (trk.normalizedChi2()<5) && std::fabs(trk.dxy(PV))<0.2 && std::fabs(trk.dz(PV))<17;
          if(trackSelection) trackSelectionMult++;
        }

      _trackSelectionMultiplicity[leptonCounter] = trackSelectionMult;
      }
    }

    if (_closeJetAngAll[leptonCounter] > 0.4) {
        _closeJetPtAll[leptonCounter] = _lPt[leptonCounter];
        //_closeJetPhiAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        //_closeJetEtaAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        //_closeJetNconstAll[leptonCounter] = 1;
        //_closeJetCSVAll[leptonCounter] = 0;
        _closeJetAngAll[leptonCounter] = 0; 
        //_ptRelAll[leptonCounter] = 0;
        _ptrel[leptonCounter] = 0;
        _ptratio[leptonCounter] = 1.;
        _closeIndex[leptonCounter] = -1;
    }   
}


void trilepton::matchCloseJet(const int leptonCounter) {
    double minDeltaR3 = 9999;
    double minDeltaR2 = 9999;
    TLorentzVector Gen1, Gen2;
    const GenParticle* mc3a = 0;
    const GenParticle* mc2a = 0;

    if(SelectedJetsAll.size() != 0) {
      Gen1.SetPtEtaPhiE(SelectedJetsAll.at(_closeIndex[leptonCounter])->pt(),SelectedJetsAll.at(_closeIndex[leptonCounter])->eta(),SelectedJetsAll.at(_closeIndex[leptonCounter])->phi(),SelectedJetsAll.at(_closeIndex[leptonCounter])->energy());
      for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
        int id = std::abs(p->pdgId());
    
        if ((id > 0 && id < 6) || (id == 21) || (id == 22)) {
         if (p->status() != 2) {
           if (fabs(p->eta()) <= 10 ) Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
           else                       continue;
         //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
           double deltaRcur = Gen1.DeltaR(Gen2);
           if (deltaRcur < minDeltaR3) {
              mc3a = &*p;
              minDeltaR3 = deltaRcur;
           }
         } else if (p->status() == 2) {
            if (fabs(p->eta()) <= 10 ) Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
            else                       continue;
            //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
            double deltaRcur = Gen1.DeltaR(Gen2);
            if (deltaRcur < minDeltaR2) {
              mc2a = &*p;
              minDeltaR2 = deltaRcur;
            }
          }
        }
      }
    } 

    if ((minDeltaR3 < 0.5 || minDeltaR2 == 9999) && mc3a!=0) {
        _closeJetPtAllMC[leptonCounter] = mc3a->pt();
        _closeJetPtAllstatus[leptonCounter] = 3;
        _partonIdMatched[leptonCounter] = mc3a->pdgId();
    } else if (mc2a!=0) {
        _closeJetPtAllMC[leptonCounter] = mc2a->pt();
        _closeJetPtAllstatus[leptonCounter] = 2;
        _partonIdMatched[leptonCounter] = mc2a->pdgId();
    } else {
        _closeJetPtAllMC[leptonCounter] = 0;
        _closeJetPtAllstatus[leptonCounter] = 0;
        _partonIdMatched[leptonCounter] = 0;
    }
}


void trilepton::fillIsoMCVars(const int leptonCounter) {
    
    if( TheGenParticles.isValid() ) {
        //if (_isolation[0] == 0) {
        //    std::cout<<"==="<<std::endl;
        //}
        _isolationMC[leptonCounter][2] = 0;
        _isolationMC[leptonCounter][3] = 0;
        for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
            //int id = std::abs(p->pdgId());
            if (p->status() == 1) {
                TLorentzVector pmc; pmc.SetPtEtaPhiM( p->pt(), p->eta(), p->phi(), p->mass() );
                TLorentzVector pm1; pm1.SetPtEtaPhiE(_lPt[leptonCounter],_lEta[leptonCounter],_lPhi[leptonCounter], _lE[leptonCounter]);
                double ang = pm1.DeltaR( pmc );
                if (ang < 0.3) {
                    _isolationMC[leptonCounter][2]+=p->pt();
                    if (fabs(p->pdgId())!= 12 && fabs(p->pdgId())!= 14 && fabs(p->pdgId())!= 16) {
                        _isolationMC[leptonCounter][3]+=p->pt();
                        //if (_isolation[0] == 0) {
                        //    std::cout<<" "<<p->pdgId()<<" "<<p->pt()<<
                        //    " "<<p->eta()<<" "<<p->phi()<<" "<<p->energy()<<" in "<<ang<<std::endl;
                        //}
                    }
                }
            }
        }
        double leptpt = _lPt[leptonCounter];
        _isolationMC[leptonCounter][2]/=leptpt;
        _isolationMC[leptonCounter][3]/=leptpt;
        _isolationMC[leptonCounter][2]-=1;
        _isolationMC[leptonCounter][3]-=1;
        //std::cout<<"*************"<<std::endl;
        //for (int i=0; i!=4; ++i)
        //    std::cout<<_isolationMC[i]<<" ";
        //std::cout<<_isolation[0]<<" "<<_closeJetPtAllMC[0]/leptpt-1<<" "<<_closeJetPtAll[0]/leptpt-1<<std::endl;
    }
}

void trilepton::fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu) {
    
    hJet_ptRaw = (jet->correctedP4("Uncorrected")).Pt();
    hJet_genPt = genpt;
    hJet_pt = jet->pt();
    hJet_phi = jet->phi();
    hJet_eta = jet->eta();
    hJet_e = jet->energy();
    
    hJet_ptLeadTrack = 0;
    
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack) hJet_ptLeadTrack = tracks[k]->pt();
    }
    
    hJet_vtx3dL = 0;
    hJet_vtx3deL = 0;
    hJet_vtxMass = 0;
    hJet_vtxPt = 0;
    
    const reco::SecondaryVertexTagInfo* scdVtx = jet->tagInfoSecondaryVertex("secondaryVertex");
    
    if (scdVtx) {
        //std::cout<<"Vertetx info: "<<scdVtx->nVertices()<<std::endl;
        if (scdVtx->nVertices()) {
            const reco::Vertex &sv1 = scdVtx->secondaryVertex(0);
            if (!sv1.isFake()) {
                Measurement1D distance1 = scdVtx->flightDistance(0, true);
                hJet_vtx3dL = distance1.value();
                hJet_vtx3deL = distance1.error();

                math::XYZTLorentzVectorD p4vtx = sv1.p4();
                hJet_vtxMass = p4vtx.M();
                hJet_vtxPt = p4vtx.Pt();
            }
        }
    }
    
    
    hJet_cef = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction();
    hJet_nconstituents = jet->getPFConstituents().size();
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta(), _runNb);
    
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    
    hJet_SoftLeptptRel = pLep.Perp(pJet.Vect());
    hJet_SoftLeptPt = mu->pt();
    hJet_SoftLeptdR = pLep.DeltaR(pJet);
    
    hJet_SoftLeptIdlooseMu = 1;
    hJet_SoftLeptId95 = 1;
}

void trilepton::fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* mu) {
    
    hJet_ptRaw = (jet->correctedP4("Uncorrected")).Pt();
    hJet_genPt = genpt;
    hJet_pt = jet->pt();
    hJet_phi = jet->phi();
    hJet_eta = jet->eta();
    hJet_e = jet->energy();
    
    hJet_ptLeadTrack = 0;
    
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack) hJet_ptLeadTrack = tracks[k]->pt();
    }
    
    hJet_vtx3dL = 0;
    hJet_vtx3deL = 0;
    hJet_vtxMass = 0;
    hJet_vtxPt = 0;
    
    const reco::SecondaryVertexTagInfo* scdVtx = jet->tagInfoSecondaryVertex("secondaryVertex");
    
    if (scdVtx) {
        //std::cout<<"Vertetx info: "<<scdVtx->nVertices()<<std::endl;
        if (scdVtx->nVertices()) {
            const reco::Vertex &sv1 = scdVtx->secondaryVertex(0);
            if (!sv1.isFake()) {
                Measurement1D distance1 = scdVtx->flightDistance(0, true);
                hJet_vtx3dL = distance1.value();
                hJet_vtx3deL = distance1.error();

                math::XYZTLorentzVectorD p4vtx = sv1.p4();
                hJet_vtxMass = p4vtx.M();
                hJet_vtxPt = p4vtx.Pt();
            }
        }
    }
    
    
    hJet_cef = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction();
    hJet_nconstituents = jet->getPFConstituents().size();
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta(), _runNb);
    
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    
    hJet_SoftLeptptRel = pLep.Perp(pJet.Vect());
    hJet_SoftLeptPt = mu->pt();
    hJet_SoftLeptdR = pLep.DeltaR(pJet);
    
    hJet_SoftLeptIdlooseMu = 1;
    hJet_SoftLeptId95 = 1;
}


DEFINE_FWK_MODULE(trilepton);

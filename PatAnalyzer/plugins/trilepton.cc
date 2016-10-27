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

// In the long run, we should get rid of these using namespace statements, this is bad coding
using namespace std;
using namespace edm;
using namespace reco;
using namespace tools;
using namespace math;
using namespace reco::tau;

trilepton::trilepton(const edm::ParameterSet & iConfig) :
  genparticleToken(                      consumes<reco::GenParticleCollection>(   iConfig.getParameter<edm::InputTag>("genPartsLabel"))),
  mvaValuesMapToken(                     consumes<edm::ValueMap<float>>(          iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
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
  _relIsoCutE(999.),//0.15
  _relIsoCutMu(999.),//0.15
  _relIsoCutEloose(999.), //0.5
  _relIsoCutMuloose(999.), //0.5
  _chargeConsistency(true),
  _minPt0(5),
  _minPt1(5.),
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
    isData              = iConfig.getUntrackedParameter<bool>("isData") ;
    SampleName          = iConfig.getUntrackedParameter<std::string>("SampleName") ;

    // eventually move this to config level
    // just a bunch of lepton triggers, might need a closer look for cleanup or additions
    triggersToSave = {"HLT_TripleMu_12_10_5", "HLT_DiMu9_Ele9_CaloIdL_TrackIdL", "HLT_Mu8_DiEle12_CaloIdL_TrackIdL", "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
                      "HLT_Mu1c7_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                      "HLT_TrkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                      "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                      "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu8_TrkIsoVVL","HLT_Mu17_TrkIsoVVL","HLT_Mu24_TrkIsoVVL","HLT_Mu34_TrkIsoVVL",
                      "HLT_Mu8","HLT_Mu17","HLT_Mu24","HLT_Mu34",
                      "HLT_Ele17_CaloIdL_TrackIdL_IsoVL","HLT_Ele23_CaloIdL_TrackIdL_IsoVL",
                      "HLT_IsoMu22", "HLT_IsoMu24", "HLT_IsoTkMu22","HLT_IsoMu18","HLT_Ele27_WPTight_Gsf","HLT_Ele27_eta2p1_WPLoose_Gsf"};
    filtersToSave  = {"Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_eeBadScFilter", "Flag_globalTightHalo2016Filter"};
}


void trilepton::beginJob()
{
    outputTree = fs->make<TTree>("trileptonTree","trileptonTree");
    Nvtx = fs->make<TH1F>("N_{vtx}", "Number of vertices;N_{vtx};events / 1"  ,    40, 0., 40.);
    
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);

    _leptonP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
        new ( (*_leptonP4)[i] ) TLorentzVector();
    }
    outputTree->Branch("_leptonP4", "TClonesArray", &_leptonP4, 32000, 0);

    // =============== GEN ================================
    _gen_leptonP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
      new ( (*_gen_leptonP4)[i] ) TLorentzVector();
    }
    outputTree->Branch("_gen_leptonP4", "TClonesArray", &_gen_leptonP4, 32000, 0);

    _gen_nuP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
      new ( (*_gen_nuP4)[i] ) TLorentzVector();
    }
    outputTree->Branch("_gen_nuP4", "TClonesArray", &_gen_nuP4, 32000, 0);

    _gen_majoP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
      new ( (*_gen_majoP4)[i] ) TLorentzVector();
    }
    outputTree->Branch("_gen_majoP4", "TClonesArray", &_gen_majoP4, 32000, 0);

     _gen_majoP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
      new ( (*_gen_majoP4)[i] ) TLorentzVector();
    }
    outputTree->Branch("_gen_majoP4", "TClonesArray", &_gen_majoP4, 32000, 0);

   
    _gen_wP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
      new ( (*_gen_wP4)[i] ) TLorentzVector();
    }
    outputTree->Branch("_gen_wP4", "TClonesArray", &_gen_wP4, 32000, 0);

     // ===============================================
    _jetP4 = new TClonesArray("TLorentzVector", 100);
    for (int i=0; i!=100; ++i) {
      new ( (*_jetP4)[i] ) TLorentzVector();
    }
    //outputTree->Branch("_jetP4", "TClonesArray", &_jetP4, 32000, 0);
    
    _jetAllP4 = new TClonesArray("TLorentzVector", 100);
    for (int i=0; i!=100; ++i) {
        new ( (*_jetAllP4)[i] ) TLorentzVector();
    }
    //outputTree->Branch("_jetAllP4", "TClonesArray", &_jetAllP4, 32000, 0);
    
    outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
    outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
    outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
    
    outputTree->Branch("_nLeptons", &_nLeptons, "_nLeptons/I");

    outputTree->Branch("_lPt", &_lPt, "_lPt[_nLeptons]/D");
    outputTree->Branch("_lEta", &_lEta, "_lEta[_nLeptons]/D");
    outputTree->Branch("_lPhi", &_lPhi, "_lPhi[_nLeptons]/D");
    outputTree->Branch("_lE", &_lE, "_lE[_nLeptons]/D");
    
    outputTree->Branch("_lPtmc", &_lPtmc, "_lPtmc[_nLeptons]/D");
    outputTree->Branch("_lEtamc", &_lEtamc, "_lEtamc[_nLeptons]/D");
    outputTree->Branch("_lPhimc", &_lPhimc, "_lPhimc[_nLeptons]/D");
    outputTree->Branch("_lEmc", &_lEmc, "_lEmc[_nLeptons]/D");
    outputTree->Branch("_lpdgmc", &_lpdgmc, "_lpdgmc[_nLeptons]/D");
    outputTree->Branch("_lchargemc", &_lchargemc, "_lchargemc[_nLeptons]/D");

    
    /*
    outputTree->Branch("_flag_lepton", &_flag_lepton, "_flag_lepton[6]/O");
    outputTree->Branch("_dR15", &_dR15 , "_dR15[6]/D");
    */
   
    outputTree->Branch("_nuPtmc", &_nuPtmc, "_nuPtmc[_nLeptons]/D");
    outputTree->Branch("_nuEtamc", &_nuEtamc, "_nuEtamc[_nLeptons]/D");
    outputTree->Branch("_nuPhimc", &_nuPhimc, "_nuPhimc[_nLeptons]/D");
    outputTree->Branch("_nuEmc", &_nuEmc, "_nuEmc[_nLeptons]/D");

    outputTree->Branch("_mtmc", &_mtmc, "_mtmc[_nLeptons]/D");
    
    outputTree->Branch("_nZboson", &_nZboson, "_nZboson/I");


    outputTree->Branch("_nEle", &_nEle, "_nEle/I");
    outputTree->Branch("_nMu", &_nMu, "_nMu/I");
    outputTree->Branch("_nTau", &_nTau, "_nTau/I");

    outputTree->Branch("_flavors", &_flavors, "_flavors[_nLeptons]/I");
    outputTree->Branch("_charges", &_charges, "_charges[_nLeptons]/D");
    outputTree->Branch("_indeces", &_indeces, "_indeces[_nLeptons]/I");
    outputTree->Branch("_isolation", &_isolation, "_isolation[_nLeptons]/D");
    outputTree->Branch("_isolationDB", &_isolationDB, "_isolationDB[_nLeptons]/D");
    //outputTree->Branch("_isolationComponents", &_isolationComponents, "_isolationComponents[4][6]/D");
    //outputTree->Branch("_isolationMC", &_isolationMC, "_isolationMC[4][6]/D");
    outputTree->Branch("_miniisolation", &_miniisolation, "_miniisolation[_nLeptons][2]/D");
    outputTree->Branch("_miniisolationCharged", &_miniisolationCharged, "_miniisolationCharged[_nLeptons][2]/D");
    outputTree->Branch("_multiisolation", &_multiisolation, "_multiisolation[_nLeptons][5]/O");
    outputTree->Branch("_ptrel", &_ptrel, "_ptrel[_nLeptons]/D");
    outputTree->Branch("_ptratio", &_ptratio, "_ptratio[_nLeptons]/D");
    outputTree->Branch("_muonSegmentComp", &_muonSegmentComp, "_muonSegmentComp[_nLeptons]/D");

    outputTree->Branch("_mll", &_mll, "_mll[3]/D");
    //outputTree->Branch("_ossf", &_ossf, "_ossf[3]/O");
    
    outputTree->Branch("_origin", &_origin, "_origin[_nLeptons]/I");
    outputTree->Branch("_originReduced", &_originReduced, "_originReduced[_nLeptons]/I");

    outputTree->Branch("_isPromptFinalState", &_isPromptFinalState, "_isPromptFinalState[_nLeptons]/O");
    outputTree->Branch("_fromHardProcessFinalState", &_fromHardProcessFinalState, "_fromHardProcessFinalState[_nLeptons]/O");
    
    outputTree->Branch("_PVchi2", &_PVchi2, "_PVchi2/D");
    outputTree->Branch("_PVerr", &_PVerr, "_PVerr[3]/D");
    
    outputTree->Branch("_ipPV", &_ipPV, "_ipPV[_nLeptons]/D");
    outputTree->Branch("_ipPVerr", &_ipPVerr, "_ipPVerr[_nLeptons]/D");
    outputTree->Branch("_ipZPV", &_ipZPV, "_ipZPV[_nLeptons]/D");
    outputTree->Branch("_ipZPVerr", &_ipZPVerr, "_ipZPVerr[_nLeptons]/D");
    
    outputTree->Branch("_ipPVmc", &_ipPVmc, "_ipPVmc[_nLeptons]/D");
    
    outputTree->Branch("_3dIP", &_3dIP, "_3dIP[_nLeptons]/D");
    outputTree->Branch("_3dIPerr", &_3dIPerr, "_3dIPerr[_nLeptons]/D");
    outputTree->Branch("_3dIPsig", &_3dIPsig, "_3dIPsig[_nLeptons]/D");
    
    
    outputTree->Branch("_mt", &_mt, "_mt[_nLeptons]/D");
    outputTree->Branch("_isloose", &_isloose, "_isloose[_nLeptons]/O");
    outputTree->Branch("_istight", &_istight, "_istight[_nLeptons]/O");
    outputTree->Branch("_istightID", &_istightID, "_istightID[_nLeptons]/O");
 
    outputTree->Branch("_trigEmulator", &_trigEmulator, "_trigEmulator[_nLeptons]/O");
    outputTree->Branch("_isotrigEmulator", &_isotrigEmulator, "_isotrigEmulator[_nLeptons]/O");
    
    outputTree->Branch("_closeJetPtAll", &_closeJetPtAll, "_closeJetPtAll[_nLeptons]/D");
    outputTree->Branch("_closeJetAngAll", &_closeJetAngAll, "_closeJetAngAll[_nLeptons]/D");
 
    outputTree->Branch("_trackSelectionMultiplicity", &_trackSelectionMultiplicity, "_trackSelectionMultiplicity[_nLeptons]/I");
    outputTree->Branch("_closeJetCSVAll", &_closeJetCSVAll, "_closeJetCSVAll[_nLeptons]/D");

    outputTree->Branch("_chargeConst", &_chargeConst, "_chargeConst[_nLeptons]/O");
    outputTree->Branch("_hitsNumber", &_hitsNumber, "_hitsNumber[_nLeptons]/I");
    outputTree->Branch("_vtxFitConversion", &_vtxFitConversion, "_vtxFitConversion[_nLeptons]/O");

    outputTree->Branch("_mvaValue", &_mvaValue, "_mvaValue[_nLeptons]/D");

    outputTree->Branch("_decayModeFinding", &_decayModeFinding, "_decayModeFinding[_nLeptons]/O");
    outputTree->Branch("_looseMVA_dR03", &_looseMVA_dR03, "_looseMVA_dR03[_nLeptons]/O");
    outputTree->Branch("_mediumMVA_dR03", &_mediumMVA_dR03, "_mediumMVA_dR03[_nLeptons]/O");

    //outputTree->Branch("_decayModeFindingOldDMs", &_decayModeFindingOldDMs, "_decayModeFindingOldDMs[_nLeptons]/O");
    outputTree->Branch("_vlooseMVAold", &_vlooseMVAold, "_vlooseMVAold[_nLeptons]/O");
    outputTree->Branch("_looseMVAold", &_looseMVAold, "_looseMVAold[_nLeptons]/O");
    outputTree->Branch("_mediumMVAold", &_mediumMVAold, "_mediumMVAold[_nLeptons]/O");
    outputTree->Branch("_tightMVAold", &_tightMVAold, "_tightMVAold[_nLeptons]/O");
    outputTree->Branch("_vtightMVAold", &_vtightMVAold, "_vtightMVAold[_nLeptons]/O");
    
    outputTree->Branch("_decayModeFindingNewDMs", &_decayModeFindingNewDMs, "_decayModeFindingNewDMs[_nLeptons]/O");
    outputTree->Branch("_vlooseMVAnew", &_vlooseMVAnew, "_vlooseMVAnew[_nLeptons]/O");
    outputTree->Branch("_looseMVAnew", &_looseMVAnew, "_looseMVAnew[_nLeptons]/O");
    outputTree->Branch("_mediumMVAnew", &_mediumMVAnew, "_mediumMVAnew[_nLeptons]/O");
    outputTree->Branch("_tightMVAnew", &_tightMVAnew, "_tightMVAnew[_nLeptons]/O");
    outputTree->Branch("_vtightMVAnew", &_vtightMVAnew, "_vtightMVAnew[_nLeptons]/O");

 
    outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
    outputTree->Branch("_n_MCTruth_PV", &_n_MCTruth_PV, "_n_MCTruth_PV/D");
    
    outputTree->Branch("_met", &_met, "_met/D");
    outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
    outputTree->Branch("HT", &HT, "HT/D");
    
    outputTree->Branch("_genmet", &_genmet, "_genmet/D");
    outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/D");
    
    outputTree->Branch("_genqpt", &_genqpt, "_genqpt/D");
    outputTree->Branch("_genqpt20", &_genqpt20, "_genqpt20/D");

    outputTree->Branch("_mompt", &_mompt, "_mompt[_nLeptons]/D");
    outputTree->Branch("_momphi", &_momphi, "_momphi[_nLeptons]/D");
    outputTree->Branch("_mometa", &_mometa, "_mometa[_nLeptons]/D");
    outputTree->Branch("_mompdg", &_mompdg, "_mompdg[_nLeptons]/I");
    
    outputTree->Branch("_n_bJets", &_n_bJets, "_n_bJets/I");
    outputTree->Branch("_n_Jets", &_n_Jets, "_n_Jets/I");
    outputTree->Branch("_bTagged", &_bTagged, "_bTagged[_n_Jets]/O");
    outputTree->Branch("_jetEta", &_jetEta, "_jetEta[_n_Jets]/D");
    outputTree->Branch("_jetPhi", &_jetPhi, "_jetPhi[_n_Jets]/D");
    outputTree->Branch("_jetPt", &_jetPt, "_jetPt[_n_Jets]/D");
    outputTree->Branch("_jetFlavour", &_jetFlavour, "_jetFlavour[_n_Jets]/I");
    outputTree->Branch("_jetE", &_jetE, "_jetE[_n_Jets]/D");
    outputTree->Branch("_csv", &_csv, "_csv[_n_Jets]/D");
    outputTree->Branch("_jetDeltaR", &_jetDeltaR, "_jetDeltaR[_n_Jets][_nLeptons]/D");
     outputTree->Branch("_clean", &_clean, "_clean[_n_Jets]/I");

    // JEC
    outputTree->Branch("_jecUnc", &_jecUnc, "_jecUnc[_n_Jets]/D");
    outputTree->Branch("_jetPtUp", &_jetPtUp, "_jetPtUp[_n_Jets]/D");
    outputTree->Branch("_jetPtDown", &_jetPtDown, "_jetPtDown[_n_Jets]/D");

    // JER
    outputTree->Branch("_matchedjetPt", &_matchedjetPt, "_matchedjetPt[_n_Jets]/D");
    outputTree->Branch("_matchedjetEta", &_matchedjetEta, "_matchedjetEta[_n_Jets]/D");
    outputTree->Branch("_matchedjetPhi", &_matchedjetPhi, "_matchedjetPhi[_n_Jets]/D");
    outputTree->Branch("_matchedjetE", &_matchedjetE, "_matchedjetE[_n_Jets]/D");
    outputTree->Branch("_matchedjetM", &_matchedjetM, "_matchedjetM[_n_Jets]/D");
    outputTree->Branch("_matchGjet", &_matchGjet, "_matchGjet[_n_Jets]/O");

    // b-tag SF
    outputTree->Branch("_btagSF", &_btagSF, "_btagSF[19][_n_Jets]/D");


    // theoretical
    outputTree->Branch("_weight", &_weight, "_weight/D");
    outputTree->Branch("_LHEweight", &_LHEweight, "_LHEweight[111]/D");
    outputTree->Branch("_LHEweightID", &_LHEweightID, "_LHEweightID[111]/D");

    outputTree->Branch("_mgluino", &_mgluino, "_mgluino/D");
    outputTree->Branch("_mchi0", &_mchi0, "_mchi0/D");
    
    outputTree->Branch("_nMajorana", &_nMajorana, "_nMajorana/I");
    outputTree->Branch("_passedMVA80", &_passedMVA80, "_passedMVA80/I");
    outputTree->Branch("_passedMVA90", &_passedMVA90, "_passedMVA90/I");
    outputTree->Branch("_findMatched", &_findMatched, "_findMatched/I");

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
      outputTree->Branch("_gen_nL", &_gen_nL, "_gen_nL/I");	    
      outputTree->Branch("_gen_lPt", &_gen_lPt, "_gen_lPt[_gen_nL]/D");
      outputTree->Branch("_gen_lE", &_gen_lE, "_gen_lE[_gen_nL]/D");
      outputTree->Branch("_gen_lEta", &_gen_lEta, "_gen_lEta[_gen_nL]/D");
      outputTree->Branch("_gen_lPhi", &_gen_lPhi, "_gen_lPhi[_gen_nL]/D");
      outputTree->Branch("_gen_lmompdg", &_gen_lmompdg, "_gen_lmompdg[_gen_nL]/I");
      outputTree->Branch("_gen_charges",&_gen_charges, "_gen_charges[_gen_nL]/D");
      outputTree->Branch("_gen_flavors",&_gen_flavors, "_gen_flavors[_gen_nL]/I");
      //Added neutrinos  ////////////////////////////////////////////////////
      outputTree->Branch("_gen_nNu", &_gen_nNu, "_gen_nNu/I");
      outputTree->Branch("_gen_nuPt", &_gen_nuPt, "_gen_nuPt[_gen_nNu]/D");
      outputTree->Branch("_gen_nuE", &_gen_nuE, "_gen_nuE[_gen_nNu]/D");
      outputTree->Branch("_gen_nuEta", &_gen_nuEta, "_gen_nuEta[_gen_nNu]/D");
      outputTree->Branch("_gen_nuPhi", &_gen_nuPhi, "_gen_nuPhi[_gen_nNu]/D");
      outputTree->Branch("_gen_numompdg", &_gen_numompdg, "_gen_numompdg[_gen_nNu]/I");
      //End of neutrino Edit////////////////////////////////////////////////////////
      //Added neutrinos Majorana  ////////////////////////////////////////////////////
      outputTree->Branch("_gen_nMajo", &_gen_nMajo, "_gen_nMajo/I");
      outputTree->Branch("_gen_majoPt", &_gen_majoPt, "_gen_majoPt[_gen_nMajo]/D");
      outputTree->Branch("_gen_majoE", &_gen_majoE, "_gen_majoE[_gen_nMajo]/D");
      outputTree->Branch("_gen_majoEta", &_gen_majoEta, "_gen_majoEta[_gen_nMajo]/D");
      outputTree->Branch("_gen_majoPhi", &_gen_majoPhi, "_gen_majoPhi[_gen_nMajo]/D");
      //Added W  ////////////////////////////////////////////////////
      outputTree->Branch("_gen_nW", &_gen_nW, "_gen_nW/I");
      outputTree->Branch("_gen_wPt", &_gen_wPt, "_gen_wPt[_gen_nW]/D");
      outputTree->Branch("_gen_wE", &_gen_wE, "_gen_wE[_gen_nW]/D");
      outputTree->Branch("_gen_wEta", &_gen_wEta, "_gen_wEta[_gen_nW]/D");
      outputTree->Branch("_gen_wPhi", &_gen_wPhi, "_gen_wPhi[_gen_nW]/D");
      outputTree->Branch("_gen_wmompdg", &_gen_wmompdg, "_gen_wmompdg[_gen_nW]/I");

   }

    
    GPM = GenParticleManager();
    
    if (isData) fMetCorrector = new OnTheFlyCorrections("Majorana/PatAnalyzer/data/", "Spring16_25nsV6_DATA", isData);
    else        fMetCorrector = new OnTheFlyCorrections("Majorana/PatAnalyzer/data/", "Spring16_25nsV6_MC",   isData);
    if (isData) _corrLevel = "L2L3Residual";
    else        _corrLevel = "L3Absolute";

    jecUnc = new JetCorrectionUncertainty(edm::FileInPath("Majorana/PatAnalyzer/data/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt").fullPath());
    
   
    //80 % eff
    looseMVA[0][0] = 0.287435;
    looseMVA[1][0] = 0.221846;
    looseMVA[2][0] = -0.303263;
    looseMVA[3][0] = 0.967083;
    looseMVA[4][0] = 0.929117;
    looseMVA[5][0] = 0.726311;
    //90 % eff
    looseMVA[0][1] = -0.083313;
    looseMVA[1][1] = -0.235222;
    looseMVA[2][1] = -0.67099;
    looseMVA[3][1] = 0.913286;
    looseMVA[4][1] = 0.805013;
    looseMVA[5][1] = 0.358969;



    //multiisolation very loose WP 
    multiConst[0][0] = 0.25;
    multiConst[0][1] = 0.67;
    multiConst[0][2] = 4.4;
    
    //multiisolation loose WP 
    multiConst[1][0] = 0.2;
    multiConst[1][1] = 0.69;
    multiConst[1][2] = 6.;
    
    //multiisolation medium WP 
    multiConst[2][0] = 0.16;
    multiConst[2][1] = 0.76;
    multiConst[2][2] = 7.2;
    
    //multiisolation tight WP 
    multiConst[3][0] = 0.12;
    multiConst[3][1] = 0.8;
    multiConst[3][2] = 7.2;
    
    //multiisolation very tight WP 
    multiConst[4][0] = 0.09;
    multiConst[4][1] = 0.84;
    multiConst[4][2] = 7.2;
    
        
    _nEventsTotal = 0;
    _nEventsFiltered = 0;
    _nEventsTotalCounted = 0;
    
    firstEvent_ = true;
}

void trilepton::endJob() {
   
    std::cout<<_nEventsTotal<<std::endl;
    std::cout<<_nEventsFiltered<<std::endl;
    
    delete fMetCorrector;
    
    
}

void trilepton::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
    using namespace edm;
    
    return;
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
    std::cout << "Available triggers:" << std::endl;
    for (unsigned int i = 0; i < trigResults->size(); ++i){
      std::cout << "  " << triggerNames.triggerName(i) << std::endl;
      for(TString triggerName : toSave){
	if(TString(triggerNames.triggerName(i)).Contains(triggerName)){
	  triggerIndices[triggerName] = i;
	}
      }
    }
    std::cout << std::endl;
  }

  // Save our triggers/flags
  for(TString triggerName : toSave){
    triggerFlags[triggerName] = trigResults->accept(triggerIndices[triggerName]);
    if(isHLT){
      triggerPrescales[triggerName] = prescales->getPrescaleForIndex(triggerIndices[triggerName]);
    }
  }
}



void trilepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{
    //bool islepton;

    _nZboson = 0;

    if(not isData) {
        //******************************************************************************************************************
        // Gen level particles                  ****************************************************************************
        //******************************************************************************************************************
        //iEvent.getByLabel("packedGenParticles", TheGenParticles);
        iEvent.getByToken(genparticleToken, TheGenParticles);
        std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW, vGenMajorana;
        if( TheGenParticles.isValid() )
        {
            GPM.SetCollection(TheGenParticles);
            GPM.Classify();
            vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),1);
            vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),1);
            vGenNPMuons = GPM.filterByStatus(GPM.getNonPromptMuons(),1);
            vGenNPElectrons = GPM.filterByStatus(GPM.getNonPromptElectrons(),1);
	    vGenMajorana = GPM.filterByStatus(GPM.getMajorana(),1);
            //std::cout<<"*************"<<std::endl;
            
            TLorentzVector Gen0;
            Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
            _genqpt = 0;
            _genqpt20 = 0;
            cout << "Particle ids in the event" << endl;
            for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ )
            {
                int id = TMath::Abs(p->pdgId());
                /*
                if(p->status() == 1){
                    const GenParticle *mom = GPM.getMother(&*p);
                    cout << id << " " << p->pt() << " " << TMath::Abs(mom->pdgId()) << endl;
                }
                */
		if(id == 9900012) _nMajorana++;
                if(id == 23)
                    _nZboson++;
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
	      if(p->status() == 0) continue;								//status 0 contains no useful information and should be skipped
	      unsigned id = TMath::Abs(p->pdgId());
	      const reco::GenStatusFlags StatusFlags;
	      //MCTruthHelper::fillGenStatusFlags(*p,StatusFlags);
	      //if( id == 11 || id ==13 || id == 15){
	      //if(p->status() == 1 || (id == 15 && p->status() == 2)){				//Tau's always decay!
	      if(id == 11 || id == 13){
		if(p->status() == 1){
		  if (nL > 29) continue;
		  _gen_lPt[nL] = p->pt();
		  _gen_lEta[nL] = p->eta();
		  _gen_lPhi[nL] = p->phi();
		  _gen_lE[nL] = p->energy();
		  const GenParticle *mom = GPM.getMother(&*p);				//Use getmother or getmotherparton here?
		  while( fabs(mom->pdgId()) == id && GPM.getMother(mom) != nullptr)
		    mom = GPM.getMother(mom);
		  if(mom != nullptr){
		    _gen_lmompdg[nL] = mom->pdgId();
		  }
		  else{
		    _gen_lmompdg[nL] = 0;  //Particles for which mompdg is 0 have no mother.
		  }
		  if(id == 11) {_gen_flavors[nL] = 0;}
		  else if(id ==13){_gen_flavors[nL] = 1;}
		  //else if(id ==15) {_gen_flavors[nL] = 2;}
		  _gen_charges[nL] = p->charge();


		  ((TLorentzVector *)_gen_leptonP4->At(nL))->SetPtEtaPhiE(p->pt(), p->eta(), p->phi(), p->energy());
		  _gen_lPt[nL] = ((TLorentzVector *)_gen_leptonP4->At(nL))->Pt();
		  _gen_lEta[nL] = ((TLorentzVector *)_gen_leptonP4->At(nL))->Eta();
		  _gen_lPhi[nL] = ((TLorentzVector *)_gen_leptonP4->At(nL))->Phi();
		  _gen_lE[nL] = ((TLorentzVector *)_gen_leptonP4->At(nL))->E();
		  ++nL;
		}
	      }	      
	      else if(id == 12 || id == 14 || id == 16){
		if(p->status() == 1){
		  if (nNu > 29) continue;
		  _gen_nuPt[nNu] = p->pt();
		  _gen_nuE[nNu] = p->energy();
		  _gen_nuEta[nNu] = p->eta();
		  _gen_nuPhi[nNu] = p->phi();
		  ((TLorentzVector *)_gen_nuP4->At(nNu))->SetPtEtaPhiE(p->pt(), p->eta(), p->phi(), p->energy());
		  //_gen_lPt[nNu] = ((TLorentzVector *)_gen_nuP4->At(nNu))->Pt();
		  //_gen_lEta[nNu] = ((TLorentzVector *)_gen_nuP4->At(nNu))->Eta();
		  //_gen_lPhi[nNu] = ((TLorentzVector *)_gen_nuP4->At(nNu))->Phi();
		  //_gen_lE[nNu] = ((TLorentzVector *)_gen_nuP4->At(nNu))->E();
		  const GenParticle *mom = GPM.getMother(&*p);				//Use getmother or getmotherparton here?
		  while( fabs(mom->pdgId()) == id && GPM.getMother(mom) != nullptr)
		    mom = GPM.getMother(mom);
		  if(mom != nullptr){
		    _gen_numompdg[nNu] = mom->pdgId();
		  }
		  else{
		    _gen_numompdg[nNu] = 0;  //Particles for which mompdg is 0 have no mother.
		  }
		 

		  
		  ++nNu;
		}
	      }
	      else if(  id == 24 ){
		//cout<<"w status:::::::::::::::::::::::::::::::::::::::: "<<p->status()<<endl;
		if(p->status() == 22 ||p->status() == 52){
		  if (nw > 29) continue;
		  _gen_wPt[nw] = p->pt();
		  _gen_wE[nw] = p->energy();
		  _gen_wEta[nw] = p->eta();
		  _gen_wPhi[nw] = p->phi();
		  if ( p->pt() == 0 &&  p->phi() == 0 && p->eta()== 0 ) continue;
		  //if (p->eta() == TMath::Power(10,10) || p->eta() == - TMath::Power(10,10)) continue;
		  //((TLorentzVector *)_gen_wP4->At(nw))->SetPtEtaPhiE(p->pt(), _gen_wEta[nw], p->phi(), p->energy());
		  //_gen_lPt[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->Pt();
		  //_gen_lEta[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->Eta();
		  //_gen_lPhi[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->Phi();
		  //_gen_lE[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->E();
		  const GenParticle *mom = GPM.getMother(&*p);				//Use getmother or getmotherparton here?
		  while( fabs(mom->pdgId()) == id && GPM.getMother(mom) != nullptr)
		    mom = GPM.getMother(mom);
		  if(mom != nullptr){
		    _gen_wmompdg[nw] = mom->pdgId();
		  }
		  else{
		    _gen_wmompdg[nw] = 0;  //Particles for which mompdg is 0 have no mother.
		  }
		  ((TLorentzVector *)_gen_wP4->At(nw))->SetPtEtaPhiE(p->pt(), p->eta(), p->phi(), p->energy());
		  //_gen_lPt[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->Pt();
		  //_gen_lEta[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->Eta();
		 //_gen_lPhi[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->Phi();
		  //_gen_lE[nw] = ((TLorentzVector *)_gen_wP4->At(nw))->E();
		  ++nw;
		}
		}




	      else if(id == 9900012){
		//cout<<"N status:::::::::::::::::::::::::::::::::::::::: "<<p->status()<<endl;
		//if(p->status() == 1){
		  cout<<"ppppppppppppppppppp========================>>>>>"<<endl;
		  _gen_majoPt[nMajo] = p->pt();
		  _gen_majoE[nMajo] = p->energy();
		  _gen_majoEta[nMajo] = p->eta();
		  _gen_majoPhi[nMajo] = p->phi();
		  ((TLorentzVector *)_gen_majoP4->At(nMajo))->SetPtEtaPhiE(p->pt(), p->eta(), p->phi(), p->energy());
		  _gen_majoPt[nMajo] = ((TLorentzVector *)_gen_majoP4->At(nMajo))->Pt();
		  _gen_majoEta[nMajo] = ((TLorentzVector *)_gen_majoP4->At(nMajo))->Eta();
		  _gen_majoPhi[nMajo] = ((TLorentzVector *)_gen_majoP4->At(nMajo))->Phi();
		  _gen_majoE[nMajo] = ((TLorentzVector *)_gen_majoP4->At(nMajo))->E();
		  ++nMajo;
		  //}
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

    Double_t weight = 1.;
    if(not isData) {
        edm::Handle<GenEventInfoProduct> pdfvariables;
        iEvent.getByToken(pdfvariablesToken, pdfvariables);
        weight = pdfvariables->weight();
        _weight=weight;
    }

    _runNb = iEvent.id().run();
    _eventNb = iEvent.id().event();
    _lumiBlock = iEvent.luminosityBlock();
   
   
    //if (_eventNb != 44802315) return;
    //if ((_eventNb != 17707) || (_runNb != 1) || (_lumiBlock != 178)) return;
    
    std::cout<<"EVENT "<<_runNb << " " << _lumiBlock << " " << _eventNb<<std::endl;
  
    
    //============ Total number of events is the sum of the events ============
    //============ in each of these luminosity blocks ============
    _nEventsTotalCounted++;
    _hCounter->Fill(0., weight);
    
   
    
    //============ Beamspot ============
    edm::Handle< reco::BeamSpot > theBeamSpot;
    //iEvent.getByLabel( IT_beamspot, theBeamSpot );
    iEvent.getByToken( IT_beamspot, theBeamSpot );
    //if( ! theBeamSpot.isValid() ) ERR( IT_beamspot ) ;
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
    //if( ! theVertices.isValid() ) ERR(IT_goodVtx ) ;
    int nvertex = theVertices->size();
    
    _n_PV = nvertex;
    Nvtx->Fill(TMath::Min(nvertex,39));
    if(! nvertex ){
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
    
    //==================================
    
    //============ PF cand ============
    edm::Handle<pat::PackedCandidateCollection> pfcands;
    iEvent.getByToken(packedPFCandidatesToken, pfcands);
    //==================================
    
    //============ Pat Muons ============
    edm::Handle< std::vector<pat::Muon> > thePatMuons;
    iEvent.getByToken( IT_muon, thePatMuons );
    //if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    //==================================
    
   
    
    //============ Pat Electrons ============
    edm::Handle< edm::View<pat::Electron> > thePatElectrons;
    iEvent.getByToken( IT_electron, thePatElectrons );
    //if( ! thePatElectrons.isValid() ) ERR( IT_electron );
    //==================================
   
    // Get MVA values and categories (optional)
    edm::Handle<edm::ValueMap<float> > mvaValues;
    iEvent.getByToken(mvaValuesMapToken,mvaValues);


    //============ Conversions ============
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    //iEvent.getByLabel("reducedEgamma","reducedConversions", theConversions);
    iEvent.getByToken(reducedEgammaToken, theConversions);
    //==================================
    
    //============ Pat Jets ============
    edm::Handle< std::vector< pat::Jet> > thePatJets;
    iEvent.getByToken(IT_jet , thePatJets );
    //if( ! thePatJets.isValid() ) ERR(IT_jet);
    //==================================
    edm::Handle<double> rhoJets;
    iEvent.getByToken(fixedGridRhoFastjetCentralNeutralToken , rhoJets);//kt6PFJets
    myRhoJets = *rhoJets;
    //==================================
    edm::Handle<double> rhoJECJets;
    iEvent.getByToken(fixedGridRhoFastjetAllToken , rhoJECJets);//kt6PFJets
    myRhoJECJets = *rhoJECJets;


    //============ Pat MET ============
    
    edm::Handle< vector<pat::MET> > ThePFMET;
    iEvent.getByToken(IT_pfmet, ThePFMET);
    //if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
    const vector<pat::MET> *pfmetcol = ThePFMET.product();
    const pat::MET *pfmet;
    pfmet = &(pfmetcol->front());
    // Old MET implementation, 2 FEB 2016
    _met = pfmet->pt();
    _met_phi = pfmet->phi();

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

    
    for( std::vector<pat::Jet>::const_iterator jet = (*thePatJets).begin(); jet != (*thePatJets).end(); jet++ ) {
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
                                                                      (&*jet)->jetArea());
        corrMetx += corr.first ;
        corrMety += corr.second;
    }
    double newmet    = sqrt(corrMetx*corrMetx + corrMety*corrMety);
    double newmetphi = atan2(corrMety, corrMetx);

    //std::cout<< "New and old met: " << newmet<<" "<<_met<<"; "<<newmetphi<<" "<<_met_phi<<" "<<std::endl;
    _met = newmet;
    _met_phi = newmetphi;
    std::cout<< "met: " <<_met<<"; "<<_met_phi<<std::endl;


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
    
    
    //std::vector<const pat::Muon* > sMu = ssbMediumMuonSelector( *thePatMuons, _minPt0, PV, _looseD0Mu);
    std::vector<const pat::Muon* > sMu = ssbLooseMuonSelector( *thePatMuons, _minPt0, PV, _looseD0Mu);
    
    //std::vector<const pat::Electron* > sEl = ssbMVAElectronSelector( *thePatElectrons, _minPt1, PV, _looseD0E, _chargeConsistency, theConversions, BS, true);

    //Taus
    edm::Handle<pat::TauCollection> PatTaus;
    iEvent.getByToken( IT_tau, PatTaus );
    std::vector<const pat::Tau* > sTau;
    for (const pat::Tau &tau : *PatTaus) {
        
        if (tau.pt() < _tauPt) continue;
        if (fabs(tau.eta()) > _tauEta) continue;
  sTau.push_back(&tau);
    }
    
    SelectedJetsAll = JetSelectorAll(*thePatJets, 5., 3.0);
    std::vector<const pat::Jet* > SelectedJets = JetSelector(*thePatJets, _jetPtCut, _jetEtaCut);

    HT = 0.;
    std::vector< const pat::Jet* > Bjets;
    
    _n_bJetsAll = 0;
    
    int n_bJetsAll30 = 0;
    
    
    _n_Jets = 0;
    _n_bJets = 0;
    
    for(unsigned int i = 0 ; i < SelectedJetsAll.size() ;i++ ){
        
        //double uncPt = (SelectedJetsAll[i]->correctedP4("Uncorrected")).Pt();
        double uncEta = (SelectedJetsAll[i]->correctedP4("Uncorrected")).Eta();
        double uncPhi = (SelectedJetsAll[i]->correctedP4("Uncorrected")).Phi();
        
        //double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[i]->jetArea(),_corrLevel);
        
        _jetEtaAll[i] = uncEta;
        _jetPhiAll[i] = uncPhi;
        _jetPtAll[i] = SelectedJetsAll[i]->pt();
        _jetEAll[i] = SelectedJetsAll[i]->energy();
        
        ((TLorentzVector *)_jetAllP4->At(i))->SetPtEtaPhiM( _jetPtAll[i], _jetEtaAll[i], _jetPhiAll[i], _jetEAll[i]);
        
        _csvAll[i] = SelectedJetsAll[i]->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
        
        if(SelectedJetsAll[i]->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") > 0.814) {
            _bTaggedAll[i] = true;
            _n_bJetsAll++;
            if (_jetPtAll[i] > _jetPtCut) {
                //bIndex[n_bJetsAll30] = i;
                n_bJetsAll30++;
            }
        } else _bTaggedAll[i] = false;
        
       
        
    }
    _n_JetsAll = SelectedJetsAll.size();
    
    int leptonCounter = 0;
    std::cout << "Muon selection started " << std::endl;
    for(unsigned int i = 0 ; i < sMu.size() ;i++ ){
        
        
        const pat::Muon *iM = sMu[i];
        //if (iM->pt() < 10) continue;
        
        //_hitsNumber[i] = 0;

        if (leptonCounter == 10) continue;

        _flavors[leptonCounter] = 1;
        _charges[leptonCounter] = iM->charge();
        _isolation[leptonCounter] = pfRelIso(iM,myRhoJECJets);
        _isolationDB[leptonCounter] = pfRelIso(iM);
        
        _miniisolation[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false, myRhoJets);
        _miniisolation[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, false, myRhoJets);

        _miniisolationCharged[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, true, myRhoJets);
        _miniisolationCharged[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, true, myRhoJets);

        _ipPV[leptonCounter] = iM->innerTrack()->dxy(PV);
        _ipPVerr[leptonCounter] = iM->innerTrack()->dxyError();
        
        _ipZPV[leptonCounter] = iM->innerTrack()->dz(PV);
        _ipZPVerr[leptonCounter] = iM->innerTrack()->dzError();
        
        
        //bool good = iM->innerTrack()->validFraction() > 0.8 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
        //bool good = iM->innerTrack()->validFraction() > 0.49 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451) && iM->isMediumMuon() && iM->isPFMuon() && (iM->isTrackerMuon() || iM->isGlobalMuon() ) && _miniisolation[leptonCounter][0] < 0.4; // short term

        //bool good = iM->innerTrack()->validFraction() > 0.8 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
        //bool good = iM->innerTrack()->validFraction() > 0.49 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451); // short term 



	//bool goodGlb = iM->isGlobalMuon() && iM->globalTrack()->normalizedChi2() < 3 && iM->combinedQuality().chi2LocalPosition < 12 && iM->combinedQuality().trkKink < 20;

	//	bool goodGlobaltest = iM->isPFMuon() && (iM->isTrackerMuon() || iM->isGlobalMuon() ) && (iM->globalTrack()->normalizedChi2() < 3) &&  iM->combinedQuality().chi2LocalPosition < 12 &&  iM->combinedQuality().trkKink < 20; 

	bool goodGlb = iM->isGlobalMuon() && iM->globalTrack()->normalizedChi2() < 3 && iM->combinedQuality().chi2LocalPosition < 12 && iM->combinedQuality().trkKink < 20;

	bool isMedium =iM->isPFMuon() && (iM->isTrackerMuon() || iM->isGlobalMuon() )&& iM->innerTrack()->validFraction() > 0.8 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);

	bool isotest=  _miniisolation[leptonCounter][0] < 0.4;


        //double _muonSegmentComp[leptonCounter] = iM->segmentCompatibility();
        _muonSegmentComp[leptonCounter] = iM->segmentCompatibility();
        
        //_istightID[leptonCounter] = good;
        
        _3dIP[leptonCounter]    = iM->dB(pat::Muon::PV3D);
        _3dIPerr[leptonCounter] = iM->edB(pat::Muon::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);
        

        //if(_miniisolation[leptonCounter][0] > 0.4) continue;
        if(_isolationDB[leptonCounter] > 1.0 ) continue;
        _isloose[leptonCounter] = goodGlb && (_3dIPsig[leptonCounter] < 4) && isMedium && isotest;
        _istight[leptonCounter] = goodGlb && (_isolationDB[i] < 0.15) && (_3dIPsig[leptonCounter] < 4)&& isMedium && isotest;


        if(!_isloose[leptonCounter]) continue;
        
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iM->pt(), iM->eta(), iM->phi(), iM->energy());
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
        std::cout<< "Muon loose selection passed -> " << _lPt[leptonCounter] << endl;

        fillCloseJetVars(leptonCounter, PV);
        
        if (not isData) {
            const GenParticle* mc = GPM.matchedMC(iM);
            if ( mc!=0 ) {
                fillMCVars(mc, leptonCounter);
                _ipPVmc[leptonCounter] = TMath::Abs(iM->innerTrack()->dxy(PVmc));
            }
            else {
                //std::cout<<"No match mu"<<std::endl;
                _origin[leptonCounter] = -1;
                _originReduced[leptonCounter] = -1;
                
                _mompt[leptonCounter] = 0;
                _momphi[leptonCounter] = 0;
                _mometa[leptonCounter] = 0;
                _mompdg[leptonCounter] = 0;
            }
            
        }
        
        
        leptonCounter++;

        
    }
    _nMu = leptonCounter;
    
    std::cout << "Electron selection started" << std::endl;
    for (size_t i = 0; i < thePatElectrons->size(); ++i){
      _passedMVA80[leptonCounter ] =0;
      _passedMVA90[leptonCounter ] =0;
      _findMatched[leptonCounter ] =-1;
      if (leptonCounter == 10) continue;

      const auto iE = thePatElectrons->ptrAt(i);
        std::cout << "Pat electron: " << iE->pt() << " " << iE->superCluster()->eta() << " " << iE->phi() << " " << std::endl;
        if (!ssbMVAElectronSelectorPassed(iE, _minPt1, PV, _looseD0E, _chargeConsistency, theConversions, BS, false)) continue;
        std::cout << "Passed electron selector" << std::endl;
        _mvaValue[leptonCounter] = (*mvaValues)[iE];

	bool passedMVA90 = false;
	if(iE->pt() <10 ){//selection 5<pt<10
	  passedMVA90 = false;
	  if (TMath::Abs(iE->eta()) < 0.8 ) {
	    passedMVA90 = _mvaValue[leptonCounter]> looseMVA[0][1];
	  } else if (TMath::Abs(iE->eta()) < 1.479 ) {
	    passedMVA90 = _mvaValue[leptonCounter]> looseMVA[1][1];
	  } else {
	    passedMVA90 = _mvaValue[leptonCounter]> looseMVA[2][1];
	  }
	}//end //selection 5<pt<10

	if (iE->pt()>=10  ){
	  passedMVA90 = false;
	  if (TMath::Abs(iE->eta()) < 0.8 ) {
	    passedMVA90 = _mvaValue[leptonCounter] > looseMVA[3][1];
	  } else if (TMath::Abs(iE->eta()) < 1.479 ) {
	    passedMVA90 = _mvaValue[leptonCounter] > looseMVA[4][1];
	  } else {
	    passedMVA90 = _mvaValue[leptonCounter] > looseMVA[5][1];
	  }
	}// end 10 pt

	bool passedMVA80 = false;
	if(iE->pt() <10 ){//selection 5<pt<10
	  passedMVA80 = false;
	  if (TMath::Abs(iE->eta()) < 0.8 ) {
	    passedMVA80 = _mvaValue[leptonCounter]> looseMVA[0][0];
	  } else if (TMath::Abs(iE->eta()) < 1.479 ) {
	    passedMVA80 = _mvaValue[leptonCounter]> looseMVA[1][0];
	  } else {
	    passedMVA80 = _mvaValue[leptonCounter]> looseMVA[2][0];
	  }
	}//end //selection 5<pt<10

	if (iE->pt()>=10){
	  passedMVA80 = false;
	  if (TMath::Abs(iE->eta()) < 0.8 ) {
	    passedMVA80 = _mvaValue[leptonCounter] > looseMVA[3][0];
	  } else if (TMath::Abs(iE->eta()) < 1.479 ) {
	    passedMVA80 = _mvaValue[leptonCounter] > looseMVA[4][0];
	  } else {
	    passedMVA80 = _mvaValue[leptonCounter] > looseMVA[5][0];
	  }
	}// end 10 pt




      
        _flavors[leptonCounter] = 0;
        _charges[leptonCounter] = iE->charge();
        _isolation[leptonCounter] = pfRelIso(&*iE, myRhoJECJets);

        cout << "Isolation value of electron: " << _isolation[leptonCounter] << endl;
        
        _ipPV[leptonCounter] = iE->gsfTrack()->dxy(PV);
        _ipZPV[leptonCounter] = iE->gsfTrack()->dz(PV);

        if(TMath::Abs(_ipPV[leptonCounter]) > 0.05) continue;
        if(TMath::Abs(_ipZPV[leptonCounter]) > 0.1) continue;
        
        _miniisolation[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, false, myRhoJets);
        _miniisolation[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.3, 10., false, false, myRhoJets);

        _miniisolationCharged[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, true, myRhoJets);
        _miniisolationCharged[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.3, 10., false, true, myRhoJets);

        _chargeConst[leptonCounter] = iE->isGsfCtfScPixChargeConsistent();
        _vtxFitConversion[leptonCounter] = ConversionTools::hasMatchedConversion(reco::GsfElectron (*&*iE), theConversions, BS);
        _hitsNumber[leptonCounter] = iE->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

        _trigEmulator[leptonCounter] = triggerEmulatorReturned(&*iE);
        _isotrigEmulator[leptonCounter] = isoTriggerEmulator(&*iE);

        _3dIP[leptonCounter]    = iE->dB(pat::Electron::PV3D);
        _3dIPerr[leptonCounter] = iE->edB(pat::Electron::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);

	_passedMVA80[leptonCounter ] = passedMVA80;
	_passedMVA90[leptonCounter ] = passedMVA90;

	
        std::cout << "Before Loose selection passed : " << iE->pt() << " " << _trigEmulator[leptonCounter] << " " << (_miniisolation[leptonCounter][0] < 0.4) << " "  << (_hitsNumber[leptonCounter] == 0)  << " " << (!_vtxFitConversion[leptonCounter]) << " " << (_3dIPsig[leptonCounter] < 4) << endl;

       
        if(_isolation[leptonCounter] > 1) continue;
        if(_3dIPsig[leptonCounter] > 4) continue;

	if (!_passedMVA90[leptonCounter]) continue;

        
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iE->pt(), iE->eta(), iE->phi(), iE->energy());
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();

       
        fillCloseJetVars(leptonCounter, PV);
        
        if (not isData) {
	  _findMatched[leptonCounter]=-1;
            const GenParticle* mc = GPM.matchedMC(&*iE);
            if ( mc!=0 ) {
                fillMCVars(mc, leptonCounter);
                //Vertex::Point PVmc = mcMom->vertex();
                _ipPVmc[leptonCounter] = TMath::Abs(iE->gsfTrack()->dxy(PVmc));
		_findMatched[leptonCounter]=1;
            }
            else {
                //std::cout<<"No match e"<<std::endl;
                _origin[leptonCounter] = -1;
                _originReduced[leptonCounter] = -1;
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

   
    if (_nEle + _nMu < 3) return;
      
    _nLeptons = leptonCounter;
   
    
    _sb = false;
   
    _n_Jets = 0;
    _n_bJets = 0;
    HT = 0;
    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){
        

        double uncPt = (SelectedJets[i]->correctedP4("Uncorrected")).Pt();
        double uncEta = (SelectedJets[i]->correctedP4("Uncorrected")).Eta();
         
        //double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),_corrLevel);
        double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJECJets, SelectedJets[i]->jetArea(),_corrLevel);
        //std::cout << "Before correction pt of jet: " << uncPt << " " << uncPt*corr << std::endl;
                                 
        if (uncPt*corr < 30) continue;
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
                double dR1 = ((TLorentzVector *)_leptonP4->At(k))->DeltaR( jt );
                //std::cout << "jet cleaning: " << dR1 << " " << ((TLorentzVector *)_leptonP4->At(k))->Pt()  << std::endl;
                clean = clean && (dR1 > 0.4) ;
            }
        }
	_clean[_n_Jets]=0;
	if (clean) _clean[_n_Jets] =1;
	if (!clean) _clean[_n_Jets] =-1;
        

        for(int j=0; j != _nLeptons; ++j){
            _jetDeltaR[_n_Jets][j] = ((TLorentzVector *)_leptonP4->At(j))->DeltaR( jt ) ;
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
        
        if (_jetPt[_n_Jets] > 30)
	  HT+= _jetPt[_n_Jets];
        _n_Jets++;
    }

    std::cout<<_runNb<<" "<<_lumiBlock<<" "<<_eventNb<<" "<<_nMu<<" "<<_nEle<<" "<<_nTau<<" "<<_n_Jets<<" "<<_n_bJets<<std::endl;
    outputTree->Fill();
}

void trilepton::fillMCVars(const GenParticle* mc, const int leptonCounter) {
    
    _lPtmc[leptonCounter] = mc->pt();
    _lEmc[leptonCounter] = mc->energy();
    _lPhimc[leptonCounter] = mc->phi();
    _lEtamc[leptonCounter] = mc->eta();
    _lpdgmc[leptonCounter] = mc->pdgId();
    _lchargemc[leptonCounter] = mc->charge();
    GPM.printInheritance(mc);
    //std::cout << "pdg of tau truth: " << mc->pdgId() << std::endl;
    
    _origin[leptonCounter] = GPM.origin(mc);
    _originReduced[leptonCounter] = GPM.originReduced(_origin[leptonCounter]);
    _isPromptFinalState[leptonCounter] = mc->isPromptFinalState();
    _fromHardProcessFinalState[leptonCounter] = mc->fromHardProcessFinalState();

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
                    /*if (_isolation == 0) {
                     TLorentzVector dauV;
                     dauV.SetPtEtaPhiE(dauts.at(counterD)->pt(), dauts.at(counterD)->eta(), dauts.at(counterD)->phi(), dauts.at(counterD)->energy());
                     double deltaR1 = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR(dauV);
                     std::cout<<counterD<<" "<<dauts.at(counterD)->pdgId()<<" "<<dauts.at(counterD)->pt()<<
                     " "<<dauts.at(counterD)->eta()<<" "<<dauts.at(counterD)->phi()<<" "<<dauts.at(counterD)->energy()<<" in "<<deltaR1<<std::endl;
                     }*/
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
        _mtmc[leptonCounter] = MT_calc(lmc, _nuPtmc[leptonCounter], _nuPhimc[leptonCounter]);
    } else {
        _mompt[leptonCounter] = 0;
        _momphi[leptonCounter] = 0;
        _mometa[leptonCounter] = 0;
        _mompdg[leptonCounter] = 0;
    }
    
    //GPM.printInheritance(&(*mc));
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

        double corr = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJECJets, SelectedJetsAll[k]->jetArea(),"L1FastJet");
        double corr2 = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJECJets, SelectedJetsAll[k]->jetArea(),_corrLevel);
        
       pJet.SetPtEtaPhiE( corr*uncorrPt, SelectedJetsAll[k]->eta(), SelectedJetsAll[k]->phi(), corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E());
       pJet-=*((TLorentzVector *)_leptonP4->At(leptonCounter));
       
       pJet*=corr2/corr;
       pJet+=*((TLorentzVector *)_leptonP4->At(leptonCounter));

              
       double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );
       if (ang < _closeJetAngAll[leptonCounter]) {
           _closeJetAngAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );
           _closeJetPtAll[leptonCounter] = pJet.Pt();
                                  
           //_closeJetEtaAll[leptonCounter] = pJet.Eta();
           //_closeJetPhiAll[leptonCounter] = pJet.Phi();
           //_closeJetEAll[leptonCounter] = pJet.E();
           //_closeJetMAll[leptonCounter] = pJet.M();
           //_closeJetNconstAll[leptonCounter] = SelectedJetsAll[k]->numberOfDaughters();
           _closeJetCSVAll[leptonCounter] = SelectedJetsAll[k]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
       
           _ptrel[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect() - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Vect());

           _ptratio[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()/(_closeJetPtAll[leptonCounter]);

           _closeIndex[leptonCounter] = k;

           int trackSelectionMult = 0;
           cout << "jet closest pt candidate: " << pJet.Pt() << endl;
           //if(_closeIndex[leptonCounter] != -1){
           for(int unsigned it=0; it!=SelectedJetsAll[_closeIndex[leptonCounter]]->numberOfDaughters(); ++it) {
             const pat::PackedCandidate * icand = dynamic_cast<const pat::PackedCandidate *> (SelectedJetsAll[_closeIndex[leptonCounter]]->daughter(it));
             const reco::Track& trk = icand->pseudoTrack();
             Double_t deta = trk.eta() - pJet.Eta();
             Double_t dphi = TVector2::Phi_mpi_pi(trk.phi() - pJet.Phi());
             bool trackSelection =  (trk.pt() > 1) && trk.charge() != 0 && (TMath::Sqrt( deta*deta+dphi*dphi ) < 0.4) && (icand->fromPV() > 1) && (trk.hitPattern().numberOfValidHits()>=8) && (trk.hitPattern().numberOfValidPixelHits()>=2) && (trk.normalizedChi2()<5) && std::fabs(trk.dxy(PV))<0.2 && std::fabs(trk.dz(PV))<17;
            if(trackSelection) trackSelectionMult++;
            /*
            if(trackSelection)
                cout << "track pt: " << trk.pt() << endl;
                */
           }
           cout << "trackMultiplicity: " << trackSelectionMult << endl;
           //}

         _trackSelectionMultiplicity[leptonCounter] = trackSelectionMult;
         }
    }

    if (_closeJetAngAll[leptonCounter] > 0.4) {
        _closeJetPtAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
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

           

    //std::cout << "Info about closest jet -> Pt: " << _closeJetPtAll[leptonCounter] << std::endl;  
                 
    /*
    for (int mi=0; mi!=5; ++mi) {
        if (_miniisolation[leptonCounter][0] < multiConst[mi][0] && ((_ptratio[leptonCounter] > multiConst[mi][1]) || (_ptrel[leptonCounter] > multiConst[mi][2])) )
            _multiisolation[leptonCounter][mi] = true;
        else
            _multiisolation[leptonCounter][mi] = false;
    }
    */
}

void trilepton::matchCloseJet(const int leptonCounter) {
    double minDeltaR3 = 9999;
    double minDeltaR2 = 9999;
    TLorentzVector Gen1, Gen2;
    const GenParticle* mc3a = 0;
    const GenParticle* mc2a = 0;

    if(SelectedJetsAll.size() != 0)
    {
      Gen1.SetPtEtaPhiE(SelectedJetsAll.at(_closeIndex[leptonCounter])->pt(),SelectedJetsAll.at(_closeIndex[leptonCounter])->eta(),SelectedJetsAll.at(_closeIndex[leptonCounter])->phi(),SelectedJetsAll.at(_closeIndex[leptonCounter])->energy());
    

        for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
            int id = TMath::Abs(p->pdgId());
        
            if ((id > 0 && id < 6) || (id == 21) || (id == 22)) {
             if (p->status() != 2) {
                if (fabs(p->eta()) <= 10 )
                    Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                else continue;
                //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
                double deltaRcur = Gen1.DeltaR(Gen2);
                if (deltaRcur < minDeltaR3) {
                    mc3a = &*p;
                    minDeltaR3 = deltaRcur;
                }
                
             } else if (p->status() == 2) {
                if (fabs(p->eta()) <= 10 )
                    Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                else continue;
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
    
    if( TheGenParticles.isValid() )
    {
        //if (_isolation[0] == 0) {
        //    std::cout<<"==="<<std::endl;
        //}
        _isolationMC[leptonCounter][2] = 0;
        _isolationMC[leptonCounter][3] = 0;
        for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
            //int id = TMath::Abs(p->pdgId());
            if (p->status() == 1) {
                TLorentzVector pmc; pmc.SetPtEtaPhiM( p->pt(), p->eta(), p->phi(), p->mass() );
                double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pmc );
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
        double leptpt = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
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
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
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
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta());
    
    
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
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
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
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta());
    
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    
    hJet_SoftLeptptRel = pLep.Perp(pJet.Vect());
    hJet_SoftLeptPt = mu->pt();
    hJet_SoftLeptdR = pLep.DeltaR(pJet);
    
    hJet_SoftLeptIdlooseMu = 1;
    hJet_SoftLeptId95 = 1;
}


DEFINE_FWK_MODULE(trilepton);

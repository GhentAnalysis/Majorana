#include "effFast.h"
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

//btagging
#include "SUSYAnalyzer/PatAnalyzer/interface/BTagCalibrationStandalone.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace tools;
using namespace math;
using namespace reco::tau;

effFast::effFast(const edm::ParameterSet & iConfig) :
//eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
//eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
//eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
//eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
genparticleToken_(consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genPartsLabel"))),
mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
pdfvariablesToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("pdfvariablesLabel"))),
IT_beamspot(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotLabel"))),
PileUpToken_(consumes<vector< PileupSummaryInfo >>(iConfig.getParameter<edm::InputTag>("slimmedAddPileupInfoLabel"))),
goodOfflinePrimaryVerticesToken_(consumes<std::vector<Vertex>>(iConfig.getParameter<edm::InputTag>("goodOfflinePrimaryVerticesLabel"))),
packedPFCandidatesToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidatesLabel"))),
IT_muon(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonLabel"))),
IT_electron(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("ElectronLabel"))),
IT_jet(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("JetLabel"))),
IT_pfmet(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("METLabel"))),
reducedEgammaToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("reducedEgammaLabel"))),
IT_tau(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("TauLabel"))),
fixedGridRhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetCentralNeutralLabel"))),
fixedGridRhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetAllLabel"))),
IT_hltresults(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTResultsLabel"))),
triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
IT_externalLHEProducer(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("exernalLHEPLabel"))),
    
//eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
//eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
//eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
_relIsoCutE(999.),//0.15
_relIsoCutMu(999.),//0.15
_relIsoCutEloose(999.), //0.5
_relIsoCutMuloose(999.), //0.5
_chargeConsistency(true),
_minPt0(5.),
_minPt1(7.),
_tightD0Mu(0.05),
_tightD0E(0.05),
_looseD0Mu(0.05), // 0.05 for sync
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
    Sample              = iConfig.getUntrackedParameter<std::string>("SampleLabel") ;
    SampleName          = iConfig.getUntrackedParameter<std::string>("SampleName") ;
    //IT_muon             = iConfig.getParameter<edm::InputTag>("MuonLabel") ;
    //IT_electron         = iConfig.getParameter<edm::InputTag>("ElectronLabel") ;
    //IT_tau              = iConfig.getParameter<edm::InputTag>("TauLabel") ;
    //IT_tauDiscriminator = iConfig.getParameter<edm::InputTag>("TauDiscriminatorLabel") ;
    //IT_jet              = iConfig.getParameter<edm::InputTag>("JetLabel");
    //IT_pfmet            = iConfig.getParameter<edm::InputTag>("METLabel")  ;
    //IT_beamspot         = iConfig.getParameter<edm::InputTag>("BeamSpotLabel");
    //IT_METFilters       = iConfig.getParameter<edm::InputTag>("METFilter");
    
    //outfile = fopen("FakeSync.txt", "w");
}


void effFast::beginJob()
{
    Nvtx           = fs->make<TH1F>("N_{vtx}"        , "Number of vertices;N_{vtx};events / 1"  ,    40, 0., 40.);
    
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);

    _hCounterSUSY = fs->make<TH2D>("hCounterSUSY", "Events counter", 80, 0., 2000., 60, 0., 1500.);

    _hCounterSUSY_WW = fs->make<TH2D>("hCounterSUSY_WW", "Events counter", 80, 0., 2000., 60, 0., 1500.);
    _hCounterSUSY_WZ = fs->make<TH2D>("hCounterSUSY_WZ", "Events counter", 80, 0., 2000., 60, 0., 1500.);
    _hCounterSUSY_ZZ = fs->make<TH2D>("hCounterSUSY_ZZ", "Events counter", 80, 0., 2000., 60, 0., 1500.);
    
    outputTree = new TTree("fakeTree","fakeTree");
    
    _leptonP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
        new ( (*_leptonP4)[i] ) TLorentzVector();
    }
    //outputTree->Branch("_leptonP4", "TClonesArray", &_leptonP4, 32000, 0);
    
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
    outputTree->Branch("_activity", &_activity, "_activity[_nLeptons]/D");
    outputTree->Branch("_miniisolation", &_miniisolation, "_miniisolation[_nLeptons][2]/D");
    outputTree->Branch("_miniisolationCharged", &_miniisolationCharged, "_miniisolationCharged[_nLeptons][2]/D");
    outputTree->Branch("_multiisolation", &_multiisolation, "_multiisolation[_nLeptons][5]/O");
    outputTree->Branch("_ptrel", &_ptrel, "_ptrel[_nLeptons]/D");
    outputTree->Branch("_ptratio", &_ptratio, "_ptratio[_nLeptons]/D");
    outputTree->Branch("_muonSegmentComp", &_muonSegmentComp, "_muonSegmentComp[_nLeptons]/D");
    outputTree->Branch("_muonDpt", &_muonDpt, "_muonDpt[_nLeptons]/D");

    outputTree->Branch("_ICHEPsoftMuonId", &_ICHEPsoftMuonId, "_ICHEPsoftMuonId[_nLeptons]/O");

    //outputTree->Branch("_mll", &_mll, "_mll[3]/D");
    //outputTree->Branch("_ossf", &_ossf, "_ossf[3]/O");


    //Triggers
    
    /*
    outputTree->Branch("_trigMu8", &_trigMu8, "_trigMu8/O");
    outputTree->Branch("_trigMu17", &_trigMu17, "_trigMu17/O");
    outputTree->Branch("_trigMu8iso", &_trigMu8iso, "_trigMu8iso/O");
    outputTree->Branch("_trigMu17iso", &_trigMu17iso, "_trigMu17iso/O");
    outputTree->Branch("_trigEle12", &_trigEle12, "_trigEle12/O");
    outputTree->Branch("_trigEle12iso", &_trigEle12iso, "_trigEle12iso/O");

    outputTree->Branch("_trigDiMuIso", &_trigDiMuIso, "_trigDiMuIso/O");
    outputTree->Branch("_trigDiMuTkIso", &_trigDiMuTkIso, "_trigDiMuTkIso/O");
    outputTree->Branch("_trigMu8Ele23Iso", &_trigMu8Ele23Iso, "_trigMu8Ele23Iso/O");
    outputTree->Branch("_trigMu23Ele12Iso", &_trigMu23Ele12Iso , "_trigMu23Ele12Iso/O");
    outputTree->Branch("_trigEle23Ele12Iso", &_trigEle23Ele12Iso , "_trigEle23Ele12Iso/O");
    
    outputTree->Branch("_trigDoubleMu8", &_trigDoubleMu8 , "_trigDoubleMu8/O");
    outputTree->Branch("_trigMu8Ele8", &_trigMu8Ele8 , "_trigMu8Ele8/O");
    outputTree->Branch("_trigDoubleEle8", &_trigDoubleEle8 , "_trigDoubleEle8/O");
    outputTree->Branch("_trigTripleMu", &_trigTripleMu , "_trigTripleMu/O");
    outputTree->Branch("_trigTripleDiMu9Ele9", &_trigTripleDiMu9Ele9 , "_trigTripleDiMu9Ele9/O");
    outputTree->Branch("_trigTripleMu8DiEle12", &_trigTripleMu8DiEle12 , "_trigTripleMu8DiEle12/O");
    outputTree->Branch("_trigTripleEle16Ele12Ele8", &_trigTripleEle16Ele12Ele8 , "_trigTripleEle16Ele12Ele8/O");
    */
    
    outputTree->Branch("_triggers2l", &_triggers2l, "_triggers2l[2][6]/O");
    outputTree->Branch("_triggers2lbkp", &_triggers2lbkp, "_triggers2lbkp[2][6]/O");
    outputTree->Branch("_triggersCS", &_triggersCS, "_triggersCS[5][5]/O");
    outputTree->Branch("_triggersCSb", &_triggersCSb, "_triggersCSb[2]/O");
    outputTree->Branch("_triggers1l", &_triggers1l, "_triggers1l[5]/O");

    outputTree->Branch("_triggers2lpresc", &_triggers2lpresc, "_triggers2lpresc[2][6]/D");
    outputTree->Branch("_triggers2lbkppresc", &_triggers2lbkppresc, "_triggers2lbkppresc[2][6]/D");
    outputTree->Branch("_triggersCSpresc", &_triggersCSpresc, "_triggersCSpresc[5][5]/D");
    outputTree->Branch("_triggersCSbpresc", &_triggersCSbpresc, "_triggersCSbpresc[2]/D");
    outputTree->Branch("_triggers1lpresc", &_triggers1lpresc, "_triggers1lpresc[5]/D");
    
    //outputTree->Branch("_index1", &_index1, "_index1/I");
    //outputTree->Branch("_index2", &_index2, "_index2/I");

    //outputTree->Branch("_sb", &_sb, "_sb/O");
    //outputTree->Branch("_doubleF", &_doubleF, "_doubleF/O");
    
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
    outputTree->Branch("_isvetoIDCutBased", &_isvetoIDCutBased, "_isvetoIDCutBased[_nLeptons]/O");
    outputTree->Branch("_islooseIDCutBased", &_islooseIDCutBased, "_islooseIDCutBased[_nLeptons]/O");
    outputTree->Branch("_ismediumIDCutBased", &_ismediumIDCutBased, "_ismediumIDCutBased[_nLeptons]/O");
    outputTree->Branch("_istightIDCutBased", &_istightIDCutBased, "_istightIDCutBased[_nLeptons]/O");


    outputTree->Branch("_trigEmulator", &_trigEmulator, "_trigEmulator[_nLeptons]/O");
    outputTree->Branch("_isotrigEmulator", &_isotrigEmulator, "_isotrigEmulator[_nLeptons]/O");
    
    outputTree->Branch("_closeJetPtAll", &_closeJetPtAll, "_closeJetPtAll[_nLeptons]/D");
    outputTree->Branch("_closeJetAngAll", &_closeJetAngAll, "_closeJetAngAll[_nLeptons]/D");
    //outputTree->Branch("_ptRel", &_ptRel, "_ptRel[6]/D");
    //outputTree->Branch("_ptRelAll", &_ptRelAll, "_ptRelAll[6]/D");
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

    /*
    outputTree->Branch("_closeJetPtAllMC", &_closeJetPtAllMC, "_closeJetPtAllMC[6]/D");
    outputTree->Branch("_closeJetPtAllstatus", &_closeJetPtAllstatus, "_closeJetPtAllstatus[6]/D");
    outputTree->Branch("_partonIdMatched", &_partonIdMatched, "_partonIdMatched[6]/I");
    outputTree->Branch("_sameParton", &_sameParton, "_sameParton[6]/O");
    
    if (_regression) {
        outputTree->Branch("_regVars", &_regVars, "_regVars[15]/D");
        
        outputTree->Branch("hJet_ptRaw", &hJet_ptRaw, "hJet_ptRaw/D");
        outputTree->Branch("hJet_genPt", &hJet_genPt, "hJet_genPt/D");
        outputTree->Branch("hJet_pt", &hJet_pt, "hJet_pt/D");
        outputTree->Branch("hJet_phi", &hJet_phi, "hJet_phi/D");
        outputTree->Branch("hJet_eta", &hJet_eta, "hJet_eta/D");
        outputTree->Branch("hJet_e", &hJet_e, "hJet_e/D");
        
        outputTree->Branch("hJet_ptLeadTrack", &hJet_ptLeadTrack, "hJet_ptLeadTrack/D");
        
        outputTree->Branch("hJet_vtx3dL", &hJet_vtx3dL, "hJet_vtx3dL/D");
        outputTree->Branch("hJet_vtx3deL", &hJet_vtx3deL, "hJet_vtx3deL/D");
        outputTree->Branch("hJet_vtxMass", &hJet_vtxMass, "hJet_vtxMass/D");
        outputTree->Branch("hJet_vtxPt", &hJet_vtxPt, "hJet_vtxPt/D");
        
        outputTree->Branch("hJet_cef", &hJet_cef, "hJet_cef/D");
        
        outputTree->Branch("hJet_nconstituents", &hJet_nconstituents, "hJet_nconstituents/D");
        outputTree->Branch("hJet_JECUnc", &hJet_JECUnc, "hJet_JECUnc/D");
        
        outputTree->Branch("hJet_SoftLeptptRel", &hJet_SoftLeptptRel, "hJet_SoftLeptptRel/D");
        outputTree->Branch("hJet_SoftLeptPt", &hJet_SoftLeptPt, "hJet_SoftLeptPt/D");
        outputTree->Branch("hJet_SoftLeptdR", &hJet_SoftLeptdR, "hJet_SoftLeptdR/D");
        
        outputTree->Branch("hJet_SoftLeptIdlooseMu", &hJet_SoftLeptIdlooseMu, "hJet_SoftLeptIdlooseMu/D");
        outputTree->Branch("hJet_SoftLeptId95", &hJet_SoftLeptId95, "hJet_SoftLeptId95/D");
    }
    */
    
    /*****    outputTree->Branch("_closeJetPt", &_closeJetPt, "_closeJetPt[3]/D");
     outputTree->Branch("_closeJetAng", &_closeJetAng, "_closeJetAng[3]/D");
     outputTree->Branch("_ptRel", &_ptRel, "_ptRel[3]/D");
     
     outputTree->Branch("_ptRelAllNotItself", &_ptRelAllNotItself, "_ptRelAllNotItself[3]/D");
     outputTree->Branch("_closeJetPtAllNotItself", &_closeJetPtAllNotItself, "_closeJetPtAllNotItself[3]/D");
     outputTree->Branch("_closeJetAngAllNotItself", &_closeJetAngAllNotItself, "_closeJetAngAllNotItself[3]/D");
     */
    /******
     outputTree->Branch("_n_bJetsAll", &_n_bJetsAll, "_n_bJetsAll/I");
     outputTree->Branch("_n_JetsAll", &_n_JetsAll, "_n_JetsAll/I");
     outputTree->Branch("_bTaggedAll", &_bTaggedAll, "_bTaggedAll[100]/O");
     outputTree->Branch("_jetEtaAll", &_jetEtaAll, "_jetEtaAll[100]/D");
     outputTree->Branch("_jetPhiAll", &_jetPhiAll, "_jetPhiAll[100]/D");
     outputTree->Branch("_jetPtAll", &_jetPtAll, "_jetPtAll[100]/D");
     outputTree->Branch("_csvAll", &_csvAll, "_csvAll[100]/D");
     */
    outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
    outputTree->Branch("_n_MCTruth_PV", &_n_MCTruth_PV, "_n_MCTruth_PV/D");
    
    outputTree->Branch("_met", &_met, "_met/D");
    outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
    outputTree->Branch("HT", &HT, "HT/D");
    
    outputTree->Branch("_genmet", &_genmet, "_genmet/D");
    outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/D");
    
    outputTree->Branch("_genqpt", &_genqpt, "_genqpt/D");

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

    
    GPM = GenParticleManager();
    
    bool isData = !(Sample=="ElectronsMC");
    if (isData)
        fMetCorrector = new OnTheFlyCorrections("Spring16_25nsV6_DATA", isData); //isData = true
    else
        fMetCorrector = new OnTheFlyCorrections("Spring16_FastSimV1_MC", isData); //isData = true
    _corrLevel = "L3Absolute";
    if (isData) _corrLevel = "L2L3Residual";

    /*
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iEventSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);

    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];

    jecUnc = new JetCorrectionUncertainty(JetCorPar);
    */
    jecUnc = new JetCorrectionUncertainty("/user/ikhvastu/CMSSW_8_0_5/src/SUSYAnalyzer/PatAnalyzer/jetfiles/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt");
    
    /*
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml");
    
    //vector<string> myManualCatWeigthsTrig;
    string the_path;
    for (unsigned i  = 0 ; i < myManualCatWeigths.size() ; i++){
        the_path = edm::FileInPath ( myManualCatWeigths[i] ).fullPath();
        myManualCatWeigthsTrig.push_back(the_path);
    }
    
    myMVATrig = new EGammaMvaEleEstimatorCSA14();
    myMVATrig->initialize("BDT",
                          EGammaMvaEleEstimatorCSA14::kNonTrigPhys14,
                          true,
                          myManualCatWeigthsTrig);
                          */
    
    looseMVA[0][0] = -0.70;
    looseMVA[1][0] = -0.83;
    looseMVA[2][0] = -0.92;
    looseMVA[0][1] = 0.87;
    looseMVA[1][1] = 0.60;
    looseMVA[2][1] = 0.17;

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
    
    /*
    calib_csvv2 = new BTagCalibration("csvv2", "/localgrid/ikhvastu/CMSSW_7_4_14/src/SUSYAnalyzer/PatAnalyzer/test/ttH_BTV_CSVv2_13TeV_2015D_20151120.csv"); // 25s version of SFs
    reader = new BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING, "iterativefit", "central");
    // JESUp
    reader_JESUp = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit","up_jes");
    // JESDown
    reader_JESDown = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_jes");
    // LFUp
    reader_LFUp = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_lf");
    // LFDown
    reader_LFDown = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_lf");
    // HFUp
    reader_HFUp = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_hf");
    // HFDown
    reader_HFDown = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_hf");
    // HFStats1Up
    reader_HFStats1Up = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_hfstats1");
    // HFStats1Down
    reader_HFStats1Down = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_hfstats1");
    // HFStats2Up
    reader_HFStats2Up = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_hfstats2");
    // HFStats2Down
    reader_HFStats2Down = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_hfstats2");
    // LFStats1Up
    reader_LFStats1Up = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_lfstats1");
    // LFStats1Down
    reader_LFStats1Down = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_lfstats1");
    // LFStats2Up
    reader_LFStats2Up = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_lfstats2");
    // LFStats2Down
    reader_LFStats2Down = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_lfstats2");
    // CFErr1Up
    reader_CFErr1Up = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_cferr1");
    // CFErr1Down
    reader_CFErr1Down = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_cferr1");
    // CFErr2Up
    reader_CFErr2Up = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "up_cferr2");
    // CFErr2Down
    reader_CFErr2Down = new  BTagCalibrationReader(calib_csvv2, BTagEntry::OP_RESHAPING,"iterativefit", "down_cferr2");
    */
}

void effFast::endJob() {
    //outputTree -> Write();
    // store nEventsTotal and nEventsFiltered in preferred way
    std::cout<<_nEventsTotal<<std::endl;
    std::cout<<_nEventsFiltered<<std::endl;
    
    delete fMetCorrector;
    
    
}

void effFast::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
    using namespace edm;

    /*
    bool changed(false);
    if(!hltConfig_.init(iRun,iSetup,"HLT",changed) ){
        edm::LogError( "CandidateTriggerObjectProducer" ) <<
        "Error! Can't initialize HLTConfigProvider";
        throw cms::Exception("HLTConfigProvider::init() returned non 0");
    }
    */

    /*
    edm::Handle<LHERunInfoProduct> run; 
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
     
    iRun.getByLabel( "externalLHEProducer", run );
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
     
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
       std::cout << iter->tag() << std::endl;
       std::vector<std::string> lines = iter->lines();
       for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
          std::cout << lines.at(iLine);
       }
    }
    */
    
    return;
}

void effFast::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{
    //bool islepton;

    _nZboson = 0;

    if (Sample=="ElectronsMC") {
        //******************************************************************************************************************
        // Gen level particles                  ****************************************************************************
        //******************************************************************************************************************
        //iEvent.getByLabel("packedGenParticles", TheGenParticles);
        iEvent.getByToken(genparticleToken_, TheGenParticles);
        std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW;
        if( TheGenParticles.isValid() )
        {
            GPM.SetCollection(TheGenParticles);
            GPM.Classify();
            vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),1);
            vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),1);
            vGenNPMuons = GPM.filterByStatus(GPM.getNonPromptMuons(),1);
            vGenNPElectrons = GPM.filterByStatus(GPM.getNonPromptElectrons(),1);
            //std::cout<<"*************"<<std::endl;
            
            TLorentzVector Gen0;
            Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
            _genqpt = 0;
            cout << "Particle ids in the event" << endl;
            for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ )
            {
                int id = TMath::Abs(p->pdgId());
                cout << id << " ";
                if(id == 23)
                    _nZboson++;
                if ( (id == 12 || id == 14 || id == 16 ) && (p->status() == 1) ) {
                    TLorentzVector Gen;
                    Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );
                    Gen0 += Gen;
                }
                if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 21 || id == 22 ) && (p->status() == 23)){
                    _genqpt += p->pt();
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
        }
    }

    if (Sample=="ElectronsMC") {

        edm::Handle<LHEEventProduct> EvtHandle ;
        iEvent.getByToken(IT_externalLHEProducer , EvtHandle ) ;

        for(int i = 0; i < 111; i++)
            _LHEweight[i] = EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
    
        for(int i = 0; i < 111; i++)
            _LHEweightID[i] = stod(EvtHandle->weights()[i].id);
    }
    
    /*
    if ( firstEvent_ ) {
     edm::Handle<LHERunInfoProduct> run; 
     typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
     
     iEvent.getByLabel( "externalLHEProducer", run );
     LHERunInfoProduct myLHERunInfoProduct = *(run.product());
     
     for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
        std::cout << iter->tag() << std::endl;
        std::vector<std::string> lines = iter->lines();
        for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
            std::cout << lines.at(iLine);
        }
     }
    }
    */

        
     //edm::Handle<LHEEventProduct> pdfvariables;
     //iEvent.getByLabel("generator", pdfvariables);
     //_weight=pdfvariables.weight();
     //


     //============= trigger ==============
    /*
     _trigMu8 = false;
     _trigMu17 = false;
     _trigMu8iso = false;
     _trigMu17iso = false;
     _trigEle12 = false;
     _trigEle12iso = false;

     _trigDiMuIso = false;
     _trigDiMuTkIso = false;
     _trigMu8Ele23Iso = false;
     _trigMu23Ele12Iso = false;
     _trigEle23Ele12Iso = false;
     _trigDoubleMu8 = false;
     _trigMu8Ele8 = false;
     _trigDoubleEle8 = false;

     _trigTripleMu = false;
     _trigTripleDiMu9Ele9 = false;
     _trigTripleMu8DiEle12 = false;
     _trigTripleEle16Ele12Ele8 = false;
     */

    /*
     edm::Handle<edm::TriggerResults> trigResults;
     iEvent.getByToken(IT_hltresults, trigResults);

     edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
     iEvent.getByToken(triggerPrescales_, triggerPrescales);

     edm::Handle<edm::TriggerResults> triggerBits;
     iEvent.getByToken(triggerBits_, triggerBits);

     for (int i=0; i!=5; ++i) {
        for (int j=0; j!=5; ++j) { _triggersCS[j][i] = 0;}
        _triggers1l[i] = 0;
     }
     for (int i=0; i!=6; ++i) {
        for (int j=0; j!=2; ++j) _triggers2l[j][i] = 0;
        for (int j=0; j!=2; ++j) _triggers2lbkp[j][i] = 0;
     }

     for (int i=0; i!=2; ++i) 
         _triggersCSb[i] = 0;
     Flag_eeBadScFilter = 0;

     if( trigResults.failedToGet() ) cout << "--- NO TRIGGER RESULTS !! ---" << endl;
     
     //==================================
     
     if( !trigResults.failedToGet() ) {
        unsigned int n_Triggers = trigResults->size();
        const edm::TriggerNames & triggerNames = iEvent.triggerNames(*trigResults);

        if ( firstEvent_ ) {
            edm::TriggerNames::Strings allTriggers( triggerNames.triggerNames() );
            std::cout << "--- Trigger Menu --- " << std::endl;
            for ( unsigned int i_Name = 0; i_Name < n_Triggers; ++i_Name ) {
                std::cout << allTriggers.at( i_Name ) << " " << triggerPrescales->getPrescaleForIndex( i_Name) << std::endl;
            }
        std::cout << "-------------------- " << std::endl;
        firstEvent_ = false;
        }
        if ( firstEvent_ ) {
            const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
            for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
                std::cout << "Trigger " << names.triggerName(i) << 
                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
                ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
                << std::endl;
            }
            firstEvent_ = false;
        }
        for(unsigned int i_Trig = 0; i_Trig < n_Triggers; ++i_Trig ) {
            if (trigResults.product()->accept(i_Trig)) {
                TString TrigPath = triggerNames.triggerName(i_Trig); 
                 double prescaleHLTL1 = 1.;
                 if (Sample!="ElectronsMC") 
                    prescaleHLTL1 = triggerPrescales->getPrescaleForIndex(i_Trig);

                for (int i=0; i!=6; ++i) {
                    for (int j=0; j!=2; ++j) {
                        if (TrigPath.Contains(_triggers2lNamesTTZ[j][i])){
                            _triggers2l[j][i] = 1;
                            _triggers2lpresc[j][i] = prescaleHLTL1;
                            //std::cout<< "Second info about prescale (trigger 2l): " << j << " " << i << " " << TrigPath<<" "<<"" <<" "<<prescaleHLTL1<<std::endl;
                        }
                    }
                    for (int j=0; j!=2; ++j) {
                        if (TrigPath.Contains(_triggers2lNamesTTZBkp[j][i])){
                            _triggers2lbkp[j][i] = 1;
                            _triggers2lbkppresc[j][i] = prescaleHLTL1;
                            //std::cout<< "Second info about prescale (trigger 2l bkp): " << j << " " << i << " " << TrigPath<<" "<<"" <<" "<<prescaleHLTL1<<std::endl;
                        }
                    }
                }
                
                for (int i=0; i!=5; ++i) {
                    for (int j=0; j!=5; ++j) {
                        if (TrigPath.Contains(_triggersCSNamesTTZ[j][i])){
                            _triggersCS[j][i] = 1;
                            _triggersCSpresc[j][i] = prescaleHLTL1;
                            //std::cout<< "Second info about prescale (trigger CS): " << j << " " << i << " " << TrigPath<<" "<<"" <<" "<<prescaleHLTL1<<std::endl;
                        }
                    }
                    if (TrigPath.Contains(_triggers1lNamesTTZ[i])){
                        _triggers1l[i] = 1;
                        _triggers1lpresc[i] = prescaleHLTL1;
                        //std::cout<< "Second info about prescale (trigger 1l): " << i << " " << TrigPath<<" "<<"" <<" "<<prescaleHLTL1<<std::endl;
                    }
                }
                for (int i=0; i!=2; ++i) {
                    if (TrigPath.Contains(_triggersCSbNamesTTZ[i])){
                        _triggersCSb[i] = 1;
                        _triggersCSbpresc[i] = prescaleHLTL1;
                        //std::cout<< "Second info about prescale (trigger csb): " << i << " " << TrigPath<<" "<<"" <<" "<<prescaleHLTL1<<std::endl;
                    }
                }
                
                if(TrigPath.Contains("Flag_eeBadScFilter"))
                    Flag_eeBadScFilter=true;

            }
        }
     }
     */
     

    realdata_ = iEvent.isRealData();

    Double_t weight = 1.;
    if (Sample=="ElectronsMC") {
        edm::Handle<GenEventInfoProduct> pdfvariables;
        iEvent.getByToken(pdfvariablesToken_, pdfvariables);
        weight = pdfvariables->weight();
        _weight=weight;
    }

    _runNb = iEvent.id().run();
    _eventNb = iEvent.id().event();
    _lumiBlock = iEvent.luminosityBlock();
   
   
    //if (_eventNb != 44802315) return;
    //if ((_eventNb != 17707) || (_runNb != 1) || (_lumiBlock != 178)) return;
    
    std::cout<<"EVENT "<<_runNb << " " << _lumiBlock << " " << _eventNb<<std::endl;
    /*
    std::cout << "Triggers: " << std::endl;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 6; j++)
            std::cout << _triggers2l[i][j] << " ";
        std::cout << std::endl;
    }
    */
    
    //============ Total number of events is the sum of the events ============
    //============ in each of these luminosity blocks ============
    _nEventsTotalCounted++;
    _hCounter->Fill(0., weight);
    
     // To run on SUSY samples
    // read T1tttt model
    /*
     float Mass1, Mass3, Mass2, Mass4, xsection;
     float xparam;
     edm::Handle< LHEEventProduct > product;
     iEvent.getByLabel("source", product);
     LHEEventProduct::comments_const_iterator c_begin = product->comments_begin();
     LHEEventProduct::comments_const_iterator c_end = product->comments_end();
     
     for( LHEEventProduct::comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
     std::cout << (*cit);
     size_t found = (*cit).find("model");
     
     if( found != std::string::npos) {
        size_t foundLength = (*cit).size();
        found = (*cit).find("T1tttt");
        //found = (*cit).find("T6ttWW");
        std::string smaller = (*cit).substr(found+6,foundLength);
        found = smaller.find("_");
        smaller = smaller.substr(found+1,smaller.size());
     
        std::istringstream iss(smaller);
        iss >> Mass1;
        iss.clear();
     
        found = smaller.find("_");
        smaller = smaller.substr(found+1,smaller.size());
        iss.str(smaller);
        iss >> Mass2;
        iss.clear();
     }
     }
     std::cout << std::endl;
     //msbottom = Mass1;
     //mchargino = Mass2;
     //mLsp = Mass3; //=const=50
     _mgluino = Mass1;
     _mchi0 = Mass2;
     std::cout<<"Masses in sample: " << Mass1<<" "<<Mass2<<" "<<Mass3<<" "<<Mass4<<" "<<xparam<<" "<<xsection<<std::endl;
     */

    // To run on SUSY samples
    // read T6ttWW model
    /*
     float Mass1, Mass3, Mass2, Mass4, xsection;
     float xparam;
     edm::Handle< LHEEventProduct > product;
     iEvent.getByLabel("source", product);
     LHEEventProduct::comments_const_iterator c_begin = product->comments_begin();
     LHEEventProduct::comments_const_iterator c_end = product->comments_end();
     
     for( LHEEventProduct::comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
     std::cout << (*cit);
     size_t found = (*cit).find("model");
     
     if( found != std::string::npos) {
        size_t foundLength = (*cit).size();
        //found = (*cit).find("T1tttt");
        found = (*cit).find("T6ttWW");
        std::string smaller = (*cit).substr(found+6,foundLength);
        found = smaller.find("_mLSP");
        smaller = smaller.substr(found+5,smaller.size());
     
        std::istringstream iss(smaller);
        iss >> Mass1;
        iss.clear();
     
        found = smaller.find("_");
        smaller = smaller.substr(found+1,smaller.size());
        iss.str(smaller);
        iss >> Mass2;

        found = smaller.find("_");
        smaller = smaller.substr(found+1,smaller.size());
        iss.str(smaller);
        iss >> Mass3;
        iss.clear();
     }
     }
     std::cout << std::endl;
     //msbottom = Mass1;
     //mchargino = Mass2;
     //mLsp = Mass3; //=const=50
     _mgluino = Mass2;
     _mchi0 = Mass3;
     std::cout<<"Masses in sample: " << Mass1<<" "<<Mass2<<" "<<Mass3<<" "<<Mass4<<" "<<xparam<<" "<<xsection<<std::endl;

    _hCounterSUSY->Fill(Mass2, Mass3, weight);
    */

    /*
     float Mass1, Mass3, Mass2, Mass4, xsection;
     float xparam;
     edm::Handle< LHEEventProduct > product;
     iEvent.getByLabel("source", product);
     LHEEventProduct::comments_const_iterator c_begin = product->comments_begin();
     LHEEventProduct::comments_const_iterator c_end = product->comments_end();
     
     for( LHEEventProduct::comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
     std::cout << (*cit);
     size_t found = (*cit).find("model");
     
     if( found != std::string::npos) {
        size_t foundLength = (*cit).size();
        found = (*cit).find("T5qqqqVV");
        //found = (*cit).find("T6ttWW");
        std::string smaller = (*cit).substr(found+8,foundLength);
        found = smaller.find("_");
        smaller = smaller.substr(found+1,smaller.size());
     
        std::istringstream iss(smaller);
        iss >> Mass1;
        iss.clear();
     
        found = smaller.find("_");
        smaller = smaller.substr(found+1,smaller.size());
        iss.str(smaller);
        iss >> Mass2;
        iss.clear();
     }
     }
     std::cout << std::endl;
     //msbottom = Mass1;
     //mchargino = Mass2;
     //mLsp = Mass3; //=const=50
     _mgluino = Mass1;
     _mchi0 = Mass2;
     std::cout<<"Masses in sample: " << Mass1<<" "<<Mass2<<" "<<Mass3<<" "<<Mass4<<" "<<xparam<<" "<<xsection<<std::endl;

    if(_nZboson == 0)
        _hCounterSUSY_WW->Fill(Mass1, Mass2, weight);
    if(_nZboson == 1)
        _hCounterSUSY_WZ->Fill(Mass1, Mass2, weight);
    if(_nZboson == 2)
        _hCounterSUSY_ZZ->Fill(Mass1, Mass2, weight);

    _hCounterSUSY->Fill(Mass1, Mass2, weight);
    */

     /*edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
     const edm::LuminosityBlock & lumi = iEvent.getLuminosityBlock();
     
     lumi.getByLabel("nEventsTotal", nEventsTotalCounter);
     _nEventsTotal += nEventsTotalCounter->value;
     _hCounter->Fill(1, nEventsTotalCounter->value);
     
     edm::Handle<edm::MergeableCounter> nEventsFilteredCounter;
     lumi.getByLabel("nEventsFiltered", nEventsFilteredCounter);
     _nEventsFiltered += nEventsFilteredCounter->value;
     _hCounter->Fill(2, nEventsFilteredCounter->value);
     _hCounter->Fill(3, double(nEventsTotalCounter->value)/double(nEventsFilteredCounter->value));
     */
    //============ Counter done ============
    
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
    if (Sample=="ElectronsMC") {
        edm::Handle<vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByToken(PileUpToken_, PupInfo);
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
    iEvent.getByToken( goodOfflinePrimaryVerticesToken_, theVertices) ;
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
    iEvent.getByToken(packedPFCandidatesToken_, pfcands);
    //==================================
    
    //============ Pat Muons ============
    edm::Handle< std::vector<pat::Muon> > thePatMuons;
    iEvent.getByToken( IT_muon, thePatMuons );
    //if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    //==================================
    
    /*
    std::cout << "Muons:" << std::endl;
    for( std::vector<pat::Muon>::const_iterator mu = (*thePatMuons).begin() ; mu != (*thePatMuons).end() ; mu++ ){
        std::cout << mu->pt() << " " << mu->eta() << " " << mu->phi() << std::endl;
    }
    */
    
    
    //============ Pat Electrons ============
    edm::Handle< edm::View<pat::Electron> > thePatElectrons;
    iEvent.getByToken( IT_electron, thePatElectrons );
    //if( ! thePatElectrons.isValid() ) ERR( IT_electron );
    //==================================
    
    /*
    std::cout << "Electrons:" << std::endl;
    for( std::vector<pat::Electron>::const_iterator el = (*thePatElectrons).begin() ; el != (*thePatElectrons).end() ; el++ ){
        std::cout << el->pt() << " " << el->eta() << " " << el->phi() << std::endl;
    }
    */

    // Get MVA values and categories (optional)
    edm::Handle<edm::ValueMap<float> > mvaValues;
    iEvent.getByToken(mvaValuesMapToken_,mvaValues);


    // Get the electron ID data from the event stream.
    // Note: this implies that the VID ID modules have been run upstream.
    // If you need more info, check with the EGM group.
    //edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
    //edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    //edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
    //iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
    //iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
    //iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);

    /*
    edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
    edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
    iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
    iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
    iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
    iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
    */
 
    //============ Conversions ============
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    //iEvent.getByLabel("reducedEgamma","reducedConversions", theConversions);
    iEvent.getByToken(reducedEgammaToken_, theConversions);
    //==================================
    
    //============ Pat Jets ============
    edm::Handle< std::vector< pat::Jet> > thePatJets;
    iEvent.getByToken(IT_jet , thePatJets );
    //if( ! thePatJets.isValid() ) ERR(IT_jet);
    //==================================
    edm::Handle<double> rhoJets;
    iEvent.getByToken(fixedGridRhoFastjetCentralNeutralToken_ , rhoJets);//kt6PFJets
    myRhoJets = *rhoJets;
    //==================================
    edm::Handle<double> rhoJECJets;
    iEvent.getByToken(fixedGridRhoFastjetAllToken_ , rhoJECJets);//kt6PFJets
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

        // previous implementation
        /*
        for(int unsigned it=0; it!=(&*jet)->numberOfDaughters(); ++it) {
            const pat::PackedCandidate * icand = dynamic_cast<const pat::PackedCandidate *> ((&*jet)->daughter(it));
            if(  fabs(icand->pdgId()) == 13 ) {
                double mupt = icand->pt();
                double muphi = icand->phi();
                for( std::vector<pat::Muon>::const_iterator mu = (*thePatMuons).begin() ; mu != (*thePatMuons).end() ; mu++ ) {
                    if ( (fabs(mu->pt() - mupt)/mupt < 0.1) && (fabs(mu->phi() - muphi) < 0.01)  ) {
                        if (mu->isGlobalMuon() || mu->isStandAloneMuon()) {
                            TLorentzVector pMu;
                            pMu.SetPtEtaPhiE(mu->pt(),mu->eta(),mu->phi(),mu->energy());
                            //pMu.SetPtEtaPhiE(icand->pt(),icand->eta(),icand->phi(),icand->energy());
                            pJet-=pMu;
                        }
                    }
                }
            }
        }
        */
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
    //std::cout<< "met: " <<_met<<"; "<<_met_phi<<std::endl;


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
        N_U_L_L // 19
    };
    
    
    //std::vector<const pat::Muon* > sMu = ssbMediumMuonSelector( *thePatMuons, _minPt0, PV, _looseD0Mu);
    std::vector<const pat::Muon* > sMu = effMuonSelector( *thePatMuons, _minPt0, PV, _looseD0Mu);
    
    //std::vector<const pat::Electron* > sEl = ssbMVAElectronSelector( *thePatElectrons, _minPt1, PV, _looseD0E, _chargeConsistency, theConversions, BS, true);

    //Taus
    edm::Handle<pat::TauCollection> PatTaus;
    iEvent.getByToken( IT_tau, PatTaus );
    std::vector<const pat::Tau* > sTau;
    for (const pat::Tau &tau : *PatTaus) {
        
        if (tau.pt() < _tauPt) continue;
        if (fabs(tau.eta()) > _tauEta) continue;
        //printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
        //       tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
        
        // old definition 24 May 2016
        /*
        if (! (tau.tauID("againstElectronLooseMVA5")
               && tau.tauID("againstMuonLoose3")
               && tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")
               && tau.tauID("decayModeFindingNewDMs"))) continue;
               */

        /*
        if (! (tau.tauID("decayModeFinding")
            && tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"))) continue;
            */

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
        
        /*if ((_jetPtAll[i] > _jetPtCut) && (fabs(_jetEtaAll[i]) < _jetEtaCut)) {
            
            _jetEta[_n_Jets] = _jetEtaAll[i];
            _jetPhi[_n_Jets] = _jetPhiAll[i];
            _jetPt[_n_Jets] = _jetPtAll[i];
            
            ((TLorentzVector *)_jetP4->At(_n_Jets))->SetPtEtaPhiM( _jetPt[_n_Jets], _jetEta[_n_Jets], _jetPhi[_n_Jets], 0 );
            
            _csv[_n_Jets] = _csvAll[i];
            
            if(_csvAll[i] > 0.679) {
                Bjets.push_back( &*SelectedJetsAll[i] );
                _bTagged[_n_Jets] = true;
                _n_bJets++;
            } else _bTagged[_n_Jets] = false;
            
            HT+= _jetPt[_n_Jets];
            _n_Jets++;
        }*/
        
        
    }
    _n_JetsAll = SelectedJetsAll.size();
    
    int leptonCounter = 0;
    //std::cout << "Muon selection started " << std::endl;
    for(unsigned int i = 0 ; i < sMu.size() ;i++ ){
        
        
        const pat::Muon *iM = sMu[i];
        //if (iM->pt() < 10) continue;
        
        //_hitsNumber[i] = 0;

        if (leptonCounter == 10) continue;

        _flavors[leptonCounter] = 1;
        _charges[leptonCounter] = iM->charge();
        _isolation[leptonCounter] = pfRelIso(iM,myRhoJets);
        _isolationDB[leptonCounter] = pfRelIso(iM);
        
        _activity[leptonCounter] = getActivityAroundLepton(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false, myRhoJets);

        _miniisolation[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false, myRhoJets);
        _miniisolation[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, false, myRhoJets);

        _miniisolationCharged[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, true, myRhoJets);
        _miniisolationCharged[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, true, myRhoJets);

        _ipPV[leptonCounter] = iM->innerTrack()->dxy(PV);
        _ipPVerr[leptonCounter] = iM->innerTrack()->dxyError();
        
        _ipZPV[leptonCounter] = iM->innerTrack()->dz(PV);
        _ipZPVerr[leptonCounter] = iM->innerTrack()->dzError();
        
        bool good_loose = iM->isPFMuon() && (iM->isGlobalMuon() || iM->isTrackerMuon());
        bool goodGlb = iM->isGlobalMuon() && iM->globalTrack()->normalizedChi2() < 3 && iM->combinedQuality().chi2LocalPosition < 12 && iM->combinedQuality().trkKink < 20;
        //bool good = iM->innerTrack()->validFraction() > 0.8 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
        bool good = iM->innerTrack()->validFraction() > 0.49 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451); // short term 

        //double _muonSegmentComp[leptonCounter] = iM->segmentCompatibility();
        _muonSegmentComp[leptonCounter] = iM->segmentCompatibility();

        _ICHEPsoftMuonId[leptonCounter] = iM->muonID("TMOneStationTight") && (iM->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 ) && (iM->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0) && (TMath::Abs(iM->innerTrack()->dxy(PV)) < 0.3) && (TMath::Abs(iM->innerTrack()->dz(PV)) < 20);

        _muonDpt[leptonCounter] = iM->innerTrack()->ptError()/iM->innerTrack()->pt();
        
        //_istightID[leptonCounter] = good;
        
        _3dIP[leptonCounter]    = iM->dB(pat::Muon::PV3D);
        _3dIPerr[leptonCounter] = iM->edB(pat::Muon::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);
        
        //std::cout<<"Muon info: Pt -> " << iM->pt() << "; Eta -> " << iM->eta()  << std::endl;
        //std::cout<<"Muon dxy -> " << _ipPV[leptonCounter] << "; dz -> " << _ipZPV[leptonCounter] << std::endl;
        //std::cout<<"Muon 3dIP -> " << _3dIPsig[leptonCounter] << std::endl;

        //_isloose[leptonCounter] = ( fabs(_ipPV[leptonCounter])<_looseD0Mu ) && (_isolation[leptonCounter]<_relIsoCutMuloose);
        //_isloose[leptonCounter] = ( fabs(_ipPV[leptonCounter])<_looseD0Mu );
        //if (_isloose[leptonCounter])
        //    _istightID[leptonCounter] = (good &&  _3dIPsig[leptonCounter]<4. );
        //else _istightID[leptonCounter] = false;

        //if(_miniisolation[leptonCounter][0] > 0.4) continue;
        //if(_isolationDB[leptonCounter] > 1.0 ) continue;
        _isloose[leptonCounter] = good_loose;
        _istight[leptonCounter] = good;


        //if(!_isloose[leptonCounter]) continue;
        //if (!_istight[leptonCounter]) continue;
        //if(_3dIPsig[leptonCounter])>4.
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iM->pt(), iM->eta(), iM->phi(), iM->energy());
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
        //std::cout<< "Muon loose selection passed -> " << _lPt[leptonCounter] << endl;

        fillCloseJetVars(leptonCounter, PV);
        
        if (Sample=="ElectronsMC") {
            //**************************************************************************************
            // MC
            //**************************************************************************************
            
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
            //fillIsoMCVars(leptonCounter);
            //matchCloseJet(leptonCounter);
            
            /*
            if (_regression) {
                std::vector <double> regVars = RegressionVars(SelectedJetsAll[_closeIndex[leptonCounter]], _closeJetPtAllMC[leptonCounter], iM);
                for (int k=0; k!=15; ++k) {
                    _regVars[k] = regVars[k];
                }
                _regVars[11] = fMetCorrector->getJECUncertainty(SelectedJetsAll[_closeIndex[leptonCounter]]->pt(),SelectedJetsAll[_closeIndex[leptonCounter]]->eta());
            }
            */
            
            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }
        
        /*
        if (_regression)
            fillRegVars(SelectedJetsAll.at(_closeIndex[leptonCounter]), _closeJetPtAllMC[leptonCounter], iM);
            */
        
        leptonCounter++;

        
    }
    //std::cout<<leptonCounter<<" "<<((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt() <<std::endl;
    _nMu = leptonCounter;
    
    //std::cout << "Electron selection started" << std::endl;
     for (size_t i = 0; i < thePatElectrons->size(); ++i){

        if (leptonCounter == 10) continue;

        const auto iE = thePatElectrons->ptrAt(i);
        //std::cout << "Pat electron: " << iE->pt() << " " << iE->superCluster()->eta() << " " << iE->phi() << " " << std::endl;
        if (!ssbMVAElectronSelectorPassed(iE, _minPt1, PV, _looseD0E, _chargeConsistency, theConversions, BS, false)) continue;
        //std::cout << "Passed electron selector" << std::endl;
        _mvaValue[leptonCounter] = (*mvaValues)[iE];

                
        _flavors[leptonCounter] = 0;
        _charges[leptonCounter] = iE->charge();
        _isolation[leptonCounter] = pfRelIso(&*iE, myRhoJets);
        _ipPV[leptonCounter] = iE->gsfTrack()->dxy(PV);
        _ipZPV[leptonCounter] = iE->gsfTrack()->dz(PV);
        
        _activity[leptonCounter] = getActivityAroundLepton(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, false, myRhoJets);
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


        double ooEmooP = fabs(1.0/iE->correctedEcalEnergy() - iE->eSuperClusterOverP()/iE->correctedEcalEnergy());
        double dEtaIn = iE->deltaEtaSuperClusterTrackAtVtx();
        double dPhiIn = iE->deltaPhiSuperClusterTrackAtVtx();
        double hOverE = iE->hcalOverEcal();
        double full5x5_sigmaIetaIeta = iE->full5x5_sigmaIetaIeta();
        
 
        if(fabs(iE->superCluster()->eta()) < 1.479){
            _isvetoIDCutBased[leptonCounter] = (ooEmooP < 0.207) && ( hOverE < 0.181) && (fabs(dEtaIn) < 0.0152) && (fabs(dPhiIn) < 0.216) && (full5x5_sigmaIetaIeta < 0.0114) && (fabs(_ipPV[leptonCounter]) < 0.0564) &&                    (fabs(_ipZPV[leptonCounter]) < 0.472) && (_hitsNumber[leptonCounter] <= 2) && (!_vtxFitConversion[leptonCounter]);
            _islooseIDCutBased[leptonCounter] = (ooEmooP < 0.102) && ( hOverE < 0.104) && (fabs(dEtaIn) < 0.0105) && (fabs(dPhiIn) < 0.115) && (full5x5_sigmaIetaIeta < 0.0103) && (fabs(_ipPV[leptonCounter]) < 0.0261) &&                   (fabs(_ipZPV[leptonCounter]) < 0.41) && (_hitsNumber[leptonCounter] <= 2) && (!_vtxFitConversion[leptonCounter]);
            _ismediumIDCutBased[leptonCounter] = (ooEmooP < 0.0174) && ( hOverE < 0.0876) && (fabs(dEtaIn) < 0.0103) && (fabs(dPhiIn) < 0.0336) && (full5x5_sigmaIetaIeta < 0.0101) && (fabs(_ipPV[leptonCounter]) < 0.0118) &&               (fabs(_ipZPV[leptonCounter]) < 0.373) && (_hitsNumber[leptonCounter] <= 2) && (!_vtxFitConversion[leptonCounter]);
            _istightIDCutBased[leptonCounter] = (ooEmooP < 0.012) && ( hOverE < 0.0597) && (fabs(dEtaIn) < 0.00926) && (fabs(dPhiIn) < 0.0336) && (full5x5_sigmaIetaIeta < 0.0101) && (fabs(_ipPV[leptonCounter]) < 0.0111) &&                (fabs(_ipZPV[leptonCounter]) < 0.0466) && (_hitsNumber[leptonCounter] <= 2) && (!_vtxFitConversion[leptonCounter]);
        }
        else{
            _isvetoIDCutBased[leptonCounter] = (ooEmooP < 0.174) && ( hOverE < 0.116) && (fabs(dEtaIn) < 0.0113) && (fabs(dPhiIn) < 0.237) && (full5x5_sigmaIetaIeta < 0.0352) && (fabs(_ipPV[leptonCounter]) < 0.222) &&                     (fabs(_ipZPV[leptonCounter]) < 0.921) && (_hitsNumber[leptonCounter] <= 3) && (!_vtxFitConversion[leptonCounter]);
            _islooseIDCutBased[leptonCounter] = (ooEmooP < 0.126) && ( hOverE < 0.0897) && (fabs(dEtaIn) < 0.00814) && (fabs(dPhiIn) < 0.182) && (full5x5_sigmaIetaIeta < 0.0301) && (fabs(_ipPV[leptonCounter]) < 0.118) &&                  (fabs(_ipZPV[leptonCounter]) < 0.822) && (_hitsNumber[leptonCounter] <= 1) && (!_vtxFitConversion[leptonCounter]);
            _ismediumIDCutBased[leptonCounter] = (ooEmooP < 0.0898) && ( hOverE < 0.0678) && (fabs(dEtaIn) < 0.00733) && (fabs(dPhiIn) < 0.114) && (full5x5_sigmaIetaIeta < 0.0283) && (fabs(_ipPV[leptonCounter]) < 0.0739) &&               (fabs(_ipZPV[leptonCounter]) < 0.602) && (_hitsNumber[leptonCounter] <= 1) && (!_vtxFitConversion[leptonCounter]);
            _istightIDCutBased[leptonCounter] = (ooEmooP < 0.00999) && ( hOverE < 0.0615) && (fabs(dEtaIn) < 0.00724) && (fabs(dPhiIn) < 0.0918) && (full5x5_sigmaIetaIeta < 0.0279) && (fabs(_ipPV[leptonCounter]) < 0.0351) &&              (fabs(_ipZPV[leptonCounter]) < 0.417) && (_hitsNumber[leptonCounter] <= 1) && (!_vtxFitConversion[leptonCounter]);
        }
        //_isloose[leptonCounter] = ( _ipPV[leptonCounter]< 0.05) && (_ipZPV[leptonCounter]< 0.1);
        //**************************************************************************************
        // Assign "clean" flag to electrons
        //**************************************************************************************
        /*
        bool Remove = false;
        for (unsigned int ii = 0 ; ii < sMu.size() ; ii++) {
            const pat::Muon *mu = sMu[ii];
            double isolation_loc = pfRelIso(mu,myRhoJets);
            if (isolation_loc < 1) {
                TLorentzVector mu1; mu1.SetPtEtaPhiE(mu->pt(),mu->eta(),mu->phi(),mu->energy());
                TLorentzVector el; el.SetPtEtaPhiE(iE->pt(),iE->eta(),iE->phi(),iE->energy());
                //std::cout << "Muon parameters: "<<mu->pt()<<" "<<mu->eta()<<" "<<mu->phi()<<std::endl;
                float dR = mu1.DeltaR(el);
                std::cout<<"Delta R between ele and muon is: "<<dR<<std::endl;
                if( dR < 0.05 ) { //for sync
                //if( dR < 0.1 ) {
                    Remove = true;
                    break;
                }
            }
        }
        if (Remove) continue; // remove 25 Jun for new MVA check
        */
        //**************************************************************************************
        //std::cout << "Loose electron is clean: " << iE->pt() << endl;

        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iE->pt(), iE->eta(), iE->phi(), iE->energy());
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();

        fillCloseJetVars(leptonCounter, PV);
        
        if (Sample=="ElectronsMC") {
            //**************************************************************************************
            // MC
            //**************************************************************************************
            
            const GenParticle* mc = GPM.matchedMC(&*iE);
            if ( mc!=0 ) {
                fillMCVars(mc, leptonCounter);
                //Vertex::Point PVmc = mcMom->vertex();
                _ipPVmc[leptonCounter] = TMath::Abs(iE->gsfTrack()->dxy(PVmc));
            }
            else {
                //std::cout<<"No match e"<<std::endl;
                _origin[leptonCounter] = -1;
                _originReduced[leptonCounter] = -1;
                _mompt[leptonCounter] = 0;
                _momphi[leptonCounter] = 0;
                _mometa[leptonCounter] = 0;
                _mompdg[leptonCounter] = 0;
            }

            //fillIsoMCVars(leptonCounter);
            //matchCloseJet(leptonCounter);
            
            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }
        
        /*
        if (_regression)
            fillRegVars(SelectedJetsAll.at(_closeIndex[leptonCounter]), _closeJetPtAllMC[leptonCounter], iE);
            */
        
        leptonCounter++;
        
    }

    _nEle = leptonCounter-_nMu;
    
    /*
    cout << "Tau selection started" << endl;

    for(unsigned int i = 0 ; i < sTau.size() ;i++ ){
        const pat::Tau *iTau = sTau[i];
        if (leptonCounter == 10) continue;

        _flavors[leptonCounter] = 2;
        _charges[leptonCounter] = iTau->charge();
        _isolation[leptonCounter] = 0;
        _ipPV[leptonCounter] = iTau->dxy();
        _ipPVsig[leptonCounter] = iTau->dxy_Sig();

        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iTau->pt(), iTau->eta(), iTau->phi(), iTau->energy());
        
        bool tau_dz_bool = TMath::Abs(Tau_dz(iTau->leadChargedHadrCand()->vertex(), *((TLorentzVector *)_leptonP4->At(leptonCounter)), PV)) < 0.2;
       
        _decayModeFinding[leptonCounter] = iTau->tauID("decayModeFinding");
        _looseMVA_dR03[leptonCounter] = iTau->tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
        _mediumMVA_dR03[leptonCounter] = iTau->tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");

        //_decayModeFindingOldDMs[leptonCounter] = iTau->tauID("decayModeFindingOldDMs");
        _vlooseMVAold[leptonCounter] = iTau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
        _looseMVAold[leptonCounter] = iTau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
        _mediumMVAold[leptonCounter] = iTau->tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
        _tightMVAold[leptonCounter] = iTau->tauID("byTightIsolationMVArun2v1DBoldDMwLT");
        _vtightMVAold[leptonCounter] = iTau->tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

        _decayModeFindingNewDMs[leptonCounter] = iTau->tauID("decayModeFindingNewDMs");
        _vlooseMVAnew[leptonCounter] = iTau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");
        _looseMVAnew[leptonCounter] = iTau->tauID("byLooseIsolationMVArun2v1DBnewDMwLT");
        _mediumMVAnew[leptonCounter] = iTau->tauID("byMediumIsolationMVArun2v1DBnewDMwLT");
        _tightMVAnew[leptonCounter] = iTau->tauID("byTightIsolationMVArun2v1DBnewDMwLT");
        _vtightMVAnew[leptonCounter] = iTau->tauID("byVTightIsolationMVArun2v1DBnewDMwLT");

        cout << "tau dz and dxy: " << TMath::Abs(Tau_dz(iTau->leadChargedHadrCand()->vertex(), *((TLorentzVector *)_leptonP4->At(leptonCounter)), PV)) << " " << TMath::Abs(_ipPV[leptonCounter]) << " " << _decayModeFinding[leptonCounter] << " " << _vlooseMVAold[leptonCounter] << endl; 
        if(!tau_dz_bool) continue;

        if (TMath::Abs(_ipPV[leptonCounter]) > 1000) continue;
    
        if(!(_decayModeFinding[leptonCounter] && _vlooseMVAold[leptonCounter])) continue; // my proposal
        //if(!(_decayModeFinding[leptonCounter] && _looseMVA_dR03[leptonCounter])) continue; // for sync purpose

        //std::cout << "Tau dz: " << TMath::Abs(Tau_dz(iTau->leadChargedHadrCand()->vertex(), *((TLorentzVector *)_leptonP4->At(leptonCounter)), PV)) << std::endl;

        //_istight[leptonCounter] = (iTau->pt() > _tauPt) && (iTau->eta() < _tauEta) && iTau->tauID("againstElectronTightMVA5") && iTau->tauID("againstMuonTight3") && iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") && iTau->tauID("decayModeFinding");
        //_istight[leptonCounter] = (iTau->pt() > _tauPt) && (iTau->eta() < _tauEta) && iTau->tauID("againstElectronLooseMVA5") && iTau->tauID("againstMuonLoose3") && iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") && iTau->tauID("decayModeFindingNewDMs") && tau_dz_bool ;

        //std::cout << "Tau info: " << (iTau->pt() > _tauPt) << " "  << (iTau->eta() < _tauEta) << " " << iTau->tauID("againstElectronTightMVA5") << " " << iTau->tauID("againstMuonTight3") << " " << iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") << " " << iTau->tauID("decayModeFinding") << " " << tau_dz_bool <<std::endl;
        //std::cout << "Tau info: " << (iTau->pt() > _tauPt) << " "  << (iTau->eta() < _tauEta) << " " << iTau->tauID("againstElectronLooseMVA5") << " " << iTau->tauID("againstMuonLoose3") << " " << iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") << " " << iTau->tauID("decayModeFindingNewDMs") << " " << tau_dz_bool <<std::endl;
        
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();

        cout << "Tau before cleaning: " << iTau->pt() << " " << iTau->eta() << " " << iTau->phi() << endl;
        
        //**************************************************************************************
        // Assign "clean" flag to taus
        //**************************************************************************************
        bool Remove = false;
        for (int itMu = 0; itMu != _nMu+_nEle; ++itMu) {
            if (_isloose[itMu]) {
                float dR = ((TLorentzVector *)_leptonP4->At(itMu))->DeltaR( *((TLorentzVector *)_leptonP4->At(leptonCounter)) );
                std::cout<< "dR between tau and light lepton: " << dR << " " << ((TLorentzVector *)_leptonP4->At(itMu))->Pt() <<std::endl;
                if( dR < 0.4 ) {
                    Remove = true;
                    break;
                }
            }
        }
        if (Remove) continue;

        cout << "Tau is clean: " << iTau->pt() << endl;

        _isloose[leptonCounter] = true;
        _istight[leptonCounter] = true;
        //**************************************************************************************
        cout << "tau pt/eta: " << _lPt[leptonCounter] << " " <<  _lEta[leptonCounter] << endl;
        
        fillCloseJetVars(leptonCounter, PV);
        
        if (Sample=="ElectronsMC" && TheGenParticles.isValid()) {
            //**************************************************************************************
            // MC
            //**************************************************************************************
            double dR15 = 10000;
            double dRother = 10000;
            const GenParticle* mc15 = 0;
            const GenParticle* mcother = 0;
            for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
                int id = TMath::Abs(p->pdgId());
                if (fabs(p->eta()) <= 10 ) {
                    TLorentzVector Gen2; Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                    double deltaRtest = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR(Gen2);
                    if (id == 15) {
                        if (deltaRtest < dR15) {
                            dR15 = deltaRtest;
                            mc15 = &*p;
                        }
                    } else {
                        if (deltaRtest < dRother) {
                            dRother = deltaRtest;
                            mcother = &*p;
                        }
                    }
                }
            }
            
            if (dR15 < 0.5 || dR15 < dRother) {
                fillMCVars(mc15, leptonCounter);
                _ipPVmc[leptonCounter] = iTau->dxy();
                _dR15[leptonCounter] = dR15;
                std::cout << "dR is: " << _dR15[leptonCounter] << std::endl;
            } else if (dRother < 10000) {
                fillMCVars(mcother, leptonCounter);
                _ipPVmc[leptonCounter] = iTau->dxy();
                _dR15[leptonCounter] = dRother;
                std::cout << "dR is: " << _dR15[leptonCounter] << std::endl;
            }
            
            //fillIsoMCVars(leptonCounter);
            //matchCloseJet(leptonCounter);

            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }

        
        leptonCounter++;
    }
*/
    
    //_nTau = leptonCounter - _nEle - _nMu;
    _nTau = 0;

    //std::cout<<leptonCounter<<" "<<_nEle + _nMu<<std::endl;
    //if (leptonCounter != 1) return;
    if (_nEle + _nMu < 1) return;

    /*
    int posCharge = 0;
    int negCharge = 0;

    for(int i = 0; i < _nEle + _nMu; i++){
        if(_charges[i] == -1)    
            negCharge++;
        if(_charges[i] == 1)    
            posCharge++;
    }

    if(posCharge < 2 && negCharge < 2)
        return;
        */

    /*
    if(SampleName=="DoubleMu"){
        if (!(_trigDiMuIso == true && _trigMu8Ele23Iso == false && _trigMu23Ele12Iso == false && _trigEle23Ele12Iso == false))
            return;
    }

    if(SampleName=="EGMu"){
        if (!(_trigDiMuIso == false && _trigMu8Ele23Iso == true && _trigMu23Ele12Iso == true && _trigEle23Ele12Iso == false))
            return;
    }

    if(SampleName=="DoubleEG"){
        if (!(_trigDiMuIso == false && _trigMu8Ele23Iso == false && _trigMu23Ele12Iso == false && _trigEle23Ele12Iso == true))
            return;
    }
    */


     
    _nLeptons = leptonCounter;
    //std::cout<<"lepton counter "<<std::endl;
    
    _sb = false;
    
    //std::cout<<((TLorentzVector*)_leptonP4->At(_index1))->Pt()<<" "<<((TLorentzVector*)_leptonP4->At(_index2))->Pt()<<std::endl;
    //std::cout<<_isolation[_index1]<<" "<<_isolation[_index2]<<std::endl;
    
    _n_Jets = 0;
    _n_bJets = 0;
    HT = 0;
    //std::cout<<"jets "<<SelectedJets.size()<<std::endl;
    //bool jetDeltaR1 = false;
    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){
        

        double uncPt = (SelectedJets[i]->correctedP4("Uncorrected")).Pt();
        double uncEta = (SelectedJets[i]->correctedP4("Uncorrected")).Eta();
         
        //double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),_corrLevel);
        double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJECJets, SelectedJets[i]->jetArea(),_corrLevel);
        //std::cout << "Before correction pt of jet: " << uncPt << " " << uncPt*corr << std::endl;
                                 
        //if (uncPt*corr < 30) continue;
        if (uncPt*corr < 25) continue;
                                         
        _jetEta[_n_Jets] = SelectedJets[i]->eta();
        _jetPhi[_n_Jets] = SelectedJets[i]->phi();
        _jetPt[_n_Jets] = uncPt*corr;
        _jetE[_n_Jets] = (SelectedJets[i]->correctedP4("Uncorrected")).E()*corr;

        //std::cout << "After correction pt of jet: " << _jetPt[_n_Jets] << std::endl;

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
        //if (!clean) continue; // Not cleaning jets with leptons 25 Jun 2016
        
        //std::cout << "Jet with pt: " << _jetPt[_n_Jets] << " is clean" <<std::endl;
        //std::cout<<"clean: "<<_jetPt[_n_Jets]<<" "<<_jetEta[_n_Jets]<<" "<<_jetPhi[_n_Jets]<<std::endl;

        for(int j=0; j != _nLeptons; ++j){
            _jetDeltaR[_n_Jets][j] = ((TLorentzVector *)_leptonP4->At(j))->DeltaR( jt ) ;
            /*
            if(_jetDeltaR[_n_Jets][j] > 1.)
                jetDeltaR1 = true;
                */
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

        /*
        if(!realdata_){
            bool matchGen=false;
            if (SelectedJets[i]->genJet()){
                matchGen=true;
                _matchedjetPt[_n_Jets] = SelectedJets[i]->genJet()->pt();
                _matchedjetEta[_n_Jets] =  SelectedJets[i]->genJet()->eta();
                _matchedjetPhi[_n_Jets] =  SelectedJets[i]->genJet()->phi();
                _matchedjetE[_n_Jets] = SelectedJets[i]->genJet()->energy();
                _matchedjetM[_n_Jets] = SelectedJets[i]->genJet()->mass();
           
            }
            else{
                _matchedjetPt[_n_Jets] = -999;
                _matchedjetEta[_n_Jets] =  -999;
                _matchedjetPhi[_n_Jets] =  -999;
                _matchedjetE[_n_Jets] = -999;
                _matchedjetM[_n_Jets] = -999;
            
           }
           _matchGjet[_n_Jets] = matchGen;
        }

        for(int ibtagSF = 0; ibtagSF < 19; ibtagSF++)
            _btagSF[ibtagSF][_n_Jets] = 0.;

        if(!realdata_){
                        
            double pt  = SelectedJets[i]->pt();
            double eta = SelectedJets[i]->eta();
            double csv = SelectedJets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            double flavor = SelectedJets[i]->hadronFlavour();
            if( csv < 0.0 ) csv = -0.05;
            if( csv > 1.0 ) csv = 1.0;
                                                                                                
            if( pt > 1000 ) pt = 999.;
            
            bool isBFlav = false;
            bool isCFlav = false;
            bool isLFlav = false;
            if( abs(flavor)==5 )      isBFlav = true;
            else if( abs(flavor)==4 ) isCFlav = true;
            else                      isLFlav = true;

            BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;
            if( isBFlav )       jf = BTagEntry::FLAV_B;
            else if( isCFlav ) jf = BTagEntry::FLAV_C;
            else                  jf = BTagEntry::FLAV_UDSG;
                                                            
            
            _btagSF[0][_n_Jets] = reader->eval(jf, eta, pt, csv);
                                                                        
            if(!isCFlav){
                _btagSF[1][_n_Jets] = reader_JESUp->eval(jf, eta, pt, csv);
                _btagSF[2][_n_Jets] =  reader_JESDown->eval(jf, eta, pt, csv);
            }

            else{
                _btagSF[1][_n_Jets] = 0.;
                _btagSF[2][_n_Jets] = 0.;
            }

            if(isBFlav){
                _btagSF[3][_n_Jets] = reader_LFUp->eval(jf, eta, pt, csv);
                _btagSF[4][_n_Jets] = reader_LFDown->eval(jf, eta, pt, csv);
                _btagSF[7][_n_Jets] = reader_HFStats1Up->eval(jf, eta, pt, csv);
                _btagSF[8][_n_Jets] = reader_HFStats1Down->eval(jf, eta, pt, csv);
                _btagSF[9][_n_Jets] = reader_HFStats2Up->eval(jf, eta, pt, csv);
                _btagSF[10][_n_Jets] = reader_HFStats2Down->eval(jf, eta, pt, csv);
            }

            if(isLFlav){
                _btagSF[5][_n_Jets] = reader_HFUp->eval(jf, eta, pt, csv);
                _btagSF[6][_n_Jets] = reader_HFDown->eval(jf, eta, pt, csv);
                _btagSF[11][_n_Jets] = reader_LFStats1Up->eval(jf, eta, pt, csv);
                _btagSF[12][_n_Jets] = reader_LFStats1Down->eval(jf, eta, pt, csv);
                _btagSF[13][_n_Jets] = reader_LFStats2Up->eval(jf, eta, pt, csv);
                _btagSF[14][_n_Jets] = reader_LFStats2Down->eval(jf, eta, pt, csv);
            }

            if(isCFlav){
                _btagSF[15][_n_Jets] = reader_CFErr1Up->eval(jf, eta, pt, csv);
                _btagSF[16][_n_Jets] = reader_CFErr1Down->eval(jf, eta, pt, csv);
                _btagSF[17][_n_Jets] = reader_CFErr2Up->eval(jf, eta, pt, csv);
                _btagSF[18][_n_Jets] = reader_CFErr2Down->eval(jf, eta, pt, csv);
            
            }

        }
    */

        if(_csv[_n_Jets] > 0.89) {
            _bTagged[_n_Jets] = true;
            _n_bJets++;
        } else _bTagged[_n_Jets] = false;
        
        if (_jetPt[_n_Jets] > 30)
            HT+= _jetPt[_n_Jets];
        _n_Jets++;
    }

    //if(!jetDeltaR1) return;
    //if(_n_Jets < 1) return;
    
    std::cout<<_runNb<<" "<<_lumiBlock<<" "<<_eventNb<<" "<<_nMu<<" "<<_nEle<<" "<<_nTau<<" "<<_n_Jets<<" "<<_n_bJets<<std::endl;
    outputTree->Fill();
}

void effFast::fillMCVars(const GenParticle* mc, const int leptonCounter) {
    
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
void effFast::fillCloseJetVars(const int leptonCounter, Vertex::Point PV) {

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
           //cout << "jet closest pt candidate: " << pJet.Pt() << endl;
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
           //cout << "trackMultiplicity: " << trackSelectionMult << endl;
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

void effFast::matchCloseJet(const int leptonCounter) {
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


void effFast::fillIsoMCVars(const int leptonCounter) {
    
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

void effFast::fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu) {
    
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

void effFast::fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* mu) {
    
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


DEFINE_FWK_MODULE(effFast);

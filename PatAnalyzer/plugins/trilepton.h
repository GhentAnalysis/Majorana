#ifndef trilepton_H
#define trilepton_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "Majorana/PatAnalyzer/interface/GenParticleManager.h"
#include "Majorana/PatAnalyzer/interface/Statistics.h"
#include "Majorana/PatAnalyzer/interface/Tools.h"
#include "Majorana/PatAnalyzer/interface/OnTheFlyCorrections.hh"
#include "Majorana/PatAnalyzer/interface/BTagCalibrationStandalone.h"


#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


//Root Classes

#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TClonesArray.h"

//Standard C++ classes
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <memory>
#include <iomanip>

const int nLeptonsMax = 30;
const int nLeptonsMax_gen = 30;
const int nMajo_gen = 20;
const int nNu_gen = 20;
const int nJetsMax = 20;

class trilepton : public edm::EDAnalyzer {
public:
    
    explicit trilepton(const edm::ParameterSet & iConfig);
    ~trilepton(){};
    
private:
    
    //virtual void analyze(edm::Event & iEvent, const edm::EventSetup & iSetup);
    virtual void beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup){};
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob(void);
    
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu);
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* el);

    void fillMCVars(const GenParticle* mc, const int leptonCounter);
    void fillCloseJetVars(const int leptonCounter, Vertex::Point PV);
    void matchCloseJet(const int leptonCounter);
    void fillIsoMCVars(const int leptonCounter);
    void getTriggerResults(const edm::Event& iEvent, bool isHLT, edm::EDGetTokenT<edm::TriggerResults> token, std::vector<TString> toSave);

    std::vector<const pat::Jet* > SelectedJetsAll;
    //edm::Handle<GenParticleCollection> TheGenParticles;
    edm::Handle<std::vector<reco::GenParticle>> TheGenParticles;
    
    
    edm::InputTag IT_tauDiscriminator;
    edm::InputTag IT_METFilters;
    
    edm::EDGetTokenT<reco::GenParticleCollection>    genparticleToken;
    edm::EDGetTokenT<edm::ValueMap<float>>           electronMvaIdMapToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>  	     electronMvaIdMap80Token;
    edm::EDGetTokenT<edm::ValueMap<bool>>  	     electronMvaIdMap90Token;
    edm::EDGetTokenT<edm::ValueMap<bool>>            electronCutBasedIdMapTightToken;
    edm::EDGetTokenT<edm::ValueMap<bool>>            electronCutBasedIdMapMediumToken;
    edm::EDGetTokenT<GenEventInfoProduct>            pdfvariablesToken;
    edm::EDGetTokenT<reco::BeamSpot>                 IT_beamspot;
    edm::EDGetTokenT<vector<PileupSummaryInfo>>      PileUpToken; 
    edm::EDGetTokenT<std::vector<Vertex>>            goodOfflinePrimaryVerticesToken;
    edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFCandidatesToken;
    edm::EDGetTokenT<std::vector<pat::Muon>>         IT_muon;
    edm::EDGetTokenT<edm::View<pat::Electron>>       IT_electron;
    edm::EDGetTokenT<std::vector<pat::Jet>>          IT_jet;
    edm::EDGetTokenT<std::vector<pat::MET>>          IT_pfmet;
    edm::EDGetTokenT<std::vector<reco::Conversion>>  reducedEgammaToken;
    edm::EDGetTokenT<std::vector<pat::Tau>>          IT_tau;
    edm::EDGetTokenT<double>                         fixedGridRhoFastjetCentralNeutralToken;
    edm::EDGetTokenT<double>                         fixedGridRhoFastjetAllToken;
    edm::EDGetTokenT<edm::TriggerResults>            triggerResultsHLTToken;
    edm::EDGetTokenT<edm::TriggerResults>            triggerResultsRECOToken;
    edm::EDGetTokenT<pat::PackedTriggerPrescales>    triggerPrescalesToken;
    edm::EDGetTokenT<LHEEventProduct>                IT_externalLHEProducer;
    bool isData;
    bool treeForFakeRate;
    std::string SampleName;


    edm::Service<TFileService> fs;
    FILE *outfile;
    
    TH1F *Nvtx;
    
    //desired output variables
    TTree* outputTree;
    
    string _corrLevel;
    

    
    double _relIsoCutE;
    double _relIsoCutMu;
    double _relIsoCutEloose;
    double _relIsoCutMuloose;
    
    bool _chargeConsistency;
  
    double _minPt0;
    double _minPt1;
    double _tightD0Mu;
    double _tightD0E;
    double _looseD0Mu;
    double _looseD0E;
    
    double _jetPtCut;
    double _jetEtaCut;
    
    double _tauPt;
    double _tauEta;
    
    bool _regression;
    
    bool firstEvent_;
    
    std::vector<std::string> myManualCatWeigths;
    vector<string> myManualCatWeigthsTrig;
    

    
    //genlevel particles
    GenParticleManager GPM;
    OnTheFlyCorrections* fMetCorrector;
    
    int _n_bJets;
    int _n_Jets;

    
    double _jetEta[nJetsMax];
    double _jetPhi[nJetsMax];
    double _jetPt[nJetsMax];
    double _jetE[nJetsMax];
    int _jetFlavour[nJetsMax];
    bool _bTagged[nJetsMax];
    double _csv[nJetsMax];
    double _jetDeltaR[nJetsMax][nLeptonsMax];

    double _weight;
    double _LHEweight[111];
    double _LHEweightID[111];
    double _mgluino;
    double _mchi0;
    
    int _n_bJetsAll;
    int _n_JetsAll;
    
    double _jetEtaAll[100];
    double _jetPhiAll[100];
    double _jetPtAll[100];
    double _jetEAll[100];
    bool _bTaggedAll[100];
    double _csvAll[100];
    int _closeIndex[nLeptonsMax];
    
    
    TClonesArray* _leptonP4;
    TClonesArray* _gen_leptonP4;
    TClonesArray* _gen_nuP4;
    TClonesArray* _gen_wP4;
    TClonesArray* _gen_majoP4;
    TClonesArray* _gen_leptonP4_var;
    TClonesArray* _gen_nuP4_var;
    TClonesArray* _gen_majoP4_var;



    TClonesArray* _jetP4;
    TClonesArray* _jetAllP4;
    
    
    int _nLeptons;
    int _nEle;
    int _nMu;
    int _nTau;
    int _nMajorana;
    
    int _eventType; //ee,mm,em
    // eWe 
    int _index1 = -1;
    int _index2 = -1;

    
    int _indeces[nLeptonsMax];
    int _flavors[nLeptonsMax];
    double _charges[nLeptonsMax];
    double _isolation[nLeptonsMax];
    double _isolation_absolute[nLeptonsMax];
    double _isolationComponents[nLeptonsMax][4];
    double _isolationMC[nLeptonsMax][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    double _miniisolation[nLeptonsMax];
    double _multiIsolationVT[nLeptonsMax];
    double _miniisolationCharged[nLeptonsMax];
    double _ptrel[nLeptonsMax];
    double _ptratio[nLeptonsMax];
    double _muonSegmentComp[nLeptonsMax];
    bool _passedCutBasedIdTight[nLeptonsMax];
    bool _passedCutBasedIdMedium[nLeptonsMax];
    bool _passedMVA80[nLeptonsMax];
    bool _passedMVA90[nLeptonsMax];
    bool _passedMVA_SUSY[nLeptonsMax][3];
    Int_t _findMatched[nLeptonsMax];
    

    double _closeJetEtaAll[nLeptonsMax];
    double _closeJetPhiAll[nLeptonsMax];
    double _closeJetEAll[nLeptonsMax];
    double _closeJetMAll[nLeptonsMax];
    double _closeJetNconstAll[nLeptonsMax];
    double _closeJetCSVAll [nLeptonsMax];


    double myRhoJets;
    double myRhoJECJets;

    int _hitsNumber[nLeptonsMax];
    bool _vtxFitConversion[nLeptonsMax];

    bool _trigEmulator[nLeptonsMax];
    bool _isotrigEmulator[nLeptonsMax];

    double _mvaValue[nLeptonsMax];

    std::map<TString, std::vector<bool>*>   leptonWorkingPoints;
    std::map<TString, std::vector<double>*> leptonConeCorrectedPt;

    int _origin[nLeptonsMax];
    int _originReduced[nLeptonsMax];
    bool _isPromptFinalState[nLeptonsMax];
    bool _fromHardProcessFinalState[nLeptonsMax];

    double _PVchi2;
    double _PVerr[3];
    double _ipPV[nLeptonsMax];
    double _ipPVsig[nLeptonsMax];
    double _ipPVerr[nLeptonsMax];
    double _ipPVmc[nLeptonsMax];
    
    double _ipZPV[nLeptonsMax];
    double _ipZPVerr[nLeptonsMax];
    
    double _3dIP[nLeptonsMax];
    double _3dIPerr[nLeptonsMax];
    double _3dIPsig[nLeptonsMax];
    
    double _mt[nLeptonsMax];
    
    double _srID[3];
    double _channel[3];
    double _mll[3];
    bool _ossf[3];
    bool _sb;
    bool _doubleF;

    double _closeJetPt[nLeptonsMax];
    double _closeJetPtAll[nLeptonsMax];

    int _trackSelectionMultiplicity[nLeptonsMax];
    
    double _closeJetAngAll[nLeptonsMax];
    double _ptRel[nLeptonsMax];
    double _ptRelAll[nLeptonsMax];
    double _closeJetPtAllMC[nLeptonsMax];
    double _closeJetPtAllstatus[nLeptonsMax];
    int _partonIdMatched[nLeptonsMax];
    bool _sameParton[nLeptonsMax];
    
    bool _isloose[nLeptonsMax];
    bool _ismedium[nLeptonsMax];
    bool _istight[nLeptonsMax];
    bool _isvetoIDCutBased[nLeptonsMax];
    bool _islooseIDCutBased[nLeptonsMax];
    bool _ismediumIDCutBased[nLeptonsMax];
    bool _istightIDCutBased[nLeptonsMax];

    bool _chargeConst[nLeptonsMax];


    double _n_MCTruth_PV;
    int _n_PV;
    
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
    
   
   
    //Generator variables
    int _gen_lmompdg[nLeptonsMax];		//Mothers of final state leptons
    int _gen_nL;
    double _gen_lPt[nLeptonsMax];		//Kinematics of final state leptons
    double _gen_lE[nLeptonsMax];
    double _gen_lPhi[nLeptonsMax];
    double _gen_lEta[nLeptonsMax];
    int _gen_flavors[nLeptonsMax];
    double _gen_charges[nLeptonsMax];

    //Generator variables
    int  _gen_nNu ;
    double _gen_nuPt[nLeptonsMax];		//Kinematics of final state leptons
    double _gen_nuE[nLeptonsMax];
    double _gen_nuPhi[nLeptonsMax];
    double _gen_nuEta[nLeptonsMax];
    int _gen_numompdg[nLeptonsMax];

    //Generator variables
    int  _gen_nW ;
    double _gen_wPt[nLeptonsMax];		//Kinematics of final state leptons
    double _gen_wE[nLeptonsMax];
    double _gen_wPhi[nLeptonsMax];
    double _gen_wEta[nLeptonsMax];
    int _gen_wmompdg[nLeptonsMax];


    //Generator variables
    int  _gen_nMajo ;
    double _gen_majoPt[nLeptonsMax];		//Kinematics of final state leptons
    double _gen_majoE[nLeptonsMax];
    double _gen_majoPhi[nLeptonsMax];
    double _gen_majoEta[nLeptonsMax];


    double _lPt[nLeptonsMax], _lEta[nLeptonsMax], _lPhi[nLeptonsMax], _lE[nLeptonsMax];
    double _ePt[nLeptonsMax], _eEta[nLeptonsMax], _ePhi[nLeptonsMax], _eE[nLeptonsMax];
    double _lPtmc[nLeptonsMax], _lEtamc[nLeptonsMax], _lPhimc[nLeptonsMax], _lEmc[nLeptonsMax], _lpdgmc[nLeptonsMax], _lchargemc[nLeptonsMax];
    double _ePtmc[nLeptonsMax], _eEtamc[nLeptonsMax], _ePhimc[nLeptonsMax], _eEmc[nLeptonsMax], _epdgmc[nLeptonsMax], _echargemc[nLeptonsMax];
    bool _flag_lepton[nLeptonsMax];
    double _dR15[nLeptonsMax];
    double _nuPtmc[nLeptonsMax], _nuEtamc[nLeptonsMax], _nuPhimc[nLeptonsMax], _nuEmc[nLeptonsMax];

    double _mtmc[nLeptonsMax];
    
    double _mompt[nLeptonsMax];
    double _momphi[nLeptonsMax];
    double _mometa[nLeptonsMax];
    int _mompdg[nLeptonsMax];

    double _met;
    double _met_phi;
    double HT;
    
    double _genmet;
    double _genmet_phi;
    
    double _genqpt;
    double _genqpt20;

    long _nEventsTotal;
    long _nEventsTotalCounted;
    long _nEventsFiltered;
    
    Vertex::Point PVmc;
    
    TH1D* _hCounter;

    
    double _regVars[15];
    double hJet_ptRaw;
    double hJet_genPt;
    double hJet_pt;
    double hJet_phi;
    double hJet_eta;
    double hJet_e;
    
    double hJet_ptLeadTrack;
    
    double hJet_vtx3dL;
    double hJet_vtx3deL;
    double hJet_vtxMass;
    double hJet_vtxPt;
    
    double hJet_cef;
    double hJet_nconstituents;
    double hJet_JECUnc;
    
    double hJet_SoftLeptptRel;
    double hJet_SoftLeptPt;
    double hJet_SoftLeptdR;
    
    double hJet_SoftLeptIdlooseMu;
    double hJet_SoftLeptId95;

   
    std::vector<TString>    filtersToSave;
    std::vector<TString>    triggersToSave;
    std::map<TString, bool> triggerFlags;
    std::map<TString, int>  triggerPrescales;
    std::map<TString, int>  triggerIndices;

    bool Flag_eeBadScFilter;
    bool Flag_HBHENoiseFilter;
    bool Flag_HBHENoiseIsoFilter;
    bool Flag_EcalDeadCellTriggerPrimitiveFilter;
    bool Flag_goodVertices;
    bool Flag_globalTightHalo2016Filter;

    //HLTConfigProvider hltConfig_;


    // JEC
    JetCorrectionUncertainty *jecUnc;

    // is Real Data
    bool realdata_;

    // JEC
    double _jecUnc[nJetsMax];
    double _jetPtUp[nJetsMax];
    double _jetPtDown[nJetsMax];

    // JER
    double _matchedjetPt[nJetsMax];
    double _matchedjetEta[nJetsMax];
    double _matchedjetPhi[nJetsMax];
    double _matchedjetE[nJetsMax];
    double _matchedjetM[nJetsMax];
    bool _matchGjet[nJetsMax];

    int _clean[nJetsMax];

    // b-tag SF
    double _btagSF[19][nJetsMax];

    BTagCalibration* calib_csvv2;
    BTagCalibrationReader* reader;
    BTagCalibrationReader* reader_JESUp;
    BTagCalibrationReader* reader_JESDown;
    BTagCalibrationReader* reader_LFUp;
    BTagCalibrationReader* reader_LFDown;
    BTagCalibrationReader* reader_HFUp;
    BTagCalibrationReader* reader_HFDown;
    BTagCalibrationReader* reader_HFStats1Up;
    BTagCalibrationReader* reader_HFStats1Down;
    BTagCalibrationReader* reader_HFStats2Up;
    BTagCalibrationReader* reader_HFStats2Down;
    BTagCalibrationReader* reader_LFStats1Up;
    BTagCalibrationReader* reader_LFStats1Down;
    BTagCalibrationReader* reader_LFStats2Up;
    BTagCalibrationReader* reader_LFStats2Down;
    BTagCalibrationReader* reader_CFErr1Up;
    BTagCalibrationReader* reader_CFErr1Down;
    BTagCalibrationReader* reader_CFErr2Up;
    BTagCalibrationReader* reader_CFErr2Down;

    int _nZboson;

};

#endif

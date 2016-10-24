#ifndef eff13_H
#define eff13_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
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


#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"
#include "DataFormats/PatCandidates/interface/Tau.h"


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

using namespace std;

const int nLeptonsMax = 6;

class eff13 : public edm::EDAnalyzer {
public:
    
    explicit eff13(const edm::ParameterSet & iConfig);
    ~eff13(){};
    
private:
    
    //virtual void analyze(edm::Event & iEvent, const edm::EventSetup & iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob(void);
    
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu);
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* el);

    void fillMCVars(const GenParticle* mc, const int leptonCounter);
    void fillCloseJetVars(const int leptonCounter);
    void matchCloseJet(const int leptonCounter);
    void fillIsoMCVars(const int leptonCounter);

    std::vector<const pat::Jet* > SelectedJetsAll;
    edm::Handle<GenParticleCollection> TheGenParticles;
    
    
    std::string Sample;
    edm::InputTag IT_muon;
    edm::InputTag IT_electron;
    edm::InputTag IT_tau;
    edm::InputTag IT_tauDiscriminator;
    edm::InputTag IT_jet;
    edm::InputTag IT_pfmet;
    edm::InputTag IT_beamspot;
    edm::InputTag IT_hltresults;
    edm::InputTag IT_METFilters;
    
    // MVA values and categories (optional)
    edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;


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
    EGammaMvaEleEstimatorCSA14* myMVATrig;
    double looseMVA[3][2]; //{{0.35, 0.20, -0.52}, {0.73, 0.57, 0.05}};//{0.8, 1.479, };
    

    //genlevel particles
    GenParticleManager GPM;
    OnTheFlyCorrections* fMetCorrector;
    
    int _n_bJets;
    int _n_Jets;
    
    double _jetEta[20];
    double _jetPhi[20];
    double _jetPt[20];
    double _jetE[20];
    bool _bTagged[20];
    double _csv[20];

    double _weight;
    
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
    TClonesArray* _jetP4;
    TClonesArray* _jetAllP4;
    
    int _nLeptons;
    int _nEle;
    int _nMu;
    int _nTau;
    
    int _eventType; //ee,mm,em
    int _index1 = -1;
    int _index2 = -1;

    
    int _indeces[nLeptonsMax];
    int _flavors[nLeptonsMax];
    double _charges[nLeptonsMax];
    double _isolation[nLeptonsMax];
    double _isolationDB[nLeptonsMax];
    double _isolationComponents[nLeptonsMax][4];
    double _isolationMC[nLeptonsMax][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    double _miniisolation[nLeptonsMax][2];
    double _activity[nLeptonsMax];
    double _miniisolationDB[nLeptonsMax];
    bool _multiisolation[nLeptonsMax][5];
    double multiConst[5][3];
    double _ptrel[nLeptonsMax];
    double _ptratio[nLeptonsMax];
    double _ptrel_unc[nLeptonsMax];

    double myRhoJets;
    double JECRhoJets;

    int _hitsNumber[nLeptonsMax];

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
    double _closeJetEtaAll[nLeptonsMax];
    double _closeJetPhiAll[nLeptonsMax];
    double _closeJetEAll[nLeptonsMax];
    double _closeJetMAll[nLeptonsMax];
    double _closeJetNconstAll[nLeptonsMax];
    double _closeJetCSVAll [nLeptonsMax];

    double _closeJetPtAll_unc[nLeptonsMax];
    
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
    bool _istightID[nLeptonsMax];
    bool _isvetoIDCutBased[nLeptonsMax];
    bool _islooseIDCutBased[nLeptonsMax];
    bool _ismediumIDCutBased[nLeptonsMax];
    bool _istightIDCutBased[nLeptonsMax];

    bool _chargeConst[nLeptonsMax];
    bool _vtxFitConversion[nLeptonsMax];

    bool _trigEmulator[nLeptonsMax];
    bool _isotrigEmulator[nLeptonsMax];

    double _mvaValue[nLeptonsMax];

    double _RhoCN;
    double _RhoAll;

    int _n_PV;
    
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
    
    
    double _lPt[nLeptonsMax], _lEta[nLeptonsMax], _lPhi[nLeptonsMax], _lE[nLeptonsMax];
    double _lPtmc[nLeptonsMax], _lEtamc[nLeptonsMax], _lPhimc[nLeptonsMax], _lEmc[nLeptonsMax], _lpdgmc[nLeptonsMax];
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

    bool _lgenParticle[nLeptonsMax];

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
};

#endif

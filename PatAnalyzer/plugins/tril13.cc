#include "tril13.h"
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

#include "PhysicsTools/CandUtils/interface/CandMatcherNew.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace tools;
using namespace math;
using namespace reco::tau;

tril13::tril13(const edm::ParameterSet & iConfig) :
mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
//eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
//eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
_relIsoCutE(999.),//0.15
_relIsoCutMu(999.),//0.15
_relIsoCutEloose(999.), //0.5
_relIsoCutMuloose(999.), //0.5
_chargeConsistency(false),
_minPt0(5.),
_minPt1(7.),
_tightD0Mu(0.05),
_tightD0E(0.05),
_looseD0Mu(0.05), // 0.05 for sync
_looseD0E(0.05), // 0.05 for sync
//_looseD0Mu(0.2),
//_looseD0E(9999999.),
//_jetPtCut(40.),
_jetPtCut(30.),
_jetEtaCut(2.5),
_tauPt(20),
_tauEta(2.3),
_regression(false)
{
    Sample              = iConfig.getUntrackedParameter<std::string>("SampleLabel") ;
    IT_muon             = iConfig.getParameter<edm::InputTag>("MuonLabel") ;
    IT_electron         = iConfig.getParameter<edm::InputTag>("ElectronLabel") ;
    IT_tau              = iConfig.getParameter<edm::InputTag>("TauLabel") ;
    //IT_tauDiscriminator = iConfig.getParameter<edm::InputTag>("TauDiscriminatorLabel") ;
    IT_jet              = iConfig.getParameter<edm::InputTag>("JetLabel");
    IT_pfmet            = iConfig.getParameter<edm::InputTag>("METLabel")  ;
    IT_beamspot         = iConfig.getParameter<edm::InputTag>("BeamSpotLabel");
    IT_hltresults       = iConfig.getParameter<edm::InputTag>("HLTResultsLabel");
    //IT_METFilters       = iConfig.getParameter<edm::InputTag>("METFilter");
    
    //outfile = fopen("FakeSync.txt", "w");
}


void tril13::beginJob()
{
    Nvtx           = fs->make<TH1F>("N_{vtx}"        , "Number of vertices;N_{vtx};events / 1"  ,    40, 0., 40.);
    
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);
    
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
    
    outputTree->Branch("_lPt", &_lPt, "_lPt[6]/D");
    outputTree->Branch("_lEta", &_lEta, "_lEta[6]/D");
    outputTree->Branch("_lPhi", &_lPhi, "_lPhi[6]/D");
    outputTree->Branch("_lE", &_lE, "_lE[6]/D");
    
    outputTree->Branch("_lPtmc", &_lPtmc, "_lPtmc[6]/D");
    outputTree->Branch("_lEtamc", &_lEtamc, "_lEtamc[6]/D");
    outputTree->Branch("_lPhimc", &_lPhimc, "_lPhimc[6]/D");
    outputTree->Branch("_lEmc", &_lEmc, "_lEmc[6]/D");
    outputTree->Branch("_lpdgmc", &_lpdgmc, "_lpdgmc[6]/D");
    outputTree->Branch("_flag_lepton", &_flag_lepton, "_flag_lepton[6]/O");
    outputTree->Branch("_dR15", &_dR15 , "_dR15[6]/D");
   
    outputTree->Branch("_nuPtmc", &_nuPtmc, "_nuPtmc[6]/D");
    outputTree->Branch("_nuEtamc", &_nuEtamc, "_nuEtamc[6]/D");
    outputTree->Branch("_nuPhimc", &_nuPhimc, "_nuPhimc[6]/D");
    outputTree->Branch("_nuEmc", &_nuEmc, "_nuEmc[6]/D");

    outputTree->Branch("_mtmc", &_mtmc, "_mtmc[6]/D");
    
    outputTree->Branch("_nLeptons", &_nLeptons, "_nLeptons/I");

    outputTree->Branch("_nEle", &_nEle, "_nEle/I");
    outputTree->Branch("_nMu", &_nMu, "_nMu/I");
    outputTree->Branch("_nTau", &_nTau, "_nTau/I");

    outputTree->Branch("_flavors", &_flavors, "_flavors[6]/I");
    outputTree->Branch("_charges", &_charges, "_charges[6]/D");
    outputTree->Branch("_indeces", &_indeces, "_indeces[6]/I");
    outputTree->Branch("_isolation", &_isolation, "_isolation[6]/D");
    outputTree->Branch("_isolationComponents", &_isolationComponents, "_isolationComponents[4][6]/D");
    outputTree->Branch("_isolationMC", &_isolationMC, "_isolationMC[4][6]/D");
    outputTree->Branch("_miniisolation", &_miniisolation, "_miniisolation[6][2]/D");
    outputTree->Branch("_multiisolation", &_multiisolation, "_multiisolation[6][5]/O");
    outputTree->Branch("_ptrel", &_ptrel, "_ptrel[6]/D");
    outputTree->Branch("_ptratio", &_ptratio, "_ptratio[6]/D");

    outputTree->Branch("_mll", &_mll, "_mll[3]/D");
    outputTree->Branch("_ossf", &_ossf, "_ossf[3]/O");

    
    //outputTree->Branch("_index1", &_index1, "_index1/I");
    //outputTree->Branch("_index2", &_index2, "_index2/I");

    //outputTree->Branch("_sb", &_sb, "_sb/O");
    //outputTree->Branch("_doubleF", &_doubleF, "_doubleF/O");
    
    outputTree->Branch("_origin", &_origin, "_origin[6]/I");
    outputTree->Branch("_originReduced", &_originReduced, "_originReduced[6]/I");
    
    outputTree->Branch("_PVchi2", &_PVchi2, "_PVchi2/D");
    outputTree->Branch("_PVerr", &_PVerr, "_PVerr[3]/D");
    
    outputTree->Branch("_ipPV", &_ipPV, "_ipPV[6]/D");
    outputTree->Branch("_ipPVerr", &_ipPVerr, "_ipPVerr[6]/D");
    outputTree->Branch("_ipZPV", &_ipZPV, "_ipZPV[6]/D");
    outputTree->Branch("_ipZPVerr", &_ipZPVerr, "_ipZPVerr[6]/D");
    
    outputTree->Branch("_ipPVmc", &_ipPVmc, "_ipPVmc[6]/D");
    
    outputTree->Branch("_3dIP", &_3dIP, "_3dIP[6]/D");
    outputTree->Branch("_3dIPerr", &_3dIPerr, "_3dIPerr[6]/D");
    outputTree->Branch("_3dIPsig", &_3dIPsig, "_3dIPsig[6]/D");
    
    
    outputTree->Branch("_mt", &_mt, "_mt[6]/D");
    outputTree->Branch("_isloose", &_isloose, "_isloose[6]/O");
    outputTree->Branch("_istight", &_istight, "_istight[6]/O");
    outputTree->Branch("_istightID", &_istightID, "_istightID[6]/O");
    outputTree->Branch("_islooseIDCutBased", &_islooseIDCutBased, "_islooseIDCutBased[6]/O");
    outputTree->Branch("_ismediumIDCutBased", &_ismediumIDCutBased, "_ismediumIDCutBased[6]/O");
    outputTree->Branch("_istightIDCutBased", &_istightIDCutBased, "_istightIDCutBased[6]/O");
    outputTree->Branch("_mvaValue", &_mvaValue, "_mvaValue[6]/D");
   
    
    outputTree->Branch("_closeJetPtAll", &_closeJetPtAll, "_closeJetPtAll[6]/D");
    outputTree->Branch("_closeJetAngAll", &_closeJetAngAll, "_closeJetAngAll[6]/D");
    outputTree->Branch("_ptRel", &_ptRel, "_ptRel[6]/D");
    outputTree->Branch("_ptRelAll", &_ptRelAll, "_ptRelAll[6]/D");

    outputTree->Branch("_hitsNumber", &_hitsNumber, "_hitsNumber[6]/I");
    
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
    
    outputTree->Branch("_met", &_met, "_met/D");
    outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
    outputTree->Branch("HT", &HT, "HT/D");
    
    outputTree->Branch("_genmet", &_genmet, "_genmet/D");
    outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/D");
    
    outputTree->Branch("_genqpt", &_genqpt, "_genqpt/D");

    outputTree->Branch("_mompt", &_mompt, "_mompt[6]/D");
    outputTree->Branch("_momphi", &_momphi, "_momphi[6]/D");
    outputTree->Branch("_mometa", &_mometa, "_mometa[6]/D");
    outputTree->Branch("_mompdg", &_mompdg, "_mompdg[6]/I");
    
    outputTree->Branch("_n_bJets", &_n_bJets, "_n_bJets/I");
    outputTree->Branch("_n_Jets", &_n_Jets, "_n_Jets/I");
    outputTree->Branch("_bTagged", &_bTagged, "_bTagged[20]/O");
    outputTree->Branch("_jetEta", &_jetEta, "_jetEta[20]/D");
    outputTree->Branch("_jetPhi", &_jetPhi, "_jetPhi[20]/D");
    outputTree->Branch("_jetPt", &_jetPt, "_jetPt[20]/D");
    outputTree->Branch("_jetM", &_jetM, "_jetM[20]/D");
    outputTree->Branch("_csv", &_csv, "_csv[20]/D");
    outputTree->Branch("_jetDeltaR", &_jetDeltaR, "_jetDeltaR[20][6]/D");
    
    outputTree->Branch("_weight", &_weight, "_weight/D");
    
    GPM = GenParticleManager();
    
    bool isData = !(Sample=="ElectronsMC");
    if (isData)
        fMetCorrector = new OnTheFlyCorrections("Summer15_25nsV2_DATA", isData); //isData = true
    else
        fMetCorrector = new OnTheFlyCorrections("Summer15_25nsV2_MC", isData); //isData = true
    _corrLevel = "L3Absolute";
    if (isData) _corrLevel = "L2L3Residual";
    
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
    
    looseMVA[0][0] = 0.35;
    looseMVA[1][0] = 0.20;
    looseMVA[2][0] = -0.52;
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
}

void tril13::endJob() {
    //outputTree -> Write();
    // store nEventsTotal and nEventsFiltered in preferred way
    std::cout<<_nEventsTotal<<std::endl;
    std::cout<<_nEventsFiltered<<std::endl;
    
    delete fMetCorrector;
    
    
}

void tril13::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{
    //bool islepton;
    if (Sample=="ElectronsMC") {
        //******************************************************************************************************************
        // Gen level particles                  ****************************************************************************
        //******************************************************************************************************************
        //iEvent.getByLabel("packedGenParticles", TheGenParticles);
        iEvent.getByLabel("prunedGenParticles", TheGenParticles);
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
            for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ )
            {
                int id = TMath::Abs(p->pdgId());
                if ( (id == 12 || id == 14 || id == 16 ) && (p->status() == 1) ) {
                    TLorentzVector Gen;
                    Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );
                    Gen0 += Gen;
                }
                if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 21 || id == 22 ) && (p->status() == 23)){
                    _genqpt += p->pt();
                }
            }
            if (Gen0.E()!=0) {
                _genmet = Gen0.Pt();
                _genmet_phi = Gen0.Phi();
            } else {
                _genmet = 0;
                _genmet_phi = 0;
            }
        }
        //******************************************************************************************************************
        //******************************************************************************************************************
        //******************************************************************************************************************
        
        //**************************************************************************************
        // MC
        //**************************************************************************************
    }
    
    /*
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
    */

    
     //edm::Handle<LHEEventProduct> pdfvariables;
     //iEvent.getByLabel("generator", pdfvariables);
     //_weight=pdfvariables.weight();

     /*
     //============= trigger ==============
     edm::Handle<TriggerResults> trigResults;
     
     iEvent.getByLabel(IT_hltresults, trigResults);
     if( trigResults.failedToGet() ) cout << "--- NO TRIGGER RESULTS !! ---" << endl;
     
     //==================================
     
     if( !trigResults.failedToGet() ) {
        unsigned int n_Triggers = trigResults->size();
        const edm::TriggerNames & triggerNames = iEvent.triggerNames(*trigResults);
        if ( firstEvent_ ) {
            edm::TriggerNames::Strings allTriggers( triggerNames.triggerNames() );
            std::cout << "--- Trigger Menu --- " << std::endl;
            for ( unsigned int i_Name = 0; i_Name < n_Triggers; ++i_Name ) {
                std::cout << allTriggers.at( i_Name ) << std::endl;
            }
        std::cout << "-------------------- " << std::endl;
        firstEvent_ = false;
        }
     }
     */

    //if (_genqpt > 100) return;
    edm::Handle<GenEventInfoProduct> pdfvariables;
    iEvent.getByLabel("generator", pdfvariables);
    Double_t weight = pdfvariables->weight();
    _weight=weight;

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
    iEvent.getByLabel( IT_beamspot, theBeamSpot );
    if( ! theBeamSpot.isValid() ) ERR( IT_beamspot ) ;
    BeamSpot::Point  BS= theBeamSpot->position();;
    //==================================
    
    //============ Primary vertices ============
    edm::InputTag IT_goodVtx = edm::InputTag("offlineSlimmedPrimaryVertices");
    edm::Handle<std::vector<Vertex> > theVertices;
    iEvent.getByLabel( "goodOfflinePrimaryVertices", theVertices) ;
    if( ! theVertices.isValid() ) ERR(IT_goodVtx ) ;
    int nvertex = theVertices->size();
    
    _n_PV = nvertex;
    Nvtx->Fill(TMath::Min(nvertex,39));
    if(! nvertex ){
        cout << "[WARNING]: No candidate primary vertices passed the quality cuts, so skipping event" << endl;
        return ;
    }
    
    Vertex::Point PV = theVertices->begin()->position();
    //std::cout<<PV.x()<<" "<<PV.y()<<" "<<PV.z()<<std::endl;
    const Vertex* PVtx = &((*theVertices)[0]);
    _PVchi2 = PVtx->chi2();
    _PVerr[0] = PVtx->xError();
    _PVerr[1] = PVtx->yError();
    _PVerr[2] = PVtx->zError();
    //==================================
    
    //============ Pat MET ============
    edm::Handle< vector<pat::MET> > ThePFMET;
    iEvent.getByLabel(IT_pfmet, ThePFMET);
    if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
    const vector<pat::MET> *pfmetcol = ThePFMET.product();
    const pat::MET *pfmet;
    pfmet = &(pfmetcol->front());
    //double rawmet = pfmet->pt();
    //double rawmetX = pfmet->px();
    //double rawmetY = pfmet->py();
    //double rawmet_phi = pfmet->phi();
    
    
    // Reco MET
    //edm::Handle< vector<reco::PFMET> > recopfMET;
    //iEvent.getByLabel("pfType1CorrectedMet", recopfMET);
    //const vector<reco::PFMET> *pfmetreco = recopfMET.product();
    //const reco::PFMET *recopfmet;
    //recopfmet = &(pfmetreco->front());
    //double pfMET= recopfmet->pt();
    //double pfMET_Phi= recopfmet->phi();
    
    //_met = recopfmet->pt();
    //_met_phi = recopfmet->phi();
    _met = pfmet->pt();
    _met_phi = pfmet->phi();
    //==================================
    
    //============ PF cand ============
    edm::Handle<pat::PackedCandidateCollection> pfcands;
    iEvent.getByLabel("packedPFCandidates", pfcands);
    //==================================
    
    //============ Pat Muons ============
    edm::Handle< std::vector<pat::Muon> > thePatMuons;
    iEvent.getByLabel( IT_muon, thePatMuons );
    if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    //==================================
    
    std::cout << "Muons:" << std::endl;
    for( std::vector<pat::Muon>::const_iterator mu = (*thePatMuons).begin() ; mu != (*thePatMuons).end() ; mu++ ){
        std::cout << mu->pt() << " " << mu->eta() << " " << mu->phi() << std::endl;
    }
    
    
    //============ Pat Electrons ============
    edm::Handle< edm::View<pat::Electron> > thePatElectrons;
    iEvent.getByLabel( IT_electron, thePatElectrons );
    if( ! thePatElectrons.isValid() ) ERR( IT_electron );
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
    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    //edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
    //iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
    iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
    //iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
 
    //============ Conversions ============
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    iEvent.getByLabel("reducedEgamma","reducedConversions", theConversions);
    //==================================
    
    
    //============ Pat Jets ============
    edm::Handle< std::vector< pat::Jet> > thePatJets;
    iEvent.getByLabel(IT_jet , thePatJets );
    if( ! thePatJets.isValid() ) ERR(IT_jet);
    //==================================
    
    edm::Handle<double> rhoJets;
    iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetCentralNeutral","") , rhoJets);//kt6PFJets
    myRhoJets = *rhoJets;
    //==================================
    
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
    std::vector<const pat::Muon* > sMu = ssbLooseMuonSelector( *thePatMuons, _minPt0, PV, _looseD0Mu);
    
    //std::vector<const pat::Electron* > sEl = ssbMVAElectronSelector( *thePatElectrons, _minPt1, PV, _looseD0E, _chargeConsistency, theConversions, BS, true);

    
    //std::cout<<sMu.size()<<" "<<sEl.size()<<std::endl;
    //if (sEl.size() + sMu.size() == 0) return;
    
    //Taus
    edm::Handle<pat::TauCollection> PatTaus;
    iEvent.getByLabel( IT_tau, PatTaus );
    std::vector<const pat::Tau* > sTau;
    for (const pat::Tau &tau : *PatTaus) {
        //
        //if (iTau.pt() < _tauPt) continue;
        //if (fabs(tau.eta()) > _tauEta) continue;
        //printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
        //       tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
        //if (! (tau.tauID("againstElectronTightMVA5")
        //       && tau.tauID("againstMuonTight3")
        //       && tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")
        //       && tau.tauID("decayModeFinding"))) continue;

        sTau.push_back(&tau);
    }
    
    SelectedJetsAll = JetSelectorAll(*thePatJets, 5., 2.6);
    std::vector<const pat::Jet* > SelectedJets = JetSelector(*thePatJets, _jetPtCut, _jetEtaCut);

    //std::cout<<sEl.size()<<" "<<sMu.size()<<" "<<sTau.size()<<std::endl;
    //if (sEl.size() + sMu.size() + sTau.size() < 3) return;
    
    
    //std::vector<const pat::Jet* > SelectedJets = JetSelector(*thePatJets, _jetPtCut, _jetEtaCut);
    //std::cout<<"Jet size "<<SelectedJetsAll.size()<<std::endl;
    //std::cout<<"Lep size "<<sMu.size()<<" "<<sEl.size()<<std::endl;
    //if (SelectedJetsAll.size() == 0) return;
    
    HT = 0.;
    std::vector< const pat::Jet* > Bjets;
    /*    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){
     _jetEta[i] = SelectedJets[i]->eta();
     _jetPhi[i] = SelectedJets[i]->phi();
     _jetPt[i] = SelectedJets[i]->pt();
     
     ((TLorentzVector *)_jetP4->At(i))->SetPtEtaPhiM( _jetPt[i], _jetEta[i], _jetPhi[i], 0 );
     
     _csv[i] = SelectedJets[i]->bDiscriminator("combinedSecondaryVertexBJetTags");
     
     if(SelectedJets[i]->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679) {
     Bjets.push_back( &*SelectedJets[i] );
     _bTagged[i] = true;
     } else _bTagged[i] = false;
     
     HT+= SelectedJets[i]->pt();
     }
     */
    _n_bJetsAll = 0;
    
    int n_bJetsAll30 = 0;
    //int bIndex[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    
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
    std::cout << " Muon selection started " << std::endl;
    for(unsigned int i = 0 ; i < sMu.size() ;i++ ){
        
        
        const pat::Muon *iM = sMu[i];
        if (iM->pt() < 10) continue;
        
        _hitsNumber[i] = 0;

        if (leptonCounter == 6) continue;

        _flavors[leptonCounter] = 1;
        _charges[leptonCounter] = iM->charge();
        _isolation[leptonCounter] = pfRelIso(iM,myRhoJets);
        
        _miniisolation[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false, myRhoJets);
        _miniisolation[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, false, myRhoJets);

        _ipPV[leptonCounter] = TMath::Abs(iM->innerTrack()->dxy(PV));
        _ipPVerr[leptonCounter] = iM->innerTrack()->dxyError();
        
        _ipZPV[leptonCounter] = iM->innerTrack()->dz(PV);
        _ipZPVerr[leptonCounter] = iM->innerTrack()->dzError();
        
        bool goodGlb = iM->isGlobalMuon() && iM->globalTrack()->normalizedChi2() < 3 && iM->combinedQuality().chi2LocalPosition < 12 && iM->combinedQuality().trkKink < 20;
        bool good = iM->innerTrack()->validFraction() >= 0.8 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
        
        //_istightID[leptonCounter] = good;
        
        _3dIP[leptonCounter]    = iM->dB(pat::Muon::PV3D);
        _3dIPerr[leptonCounter] = iM->edB(pat::Muon::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);
        
        //std::cout<<i<<" "<<_isolation[leptonCounter]<<" "<<fabs(_ipPV[leptonCounter])<<" "<< _3dIPsig[leptonCounter]<<std::endl;

        //_isloose[leptonCounter] = ( fabs(_ipPV[leptonCounter])<_looseD0Mu ) && (_isolation[leptonCounter]<_relIsoCutMuloose);
        _isloose[leptonCounter] = ( fabs(_ipPV[leptonCounter])<_looseD0Mu );
        if (_isloose[leptonCounter])
            _istightID[leptonCounter] = (good &&  _3dIPsig[leptonCounter]<4. );
        else _istightID[leptonCounter] = false;
        //std::cout<<_istight[leptonCounter]<<" "<<_isloose[leptonCounter]<<std::endl;

        std::cout << iM->pt() << " " << iM->eta() << " " << iM->phi() << std::endl;

        if (!_istightID[leptonCounter]) continue;
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iM->pt(), iM->eta(), iM->phi(), iM->energy());
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
        
        fillCloseJetVars(leptonCounter);
        
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
            fillIsoMCVars(leptonCounter);
            matchCloseJet(leptonCounter);
            
            if (_regression) {
                std::vector <double> regVars = RegressionVars(SelectedJetsAll[_closeIndex[leptonCounter]], _closeJetPtAllMC[leptonCounter], iM);
                for (int k=0; k!=15; ++k) {
                    _regVars[k] = regVars[k];
                }
                _regVars[11] = fMetCorrector->getJECUncertainty(SelectedJetsAll[_closeIndex[leptonCounter]]->pt(),SelectedJetsAll[_closeIndex[leptonCounter]]->eta());
            }
            
            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }
        
        if (_regression)
            fillRegVars(SelectedJetsAll.at(_closeIndex[leptonCounter]), _closeJetPtAllMC[leptonCounter], iM);
        
        leptonCounter++;

        
    }
    //std::cout<<leptonCounter<<" "<<((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt() <<std::endl;
    _nMu = leptonCounter;
    
    std::cout << "Electron selection started" << std::endl;
     for (size_t i = 0; i < thePatElectrons->size(); ++i){

        const auto iE = thePatElectrons->ptrAt(i);
        std::cout << "Pat electron: " << iE->pt() << " " << iE->superCluster()->eta() << " " << iE->phi() << " " << std::endl;
        if (!ssbMVAElectronSelectorPassed(iE, _minPt1, PV, _looseD0E, _chargeConsistency, theConversions, BS, false)) continue;
        
        //double mvaValueE = myMVATrig->mvaValue(*iE,false);
        
        //std::cout<<"MVA value "<<mvaValueE<<"; eta "<< iE->eta()<<std::endl;
        
        /*
        if (TMath::Abs(iE->eta()) < 0.8 ) {
            passed = mvaValueE > looseMVA[0][0];
        } else if (TMath::Abs(iE->eta()) < 1.479 ) {
            passed = mvaValueE > looseMVA[1][0];
        } else {
            passed = mvaValueE > looseMVA[2][0];
        }
        */
        
        //std::cout << "mvaValueE:\t" << mvaValueE << std::endl;
        //std::cout << "MVA Passed:\t" << passed << std::endl;

        //if (!passed) continue;
        
        if (leptonCounter == 6) continue;
        
        _flavors[leptonCounter] = 0;
        _charges[leptonCounter] = iE->charge();
        _isolation[leptonCounter] = pfRelIso(&*iE, myRhoJets);
        _ipPV[leptonCounter] = TMath::Abs(iE->gsfTrack()->dxy(PV));
        
        _miniisolation[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, false, myRhoJets);
        _miniisolation[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.3, 10., false, false, myRhoJets);
        
        //_isloose[leptonCounter] = ( fabs(_ipPV[leptonCounter])<_looseD0E ) && (_miniisolation[leptonCounter][0]<_relIsoCutEloose);
        //**************************************************************************************
        // Assign "clean" flag to electrons
        //**************************************************************************************
        /*
        bool Remove = false;
        for (unsigned int ii = 0 ; ii < sMu.size() ; ii++) {
            const pat::Muon *mu = sMu[ii];
            double iso = pfRelIso(mu,myRhoJets);
            if (iso < 0.5) {
                TLorentzVector mu1; mu1.SetPtEtaPhiE(mu->pt(),mu->eta(),mu->phi(),mu->energy());
                TLorentzVector el; el.SetPtEtaPhiE(iE->pt(),iE->eta(),iE->phi(),iE->energy());
                //std::cout << "Muon parameters: "<<mu->pt()<<" "<<mu->eta()<<" "<<mu->phi()<<std::endl;
                float dR = mu1.DeltaR(el);
                //std::cout<<"Delta R is: "<<dR<<std::endl;
                if( dR < 0.05 ) { //for sync
                //if( dR < 0.1 ) {
                    Remove = true;
                    break;
                }
            }
        }
        if (Remove) continue;
        */
        //**************************************************************************************

        _3dIP[leptonCounter]    = iE->dB(pat::Electron::PV3D);
        _3dIPerr[leptonCounter] = iE->edB(pat::Electron::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);
        
        //std::cout << "Passed loose:\t" << _isloose[leptonCounter] << std::endl;
        //std::cout << "Passed 3dIPsig:\t" << (_3dIPsig[leptonCounter] < 4) << std::endl;
        std::cout << "Electron isolation: " << _isolation[leptonCounter] << std::endl;

        /*
        bool passedMVA = false;
        _mvaValue[leptonCounter] = (*mvaValues)[iE];
        if (TMath::Abs(iE->eta()) < 0.8 ) {
            passedMVA = _mvaValue[leptonCounter]> looseMVA[0][1];
        } else if (TMath::Abs(iE->eta()) < 1.479 ) {
            passedMVA = _mvaValue[leptonCounter]> looseMVA[1][1];
        } else {
            passedMVA = _mvaValue[leptonCounter]> looseMVA[2][1];
        }
        */
        
        //bool isPassTight = (*tight_id_decisions)[iE];
        //bool isPassLoose  = (*loose_id_decisions)[iE];
        bool isPassMedium = (*medium_id_decisions)[iE];
        _istightID[leptonCounter] = isPassMedium ;
        /*
        _islooseIDCutBased[leptonCounter] = isPassLoose;
        _ismediumIDCutBased[leptonCounter] = isPassMedium;
        _istightIDCutBased[leptonCounter] = isPassTight;
        */

        //if (_isloose[leptonCounter])
        _istight[leptonCounter] = isPassMedium && _3dIPsig[leptonCounter] < 4;
        //else _istight[leptonCounter] = false;
        std::cout<<"Electron's MVA: " << _istight[leptonCounter]<<" "<<(*mvaValues)[iE]<<std::endl;
        
        if (!_istight[leptonCounter]) continue;

        std::cout << iE->pt() << " " << iE->eta() << " " << iE->phi() << std::endl;
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iE->pt(), iE->eta(), iE->phi(), iE->energy());
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
        
        fillCloseJetVars(leptonCounter);
        
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

            fillIsoMCVars(leptonCounter);
            matchCloseJet(leptonCounter);
            
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
    for(unsigned int i = 0 ; i < sTau.size() ;i++ ){
        const pat::Tau *iTau = sTau[i];
        if (leptonCounter == 6) continue;

        _flavors[leptonCounter] = 2;
        _charges[leptonCounter] = iTau->charge();
        _isolation[leptonCounter] = 0;
        _ipPV[leptonCounter] = iTau->dxy();
        _ipPVsig[leptonCounter] = iTau->dxy_Sig();
        
        //bool tau_dz_bool = TMath::Abs(Tau_dz(iTau->leadChargedHadrCand()->vertex(), *((TLorentzVector *)_leptonP4->At(leptonCounter)), PV)) < 0.2;

        //std::cout << "Tau dz: " << TMath::Abs(Tau_dz(iTau->leadChargedHadrCand()->vertex(), *((TLorentzVector *)_leptonP4->At(leptonCounter)), PV)) << std::endl;

        _istight[leptonCounter] = (iTau->pt() > _tauPt) && (iTau->eta() < _tauEta) && iTau->tauID("againstElectronTightMVA5") && iTau->tauID("againstMuonTight3") && iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") && iTau->tauID("decayModeFinding");
        //_istight[leptonCounter] = (iTau->pt() > _tauPt) && (iTau->eta() < _tauEta) && iTau->tauID("againstElectronLooseMVA5") && iTau->tauID("againstMuonLoose3") && iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") && iTau->tauID("decayModeFindingNewDMs") && tau_dz_bool ;

        //std::cout << "Tau info: " << (iTau->pt() > _tauPt) << " "  << (iTau->eta() < _tauEta) << " " << iTau->tauID("againstElectronTightMVA5") << " " << iTau->tauID("againstMuonTight3") << " " << iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") << " " << iTau->tauID("decayModeFinding") << " " << tau_dz_bool <<std::endl;
        //std::cout << "Tau info: " << (iTau->pt() > _tauPt) << " "  << (iTau->eta() < _tauEta) << " " << iTau->tauID("againstElectronLooseMVA5") << " " << iTau->tauID("againstMuonLoose3") << " " << iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") << " " << iTau->tauID("decayModeFindingNewDMs") << " " << tau_dz_bool <<std::endl;
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iTau->pt(), iTau->eta(), iTau->phi(), iTau->energy());
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();

        */
        
        //**************************************************************************************
        // Assign "clean" flag to taus
        //**************************************************************************************
        /*
        bool Remove = false;
        for (int itMu = 0; itMu != _nMu+_nEle; ++itMu) {
            if (_istight[itMu]) {
                float dR = ((TLorentzVector *)_leptonP4->At(itMu))->DeltaR( *((TLorentzVector *)_leptonP4->At(leptonCounter)) );
                //std::cout<<dR<<std::endl;
                if( dR < 0.3 ) {
                    Remove = true;
                    break;
                }
            }
        }
        if (Remove) continue;
        */
        //**************************************************************************************
        
        
    /*
        fillCloseJetVars(leptonCounter);
        
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
                //std::cout << "dR is: " << _dR15[leptonCounter] << std::endl;
            } else if (dRother < 10000) {
                fillMCVars(mcother, leptonCounter);
                _ipPVmc[leptonCounter] = iTau->dxy();
                _dR15[leptonCounter] = dRother;
                //std::cout << "dR is: " << _dR15[leptonCounter] << std::endl;
            }
            
            fillIsoMCVars(leptonCounter);
            matchCloseJet(leptonCounter);

            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }

        
        leptonCounter++;
    }
    
    _nTau = leptonCounter - _nEle - _nMu;
    */

    //std::cout<<leptonCounter<<" "<<_nEle + _nMu<<std::endl;
    if (leptonCounter < 3) return;
    //if (_nEle + _nMu == 0) return;
    
    _nLeptons = leptonCounter;
    //std::cout<<"lepton counter "<<std::endl;
    
    _sb = false;
    
    //std::cout<<((TLorentzVector*)_leptonP4->At(_index1))->Pt()<<" "<<((TLorentzVector*)_leptonP4->At(_index2))->Pt()<<std::endl;
    //std::cout<<_isolation[_index1]<<" "<<_isolation[_index2]<<std::endl;
    
    _n_Jets = 0;
    _n_bJets = 0;
    HT = 0;
    //std::cout<<"jets "<<SelectedJets.size()<<std::endl;
    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){
        _jetEta[_n_Jets] = SelectedJets[i]->eta();
        _jetPhi[_n_Jets] = SelectedJets[i]->phi();
        _jetPt[_n_Jets] = SelectedJets[i]->pt();
        _jetM[_n_Jets] = SelectedJets[i]->mass();
        
        TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[_n_Jets],_jetEta[_n_Jets],_jetPhi[_n_Jets], _jetM[_n_Jets]);
        //std::cout<<_jetPt[_n_Jets]<<" "<<_jetEta[_n_Jets]<<" "<<_jetPhi[_n_Jets]<<std::endl;
        /*
        TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[_n_Jets],_jetEta[_n_Jets],_jetPhi[_n_Jets],0);
        bool clean = true;
        for (int k=0; k!=_nLeptons; ++k) {
            double dR1 = ((TLorentzVector *)_leptonP4->At(k))->DeltaR( jt );
            clean = clean && (dR1 > 0.4) ;
        }
        if (!clean) continue;
        */
        //std::cout<<"clean: "<<_jetPt[_n_Jets]<<" "<<_jetEta[_n_Jets]<<" "<<_jetPhi[_n_Jets]<<std::endl;

        for(int j=0; j != _nLeptons; ++j){
            _jetDeltaR[_n_Jets][j] = ((TLorentzVector *)_leptonP4->At(j))->DeltaR( jt ) ;
        }

        _csv[_n_Jets] = SelectedJets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        
        if(_csv[_n_Jets] > 0.814) {
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

void tril13::fillMCVars(const GenParticle* mc, const int leptonCounter) {
    
    _lPtmc[leptonCounter] = mc->pt();
    _lEmc[leptonCounter] = mc->energy();
    _lPhimc[leptonCounter] = mc->phi();
    _lEtamc[leptonCounter] = mc->eta();
    _lpdgmc[leptonCounter] = mc->pdgId();
    //std::cout << "pdg of tau truth: " << mc->pdgId() << std::endl;
    
    _origin[leptonCounter] = GPM.origin(mc);
    _originReduced[leptonCounter] = GPM.originReduced(_origin[leptonCounter]);
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

void tril13::fillCloseJetVars(const int leptonCounter) {
    _closeJetPtAll[leptonCounter] = 0;
    _closeJetAngAll[leptonCounter] = 10000;
    _ptRelAll[leptonCounter] = 0;
    
    _closeIndex[leptonCounter] = 0;
    TLorentzVector pJet;
    //std::cout << "Selected Jets Number:\t"<< SelectedJetsAll.size() << std::endl;
    for(unsigned int k = 0 ; k < SelectedJetsAll.size() ;k++ ){
        double uncPt = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Pt();
        double uncE = (SelectedJetsAll[k]->correctedP4("Uncorrected")).E();
        double corr = fMetCorrector->getJetCorrectionRawPt( uncPt, _jetEtaAll[k], myRhoJets, SelectedJetsAll[k]->jetArea(),"L1FastJet");
        double corr2 = fMetCorrector->getJetCorrectionRawPt( uncPt, _jetEtaAll[k], myRhoJets, SelectedJetsAll[k]->jetArea(),"L3Absolute"); // - for MC 
        _jetPtAll[k] = (uncPt* corr - _lPt[leptonCounter])* corr2/corr + _lPt[leptonCounter];
        _jetEAll[k] = (uncE* corr - _lE[leptonCounter])* corr2/corr + _lE[leptonCounter];

        pJet.SetPtEtaPhiM( _jetPtAll[k], _jetEtaAll[k], _jetPhiAll[k], _jetEAll[k]);
        double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );
        if (ang < _closeJetAngAll[leptonCounter]) {
            //std::cout << "This is true!" << std::endl;
            _closeJetAngAll[leptonCounter] = ang;
            _closeJetPtAll[leptonCounter] = _jetPtAll[k];
            _ptrel[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect() - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Vect());
            _ptratio[leptonCounter] = _lPt[leptonCounter] / _closeJetPtAll[leptonCounter];
            _ptRelAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect());
            
            _closeIndex[leptonCounter] = k;
        }
    }
    for (int mi=0; mi!=5; ++mi) {
        if (_miniisolation[leptonCounter][0] < multiConst[mi][0] && (_ptratio[leptonCounter] > multiConst[mi][1] || _ptrel[leptonCounter] > multiConst[mi][2]) )
            _multiisolation[leptonCounter][mi] = true;
        else
            _multiisolation[leptonCounter][mi] = false;
    }
}

void tril13::matchCloseJet(const int leptonCounter) {
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


void tril13::fillIsoMCVars(const int leptonCounter) {
    
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

void tril13::fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu) {
    
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

void tril13::fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* mu) {
    
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


DEFINE_FWK_MODULE(tril13);

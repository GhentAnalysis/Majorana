#ifndef Tools_H
#define Tools_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

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

namespace tools {
    
    struct event_info
    {
        std::vector<int> el,mu;
        
        float met;
        float ht;
        int channel;
        float weight;
        ULong64_t event;
        int  lumi;
        uint run;
        TString comment;
        event_info():  met(-1.), ht(-1.), channel(0), weight(1.), event(0), lumi(0),
        run(0), comment("") {}
        
    } ;
    
    double getActivityAroundLepton(edm::Handle<pat::PackedCandidateCollection> pfcands,
                                 const reco::Candidate* ptcl,
                                 double r_iso_min, double r_iso_max, double kt_scale,
                                 bool use_pfweight, bool charged_only,
                                 double corrTerm);
    
    double getActivityAroundLeptonDB(edm::Handle<pat::PackedCandidateCollection> pfcands,
                                 const reco::Candidate* ptcl,
                                 double r_iso_min, double r_iso_max, double kt_scale,
                                 bool use_pfweight, bool charged_only,
                                 double corrTerm);
    
    double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                                 const reco::Candidate* ptcl,
                                 double r_iso_min, double r_iso_max, double kt_scale,
                                 bool use_pfweight, bool charged_only,
                                 double corrTerm);
    
    double getPFIsolationDB(edm::Handle<pat::PackedCandidateCollection> pfcands,
                                 const reco::Candidate* ptcl,
                                 double r_iso_min, double r_iso_max, double kt_scale,
                                 bool use_pfweight, bool charged_only,
                                 double corrTerm);
    
    std::vector<const pat::Muon* > effMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                               double v_muon_pt,
                                                               reco::Vertex::Point PV,
                                                               double v_muon_d0);


    bool effMVAElectronSelectorPassed(const edm::Ptr<pat::Electron>  el,
                                                                          double v_electron_pt,
                                                                          reco::Vertex::Point PV,
                                                                          double v_electron_d0,
                                                                          bool bool_electron_chargeConsistency,
                                                                          edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                                          reco::BeamSpot::Point BS,
                                                                          bool tight);
    
    
    double pfRelIso(const pat::Muon *mu);
    double pfRelIso(const pat::Electron *el, double myRho);
    double pfAbsIso(const pat::Muon *mu);
    double pfAbsIso(const pat::Electron *el, double myRho);
    double pfRelIso(const pat::Muon *mu, double myRho);
    double pfAbsIso(const pat::Muon *mu, double myRho);

    bool triggerEmulator(const pat::Electron *iE);
    bool triggerEmulatorReturned(const pat::Electron *iE);
    bool isoTriggerEmulator(const pat::Electron *iE);
    
    float dEtaInSeed(const pat::Electron* ele);
    bool isLooseCutBasedElectronWithoutIsolation(const pat::Electron* ele);
    bool isTightCutBasedElectronWithoutIsolation(const pat::Electron* ele);
    
    bool cleanUp(const pat::Muon* testMu,
                       std::vector<const pat::Muon* > allMu,
                       const char* cutName);
    
    bool cleanUp(const pat::Electron* testMu,
                       std::vector<const pat::Electron* > allMu,
                 const char* cutName);


    std::vector<const pat::Jet* > EMJetSelector(const std::vector<pat::Jet>  & thePatJets,
                                              double  value_jet_et,
                                              double  value_jet_eta);
    
    
    std::vector<const pat::Jet* > JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                              double  value_jet_et,
                                              double  value_jet_eta);
    
    
    std::vector<const pat::Jet* > JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                              double  value_jet_et,
                                              double  value_jet_eta,
                                              std::vector<const pat::Electron*> vElectrons,
                                              std::vector<const pat::Muon*> vMuons);

    
    std::vector<const pat::Jet* > JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                              double  value_jet_et,
                                              double  value_jet_eta,
                                              std::vector<const pat::Electron*> vElectrons,
                                              std::vector<const pat::Muon*> vMuons,
                                              std::map<const reco::PFTau*, int> SelectedTaus);

    std::vector<const pat::Photon* > PhotonSelector(const std::vector<pat::Photon>  & thePatPhotons,
                                              double  value_photon_et,
                                              double  value_photon_eta,
                                              std::vector<const pat::Electron*> vElectrons,
                                              std::vector<const pat::Muon*> vMuons,
                                              std::map<const reco::PFTau*, int> SelectedTaus);
    
    
    
    std::map<const reco::PFTau*, int >  TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                 double v_tau_pt,
                                                 double v_tau_eta,
                                                 edm::Handle<reco::PFTauDiscriminator> & electron,
                                                 edm::Handle<reco::PFTauDiscriminator> & muon,
                                                 edm::Handle<reco::PFTauDiscriminator> & iso,
                                                 edm::Handle<reco::PFTauDiscriminator> & decay);

    std::map<const reco::PFTau*, int >  TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                    double v_tau_pt,
                                                    double v_tau_eta,
                                                    edm::Handle<reco::PFTauDiscriminator> & electron,
                                                    edm::Handle<reco::PFTauDiscriminator> & muon,
                                                    edm::Handle<reco::PFTauDiscriminator> & iso,
                                                    edm::Handle<reco::PFTauDiscriminator> & decay,
                                                    std::vector<const pat::Muon*> & thePatMuons,
                                                    const std::vector<const pat::Electron*>  & thePatElectrons);
    
    double MT_calc(TLorentzVector Vect, double MET, double MET_Phi);
    double Tau_dz(ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> >, TLorentzVector, reco::Vertex::Point);
    
    double Mll_calc(TLorentzVector Vect1, TLorentzVector Vect2);
    int srID(double met, double mt, double mll, double channel);
    
    double JER (double eta);
    float quadsum(float a, float b);
    float smear_pt_res(float pt, float genpt, float eta);
    double evalEt(double pt, double eta, double phi, double e);
    double evalMt(double pt, double eta, double phi, double e);
    
    
    std::vector<double> RegressionVars(const pat::Jet *jet, float genpt, const pat::Muon* mu);

    bool passMultiIsolation(TString level, double mini_iso, double jetPtRatio, double jetPtRel);
    struct effAreaForRange{
      double min;
      double max;
      double effArea;
    };


    void readEffAreas(std::string fileName, int pdgId);
    double getEffArea(int pdgId, double eta);
}

#endif

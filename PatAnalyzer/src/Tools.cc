#include "Majorana/PatAnalyzer/interface/Tools.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

double tools::getActivityAroundLepton(edm::Handle<pat::PackedCandidateCollection> pfcands,
                      const reco::Candidate* ptcl,
                      double r_iso_min, double r_iso_max, double kt_scale,
                      bool use_pfweight, bool charged_only,
                      double myRho) {
    
    if (ptcl->pt()<5.) return 99999.;

    double CorrectedTerm=0.0;

    if(ptcl->isMuon()) {

        double  Aeff[ 5 ] = { 0.0735, 0.0619, 0.0465, 0.0433, 0.0577};
        double etaMuCd = ptcl->eta();
    
        if( TMath::Abs( etaMuCd ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
        else if( TMath::Abs( etaMuCd ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
        else if( TMath::Abs( etaMuCd ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
        else if( TMath::Abs( etaMuCd ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
        else  CorrectedTerm = myRho * Aeff[ 4 ];
    }

    else if(ptcl->isElectron()){

        double  Aeff[ 7 ] = { 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687};
        double etaElCd = dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta();
        if( TMath::Abs( etaElCd ) < 1.0                                            )   CorrectedTerm = myRho * Aeff[ 0 ];
            else if( TMath::Abs( etaElCd ) > 1.0   && TMath::Abs( etaElCd ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
            else if( TMath::Abs( etaElCd ) > 1.479 && TMath::Abs( etaElCd ) < 2.0    )   CorrectedTerm = myRho * Aeff[ 2 ];
            else if( TMath::Abs( etaElCd ) > 2.0   && TMath::Abs( etaElCd ) < 2.2    )   CorrectedTerm = myRho * Aeff[ 3 ];
            else if( TMath::Abs( etaElCd ) > 2.2   && TMath::Abs( etaElCd ) < 2.3    )   CorrectedTerm = myRho * Aeff[ 4 ];
            else if( TMath::Abs( etaElCd ) > 2.3   && TMath::Abs( etaElCd ) < 2.4    )   CorrectedTerm = myRho * Aeff[ 5 ];
            else CorrectedTerm = myRho * Aeff[6];

    }
    

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
        if (fabs(dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
        //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = std::max(r_iso_min,std::min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
        if (abs(pfc.pdgId())<7) continue;
        
        double dr = deltaR(pfc, *ptcl);
        double drMax = 0.4;
        if (dr < r_iso) continue;
        if (dr > drMax) continue;
        
        //////////////////  NEUTRALS  /////////////////////////
        if (pfc.charge()==0){
            if (pfc.pt()>ptThresh) {
                double wpf(1.);
                if (use_pfweight){
                    double wpv(0.), wpu(0.);
                    for (const pat::PackedCandidate &jpfc : *pfcands) {
                        double jdr = deltaR(pfc, jpfc);
                        if (pfc.charge()!=0 || jdr<0.00001) continue;
                        double jpt = jpfc.pt();
                        if (pfc.fromPV()>1) wpv *= jpt/jdr;
                        else wpu *= jpt/jdr;
                    }
                    wpv = log(wpv);
                    wpu = log(wpu);
                    wpf = wpv/(wpv+wpu);
                }
                /////////// PHOTONS ////////////
                if (abs(pfc.pdgId())==22) {
                    if(dr < deadcone_ph) continue;
                    iso_ph += wpf*pfc.pt();
                    /////////// NEUTRAL HADRONS ////////////
                } else if (abs(pfc.pdgId())==130) {
                    if(dr < deadcone_nh) continue;
                    iso_nh += wpf*pfc.pt();
                }
            }
            //////////////////  CHARGED from PV  /////////////////////////
        } else if (pfc.fromPV()>1){
            if (abs(pfc.pdgId())==211) {
                if(dr < deadcone_ch) continue;
                iso_ch += pfc.pt();
            }
            //////////////////  CHARGED from PU  /////////////////////////
        } else {
            if (pfc.pt()>ptThresh){
                if(dr < deadcone_pu) continue;
                iso_pu += pfc.pt();
            }
        }
    }
    double iso(0.);
    if (charged_only){
        iso = iso_ch;
    } else {
        iso = iso_ph + iso_nh;
        iso -= CorrectedTerm/(0.3 * 0.3 / (0.4 * 0.4 - r_iso * r_iso));
        if (iso>0) 
            iso += iso_ch;
        else 
            iso = iso_ch;
    }
    iso = iso/ptcl->pt();
    
    return iso;
}

double tools::getActivityAroundLeptonDB(edm::Handle<pat::PackedCandidateCollection> pfcands,
                      const reco::Candidate* ptcl,
                      double r_iso_min, double r_iso_max, double kt_scale,
                      bool use_pfweight, bool charged_only,
                      double myRho) {
    
    if (ptcl->pt()<5.) return 99999.;

    /*
    double CorrectedTerm=0.0;

    if(ptcl->isMuon()) {

        double  Aeff[ 5 ] = { 0.0735, 0.0619, 0.0465, 0.0433, 0.0577};
    
        if( TMath::Abs( ptcl->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
        else if( TMath::Abs( ptcl->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
        else if( TMath::Abs( ptcl->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
        else if( TMath::Abs( ptcl->eta() ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
        else  CorrectedTerm = myRho * Aeff[ 4 ];
    }

    else if(ptcl->isElectron()){

        double  Aeff[ 7 ] = { 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687};
        if( TMath::Abs( ptcl->eta() ) < 1.0                                            )   CorrectedTerm = myRho * Aeff[ 0 ];
            else if( TMath::Abs( ptcl->eta() ) > 1.0   && TMath::Abs( ptcl->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
            else if( TMath::Abs( ptcl->eta() ) > 1.479 && TMath::Abs( ptcl->eta() ) < 2.0    )   CorrectedTerm = myRho * Aeff[ 2 ];
            else if( TMath::Abs( ptcl->eta() ) > 2.0   && TMath::Abs( ptcl->eta() ) < 2.2    )   CorrectedTerm = myRho * Aeff[ 3 ];
            else if( TMath::Abs( ptcl->eta() ) > 2.2   && TMath::Abs( ptcl->eta() ) < 2.3    )   CorrectedTerm = myRho * Aeff[ 4 ];
            else if( TMath::Abs( ptcl->eta() ) > 2.3   && TMath::Abs( ptcl->eta() ) < 2.4    )   CorrectedTerm = myRho * Aeff[ 5 ];
            else CorrectedTerm = myRho * Aeff[6];

    }
    */
    

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
        if (fabs(dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
        //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = std::max(r_iso_min,std::min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
        if (abs(pfc.pdgId())<7) continue;
        
        double dr = deltaR(pfc, *ptcl);
        double drMax = 0.4;
        if (dr < r_iso) continue;
        if (dr > drMax) continue;
        
        //////////////////  NEUTRALS  /////////////////////////
        if (pfc.charge()==0){
            if (pfc.pt()>ptThresh) {
                double wpf(1.);
                if (use_pfweight){
                    double wpv(0.), wpu(0.);
                    for (const pat::PackedCandidate &jpfc : *pfcands) {
                        double jdr = deltaR(pfc, jpfc);
                        if (pfc.charge()!=0 || jdr<0.00001) continue;
                        double jpt = jpfc.pt();
                        if (pfc.fromPV()>1) wpv *= jpt/jdr;
                        else wpu *= jpt/jdr;
                    }
                    wpv = log(wpv);
                    wpu = log(wpu);
                    wpf = wpv/(wpv+wpu);
                }
                /////////// PHOTONS ////////////
                if (abs(pfc.pdgId())==22) {
                    if(dr < deadcone_ph) continue;
                    iso_ph += wpf*pfc.pt();
                    /////////// NEUTRAL HADRONS ////////////
                } else if (abs(pfc.pdgId())==130) {
                    if(dr < deadcone_nh) continue;
                    iso_nh += wpf*pfc.pt();
                }
            }
            //////////////////  CHARGED from PV  /////////////////////////
        } else if (pfc.fromPV()>1){
            if (abs(pfc.pdgId())==211) {
                if(dr < deadcone_ch) continue;
                iso_ch += pfc.pt();
            }
            //////////////////  CHARGED from PU  /////////////////////////
        } else {
            if (pfc.pt()>ptThresh){
                if(dr < deadcone_pu) continue;
                iso_pu += pfc.pt();
            }
        }
    }
    double iso(0.);
    if (charged_only){
        iso = iso_ch;
    } else {
        iso = iso_ph + iso_nh;
        iso -= 0.5 * iso_pu;
        if (iso>0) 
            iso += iso_ch;
        else 
            iso = iso_ch;
    }
    iso = iso/ptcl->pt();
    
    return iso;
}




double tools::getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands, const reco::Candidate* ptcl,
                             double r_iso_min, double r_iso_max, double kt_scale, bool use_pfweight, bool charged_only, double myRho) {
    
    if(ptcl->pt()<5.) return 99999.;

    double CorrectedTerm=0.0;

    if(ptcl->isMuon()){

        double  Aeff[ 5 ] = { 0.0735, 0.0619, 0.0465, 0.0433, 0.0577};
        double etaMuCd = ptcl->eta();
    
        if( TMath::Abs( etaMuCd ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
        else if( TMath::Abs( etaMuCd ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
        else if( TMath::Abs( etaMuCd ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
        else if( TMath::Abs( etaMuCd ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
        else  CorrectedTerm = myRho * Aeff[ 4 ];
    }

    else if(ptcl->isElectron()){

        double  Aeff[ 7 ] = { 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687};
        double etaElCd = dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta();
        if( TMath::Abs( etaElCd ) < 1.0                                            )   CorrectedTerm = myRho * Aeff[ 0 ];
            else if( TMath::Abs( etaElCd ) > 1.0   && TMath::Abs( etaElCd ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
            else if( TMath::Abs( etaElCd ) > 1.479 && TMath::Abs( etaElCd ) < 2.0    )   CorrectedTerm = myRho * Aeff[ 2 ];
            else if( TMath::Abs( etaElCd ) > 2.0   && TMath::Abs( etaElCd ) < 2.2    )   CorrectedTerm = myRho * Aeff[ 3 ];
            else if( TMath::Abs( etaElCd ) > 2.2   && TMath::Abs( etaElCd ) < 2.3    )   CorrectedTerm = myRho * Aeff[ 4 ];
            else if( TMath::Abs( etaElCd ) > 2.3   && TMath::Abs( etaElCd ) < 2.4    )   CorrectedTerm = myRho * Aeff[ 5 ];
            else CorrectedTerm = myRho * Aeff[6];

    }
 
    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
        if (fabs(dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
        //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = std::max(r_iso_min,std::min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
        if (abs(pfc.pdgId())<7) continue;
        
        double dr = deltaR(pfc, *ptcl);
        if (dr > r_iso) continue;
        
        //////////////////  NEUTRALS  /////////////////////////
        if (pfc.charge()==0){
            if (pfc.pt()>ptThresh) {
                double wpf(1.);
                if (use_pfweight){
                    double wpv(0.), wpu(0.);
                    for (const pat::PackedCandidate &jpfc : *pfcands) {
                        double jdr = deltaR(pfc, jpfc);
                        if (pfc.charge()!=0 || jdr<0.00001) continue;
                        double jpt = jpfc.pt();
                        if (pfc.fromPV()>1) wpv *= jpt/jdr;
                        else wpu *= jpt/jdr;
                    }
                    wpv = log(wpv);
                    wpu = log(wpu);
                    wpf = wpv/(wpv+wpu);
                }
                /////////// PHOTONS ////////////
                if (abs(pfc.pdgId())==22) {
                    if(dr < deadcone_ph) continue;
                    iso_ph += wpf*pfc.pt();
                    /////////// NEUTRAL HADRONS ////////////
                } else if (abs(pfc.pdgId())==130) {
                    if(dr < deadcone_nh) continue;
                    iso_nh += wpf*pfc.pt();
                }
            }
            //////////////////  CHARGED from PV  /////////////////////////
        } else if (pfc.fromPV()>1){
            if (abs(pfc.pdgId())==211) {
                if(dr < deadcone_ch) continue;
                iso_ch += pfc.pt();
            }
            //////////////////  CHARGED from PU  /////////////////////////
        } else {
            if (pfc.pt()>ptThresh){
                if(dr < deadcone_pu) continue;
                iso_pu += pfc.pt();
            }
        }
    }
    double iso(0.);
    if (charged_only){
        iso = iso_ch;
    } else {
        iso = iso_ph + iso_nh;
        iso -= CorrectedTerm/(0.3 * 0.3 / (r_iso * r_iso));
        if (iso>0) 
            iso += iso_ch;
        else 
            iso = iso_ch;
    }
    iso = iso/ptcl->pt();
    
    return iso;
}

double tools::getPFIsolationDB(edm::Handle<pat::PackedCandidateCollection> pfcands,
                      const reco::Candidate* ptcl,
                      double r_iso_min, double r_iso_max, double kt_scale,
                      bool use_pfweight, bool charged_only,
                      double myRho) {
    
    if (ptcl->pt()<5.) return 99999.;

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
        if (fabs(dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
        //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = std::max(r_iso_min,std::min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
        if (abs(pfc.pdgId())<7) continue;
        
        double dr = deltaR(pfc, *ptcl);
        if (dr > r_iso) continue;
        
        //////////////////  NEUTRALS  /////////////////////////
        if (pfc.charge()==0){
            if (pfc.pt()>ptThresh) {
                double wpf(1.);
                if (use_pfweight){
                    double wpv(0.), wpu(0.);
                    for (const pat::PackedCandidate &jpfc : *pfcands) {
                        double jdr = deltaR(pfc, jpfc);
                        if (pfc.charge()!=0 || jdr<0.00001) continue;
                        double jpt = jpfc.pt();
                        if (pfc.fromPV()>1) wpv *= jpt/jdr;
                        else wpu *= jpt/jdr;
                    }
                    wpv = log(wpv);
                    wpu = log(wpu);
                    wpf = wpv/(wpv+wpu);
                }
                /////////// PHOTONS ////////////
                if (abs(pfc.pdgId())==22) {
                    if(dr < deadcone_ph) continue;
                    iso_ph += wpf*pfc.pt();
                    /////////// NEUTRAL HADRONS ////////////
                } else if (abs(pfc.pdgId())==130) {
                    if(dr < deadcone_nh) continue;
                    iso_nh += wpf*pfc.pt();
                }
            }
            //////////////////  CHARGED from PV  /////////////////////////
        } else if (pfc.fromPV()>1){
            if (abs(pfc.pdgId())==211) {
                if(dr < deadcone_ch) continue;
                iso_ch += pfc.pt();
            }
            //////////////////  CHARGED from PU  /////////////////////////
        } else {
            if (pfc.pt()>ptThresh){
                if(dr < deadcone_pu) continue;
                iso_pu += pfc.pt();
            }
        }
    }
    double iso(0.);
    if (charged_only){
        iso = iso_ch;
    } else {
        iso = iso_ph + iso_nh;
        iso -= 0.5 * iso_pu;
        if (iso>0) 
            iso += iso_ch;
        else 
            iso = iso_ch;
    }
    iso = iso/ptcl->pt();
    
    return iso;
}

std::vector<const pat::Muon* > tools::effMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                            double v_muon_pt,
                                                            reco::Vertex::Point PV,
                                                            double v_muon_d0)
{
    double v_muon_eta = 2.4;
    //double v_muon_dz = 0.1; // 0.1 for sync
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
        //std::cout << mu->pt() << "\t" << TMath::Abs( mu->eta() ) << "\t" << mu->isGlobalMuon() << "\t" << mu->isPFMuon() << "\t" << mu->isTrackerMuon() << std::endl;
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        
        /*
        if ( !(mu->isGlobalMuon()  || mu->isTrackerMuon()) ) continue;
        if ( !mu->isPFMuon() ) continue;
        */
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        //std::cout << innerTrack.isNull() << "\t" << TMath::Abs(innerTrack->dxy(PV)) << "\t" << TMath::Abs(innerTrack->dz(PV)) << std::endl;
        if( innerTrack.isNull() ) continue;
        
        //if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        //if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}




//Muon pfRelIso
double tools::pfRelIso(const pat::Muon *mu)
{
	double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR03().sumPhotonEt;
    double beta = mu->pfIsolationR03().sumPUPt;
    double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt();
    return pfRelIsoMu;
}
// Muon absolute isolation
double tools::pfAbsIso(const pat::Muon *mu)
{
 	//******************** absolute isolation 
	double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR03().sumPhotonEt;
    double beta = mu->pfIsolationR03().sumPUPt;
    double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) );
    return pfRelIsoMu;
}


std::map<int, std::vector<tools::effAreaForRange>*> effAreas; 


void tools::readEffAreas(std::string fileName, int pdgId){
  effAreas[pdgId] = new std::vector<tools::effAreaForRange>();
  std::ifstream file(fileName);
  std::string line;
  while(std::getline(file, line)){
    if(line.find('#') == 0) continue;
    std::stringstream linestream(line);
    tools::effAreaForRange ea;
    linestream >> ea.min >> ea.max >> ea.effArea;
    effAreas[pdgId]->push_back(ea);
  } 
}

double tools::getEffArea(int pdgId, double eta){
  for(auto ea : *effAreas[pdgId]){
    if(ea.min < std::abs(eta) and ea.max > std::abs(eta)) return ea.effArea;
  }
  return 0;
}


double tools::pfRelIso(const pat::Electron *iE, double myRho){
    return tools::pfAbsIso(iE, myRho)/iE->pt();
}

double tools::pfAbsIso(const pat::Electron *iE, double myRho){
    if(!effAreas[11]) readEffAreas(edm::FileInPath("Majorana/PatAnalyzer/src/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt").fullPath(), 11);

    double CorrectedTerm = myRho*getEffArea(11, iE->superCluster()->eta());
    return (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - CorrectedTerm ) );
}


double tools::pfRelIso(const pat::Muon *mu, double myRho)
{
    //double  Aeff[ 5 ] = { 0.0913, 0.0765, 0.0546, 0.0728, 0.1177};
    double  Aeff[ 5 ] = { 0.0735, 0.0619, 0.0465, 0.0433, 0.0577};
    double CorrectedTerm=0.0;
    
    if( TMath::Abs( mu->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( mu->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( mu->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( mu->eta() ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
    else  CorrectedTerm = myRho * Aeff[ 4 ];
    
    double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR03().sumPhotonEt;
    //double beta = mu->pfIsolationR03().sumPUPt;
    
    double pfRelIsoE = (chargedHadronIso + TMath::Max(0.0, neutralHadronIso + photonIso - CorrectedTerm ) ) /mu->pt() ;
    
    return pfRelIsoE;
}



bool tools::triggerEmulator(const pat::Electron *iE) {

    bool passed = true;
    //if( TMath::Abs(1.0/iE->ecalEnergy() - iE->eSuperClusterOverP()/iE->ecalEnergy()) > 0.01 ) passed = false;
    if( TMath::Abs(1.0/iE->correctedEcalEnergy() - iE->eSuperClusterOverP()/iE->correctedEcalEnergy()) > 0.01 ) passed = false;
                
    else if( TMath::Abs(iE->superCluster()->eta()) < 1.479  ) {
         if( TMath::Abs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.04  ) passed = false;
         else if( TMath::Abs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) passed = false;
         //else if( TMath::Abs(iE->full5x5_sigmaIetaIeta()) > 0.011 ) passed = false;
         else if( iE->full5x5_sigmaIetaIeta() > 0.011 ) passed = false;
         //else if( TMath::Abs(iE->hadronicOverEm())  > 0.08  ) passed = false;
         //else if( TMath::Abs(iE->hcalOverEcal())  > 0.08  ) passed = false;
         else if( iE->hcalOverEcal()  > 0.08  ) passed = false;
    }

    else {
         if( TMath::Abs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.08 ) passed = false;
         else if( TMath::Abs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) passed = false;
         //else if( TMath::Abs(iE->full5x5_sigmaIetaIeta()) > 0.031 ) passed = false;
         else if( iE->full5x5_sigmaIetaIeta() > 0.031 ) passed = false;
         //else if( TMath::Abs(iE->hadronicOverEm()) > 0.08 ) passed = false;
         //else if( TMath::Abs(iE->hcalOverEcal())  > 0.08  ) passed = false;
         else if( iE->hcalOverEcal()  > 0.08  ) passed = false;
    }
                        
    return passed;
}

bool tools::triggerEmulatorReturned(const pat::Electron *iE) {

    bool passed = true;
    double etasc = iE->superCluster()->eta();
    
    if (iE->hadronicOverEm()>=(0.10-0.03*(fabs(etasc)>1.479))) passed = false;
    //if (iE->hcalOverEcal()>=(0.10-0.03*(fabs(etasc)>1.479))) passed = false;
    if (fabs(iE->deltaEtaSuperClusterTrackAtVtx())>=(0.01-0.002*(fabs(etasc)>1.479))) passed = false;
    if (fabs(iE->deltaPhiSuperClusterTrackAtVtx())>=(0.04+0.03*(fabs(etasc)>1.479))) passed = false;
    double eInvMinusPInv = iE->ecalEnergy()>0. ? (1.0/iE->ecalEnergy() - iE->eSuperClusterOverP()/iE->ecalEnergy()) : 9e9;
    if (eInvMinusPInv<=-0.05) passed = false;
    if (eInvMinusPInv>=(0.01-0.005*(fabs(etasc)>1.479))) passed = false;
    if (iE->full5x5_sigmaIetaIeta()>=(0.011+0.019*(fabs(etasc)>1.479))) passed = false;
    
    return passed;
}


bool tools::isoTriggerEmulator(const pat::Electron *iE) {

    bool passed = true;
            
    //if (iE->pfIsolationVariables().sumChargedHadronPt/iE->pt() > 0.2) passed = false;
    if (iE->dr03TkSumPt()/iE->pt() > 0.2) passed = false;
    else if (iE->hcalPFClusterIso()/iE->pt() > 0.25) passed = false;
    else if (iE->ecalPFClusterIso()/iE->pt() > 0.45) passed = false;
                        
     return passed;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float tools::dEtaInSeed(const pat::Electron* ele){
  if(ele->superCluster().isNonnull() and ele->superCluster()->seed().isNonnull()) return ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta();
  else                                                                            return std::numeric_limits<float>::max();
}

bool tools::isLooseCutBasedElectronWithoutIsolation(const pat::Electron* ele){
    if(not (ele->isEB() or ele->isEE())) return false;

    float eInvMinusPInv = std::abs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();

    if(ele->full5x5_sigmaIetaIeta()               >= (ele->isEB() ? 0.11    : 0.0314 ))       return false;
    if(abs(dEtaInSeed(ele))                       >= (ele->isEB() ? 0.00477 : 0.00868))       return false;
    if(abs(ele->deltaPhiSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.222   : 0.213  ))       return false;
    if(ele->hadronicOverEm()                      >= (ele->isEB() ? 0.298   : 0.101  ))       return false;
    if(eInvMinusPInv                              >= (ele->isEB() ? 0.241   : 0.14   ))       return false;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1) return false;
    if(not ele->passConversionVeto())                                                         return false;
    return true;
}

// It has been checked that (this + isolation cut) == cutBasedElectronID-Summer16-80X-V1-tight 
bool tools::isTightCutBasedElectronWithoutIsolation(const pat::Electron* ele){
    if(not (ele->isEB() or ele->isEE())) return false;

    float eInvMinusPInv = std::abs(1.0 - ele->eSuperClusterOverP())/ele->ecalEnergy();
    if(ele->full5x5_sigmaIetaIeta()               >= (ele->isEB() ? 0.00998 : 0.0292))        return false;
    if(abs(dEtaInSeed(ele))                       >= (ele->isEB() ? 0.00308 : 0.00605))       return false;
    if(abs(ele->deltaPhiSuperClusterTrackAtVtx()) >= (ele->isEB() ? 0.0816  : 0.0394))        return false;
    if(ele->hadronicOverEm()                      >= (ele->isEB() ? 0.0414  : 0.0641))        return false;
    if(eInvMinusPInv                              >= (ele->isEB() ? 0.0129  : 0.0129))        return false;
    if(ele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >  1) return false;
    if(not ele->passConversionVeto())                                                         return false;
    return true;
}

//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tools::Tau_dz(ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > vtx, TLorentzVector p4, reco::Vertex::Point PV){

    /*
    std::cout << "Check Tau_dz function" << std::endl;
    std::cout << vtx.X() << " " << vtx.Y() << " "<< vtx.Z() << std::endl;
    std::cout << p4.Px() << " " << p4.Py() << " "<< p4.Pz()  << " "<< p4.Pt() << std::endl;
    std::cout << PV.x() << " " << PV.y() << " "<< PV.z() << std::endl;
    std::cout << "Check Tau_dz function Happy ending" << std::endl;
    */

    double dz = (vtx.Z() - PV.z()) - ((vtx.X() - PV.x())*p4.Px()+(vtx.Y()-PV.y())*p4.Py())/ p4.Pt() *  p4.Pz()/ p4.Pt();
    return dz;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool bool_jet_id= true;
    

    std::vector< const pat::Jet* > vJets;

    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        //std::cout << "Jet Pt: " << jet->pt() << std::endl;
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;

        double eta = jet->eta();

        double NHF = jet->neutralHadronEnergyFraction();
        double NEMF = jet->neutralEmEnergyFraction();
        double CHF = jet->chargedHadronEnergyFraction();
        //double MUF = jet->muonEnergyFraction();
        double CEMF = jet->chargedEmEnergyFraction();
        double NumConst = jet->chargedMultiplicity()+jet->neutralMultiplicity();
        double CHM = jet->chargedMultiplicity(); 
        double NumNeutralParticles = jet->neutralMultiplicity(); 

        bool looseJetID;

        if( bool_jet_id )
	    {

            //std::cout << "Jets ingridients: " << jet->pt() << " " << jet->eta() << " " << NHF << " " << NEMF << " " << CHF << " " << MUF << " " << CEMF << " " << NumConst << " " << CHM << " " << NumNeutralParticles << std::endl;
            
            /*
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() + jet->photonMultiplicity() ) < 2 ) continue;
            if( ( jet->chargedMultiplicity() + jet->neutralMultiplicity() ) < 2 ) continue;
            if(  jet->muonEnergyFraction() >=0.8 ) continue;


            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
            */
            if(abs(eta) < 3.0)
                looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4); 
            else
                looseJetID = (NEMF<0.90 && NumNeutralParticles>10);

	    }

        //std::cout << "Jet Pt: " << jet->pt() << " Jet Eta: " << jet->eta() << std::endl;  
        if (looseJetID)
            vJets.push_back( &*jet );
        
        /*unsigned int nConst = jet->getPFConstituents().size();
        std::cout<<"Number of constituents "<<nConst<<std::endl;
        for (unsigned int i=0; i!=nConst; ++i) {
            std::cout<<jet->getPFConstituent(i)->reco::LeafCandidate::vz()<<std::endl;
        }*/
        
        
    }
    return vJets;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    //double  v_jet_leptonVetoDR = -1.;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
	    
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons,
                                                 std::map<const reco::PFTau*, int> SelectedTaus)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }

        for(std::map<const reco::PFTau*, int >::iterator it = SelectedTaus.begin() ; it != SelectedTaus.end() ;it++ ){
            
            const reco::PFTau *itau = it->first;

            float dphi = TMath::ACos( TMath::Cos( itau->phi()-jet->phi() ) );
            float deta = itau->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}


bool tools::passMultiIsolation(TString level, double mini_iso, double jetPtRatio, double jetPtRel){
    if(level == "VL") return mini_iso < 0.25 && (jetPtRatio > 0.67 || jetPtRel > 4.4);
    if(level == "L")  return mini_iso < 0.20 && (jetPtRatio > 0.69 || jetPtRel > 6.0);
    if(level == "M")  return mini_iso < 0.16 && (jetPtRatio > 0.76 || jetPtRel > 7.2);
    if(level == "T")  return mini_iso < 0.12 && (jetPtRatio > 0.80 || jetPtRel > 7.2);
    if(level == "VT") return mini_iso < 0.09 && (jetPtRatio > 0.84 || jetPtRel > 7.2);
    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Photon* > tools::PhotonSelector(const std::vector<pat::Photon>  & thePatPhotons,
                                                 double  v_photon_pt,
                                                 double  v_photon_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons,
                                                 std::map<const reco::PFTau*, int> SelectedTaus)
{
    bool    bool_photon_id= true;
    double  v_photon_leptonVetoDR=0.2;
    
    std::vector< const pat::Photon* > vJets;
    
    for( std::vector<pat::Photon>::const_iterator jet = thePatPhotons.begin(); jet != thePatPhotons.end(); jet++ )
	{
        //std::cout<<jet->pt()<<" "<<jet->eta()<<" "<<jet->hadTowOverEm()<<" "<<jet->sigmaIetaIeta()<<std::endl;
        if( jet->pt() < v_photon_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_photon_eta) continue;
        if( bool_photon_id )
	    {
            if( jet->hadTowOverEm() >= 0.05 ) continue;
            if( TMath::Abs(jet->sigmaIetaIeta()) > 0.012 ) continue;
            if( TMath::Abs( jet->eta() ) < 1.479 )
            {
                if( jet->hadTowOverEm() >= 0.05 ) continue;
                if( TMath::Abs(jet->sigmaIetaIeta()) > 0.034 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_photon_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_photon_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        for(std::map<const reco::PFTau*, int >::iterator it = SelectedTaus.begin() ; it != SelectedTaus.end() ;it++ ){
            
            const reco::PFTau *itau = it->first;
            
            float dphi = TMath::ACos( TMath::Cos( itau->phi()-jet->phi() ) );
            float deta = itau->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_photon_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}


double tools::MT_calc(TLorentzVector Vect, double MET, double MET_Phi){
    
    double MT=sqrt(2* Vect.Pt() * MET * ( 1 - (TMath::Cos(Vect.Phi() - MET_Phi )) ) );
    
    return MT;
}

double tools::Mll_calc(TLorentzVector Vect1, TLorentzVector Vect2){
    return (Vect1 + Vect2).Mag();
}

//ID=#MET+5*(#MT+3*(#Mll+3*#category))
int tools::srID(double met, double mt, double mll, double channel) {
    if ((channel > 0) && (channel!=4))
        return TMath::Min(int(met/50),4) + 5*(int(mt>120)+int(mt>160) + 3*(2*int(mll>100) + 3*channel));
    else 
        return TMath::Min(int(met/50),4) + 5*(int(mt>120)+int(mt>160) + 3*(int(mll>75) + int(mll>105) + 3*channel));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//*****************************************************************************************************************************
//**** Tau Selector ***********************************************************************************************************
//*****************************************************************************************************************************
std::map<const reco::PFTau*, int > tools::TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                    double v_tau_pt,
                                                    double v_tau_eta,
                                                    edm::Handle<reco::PFTauDiscriminator> & electron,
                                                    edm::Handle<reco::PFTauDiscriminator> & muon,
                                                    edm::Handle<reco::PFTauDiscriminator> & iso,
                                                    edm::Handle<reco::PFTauDiscriminator> & decay){
    std::map<const reco::PFTau*, int> vTaus;
    
    for( unsigned  i=0; i<PFTaus->size(); i++ ) {
        
        //std::cout<<"Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
        //std::cout<<"pt threshold"<<std::endl;
        
        if((*PFTaus)[i].pt()<v_tau_pt) continue;
        
        //std::cout<<"eta cut"<<std::endl;
        
        if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
        
        reco::PFTauRef tauCandidate(PFTaus, i);
        //std::cout<<"electron discr."<<std::endl;
        if( (*electron)[tauCandidate] < 0.5 ) continue;
        
        //std::cout<<"electron discr."<<std::endl;
        if( (*muon)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"muon discr."<<std::endl;
        if( (*iso)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"iso discr."<<std::endl;
        if( (*decay)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"tau pass"<<std::endl;
        vTaus.insert(std::pair<const reco::PFTau*, int >(&((*PFTaus)[i]), i));
    }
    return vTaus;
}


std::map<const reco::PFTau*, int > tools::TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                      double v_tau_pt,
                                                      double v_tau_eta,
                                                      edm::Handle<reco::PFTauDiscriminator> & electron,
                                                      edm::Handle<reco::PFTauDiscriminator> & muon,
                                                      edm::Handle<reco::PFTauDiscriminator> & iso,
                                                      edm::Handle<reco::PFTauDiscriminator> & decay,
                                                      std::vector<const pat::Muon*> & thePatMuons,
                                                      const std::vector<const pat::Electron*>  & thePatElectrons){
    std::map<const reco::PFTau*, int> vTaus;
    
    for( unsigned  i=0; i<PFTaus->size(); i++ ) {
        
        //std::cout<<"Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
        //std::cout<<"pt threshold"<<std::endl;
        
        if((*PFTaus)[i].pt()<v_tau_pt) continue;
        
        //std::cout<<"eta cut"<<std::endl;
        
        if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
        
        reco::PFTauRef tauCandidate(PFTaus, i);
        //std::cout<<"electron discr."<<std::endl;
        if( (*electron)[tauCandidate] < 0.5 ) continue;
        
        //std::cout<<"electron discr."<<std::endl;
        if( (*muon)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"muon discr."<<std::endl;
        if( (*iso)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"iso discr."<<std::endl;
        if( (*decay)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"tau pass"<<std::endl;
        
        bool Remove = false;
        TLorentzVector Tau; Tau.SetPtEtaPhiE( (*PFTaus)[i].pt(), (*PFTaus)[i].eta(), (*PFTaus)[i].phi(), (*PFTaus)[i].energy() );
        for( std::vector<const pat::Muon *>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
            TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
            float dR = Mu.DeltaR( Tau );
            if( dR < 0.1 ) {
                Remove = true;
                break;
            }
        }
        if (!Remove) {
            for( std::vector<const pat::Electron *>::const_iterator mu = thePatElectrons.begin() ; mu != thePatElectrons.end() ; mu++ ) {
                TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
                float dR = Mu.DeltaR( Tau );
                if( dR < 0.1 ) {
                    Remove = true;
                    break;
                }
            }
        }
        if (!Remove)
            vTaus.insert(std::pair<const reco::PFTau*, int >(&((*PFTaus)[i]), i));
    }
    return vTaus;
}


bool tools::cleanUp(const pat::Muon* testMu,
                   std::vector<const pat::Muon* > allMu,
                   const char* cutName)
{
    
    TLorentzVector TLVec; TLVec.SetPtEtaPhiE( testMu->pt(), testMu->eta(), testMu->phi(), testMu->energy() );
    TString CutName = cutName;

    for (unsigned int i=0; i!=allMu.size(); ++i) {
        if (testMu->charge() + allMu[i]->charge() != 0) continue;
    
        TLorentzVector P1; P1.SetPtEtaPhiE( allMu[i]->pt(), allMu[i]->eta(), allMu[i]->phi(), allMu[i]->energy() );
        
        float Mass = ( P1 + TLVec ).M();
        if( CutName == "Z" ){
            if( Mass < 106. && Mass > 76. )  {
                return true;
            }
        }
        else if( CutName == "gammastar" ){
            if( Mass < 12. )  {
                return true;
            }
        }
    }
    
    return false;
}

bool tools::cleanUp(const pat::Electron* testMu,
                   std::vector<const pat::Electron* > allMu,
                   const char* cutName)
{
    
    TLorentzVector TLVec; TLVec.SetPtEtaPhiE( testMu->pt(), testMu->eta(), testMu->phi(), testMu->energy() );
    TString CutName = cutName;

    for (unsigned int i=0; i!=allMu.size(); ++i) {
        if (testMu->charge() + allMu[i]->charge() != 0) continue;
        
        TLorentzVector P1; P1.SetPtEtaPhiE( allMu[i]->pt(), allMu[i]->eta(), allMu[i]->phi(), allMu[i]->energy() );
        
        float Mass = ( P1 + TLVec ).M();
        if( CutName == "Z" ){
            if( Mass < 106. && Mass > 76. )  {
                return true;
            }
        }
        else if( CutName == "gammastar" ){
            if( Mass < 12. )  {
                return true;
            }
        }
    }
    
    return false;
}



double tools::JER (double eta) {

    double feta = fabs(eta);
    double scf = 1;
    if (feta < 0.5)
        scf = 1.052;
    else if (feta < 1.1)
        scf = 1.057;
    else if (feta < 1.7)
        scf = 1.096;
    else if (feta < 2.3)
        scf = 1.134;
    else scf = 1.288;
    
    return scf;
}

float tools::quadsum(float a, float b) {
    return sqrt(a*a+b*b);
}

//______________________________________________________________________________
#define JERSUMMER11
float tools::smear_pt_res(float pt, float genpt, float eta)
{
    eta = fabs(eta);
    if (genpt>15. && (fabs(pt - genpt) / pt)<0.5) {  // limit the effect to the core
        double res    = 1.0;
        //double resErr = 0.0;
#ifdef JERSUMMER11
        // from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        if (eta <= 0.5) {
            res    = 1.052;
            //resErr = quadsum(0.012, 0.062);
        } else if (0.5 < eta && eta <= 1.1) {
            res    = 1.057;
            //resErr = quadsum(0.012, 0.062);
        } else if (1.1 < eta && eta <= 1.7) {
            res    = 1.096;
            //resErr = quadsum(0.017, 0.063);
        } else if (1.7 < eta && eta <= 2.3) {
            res    = 1.134;
            //resErr = quadsum(0.035, 0.087);
        } else {
            res    = 1.288;
            //resErr = quadsum(0.127, 0.155);
        }
#else
        // from VHbb analysis
        if (eta <= 1.1) {
            res    = 1.05;
            //resErr = 0.05;
        } else if (1.1 < eta && eta <= 2.5) {
            res    = 1.10;
            //resErr = 0.10;
        } else {
            res    = 1.30;
            //resErr = 0.20;
        }
#endif
        float deltapt = (pt - genpt) * res;
        return TMath::Max(float(0.), genpt + deltapt);
    }
    return pt;
}

//______________________________________________________________________________
#include "TLorentzVector.h"
double tools::evalEt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Et();
}

double tools::evalMt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Mt();
}

std::vector<double> tools::RegressionVars(const pat::Jet *jet, float genpt, const pat::Muon* mu) {
    std::vector<double> vars;
    
    vars.push_back(smear_pt_res((jet->correctedP4("Uncorrected")).Pt(), genpt, jet->eta()));
    vars.push_back(jet->pt());
    vars.push_back(evalEt(jet->pt(), jet->eta(), jet->phi(), jet->energy()));
    vars.push_back(evalMt(jet->pt(), jet->eta(), jet->phi(), jet->energy()));
    
    double hJet_ptLeadTrack = 0;
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
    }
    vars.push_back(hJet_ptLeadTrack);
    
    double hJet_vtx3dL = 0;
    double hJet_vtx3deL = 0;
    double hJet_vtxMass = 0;
    double hJet_vtxPt = 0;
    
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
    } //else
      //  std::cout<<"no info"<<std::endl;
    
    vars.push_back(TMath::Max(0.,hJet_vtx3dL));
    vars.push_back(TMath::Max(0.,hJet_vtx3deL));
    vars.push_back(TMath::Max(0.,hJet_vtxMass));
    vars.push_back(TMath::Max(0.,hJet_vtxPt));
    
    vars.push_back(jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction());
    //vars.push_back(jet->getPFConstituents().size());
    vars.push_back(jet->numberOfDaughters());
    
    vars.push_back(0.); //JECUnc in the main file
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    vars.push_back(pLep.Perp(pJet.Vect()));
    vars.push_back(mu->pt());
    vars.push_back(pLep.DeltaR(pJet));
    
    /*values[0] = "breg_rawptJER := smear_pt_res(hJet_ptRaw, hJet_genPt, hJet_eta)";
    values[1] = "breg_pt := hJet_pt";
    values[2] = "breg_et := evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[3] = "breg_mt := evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[4] = "breg_leadtrackpt := hJet_ptLeadTrack";
    values[5] = "breg_vtx3dL := max(0,hJet_vtx3dL)";
    values[6] = "breg_vtx3deL := max(0,hJet_vtx3deL)";
    values[7] = "breg_vtxMass := max(0,hJet_vtxMass)";
    values[8] = "breg_vtxPt := max(0,hJet_vtxPt)";
    values[9] = "breg_cef := hJet_cef";
    values[10] = "breg_ntot := hJet_nconstituents";
    values[11] = "breg_eJEC := hJet_JECUnc";
    values[12] = "breg_softlepptrel := max(0,hJet_SoftLeptptRel*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[13] = "breg_softleppt := max(0,hJet_SoftLeptPt*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[14] = "breg_softlepdR := max(0,hJet_SoftLeptdR*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[15] = "breg_evt_rho25 := rho25";*/
    
    return vars;
    
}


#include "Majorana/PatAnalyzer/interface/Tools.h"

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
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
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
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
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

double tools::getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
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
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
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
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
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
    // change to R04 for Deniz analysis
    double chargedHadronIso = mu->pfIsolationR04().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR04().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR04().sumPhotonEt;
    double beta = mu->pfIsolationR04().sumPUPt;
    double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
    return pfRelIsoMu;
}

double tools::pfRelIso(const pat::Electron *iE, double myRho)
{
    //double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13  };
    //double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14  };
    //double  Aeff[ 5 ] = { 0.1013, 0.0988, 0.0572, 0.0842, 0.1530};
    double  Aeff[ 7 ] = { 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687};

    double CorrectedTerm=0.0;
    if( TMath::Abs( iE->superCluster()->eta() ) < 1.0                                            )   CorrectedTerm = myRho * Aeff[ 0 ];
     else if( TMath::Abs( iE->superCluster()->eta() ) > 1.0   && TMath::Abs( iE->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
     else if( TMath::Abs( iE->superCluster()->eta() ) > 1.479 && TMath::Abs( iE->superCluster()->eta() ) < 2.0    )   CorrectedTerm = myRho * Aeff[ 2 ];
     else if( TMath::Abs( iE->superCluster()->eta() ) > 2.0   && TMath::Abs( iE->superCluster()->eta() ) < 2.2    )   CorrectedTerm = myRho * Aeff[ 3 ];
     else if( TMath::Abs( iE->superCluster()->eta() ) > 2.2   && TMath::Abs( iE->superCluster()->eta() ) < 2.3    )   CorrectedTerm = myRho * Aeff[ 4 ];
     else if( TMath::Abs( iE->superCluster()->eta() ) > 2.3   && TMath::Abs( iE->superCluster()->eta() ) < 2.4    )   CorrectedTerm = myRho * Aeff[ 5 ];
     else  CorrectedTerm = myRho * Aeff[ 6 ];
    
    /*if( TMath::Abs( iE->superCluster()->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
    else  CorrectedTerm = myRho * Aeff[ 4 ];*/

    /*
    if( TMath::Abs( iE->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( iE->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( iE->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( iE->eta() ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
    else  CorrectedTerm = myRho * Aeff[ 4 ];
    */

    //double pfRelIsoE = (iE->chargedHadronIso() + TMath::Max(0.0, iE->neutralHadronIso() + iE->photonIso() - CorrectedTerm ) ) /iE->pt() ;
     //std::cout << iE->pfIsolationVariables().sumChargedHadronPt << " " << iE->pfIsolationVariables().sumNeutralHadronEt << " " <<  iE->pfIsolationVariables().sumPhotonEt << " " << CorrectedTerm << " " << iE->pt() << std::endl ;

    // Iso 0.3, 22 Jan 2016
    cout << "Info about iso: " << iE->pfIsolationVariables().sumChargedHadronPt << " " << iE->pfIsolationVariables().sumNeutralHadronEt << " " << iE->pfIsolationVariables().sumPhotonEt << " " << CorrectedTerm << endl; 
    double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - CorrectedTerm ) ) /iE->pt() ;

    return pfRelIsoE;

     /* Iso 0.4, 22 Jan 2016
    double chargedHadronIso = iE->chargedHadronIso();
    double neutralHadronIso = iE->neutralHadronIso();
    double photonIso = iE->photonIso();
    //double beta = mu->pfIsolationR03().sumPUPt;
    
    double pfRelIsoE = (chargedHadronIso + TMath::Max(0.0, neutralHadronIso + photonIso - CorrectedTerm/(0.3 * 0.3 / (0.4 * 0.4)) ) ) / iE->pt() ;
    
    return pfRelIsoE;
    */
}

/*
double tools::pfRelIso(const edm::Ptr<reco::GsfElectron> iE, double myRho)
{
    //double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13  };
    //double  Aeff[ 5 ] = { 0.1013, 0.0988, 0.0572, 0.0842, 0.1530};

    double  Aeff[ 7 ] = { 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687};
    double CorrectedTerm=0.0;
    if( TMath::Abs( iE->eta() ) < 1.0                                            )   CorrectedTerm = myRho * Aeff[ 0 ];
     else if( TMath::Abs( iE->eta() ) > 1.0   && TMath::Abs( iE->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
     else if( TMath::Abs( iE->eta() ) > 1.479 && TMath::Abs( iE->eta() ) < 2.0    )   CorrectedTerm = myRho * Aeff[ 2 ];
     else if( TMath::Abs( iE->eta() ) > 2.0   && TMath::Abs( iE->eta() ) < 2.2    )   CorrectedTerm = myRho * Aeff[ 3 ];
     else if( TMath::Abs( iE->eta() ) > 2.2   && TMath::Abs( iE->eta() ) < 2.3    )   CorrectedTerm = myRho * Aeff[ 4 ];
     else if( TMath::Abs( iE->eta() ) > 2.3   && TMath::Abs( iE->eta() ) < 2.4    )   CorrectedTerm = myRho * Aeff[ 5 ];
    else CorrectedTerm = myRho * Aeff[6];
    
    if( TMath::Abs( iE->superCluster()->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
    else  CorrectedTerm = myRho * Aeff[ 4 ];


    if( TMath::Abs( iE->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( iE->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( iE->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( iE->eta() ) < 2.2  )   CorrectedTerm = myRho * Aeff[ 3 ];
    else  CorrectedTerm = myRho * Aeff[ 4 ];

    //double pfRelIsoE = (iE->chargedHadronIso() + TMath::Max(0.0, iE->neutralHadronIso() + iE->photonIso() - CorrectedTerm ) ) /iE->pt() ;
    double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - CorrectedTerm ) ) /iE->pt() ;

    return pfRelIsoE;
}
*/

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

//Muon Selector

std::vector<const pat::Muon* > tools::fakeMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    //double v_muon_dz = 0.2;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        //if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

std::vector<const pat::Muon* > tools::ewkMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    //double v_muon_dz = 0.2; <2012
    double v_muon_dz = 0.5;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
        
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;

        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<const pat::Electron* > tools::ewkElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    double v_electron_eta=2.4;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        } 
     /*   std::cout<<TMath::Abs(el->eta())<<std::endl;
        std::cout<<TMath::Abs(el->superCluster()->eta())<<std::endl;
        std::cout<<el->pt()<<std::endl;
        std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
        std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;
        std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;
        std::cout<<"Conversion  "<<ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS)<<std::endl;
        std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
        
        std::cout<<"********"<<std::endl;*/
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        //if (!el->trackerDrivenSeed() ) continue;

        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;

        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
    
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
                
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;

        if( TMath::Abs( el->eta()) < 1.5  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else if( TMath::Abs( el->eta() ) < 2.4 )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<const pat::Electron* > tools::fakeElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.1;
    double v_electron_eta=2.5;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        /*   std::cout<<TMath::Abs(el->eta())<<std::endl;
         std::cout<<TMath::Abs(el->superCluster()->eta())<<std::endl;
         std::cout<<el->pt()<<std::endl;
         std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
         std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;
         std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;
         std::cout<<"Conversion  "<<ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS)<<std::endl;
         std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
         
         std::cout<<"********"<<std::endl;*/
        
        if( el->pt() < v_electron_pt ) continue;
        
        //std::cout<<TMath::Abs(el->eta())<<std::endl;
        
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;

        //std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
        //std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;

        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;

        //std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;

        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        //std::cout<<vtxFitConversion<<std::endl;
        if( vtxFitConversion )  continue;
        //std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
        
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}


std::vector<const pat::Electron* > tools::phys14LooseElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                               double v_electron_pt,
                                                               reco::Vertex::Point PV,
                                                               double v_electron_d0,
                                                               bool bool_electron_chargeConsistency,
                                                               edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                               reco::BeamSpot::Point BS)
{
    //double v_electron_dz = 0.54342;
    double v_electron_eta=2.5;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        if( el->pt() < v_electron_pt ) continue;

        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;

        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.072624  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.012442 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.010557 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.121476  ) continue;  //recommended is 0.12 but HLT applies 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.221803 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.173670  ) continue;
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.145129 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.010654 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.032602 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.131862 ) continue;   /// at the HLT 0.075  recommended is 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.142283 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.198444  ) continue;
        }
        
        vElectrons.push_back(&*el );
    }
    return vElectrons;
}

std::vector<const pat::Electron* > tools::csa14MediumElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                                     double v_electron_pt,
                                                                     reco::Vertex::Point PV,
                                                                     double v_electron_d0,
                                                                     bool bool_electron_chargeConsistency,
                                                                     edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                                     reco::BeamSpot::Point BS)
{
    //double v_electron_dz = 0.54342;
    double v_electron_eta=2.5;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        if( el->pt() < v_electron_pt ) continue;
        
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        
        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.0323  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.0106 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.0107 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.067  ) continue;  //recommended is 0.12 but HLT applies 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.1043 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.22310  ) continue;
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.0455 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.0108 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.0318 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.097 ) continue;   /// at the HLT 0.075  recommended is 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.1201 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.7523  ) continue;
        }
        
        vElectrons.push_back(&*el );
    }
    return vElectrons;
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

std::vector<const pat::Jet* > tools::JetSelectorAll(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool    bool_jet_id = false;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{

        
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() + jet->photonMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
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
                looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4); 
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


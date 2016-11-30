#include "Majorana/PatAnalyzer/interface/Tools.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"


/*
 * Stuff for reading in and retrieving the effective areas
 */
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


/*
 * mini isolation
 */
double tools::getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                        const reco::Candidate* ptcl,
                        double r_iso_min, double r_iso_max, double kt_scale, double rho,
                        bool charged_only, bool deltaBeta){

    if (ptcl->pt()<5.) return 99999.;

    double absEta = 0;
    if(ptcl->isElectron())  absEta = abs(dynamic_cast<const pat::Electron*>(ptcl)->superCluster()->eta());
    else if(ptcl->isMuon()) absEta = abs(ptcl->eta());


    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron() and absEta>1.479){ deadcone_ch = 0.015;  deadcone_pu = 0.015; deadcone_ph = 0.08; deadcone_nh = 0;}
    else if(ptcl->isMuon())                { deadcone_ch = 0.0001; deadcone_pu = 0.01;  deadcone_ph = 0.01; deadcone_nh = 0.01;}

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh = ptcl->isElectron()? 0. : 0.5;

    double max_pt = kt_scale/r_iso_min;
    double min_pt = kt_scale/r_iso_max;
    double r_iso  = kt_scale/std::max(std::min(ptcl->pt(), max_pt), min_pt);

    for(const pat::PackedCandidate &pfc : *pfcands){
      if(abs(pfc.pdgId())<7) continue;

      double dr = deltaR(pfc, *ptcl);
      if(dr > r_iso) continue;

      if(pfc.charge()==0){								// Neutral
        if(pfc.pt()>ptThresh){
          if(abs(pfc.pdgId())==22 and dr > deadcone_ph)        iso_ph += pfc.pt();	// Photons
          else if (abs(pfc.pdgId())==130 and dr > deadcone_nh) iso_nh += pfc.pt();	// Neutral hadrons
        }
      } else if (pfc.fromPV()>1){
        if(abs(pfc.pdgId())==211 and dr > deadcone_ch) iso_ch += pfc.pt();		// Charged from PV
      } else if(pfc.pt()>ptThresh and dr > deadcone_pu) iso_pu += pfc.pt();		// Charged from PU
    }

    double CorrectedTerm = 0;
    if(ptcl->isMuon())          CorrectedTerm = rho*getEffArea(13, ptcl->eta());
    else if(ptcl->isElectron()) CorrectedTerm = rho*getEffArea(11, dynamic_cast<const pat::Electron*>(ptcl)->superCluster()->eta());

    double iso;
    if(charged_only)   iso = iso_ch;
    else if(deltaBeta) iso = iso_ch + std::max(0., iso_ph + iso_nh - 0.5*iso_pu);
    else               iso = iso_ch + std::max(0., iso_ph + iso_nh - CorrectedTerm*(r_iso*r_iso)/(0.3*0.3));

    return iso/ptcl->pt();
}



/*
 * pfRelIso and pfAbsIso
 */

double tools::pfRelIso(const pat::Muon *mu){
    return tools::pfAbsIso(mu)/mu->pt();
}

double tools::pfAbsIso(const pat::Muon *mu){
    double beta = mu->pfIsolationR03().sumPUPt;
    return (mu->pfIsolationR03().sumChargedHadronPt + std::max(0.0, mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - 0.5*beta));
}

double tools::pfRelIso(const pat::Electron *iE, double myRho){
    return tools::pfAbsIso(iE, myRho)/iE->pt();
}

double tools::pfAbsIso(const pat::Electron *iE, double myRho){
    double CorrectedTerm = myRho*getEffArea(11, iE->superCluster()->eta());
    return (iE->pfIsolationVariables().sumChargedHadronPt + std::max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - CorrectedTerm));
}

double tools::pfRelIso(const pat::Muon *mu, double myRho){
    return tools::pfAbsIso(mu, myRho)/mu->pt();
}

double tools::pfAbsIso(const pat::Muon *mu, double myRho){
    double CorrectedTerm = myRho*getEffArea(13, mu->eta());
    return (mu->pfIsolationR03().sumChargedHadronPt + std::max(0.0, mu->pfIsolationR03().sumNeutralHadronEt + mu->pfIsolationR03().sumPhotonEt - CorrectedTerm));
}


/*
 * Trigger emulation: NOT YET CHECKED/VERIFIED, probably not up to date
 */
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


/*
 * Cut based electron id
 */
float tools::dEtaInSeed(const pat::Electron* ele){
  if(ele->superCluster().isNonnull() and ele->superCluster()->seed().isNonnull()) return ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta();
  else                                                                            return std::numeric_limits<float>::max();
}

bool tools::passed_loose_MVA_FR(const pat::Electron* iE, double mvaValue){
	
	bool passedMVA_loose = false;
	if(iE->pt() > 10 && iE->pt() < 20){//selection 5<pt<10
		passedMVA_loose = false;
		if (TMath::Abs(iE->eta()) < 0.8 ) {
		    passedMVA_loose = mvaValue > -0.86;
		} else if (TMath::Abs(iE->eta()) < 1.479 ) {
		    passedMVA_loose = mvaValue > -0.85;
		} else {
		    passedMVA_loose = mvaValue > -0.81;
		}
	}//end //selection 5<pt<10

	if (iE->pt()>=20  ){
		passedMVA_loose = false;
		if (TMath::Abs(iE->eta()) < 0.8 ) {
		    passedMVA_loose = mvaValue > -0.96;
		} else if (TMath::Abs(iE->eta()) < 1.479 ) {
		    passedMVA_loose = mvaValue > -0.96;
		} else {
		    passedMVA_loose = mvaValue > -0.95;
		}
	}// end 10 pt
	return passedMVA_loose;
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


/*
 * jet selector: NOT YET CHECKED
 */
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


/*
 * Definitions for multi-isolation and cone-corrected lepton pt's
 */
std::map<TString, std::vector<double>> multiIso;
void tools::initMultiIsoConstants(){
  multiIso["VL"] = {0.25, 0.67, 4.4};
  multiIso["L"]  = {0.20, 0.69, 6.0};
  multiIso["M"]  = {0.16, 0.76, 7.2};
  multiIso["T"]  = {0.12, 0.80, 7.2};
  multiIso["VT"] = {0.09, 0.84, 7.2};
}

bool tools::passMultiIsolation(TString level, double miniIso, double jetPtRatio, double jetPtRel){
  return miniIso < multiIso[level][0] and (jetPtRatio > multiIso[level][1] or jetPtRel > multiIso[level][2]);
}

double tools::leptonConeCorrectedPt(double pt, TString level, double miniIso, double jetPtRatio, double jetPtRel){
  if (jetPtRel > multiIso[level][2]) return pt*(1+std::max(0., miniIso - multiIso[level][0]));
  else                               return pt*(std::max(1., multiIso[level][1]/jetPtRatio));
}



/*
 * Some helper functions
 */
double tools::MT_calc(TLorentzVector vect, double met, double met_phi){
    return sqrt(2*vect.Pt()*met*(1 - (TMath::Cos(vect.Phi() - met_phi))));
}

double tools::Mll_calc(TLorentzVector v1, TLorentzVector v2){
    return (v1 + v2).Mag();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




// Tau functions below not yet used/checked
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
            vTaus.insert(std::pair<const reco::PFTau*,   int >(&((*PFTaus)[i]), i));
    }
    return vTaus;
}

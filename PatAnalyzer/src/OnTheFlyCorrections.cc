#include "Majorana/PatAnalyzer/interface/OnTheFlyCorrections.hh"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <math.h>


OnTheFlyCorrections::OnTheFlyCorrections(std::string path, std::string gt, bool isdata){
        std::vector<std::string> runs;
        if(isdata) runs = {"BCD_DATA", "E_DATA", "F_DATA", "p2_DATA"};
        else       runs = {"_MC"};

        for(std::string run : runs){
	  jetUncertainties[run] = new JetCorrectionUncertainty(edm::FileInPath(path+gt+run+"_Uncertainty_AK4PFchs.txt").fullPath());
          std::vector<JetCorrectorParameters> jcParam;
          jcParam.push_back(JetCorrectorParameters(edm::FileInPath(path+gt+run+"_L1FastJet_AK4PFchs.txt").fullPath()));
          jcParam.push_back(JetCorrectorParameters(edm::FileInPath(path+gt+run+"_L2Relative_AK4PFchs.txt").fullPath()));
          jcParam.push_back(JetCorrectorParameters(edm::FileInPath(path+gt+run+"_L3Absolute_AK4PFchs.txt").fullPath()));
          if(isdata) jcParam.push_back(JetCorrectorParameters(edm::FileInPath(path+gt+run+"_L2L3Residual_AK4PFchs.txt").fullPath()));
	  jetCorrectors[run] = new FactorizedJetCorrector(jcParam);
        }

        fIsData = isdata;
}

std::string getRunName(int runNumber){
  if(runNumber < 1)           return "_MC";
  else if(runNumber < 276811) return "BCD_DATA";
  else if(runNumber < 277420) return "E_DATA";
  else if(runNumber < 278801) return "E_DATA";
  else                        return "p2_DATA";

}


float OnTheFlyCorrections::getJECUncertainty(float pt, float eta, int runnumber){
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;
	jetUncertainties[fIsData ? getRunName(runnumber) : "_MC"]->setJetPt(pt);
	jetUncertainties[fIsData ? getRunName(runnumber) : "_MC"]->setJetEta(eta);
	return jetUncertainties[fIsData ? getRunName(runnumber) : "_MC"]->getUncertainty(true);
}

float OnTheFlyCorrections::getJetCorrection(float pt, float corr, float eta, float rho, float area, std::string level = "L3Absolute"){
	// sets the pT back to raw and returns the raw-pT correction factor
	return getJetCorrectionRawPt(pt/corr, eta, rho, area, level);
}

float OnTheFlyCorrections::getJetCorrectionRawPt(float rawpt, float eta, float rho, float area, std::string level = "L3Absolute", int runnumber){
	// slighly redundant considering we have what we have below, but I think that's what frederic was thinking about
        FactorizedJetCorrector* jetCorrector = jetCorrectors[fIsData ? getRunName(runnumber) : "_MC"];
	jetCorrector->setJetEta(eta);
	jetCorrector->setRho(rho);
	jetCorrector->setJetA(area);
	jetCorrector->setJetPt(rawpt); // new raw-pT here...! this is called with fTR->JPt[jetindex]/fTR->JEcorr[jetindex]; in the useranalysisbase.
	std::vector< float > corrections = jetCorrector->getSubCorrections();

	if (level == "L1FastJet")    return corrections.front();
	if (level == "L2Relative")   return corrections[1];
	if (level == "L2L3Residual") return corrections.back();
	return corrections[2];
}



std::pair<float,float> OnTheFlyCorrections::getCorrections(float rawpt, float raweta, float rawnomupt, float phi, float emf, float rho, float area, int runnumber) {
  std::pair<float, float> corr = std::pair<float, float>(0, 0);
  FactorizedJetCorrector* jetCorrector = jetCorrectors[fIsData ? getRunName(runnumber) : "_MC"];
  jetCorrector->setJetEta(raweta);
  jetCorrector->setRho(rho);
  jetCorrector->setJetA(area);
  jetCorrector->setJetPt(rawpt);
  
  std::vector< float > corrections = jetCorrector->getSubCorrections();
  

  float l1corrpt   = rawnomupt*corrections.front(); // l1fastjet corrections were pushed pack first
  float fullcorrpt = rawnomupt*corrections.back();  // full corrections are the last in the vector


  // the corretions for the MET are the difference between l1fastjet and the full corrections on the jet!
  if(emf > 0.9 or fullcorrpt < 15.) return corr;       // skip jets with EMF > 0.9
  
  corr.first  = getPx(l1corrpt - fullcorrpt, phi); // fill the px of the correction in the pair.first
  corr.second = getPy(l1corrpt - fullcorrpt, phi); // and the py in the pair.second
  
  return corr;
}

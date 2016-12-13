#include <vector>
#include <iostream>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

class OnTheFlyCorrections {
	public:
		OnTheFlyCorrections (std::string path, std::string gt, bool isData);
		~OnTheFlyCorrections(){};

		bool fIsData;
		std::map<std::string, FactorizedJetCorrector*>   jetCorrectors;
		std::map<std::string, JetCorrectionUncertainty*> jetUncertainties;

		std::pair< float, float > getCorrections( float rawpt, float raweta, float rawnomupt, float phi, float emf, float rho, float area, int runnumber = -1); // for on the fly corrections
		float getJetCorrection     (float pt, float corr, float eta, float rho, float area, std::string level);     // this function returns, for a given jet the correction factor
		float getJetCorrectionRawPt(float pt,             float eta, float rho, float area, std::string level, int runnumber = -1 );     // same as above, for people who want to call it with the raw-pt already
		float getJECUncertainty(float, float, int runnumber = -1);

		float getPx(float pt, float phi){ return pt*cos(phi); };
		float getPy(float pt, float phi){ return pt*sin(phi); };
};


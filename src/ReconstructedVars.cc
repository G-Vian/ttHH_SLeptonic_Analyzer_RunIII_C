//#include "TTH/CommonClassifier/interface/ReconstructedVars.h"
#include "ReconstructedVars.h"

using namespace std;

ReconstructedVars::ReconstructedVars(bool ttH, bool ttZ, bool Higgs, bool Z, bool ZH):
    reconstruct_ttH(ttH), reconstruct_ttZ(ttZ), reconstruct_Higgs(Higgs), reconstruct_Z(Z), reconstruct_ZH(ZH),
    chi2reco_ttH("RecoTTH_", "Higgs", 116.3, 17.6),
    chi2reco_ttZ("RecoTTZ_", "Z", 86.8, 13.5),
    chi2reco_Higgs("RecoHiggs_", "Higgs", 116.3, 17.6),
    chi2reco_Z("RecoZ_", "Z", 86.8, 13.5),
    chi2reco_ZH("RecoZH_", "ZH", 116.3, 17.6)
{
    // ==================================================
    //init all variables potentially used
    if(reconstruct_ttH)
	{
	    // angles of top system
	    std::map<std::string, double> topSystemAngles = chi2reco_ttH.GetVariableMap_TopSystemAngles();
	    for (auto it = topSystemAngles.begin(); it != topSystemAngles.end(); it++)
		variableMap[it->first] = -999.;

	    // objects of top system
	    std::map<std::string, double> topSystemObjects = chi2reco_ttH.GetVariableMap_TopSystemObjects();
	    for (auto it = topSystemObjects.begin(); it != topSystemObjects.end(); it++)
		variableMap[it->first] = -999.;

	    // angles of higgs system
	    std::map<std::string, double> higgsSystemAngles = chi2reco_ttH.GetVariableMap_BosonSystemAngles();
	    for (auto it = higgsSystemAngles.begin(); it != higgsSystemAngles.end(); it++)
		variableMap[it->first] = -999.;

	    // objects of higgs system
	    std::map<std::string, double> higgsSystemObjects = chi2reco_ttH.GetVariableMap_BosonSystemObjects();
	    for (auto it = higgsSystemObjects.begin(); it != higgsSystemObjects.end(); it++)
		variableMap[it->first] = -999.;
	}

    
    if(reconstruct_ttZ)
	{
	    // angles of top system
	    std::map<std::string, double> topSystemAngles = chi2reco_ttZ.GetVariableMap_TopSystemAngles();
	    for (auto it = topSystemAngles.begin(); it != topSystemAngles.end(); it++)
		variableMap[it->first] = -999.;

	    // objects of top system
	    std::map<std::string, double> topSystemObjects = chi2reco_ttZ.GetVariableMap_TopSystemObjects();
	    for (auto it = topSystemObjects.begin(); it != topSystemObjects.end(); it++)
		variableMap[it->first] = -999.;

	    // angles of higgs system
	    std::map<std::string, double> ZSystemAngles = chi2reco_ttZ.GetVariableMap_BosonSystemAngles();
	    for (auto it = ZSystemAngles.begin(); it != ZSystemAngles.end(); it++)
		variableMap[it->first] = -999.;

	    // objects of higgs system
	    std::map<std::string, double> ZSystemObjects = chi2reco_ttZ.GetVariableMap_BosonSystemObjects();
	    for (auto it = ZSystemObjects.begin(); it != ZSystemObjects.end(); it++)
		variableMap[it->first] = -999.;
	}
    
    if(reconstruct_Higgs)
	{
	    // objetcs of higgs 
	    std::map<std::string, double> higgsOnlyObjects = chi2reco_Higgs.GetVariableMap_BosonOnlyObjects();
	    for (auto it = higgsOnlyObjects.begin(); it != higgsOnlyObjects.end(); it++)
		variableMap[it->first] = -999.;
        
	    // angles of higgs 
	    std::map<std::string, double> higgsOnlyAngles = chi2reco_Higgs.GetVariableMap_BosonOnlyAngles();
	    for (auto it = higgsOnlyAngles.begin(); it != higgsOnlyAngles.end(); it++)
		variableMap[it->first] = -999.;
	}
    
    if(reconstruct_Z)
	{
	    // objects of Z 
	    std::map<std::string, double> ZOnlyObjects = chi2reco_Z.GetVariableMap_BosonOnlyObjects();
	    for (auto it = ZOnlyObjects.begin(); it != ZOnlyObjects.end(); it++)
		variableMap[it->first] = -999.;
        
	    // angles of Z
	    std::map<std::string, double> ZOnlyAngles = chi2reco_Z.GetVariableMap_BosonOnlyAngles();
	    for (auto it = ZOnlyAngles.begin(); it != ZOnlyAngles.end(); it++)
		variableMap[it->first] = -999.;
	}
    
    if(reconstruct_ZH)
	{
	    // objects of Z
	    std::map<std::string, double> ZHOnlyObjects = chi2reco_ZH.GetVariableMap_BosonOnlyObjects();
	    for (auto it = ZHOnlyObjects.begin(); it != ZHOnlyObjects.end(); it++)
		variableMap[it->first] = -999.;
        
	    // angles of Z
	    std::map<std::string, double> ZHOnlyAngles = chi2reco_ZH.GetVariableMap_BosonOnlyAngles();
	    for (auto it = ZHOnlyAngles.begin(); it != ZHOnlyAngles.end(); it++)
		variableMap[it->first] = -999.;
	}
}

ReconstructedVars::~ReconstructedVars() {}


void ReconstructedVars::ResetVariableMap()
{
    for (auto it = variableMap.begin(); it != variableMap.end(); it++)
        it->second = -999.;
}

void ReconstructedVars::SetWP(double WP){
    btagMcut = WP;
}

void ReconstructedVars::SetLooseWP(double LooseWP){
    btagLooseMcut = LooseWP;
}




std::map<std::string, double> ReconstructedVars::GetReconstructedVars(
								      const std::vector<TLorentzVector> &selectedLeptonP4,
								      const std::vector<TLorentzVector> &selectedJetP4,
								      const std::vector<double> &selectedJetCSV,
								      const TLorentzVector &metP4)
{
    //check if CSVWP is set properly
    if(btagMcut<0){
	std::cerr << "Please set CSV Wp properly\n" << std::endl;
        throw std::exception();
    }
    
    //check if CSVWP is set properly
    if(btagLooseMcut<0){
	std::cerr << "Please set Loose CSV Wp properly\n" << std::endl;
        throw std::exception();
    }
    // Reset all map entries to -999 so that noting is left over from the last event
    ResetVariableMap();

    // get number of jets
    int nJets = selectedJetP4.size();

    // Start Reconstruction
    if(reconstruct_ttH)
	{
	    chi2reco_ttH.ReconstructTTXSystem(selectedLeptonP4[0], selectedJetP4, selectedJetCSV, metP4, btagMcut);

	    ///////////// fill map /////////////
	    // angles of top system
	    std::map<std::string, double> topSystemAngles = chi2reco_ttH.FillVariableMap_TopSystemAngles();
	    for (auto it = topSystemAngles.begin(); it != topSystemAngles.end(); it++)
		variableMap[it->first] = it->second;

	    // objects of top system
	    std::map<std::string, double> topSystemObjects = chi2reco_ttH.FillVariableMap_TopSystemObjects();
	    for (auto it = topSystemObjects.begin(); it != topSystemObjects.end(); it++)
		variableMap[it->first] = it->second;

	    // only fill higgs system variables if more/equal 6 jets
	    if (nJets >= 6)
		{
		    // angles of higgs system
		    std::map<std::string, double> higgsSystemAngles = chi2reco_ttH.FillVariableMap_BosonSystemAngles();
		    for (auto it = higgsSystemAngles.begin(); it != higgsSystemAngles.end(); it++)
			variableMap[it->first] = it->second;

		    // objects of higgs system
		    std::map<std::string, double> higgsSystemObjects = chi2reco_ttH.FillVariableMap_BosonSystemObjects();
		    for (auto it = higgsSystemObjects.begin(); it != higgsSystemObjects.end(); it++)
			variableMap[it->first] = it->second;
		}
	}

    if(reconstruct_ttZ)
	{
	    chi2reco_ttZ.ReconstructTTXSystem(selectedLeptonP4[0], selectedJetP4, selectedJetCSV, metP4, btagMcut);

	    ///////////// fill map /////////////
	    // angles of top system
	    std::map<std::string, double> topSystemAngles = chi2reco_ttZ.FillVariableMap_TopSystemAngles();
	    for (auto it = topSystemAngles.begin(); it != topSystemAngles.end(); it++)
		variableMap[it->first] = it->second;

	    // objects of top system
	    std::map<std::string, double> topSystemObjects = chi2reco_ttZ.FillVariableMap_TopSystemObjects();
	    for (auto it = topSystemObjects.begin(); it != topSystemObjects.end(); it++)
		variableMap[it->first] = it->second;

	    // only fill Z system variables if more/equal 6 jets
	    if (nJets >= 6)
		{
		    // angles of higgs system
		    std::map<std::string, double> ZSystemAngles = chi2reco_ttZ.FillVariableMap_BosonSystemAngles();
		    for (auto it = ZSystemAngles.begin(); it != ZSystemAngles.end(); it++)
			variableMap[it->first] = it->second;

		    // objects of Z system
		    std::map<std::string, double> ZSystemObjects = chi2reco_ttZ.FillVariableMap_BosonSystemObjects();
		    for (auto it = ZSystemObjects.begin(); it != ZSystemObjects.end(); it++)
			variableMap[it->first] = it->second;
		}
	}
    
    
    //needs to be done after system reco because Init will probably not be called again
    if(reconstruct_Higgs)
	{
	    chi2reco_Higgs.ReconstructBosonOnly(selectedJetP4, selectedJetCSV, btagMcut, btagLooseMcut);

	    ///////////// fill map /////////////

	    // only fill boson variables if more/equal 4 jets
	    if (nJets >= 4)
		{
		    // objects of higgs system
		    std::map<std::string, double> HiggsOnlyObjects = chi2reco_Higgs.FillVariableMap_BosonOnlyObjects();
		    for (auto it = HiggsOnlyObjects.begin(); it != HiggsOnlyObjects.end(); it++)
			variableMap[it->first] = it->second;
            
		    // objects of higgs system
		    std::map<std::string, double> HiggsOnlyAngles = chi2reco_Higgs.FillVariableMap_BosonOnlyAngles();
		    for (auto it = HiggsOnlyAngles.begin(); it != HiggsOnlyAngles.end(); it++)
			variableMap[it->first] = it->second;
		}
	}
    
    if(reconstruct_Z)
	{
	    chi2reco_Z.ReconstructBosonOnly(selectedJetP4, selectedJetCSV, btagMcut, btagLooseMcut);

	    ///////////// fill map /////////////

	    // only fill boson variables if more/equal 4 jets
	    if (nJets >= 4)
		{
		    // objects of higgs system
		    std::map<std::string, double> ZOnlyObjects = chi2reco_Z.FillVariableMap_BosonOnlyObjects();
		    for (auto it = ZOnlyObjects.begin(); it != ZOnlyObjects.end(); it++)
			variableMap[it->first] = it->second;
            
		    // objects of higgs system
		    std::map<std::string, double> ZOnlyAngles = chi2reco_Z.FillVariableMap_BosonOnlyAngles();
		    for (auto it = ZOnlyAngles.begin(); it != ZOnlyAngles.end(); it++)
			variableMap[it->first] = it->second;
		}
	}
    
    if(reconstruct_ZH)
	{
	    chi2reco_ZH.ReconstructZH(selectedJetP4, selectedJetCSV, btagMcut, btagLooseMcut);

	    ///////////// fill map /////////////

	    // only fill boson variables if more/equal 4 jets
	    if (nJets >= 4)
		{
		    // objects of higgs system
		    std::map<std::string, double> ZHOnlyObjects = chi2reco_ZH.FillVariableMap_BosonOnlyObjects();
		    for (auto it = ZHOnlyObjects.begin(); it != ZHOnlyObjects.end(); it++)
			variableMap[it->first] = it->second;
            
		    // objects of higgs system
		    std::map<std::string, double> ZHOnlyAngles = chi2reco_ZH.FillVariableMap_BosonOnlyAngles();
		    for (auto it = ZHOnlyAngles.begin(); it != ZHOnlyAngles.end(); it++)
			variableMap[it->first] = it->second;
		}
	}
    
    
    
    return variableMap;
}


std::map<std::string, double> ReconstructedVars::GetVariables()
{
    ResetVariableMap();
    return variableMap;
}

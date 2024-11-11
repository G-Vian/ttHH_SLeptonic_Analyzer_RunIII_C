//#include "TTH/CommonClassifier/interface/HypothesisCombinatorics.h"
#include "HypothesisCombinatorics.h"

HypothesisCombinatorics::HypothesisCombinatorics(const std::string& optional_varstring){
    if(!optional_varstring.empty()) FillVariableNameList(optional_varstring, optional_varlist);
    //     TMVA::PyMethodBase::PyInitialize();
}

HypothesisCombinatorics::~HypothesisCombinatorics()
{
    //     TMVA::PyMethodBase::PyFinalize();
}


void HypothesisCombinatorics::FillVariableNameList(const TString& variables, std::vector<std::string>& vec_to_fill){
    TString tosplit = "";
    if(variables.EndsWith(".csv")){
	ifstream in(variables.Data());
	if(!in.is_open()) return;
	std::cout << "DEBUG: reading from file " << variables.Data() << std::endl;
	std::string line;
	while(std::getline(in, line)){
	    if(!tosplit.EqualTo("")) tosplit += "," + line;
	    else tosplit = line;
	}
    }
    else tosplit = variables;
    TObjArray* tokens = tosplit.Tokenize(",");
    for(int i = 0; i<tokens->GetEntries(); i++)
	{
	    vec_to_fill.push_back( ((TObjString *)(tokens->At(i)))->String().Data() );
	}
}

void HypothesisCombinatorics::SetBtagWP(const double& in_btagWP){
    mvars->SetWP(in_btagWP);
}

void HypothesisCombinatorics::SetBtagLooseWP(const double& in_btagLooseWP){
    mvars->SetLooseWP(in_btagLooseWP);
}

std::vector<std::string> HypothesisCombinatorics::GetVariableNames() const{
    std::vector<std::string> outvector;
  
    for(auto& it : mvars->GetVariables()){
	outvector.push_back(it.first);
    }

    outvector.push_back(bdtoutput_name);
    return outvector;
}

std::map<std::string, float> HypothesisCombinatorics::GetBestPermutation(const std::vector<TLorentzVector> &selectedLeptonP4,
                                                                         const std::vector<TLorentzVector> &selectedJetP4,
                                                                         const std::vector<double> &selectedJetCSV,
                                                                         const TLorentzVector &metP4)
{
    Float_t best_bdtout = -1;
    std::vector<int> permutation;
    std::vector<int> best_reco_idx;
    std::map<std::string, float> varmap;
    
    //TODO jet cleaning
    
    unsigned int NJets = selectedJetP4.size();
    std::cout << "NJets: "<< NJets << "\n";
    std::cout << "minJets: "<< minJets << "\n";
    
    std::cout << "BDT: "<< bdt_name << "\n";
    
    int i = 0;
    int j = 0;
    if(NJets < minJets)
	{
	    // std::cout << "WARNING: found not enough jets! (found: " << NJets << ", required: " << minJets << ")\n";
	    mvars->ResetVariableMap();
	    i++;
	    std::cout << "Reset variables for lacking enough jets " << i << "\n";
	}
    else
	{
	    // iterate permutations to find best
	    for(unsigned int permutation_idx=0; permutation_idx < getPermutator()->get_Npermutations(NJets); permutation_idx++)
		{
		    getPermutator()->get_permutation(&permutation, NJets, permutation_idx);
		    if(mvars->SkipEvent(selectedJetP4, selectedJetCSV, permutation)) continue;
            
		    mvars->FillMVAvarMap(selectedLeptonP4, selectedJetP4, selectedJetCSV, metP4, permutation);
            
		    // eval BDT
		    Float_t recbdtout = reader_th->EvaluateMVA(bdt_name);
            
		    if(recbdtout > best_bdtout)
			{
			    best_bdtout   = recbdtout;
			    best_reco_idx = permutation;
			}
		}
	    if(best_reco_idx.size() == 0){
		j++;
		std::cout << "WARNING: found no list of best permutation indices! Resetting vars" << j << "\n";
		mvars->ResetVariableMap();
            
	    }
	    else
		{
		    mvars->FillMVAvarMap(selectedLeptonP4, selectedJetP4, selectedJetCSV, metP4, best_reco_idx);
		}
	}
    
    varmap = mvars->GetVariables();
    varmap[bdtoutput_name.c_str()] = best_bdtout;
    
    return varmap;
}

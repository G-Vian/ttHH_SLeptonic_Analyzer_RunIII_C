#ifndef _AngularVariables_h
#define _AngularVariables_h

#include <vector>
#include <map>
#include "TVector.h"
#include "TLorentzVector.h"

//using namespace std;

class AngularVariables{
    
    public:
        AngularVariables();
        virtual ~AngularVariables();
        void ReconstructTopSystem(const TLorentzVector& lepton, const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, TLorentzVector met, double CSVMwp);
        std::map<TString,double> GetAngularVariables();
        
    private:
        //double CSVLwp, CSVMwp, CSVTwp;
        TLorentzVector hadtop;
        TLorentzVector leptop;
        TLorentzVector Lepton;
        TLorentzVector hadb;
        TLorentzVector lepb;
        TLorentzVector hadw;
        TLorentzVector lepw;
        TLorentzVector mET;

        // mass templates
        const double top_mass_had = 165;
        const double top_mass_lep = 168;
        const double W_mass = 80.0;
 
        const double sigma_hadW   = 10;
        const double sigma_hadTop = 17;
        const double sigma_lepTop = 26;
        
        
    
    
};

#endif

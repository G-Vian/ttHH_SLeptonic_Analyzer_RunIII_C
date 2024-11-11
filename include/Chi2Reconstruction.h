#ifndef _Chi2Reconstruction_h
#define _Chi2Reconstruction_h

#include <vector>
#include <map>
#include "TVector.h"
#include "TLorentzVector.h"

//using namespace std;

class Chi2Reconstruction{
    
    public:
        Chi2Reconstruction(std::string tag, std::string bosonName, double templateMass_Boson, double templateSigma_Boson);

        virtual ~Chi2Reconstruction();

        void InitReconstruction(const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, double CSVMwp);
        void ReconstructNeutrino(const TLorentzVector& lepton, TLorentzVector met);

        void ReconstructTTXSystem(const TLorentzVector& lepton, const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, TLorentzVector met, double CSVMwp);
        
        void ReconstructBosonOnly(const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, double CSVMwp, double LooseCSVMwp);
    
        void ReconstructZH(const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, double CSVMwp, double LooseCSVMwp);

        double GetDeltaPhi(double phi1, double phi2);

        // top system 
        std::map<std::string,double> GetVariableMap_TopSystemAngles();
        std::map<std::string,double> GetVariableMap_TopSystemObjects();
        std::map<std::string,double> FillVariableMap_TopSystemAngles();
        std::map<std::string,double> FillVariableMap_TopSystemObjects();
        
        // higgs system
        std::map<std::string,double> GetVariableMap_BosonSystemAngles();
        std::map<std::string,double> GetVariableMap_BosonSystemObjects();
        std::map<std::string,double> FillVariableMap_BosonSystemAngles();
        std::map<std::string,double> FillVariableMap_BosonSystemObjects();
        
        //boson Only
        std::map<std::string, double> GetVariableMap_BosonOnlyAngles();
        std::map<std::string, double> GetVariableMap_BosonOnlyObjects();
        std::map<std::string, double> FillVariableMap_BosonOnlyAngles();
        std::map<std::string, double> FillVariableMap_BosonOnlyObjects();
        //HERE ARE THE FILL FUNCTIONS MISSING
        

    private:
        // tags
        std::string tag;
        std::string bosonName;


        // variables for reconstruction
        double btagCut;
        double btagLooseCut;
        double metPz[2];
        double chi;

        int nJets;
        int nBtags;
        int nLooseBtags;
        int nOtherJets; 

        bool reconstructBoson;
        int allowedMistags;

        // best candidates
        TLorentzVector best_Boson;
        TLorentzVector best_B1;
        TLorentzVector best_B2;
        TLorentzVector best_BosonOnly1;
        TLorentzVector best_BosonOnly2;
        TLorentzVector best_B1Only;
        TLorentzVector best_B2Only;
        TLorentzVector best_B3Only;
        TLorentzVector best_B4Only;

        TLorentzVector best_topLep;
        TLorentzVector best_topHad;
        TLorentzVector best_bHad;
        TLorentzVector best_bLep;
        TLorentzVector best_wHad;
        TLorentzVector best_wLep;
        TLorentzVector best_lepton;
        TLorentzVector best_mET;
        
        // best chi2 values
        double minChi = 1000000.;
        double best_chi2_topHad = 10000.;
        double best_chi2_topLep = 10000.;
        double best_chi2_wHad = 10000.;
        double best_chi2_Boson = 10000.;
        double best_chi2_BosonOnly = 10000.;

        // mass templates
        const double templateMass_W = 83.7;
        const double templateMass_topHad = 173.1;
        const double templateMass_topLep = 168; // not updated
 
        const double templateSigma_W = 11.9;
        const double templateSigma_topHad = 18.5;
        const double templateSigma_topLep = 26; // not updated
    
        const double templateMass_Z = 86.8;
        const double templateSigma_Z = 13.5;
        
        // template masses determined by init
        double templateMass_Boson;
        double templateSigma_Boson;
};

#endif

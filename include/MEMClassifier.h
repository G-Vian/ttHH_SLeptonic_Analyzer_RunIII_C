
#ifndef MEMCLASSIFIER_H
#define MEMCLASSIFIER_H

#include "Integrand.h"
#include "JetLikelihood.h"
#include "TFile.h"
#include "TH3D.h"
#include "TLorentzVector.h"

class MEMResult {
public:
    double blr_4b;
    double blr_2b;

    //likelihood ratio
    double p;

    //individual signal and background probabilities
    double p_sig;
    double p_bkg;

    //Integration uncertainties of the probabilities
    double p_err_sig;
    double p_err_bkg;

    //Number of permutations per hypothesis
    double n_perm_sig;
    double n_perm_bkg;

    vector<double> p_variated;


};

class MEMClassifier {

public:

    // Jet type: Needed to decide:
    //   - which jets to put into which slot of the MEM
    //   - which transfer function to use

    // Hypothesis
    enum Hypothesis {
        HYPO_DEBUG,         // Don't run the MEM. Print debug info
        SL_0W2H2T,     // Default SL MEM: Integrate over light jets
        SL_1W2H2T,     // Default SL MEM: Integrate over light jets
        SL_2W2H2T,     // Fully reconstructed hypothesis
        SL_2W2H2T_SJ,  // Boosted SL MEM: 2 light jets + bjets
        DL_0W2H2T,     // Default DL MEM
    };

    enum JetType {
        RESOLVED,
        BOOSTED_LIGHT,
        BOOSTED_B
    };

    enum Systematic {
        NOMINAL,
        JESUP,
        JESDOWN,
        JERUP,
        JERDOWN
    };

    //BTag Medium Working Point
    float btagWP;

    //The constructor loads the transfer functions and b-tag PDF-s
    MEMClassifier(int verbosity, const char* _btag_prefix, const char* _era);
    MEMClassifier();
    ~MEMClassifier();

    // Call this method to return the MEM output, provide all necessary inputs.
    // Jet CSV and type should be sorted the same way as jet p4.
    // We could also write a class to contain the jet CSV, type, and p4 information
    std::vector<MEMResult> GetOutput(
        const Hypothesis hypo,
        const std::vector<TLorentzVector>& selectedLeptonP4,
        const std::vector<double>& selectedLeptonCharge,
        const std::vector<vector<TLorentzVector>>& selectedJetP4,
        const std::vector<vector<double>>& selectedJetCSV,
        const std::vector<vector<JetType>>& selectedJetType,
        const std::vector<vector<vector<double>>>& selectedJetVariation,
        const std::vector<TLorentzVector>& looseSelectedJetP4,
        const std::vector<double>& looseSelectedJetCSV,
        TLorentzVector& metP4,
        const std::vector<bool>& changes_jet_category,
        int ncalls=-1
    );

    //Default method to get the MEM output
    std::vector<MEMResult> GetOutput(
        const std::vector<TLorentzVector>& selectedLeptonP4,
        const std::vector<double>& selectedLeptonCharge,
        const std::vector<vector<TLorentzVector>>& selectedJetP4,
        const std::vector<vector<double>>& selectedJetCSV,
        const std::vector<vector<JetType>>& selectedJetType,
        const std::vector<vector<vector<double>>>& selectedJetVariation,
        TLorentzVector& metP4,
        const std::vector<bool>& changes_jet_category
    );

    void setup_mem(
        const Hypothesis hypo,
        const std::vector<TLorentzVector>& selectedLeptonP4,
        const std::vector<double>& selectedLeptonCharge,
        const std::vector<TLorentzVector>& selectedJetP4,
        const std::vector<double>& selectedJetCSV,
        const std::vector<JetType>& selectedJetType,
        const std::vector<vector<double>>& selectedJetVariation,
        const std::vector<TLorentzVector>& looseSelectedJetP4,
        const std::vector<double>& looseSelectedJetCSV,
        TLorentzVector& metP4,
        std::vector<MEM::Object*>& objs,
        MEMResult& res,
        const std::vector<bool>& changes_jet_category
    );

    void setup_mem_impl(
        const std::vector<TLorentzVector>& selectedLeptonP4,
        const std::vector<double>& selectedLeptonCharge,
        const std::vector<TLorentzVector>& selectedJetP4,
        const std::vector<double>& selectedJetCSV,
        const std::vector<JetType>& selectedJetType,
        const std::vector<vector<double>>& selectedJetVariation,
        const std::vector<TLorentzVector>& looseSelectedJetP4,
        const std::vector<double>& looseSelectedJetCSV,
        TLorentzVector& metP4,
        std::vector<MEM::Object*>& objs,
        MEMResult& res,
        const std::vector<bool>& changes_jet_category
    );

    void setup_mem_sl_2w2h2t_sj(
        const std::vector<TLorentzVector>& selectedLeptonP4,
        const std::vector<double>& selectedLeptonCharge,
        const std::vector<TLorentzVector>& selectedJetP4,
        const std::vector<double>& selectedJetCSV,
        const std::vector<JetType>& selectedJetType,
        const std::vector<vector<double>>& selectedJetVariation,
        const std::vector<TLorentzVector>& looseSelectedJetP4,
        const std::vector<double>& looseSelectedJetCSV,
        TLorentzVector& metP4,
        std::vector<MEM::Object*>& objs,
        MEMResult& res,
        const std::vector<bool>& changes_jet_category
    );

    // returns the category of the last evaluated Event
    std::string GetCategoryOfLastEvaluation() const;

    double GetBTagLikelihoodRatio(
        const std::vector<TLorentzVector>& selectedJetP4,
        const std::vector<double>& selectedJetCSV,
        std::vector<unsigned int>& out_best_perm,
        double& out_P_4b,
        double& out_P_2b
    );

  private:
    //Holds the transfer functions
    TFile* transfers;
    //holds the b-tag PDF-s
    TFile* btagfile;

    //This is the MEM workhorse from the ETH side
    MEM::Integrand* integrand;
    MEM::MEMConfig cfg;

    MEM::JetLikelihood* blr;

    //Convenience functions to construct MEM input objects
    MEM::Object* make_jet(double pt, double eta, double phi, double mass, double istagged, double csv, bool is_subjet) const;
    MEM::Object* make_lepton(double pt, double eta, double phi, double mass, double charge) const;

    // Returns the transfer function corresponding to a jet flavour and eta
    TF1* getTransferFunction(const char* flavour, double eta) const;
    double GetJetBProbability(const char* prefix, const char* flavour, double pt, double eta, double bdisc);
    MEM::JetProbability GetJetBProbabilities(const TLorentzVector& p4, double bdisc);
    
    TH3D* GetBTagPDF(const char* prefix, const char* flavour);

    long unsigned int numMaxJets = 9;
    long unsigned int numMaxJetsBLR = 9;

    float btagcut = 0.8484;

    const char* btag_prefix;

    const char* era;
};

#endif

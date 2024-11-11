#ifndef TTH_ReconstructedVars_H
#define TTH_ReconstructedVars_H
#include <vector>
#include <map>
#include <math.h>
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include "TTH/CommonClassifier/interface/Chi2Reconstruction.h"

// class to evaluate lepton plus jets BDT set
class ReconstructedVars
{

  public:
    ReconstructedVars(bool reconstruct_ttH, bool reconstruct_ttZ, bool Higgs, bool Z, bool ZH);
    ~ReconstructedVars();

    // Call this method to return the Reconstruced vars, provide all necessary inputs. Jet CSV should be sorted the same way as jet p4.
    // We could also write a class to contain the jet CSV and p4 information
    std::map<std::string, double> GetReconstructedVars(
        const std::vector<TLorentzVector> &selectedLeptonP4,
        const std::vector<TLorentzVector> &selectedJetP4,
        const std::vector<double> &selectedJetCSV,
        const TLorentzVector &metP4);

    // return the variable names with default initialization
    std::map<std::string, double> GetVariables();
    void SetWP(double WP);
    void SetLooseWP(double LooseWP);

  private:
    void ResetVariableMap();

    std::string category;
    std::map<std::string, double> variableMap;
    double btagMcut = -99;
    double btagLooseMcut = -99;

    bool reconstruct_ttH;
    bool reconstruct_ttZ;
    bool reconstruct_Higgs;
    bool reconstruct_Z;
    bool reconstruct_ZH;

    Chi2Reconstruction chi2reco_ttH;
    Chi2Reconstruction chi2reco_ttZ;
    Chi2Reconstruction chi2reco_Higgs;
    Chi2Reconstruction chi2reco_Z;
    Chi2Reconstruction chi2reco_ZH;
};

#endif

#ifndef TTH_MVAVARS_H
#define TTH_MVAVARS_H
#include <vector>
#include <map>
#include <math.h>
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include "CommonBDTvars.h"
#include "MEMClassifier.h"
#include "AngularVariables.h"
#include "MVAvarsBase.h"

// class to evaluate lepton plus jets BDT set
class MVAvars : public MVAvarsBase
{

  public:
    MVAvars(const char* era="2018");
    ~MVAvars();

    void FillMVAvarMap(const std::vector<TLorentzVector> &selectedLeptonP4,
                      const std::vector<TLorentzVector> &selectedJetP4,
                      const std::vector<double> &selectedJetCSV,
                      const TLorentzVector &metP4);

  private:
    CommonBDTvars bdtvar;
    MEMClassifier mem;
    // MEMClassifier(int verbosity, const char* _btag_prefix, const char* _era);

};

#endif

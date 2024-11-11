#ifndef TTH_THHHYPOTHESISCOMBINATORICS_H
#define TTH_THHHYPOTHESISCOMBINATORICS_H

#include <vector>
#include <map>
#include "TLorentzVector.h"
#include "TTH/CommonClassifier/interface/HypothesisCombinatorics.h"
#include "TTH/CommonClassifier/interface/MVAvarsJABDTthh.h"


class thhHypothesisCombinatorics : public HypothesisCombinatorics
{
    public:
        thhHypothesisCombinatorics(const std::string& weightpath, const std::string& optional_varstring);
        ~thhHypothesisCombinatorics();

        PermutationManager* getPermutator(){return &permutator;}

    private:
        static const unsigned int minJets = 6;
        static PermutationManager permutator;
};

#endif

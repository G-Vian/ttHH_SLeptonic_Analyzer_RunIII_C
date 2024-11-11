////////////////////////////////////////////////////////////////////////////////
//
// PermutationManager
// ---------------
//
//            13/05/2019 Marco Link
////////////////////////////////////////////////////////////////////////////////


#include "PermutationManager.h"


// #include <cassert>
#include <iostream>
#include <fstream>
// #include <sstream>


using namespace std;


////////////////////////////////////////////////////////////////////////////////
// construction / destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
PermutationManager::PermutationManager(unsigned int N_idx, unsigned int Jets_Max) : Nidx(N_idx), JetsMax(Jets_Max)
{
    if(Nidx > JetsMax)
    {
        cout << "error: number of indizes > max number of jets" << endl;
    }

    for(unsigned int Njets=Nidx; Njets<=JetsMax; Njets++)
    {
        std::vector<std::vector<int>> temp;

        unsigned int comb[Nidx];
        for(unsigned int idx=0; idx<Nidx; idx++)
        {
            comb[idx] = Njets-1;
        }

        unsigned int sum = 999;
        bool reduce;
        bool skip;

        // fill permutations
        while(sum > 0)
        {
            std::vector<int> combination;
            reduce = true;
            skip = false;
            sum = 0;
            for(unsigned int idx=0; idx<Nidx; idx++)
            {
                combination.push_back(comb[idx]);

                // skip combinations with double indizes
                for(unsigned int idx2=0; idx2<idx; idx2++)
                {
                    if(comb[idx] == comb[idx2]) skip=true;
                }
            }
            if(!skip) temp.push_back(combination);

            for(unsigned int idx=0; idx<Nidx; idx++)
            {
                // reduce by 1
                if(reduce && comb[idx]>0)
                {
                    comb[idx] = comb[idx] - 1;
                    reduce = false;
                }
                else if (reduce)
                {
                    comb[idx] = Njets-1;
                }

                sum += comb[idx];
            }
        }
        permutations.push_back(temp);
    }
}


//______________________________________________________________________________
PermutationManager::~PermutationManager()
{

}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////
unsigned int PermutationManager::get_Npermutations(unsigned int Njets)
{
    return permutations.at(getIndex(Njets)).size();
}


void PermutationManager::get_permutation(std::vector<int>* permutation, unsigned int Njets, unsigned int i)
{
    if(i >= get_Npermutations(Njets)) cout << "error: permutation index out of range!" << endl;

    *permutation = permutations[getIndex(Njets)].at(i);
}


void PermutationManager::show()
{
    cout << "Permutationmanager for " << Nidx << " indizes" << endl;
    cout << "NJets (Index),     theo,       generated permutations" << endl;
    for(unsigned int i=Nidx; i<=JetsMax; i++)
    {
        cout << setw(2) << i << "(" << setw(2) << getIndex(i) << "), " << setw(15) << TMath::Factorial(i)/TMath::Factorial(i - Nidx) << ", " << setw(15) << get_Npermutations(i) << endl;
    }
    cout << endl;
}


void PermutationManager::show(unsigned int Njets)
{
    std::vector<std::vector<int>> combinations = permutations.at(getIndex(Njets));

    cout << "permutations for " << Nidx << " indizes and " << Njets << " jets" << endl;

    for(unsigned int n=0; n<get_Npermutations(Njets); n++)
    {
        for(unsigned int idx=0; idx<Nidx; idx++)
        {
            cout << combinations.at(n).at(idx) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}


unsigned int PermutationManager::getIndex(unsigned int Njets)
{
    if(Njets > JetsMax)
    {
        cout << "error: permutation not generated for " << Njets << " increase maximum jet number!" << endl;
        return 0;
    }
    else if(Njets < Nidx)
    {
        cout << "error: jet number " << Njets << " smaller than number of indizes! Automagically fixing problem ... NOT" << endl;
        return 0;
    }

    return Njets - Nidx;
}

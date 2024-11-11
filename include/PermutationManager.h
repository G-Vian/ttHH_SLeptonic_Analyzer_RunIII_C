#ifndef THANALYSIS_PERMUTATIONMANAGER_H
#define THANALYSIS_PERMUTATIONMANAGER_H

#include <iostream>
#include <iomanip>
#include <vector>
#include "TMath.h"



using namespace std;

class PermutationManager
{
  //
  // construction / destruction
  //
public:
  PermutationManager(unsigned int N_idx, unsigned int Jets_Max);
  virtual ~PermutationManager();


public:
    unsigned int get_Npermutations(unsigned int Njets);
    void get_permutation(std::vector<int>* permutation, unsigned int Njets, unsigned int i);

    void show();
    void show(unsigned int Njets);




private:
    const unsigned int Nidx;
    const unsigned int JetsMax;
    std::vector<std::vector<std::vector<int>>> permutations;

    unsigned int getIndex(unsigned int Njets);
};

#endif

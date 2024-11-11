#include "../include/Chi2Reconstruction.h"
#include <iostream>
using namespace std;

Chi2Reconstruction::Chi2Reconstruction(std::string t, std::string n, double mBoson, double sBoson):
    tag(t), bosonName(n), templateMass_Boson(mBoson), templateSigma_Boson(sBoson) {}

Chi2Reconstruction::~Chi2Reconstruction() {}

void Chi2Reconstruction::InitReconstruction(const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, double CSVMwp) {
    
    // init best candidates empty
    best_Boson.SetPxPyPzE(0.,0.,0.,0.);

    best_B1.SetPxPyPzE(0.,0.,0.,0.);
    best_B2.SetPxPyPzE(0.,0.,0.,0.);

    best_topLep.SetPxPyPzE(0.,0.,0.,0.);
    best_topHad.SetPxPyPzE(0.,0.,0.,0.);

    best_bHad.SetPxPyPzE(0.,0.,0.,0.);
    best_bLep.SetPxPyPzE(0.,0.,0.,0.);

    best_wHad.SetPxPyPzE(0.,0.,0.,0.);
    best_wLep.SetPxPyPzE(0.,0.,0.,0.);

    best_lepton.SetPxPyPzE(0.,0.,0.,0.);
    best_mET.SetPxPyPzE(0.,0.,0.,0.);
    
    /*
    best_BosonOnly.SetPxPyPzE(0.,0.,0.,0.); //vars for boson only reconstruction
    best_B1Only.SetPxPyPzE(0.,0.,0.,0.);
    best_B2Only.SetPxPyPzE(0.,0.,0.,0.); */

    // z components of MET
    metPz[0] = 0.;
    metPz[1] = 0.;
    
    // chi2 default value
    chi = 1e6;
    minChi = 1e5;

    // get number of tags and btags
    btagCut = CSVMwp;
    nJets = int(jets.size());
    int nBtags = 0;
    for (int i = 0; i < nJets; i++)
        if (csvs[i] > btagCut) nBtags++;
    int nUntags = nJets-nBtags;
        
    // only try to reconstruct Higgs/Z if more than six jets found
    reconstructBoson = false;
    if (nJets >= 6)
        reconstructBoson = true;

    // figure out number of allowed mistags
    allowedMistags = 0;
    // allow one mistag if only 4 jets
    if(nJets==4) allowedMistags += 1;

    // allow as many mistags as needed to account for hadronic W
    if(nUntags==0) allowedMistags += 2;
    if(nUntags==1) allowedMistags += 1;

    // allow as many mistags as needed to account for 4 bs
    if(reconstructBoson)
    {
        if(nBtags==2) allowedMistags += 2;
        if(nBtags==3) allowedMistags += 1;
    }

    //std::cout << std::endl;
    //std::cout << "jets|tags " << nJets << "|" << nBtags << std::endl;
    //std::cout << "allowedMistags: " << allowedMistags << std::endl;

}



///////////////////////////// Neutrino /////////////////////////////////////////////////////////
void Chi2Reconstruction::ReconstructNeutrino(const TLorentzVector& lepton, TLorentzVector met) 
{
    
    double energyLep = lepton.E();
    double a = (templateMass_W*templateMass_W/(2.0*energyLep)) + (lepton.Px()*met.Px() + lepton.Py()*met.Py())/energyLep;

    double radical = (2.0*lepton.Pz()*a/energyLep)*(2.0*lepton.Pz()*a/energyLep);
    radical = radical - 4.0*(1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep))*(met.Px()*met.Px() + met.Py()*met.Py()- a*a);
    // set radical to zero if smaller than zero
    if (radical < 0.0)
        radical = 0;

    metPz[0] = (lepton.Pz()*a/energyLep) + 0.5*sqrt(radical);
    metPz[0] = metPz[0] / (1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep));
    metPz[1] = (lepton.Pz()*a/energyLep) - 0.5*sqrt(radical);
    metPz[1] = metPz[1] / (1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep));


}

void Chi2Reconstruction::ReconstructBosonOnly(const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, double CSVMwp, double LooseCSVMwp)
{    
    
    // chi2 default value (class vars)
    //chi = 1e6;
    //minChi = 1e5;
    
    //set csv cut and get jet number
    btagCut = CSVMwp;
    btagLooseCut = LooseCSVMwp;
    nJets = int(jets.size());
    
    // get number of loose btags and btags
    nBtags = 0;
    nLooseBtags = 0;
    nOtherJets = 0;
    std::vector<int> BJet_order;
    BJet_order.resize(nJets);
    std::vector<int> LooseBJet_order;
    LooseBJet_order.resize(nJets);
    std::vector<int> OtherJet_order;
    OtherJet_order.resize(nJets);
    for (int i = 0; i < nJets; i++){
        if (csvs[i] > btagCut) {
            nBtags++;
            BJet_order[nBtags - 1] = i;
        }else if (csvs[i] > btagLooseCut && csvs[i] <= btagCut) {
            nLooseBtags++;
            LooseBJet_order[nLooseBtags - 1] = i;
        }else{
            nOtherJets++;
            OtherJet_order[nOtherJets - 1] = i;
        }
    }

    
    best_BosonOnly1.SetPxPyPzE(0.,0.,0.,0.); //vars for boson only reconstruction
    best_BosonOnly2.SetPxPyPzE(0.,0.,0.,0.); //vars for boson only reconstruction
    best_B1Only.SetPxPyPzE(0.,0.,0.,0.);
    best_B2Only.SetPxPyPzE(0.,0.,0.,0.);
    best_B3Only.SetPxPyPzE(0.,0.,0.,0.);
    best_B4Only.SetPxPyPzE(0.,0.,0.,0.);
    
    //init LVector that is to be used
    TLorentzVector boson1;
    TLorentzVector boson2;
    TLorentzVector current_boson1;
    TLorentzVector current_boson2;
    TLorentzVector B1;
    TLorentzVector B2;
    TLorentzVector B3;
    TLorentzVector B4;
    
    double current_chi2_Boson = 1e5;
    double chi2_Boson = 1e5;
    //first loop for finding the boson
    if (nBtags > 3){
        for(int i=0; i<nBtags; i++) {
        //second loop
            for(int j=0; j<nBtags; j++) {
                if(i==j) continue;
                // third loop
                for(int k=0; k<nBtags; k++) {
                    if(k==j) continue;
                    if(k==i) continue;
                    //forth loop
                    for(int l=0; l<nBtags; l++) {
                        if(l==k) continue;
                        if(l==j) continue;
                        if(l==i) continue;
                        //add jets to obatin boson
                        current_boson1 = jets.at(BJet_order[i])+jets.at(BJet_order[j]);
                        current_boson2 = jets.at(BJet_order[k])+jets.at(BJet_order[l]);
            
                        //calculate corresponding chi2
                        current_chi2_Boson = pow((current_boson1.M()-templateMass_Boson)/templateSigma_Boson, 2) + pow((current_boson2.M()-templateMass_Boson)/templateSigma_Boson, 2);
            
                        // save best solution
                        if(current_chi2_Boson < chi2_Boson)
                        {
                            chi2_Boson = current_chi2_Boson;
                            boson1 = current_boson1;
                            boson2 = current_boson2;
                            B1 = jets.at(BJet_order[i]);
                            B2 = jets.at(BJet_order[j]);
                            B3 = jets.at(BJet_order[k]);
                            B4 = jets.at(BJet_order[l]);
                        }
                    }
                }
            }
        }
       
        best_BosonOnly1 = boson1;
        best_BosonOnly2 = boson2;
        best_B1Only = B1;
        best_B2Only = B2;
        best_B3Only = B3;
        best_B4Only = B4;
        
        best_chi2_BosonOnly = chi2_Boson;
    }else if (nBtags == 3 && nLooseBtags > 0){
        for(int i=0; i<nBtags; i++) {
        //second loop
            for(int j=0; j<nBtags; j++) {
                if(i==j) continue;
                // third loop
                for(int k=0; k<nBtags; k++) {
                    if(k==j) continue;
                    if(k==i) continue;
                    //forth loop
                    for(int l=0; l<nLooseBtags; l++) {
                        //add jets to obatin boson
                        current_boson1 = jets.at(BJet_order[i])+jets.at(BJet_order[j]);
                        current_boson2 = jets.at(BJet_order[k])+jets.at(LooseBJet_order[l]);
            
                        //calculate corresponding chi2
                        current_chi2_Boson = pow((current_boson1.M()-templateMass_Boson)/templateSigma_Boson, 2) + pow((current_boson2.M()-templateMass_Boson)/templateSigma_Boson, 2);
            
                        // save best solution
                        if(current_chi2_Boson < chi2_Boson)
                        {
                            chi2_Boson = current_chi2_Boson;
                            boson1 = current_boson1;
                            boson2 = current_boson2;
                            B1 = jets.at(BJet_order[i]);
                            B2 = jets.at(BJet_order[j]);
                            B3 = jets.at(BJet_order[k]);
                            B4 = jets.at(LooseBJet_order[l]);
                        }
                    }
                }
            }
        }
        
        best_BosonOnly1 = boson1;
        best_BosonOnly2 = boson2;
        best_B1Only = B1;
        best_B2Only = B2;
        best_B3Only = B3;
        best_B4Only = B4;
    
        best_chi2_BosonOnly = chi2_Boson;
        
    }else if (nBtags == 3 && nLooseBtags == 0 && nOtherJets > 0){
        for(int i=0; i<nBtags; i++) {
        //second loop
            for(int j=0; j<nBtags; j++) {
                if(i==j) continue;
                // third loop
                for(int k=0; k<nBtags; k++) {
                    if(k==j) continue;
                    if(k==i) continue;
                    //forth loop
                    for(int l=0; l<nOtherJets; l++) {
                        //add jets to obatin boson
                        current_boson1 = jets.at(BJet_order[i])+jets.at(BJet_order[j]);
                        current_boson2 = jets.at(BJet_order[k])+jets.at(OtherJet_order[l]);
            
                        //calculate corresponding chi2
                        current_chi2_Boson = pow((current_boson1.M()-templateMass_Boson)/templateSigma_Boson, 2) + pow((current_boson2.M()-templateMass_Boson)/templateSigma_Boson, 2);
            
                        // save best solution
                        if(current_chi2_Boson < chi2_Boson)
                        {
                            chi2_Boson = current_chi2_Boson;
                            boson1 = current_boson1;
                            boson2 = current_boson2;
                            B1 = jets.at(BJet_order[i]);
                            B2 = jets.at(BJet_order[j]);
                            B3 = jets.at(BJet_order[k]);
                            B4 = jets.at(OtherJet_order[l]);
                        }
                    }
                }
            }
        }
        
        best_BosonOnly1 = boson1;
        best_BosonOnly2 = boson2;
        best_B1Only = B1;
        best_B2Only = B2;
        best_B3Only = B3;
        best_B4Only = B4;
    
        best_chi2_BosonOnly = chi2_Boson;
    }
}

void Chi2Reconstruction::ReconstructZH(const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, double CSVMwp, double LooseCSVMwp)
{
    
    // chi2 default value (class vars)
    //chi = 1e6;
    //minChi = 1e5;
    
    //set csv cut and get jet number
    btagCut = CSVMwp;
    btagLooseCut = LooseCSVMwp;
    nJets = int(jets.size());
    
    // get number of loose btags and btags
    nBtags = 0;
    nLooseBtags = 0;
    nOtherJets = 0;
    std::vector<int> BJet_order;
    BJet_order.resize(nJets);
    std::vector<int> LooseBJet_order;
    LooseBJet_order.resize(nJets);
    std::vector<int> OtherJet_order;
    OtherJet_order.resize(nJets);
    for (int i = 0; i < nJets; i++){
        if (csvs[i] > btagCut) {
            nBtags++;
            BJet_order[nBtags - 1] = i;
        }else if (csvs[i] > btagLooseCut && csvs[i] <= btagCut) {
            nLooseBtags++;
            LooseBJet_order[nLooseBtags - 1] = i;
        }else{
            nOtherJets++;
            OtherJet_order[nOtherJets - 1] = i;
        }
    }

    
    best_BosonOnly1.SetPxPyPzE(0.,0.,0.,0.); //vars for boson only reconstruction
    best_BosonOnly2.SetPxPyPzE(0.,0.,0.,0.); //vars for boson only reconstruction
    best_B1Only.SetPxPyPzE(0.,0.,0.,0.);
    best_B2Only.SetPxPyPzE(0.,0.,0.,0.);
    best_B3Only.SetPxPyPzE(0.,0.,0.,0.);
    best_B4Only.SetPxPyPzE(0.,0.,0.,0.);
    
    //init LVector that is to be used
    TLorentzVector boson1;
    TLorentzVector boson2;
    TLorentzVector current_boson1;
    TLorentzVector current_boson2;
    TLorentzVector B1;
    TLorentzVector B2;
    TLorentzVector B3;
    TLorentzVector B4;
    
    double current_chi2_Boson = 1e5;
    double chi2_Boson = 1e5;
    //first loop for finding the boson
    if (nBtags > 3){
        for(int i=0; i<nBtags; i++) {
        //second loop
            for(int j=0; j<nBtags; j++) {
                if(i==j) continue;
                // third loop
                for(int k=0; k<nBtags; k++) {
                    if(k==j) continue;
                    if(k==i) continue;
                    //forth loop
                    for(int l=0; l<nBtags; l++) {
                        if(l==k) continue;
                        if(l==j) continue;
                        if(l==i) continue;
                        //add jets to obatin boson
                        current_boson1 = jets.at(BJet_order[i])+jets.at(BJet_order[j]);
                        current_boson2 = jets.at(BJet_order[k])+jets.at(BJet_order[l]);
            
                        //calculate corresponding chi2
                        current_chi2_Boson = pow((current_boson1.M()-templateMass_Z)/templateSigma_Z, 2) + pow((current_boson2.M()-templateMass_Boson)/templateSigma_Boson, 2);
            
                        // save best solution
                        if(current_chi2_Boson < chi2_Boson)
                        {
                            chi2_Boson = current_chi2_Boson;
                            boson1 = current_boson1;
                            boson2 = current_boson2;
                            B1 = jets.at(BJet_order[i]);
                            B2 = jets.at(BJet_order[j]);
                            B3 = jets.at(BJet_order[k]);
                            B4 = jets.at(BJet_order[l]);
                        }
                    }
                }
            }
        }
       
        best_BosonOnly1 = boson1;
        best_BosonOnly2 = boson2;
        best_B1Only = B1;
        best_B2Only = B2;
        best_B3Only = B3;
        best_B4Only = B4;
        
        best_chi2_BosonOnly = chi2_Boson;
    }else if (nBtags == 3 && nLooseBtags > 0){
        for(int i=0; i<nBtags; i++) {
        //second loop
            for(int j=0; j<nBtags; j++) {
                if(i==j) continue;
                // third loop
                for(int k=0; k<nBtags; k++) {
                    if(k==j) continue;
                    if(k==i) continue;
                    //forth loop
                    for(int l=0; l<nLooseBtags; l++) {
                        //add jets to obatin boson
                        current_boson1 = jets.at(BJet_order[i])+jets.at(BJet_order[j]);
                        current_boson2 = jets.at(BJet_order[k])+jets.at(LooseBJet_order[l]);
                      
            
                        //calculate corresponding chi2
                        current_chi2_Boson = pow((current_boson1.M()-templateMass_Z)/templateSigma_Z, 2) + pow((current_boson2.M()-templateMass_Boson)/templateSigma_Boson, 2);
            
                        // save best solution
                        if(current_chi2_Boson < chi2_Boson)
                        {
                            chi2_Boson = current_chi2_Boson;
                            boson1 = current_boson1;
                            boson2 = current_boson2;
                            B1 = jets.at(BJet_order[i]);
                            B2 = jets.at(BJet_order[j]);
                            B3 = jets.at(BJet_order[k]);
                            B4 = jets.at(LooseBJet_order[l]);
                        }
                    }
                }
            }
        }
        
        best_BosonOnly1 = boson1;
        best_BosonOnly2 = boson2;
        best_B1Only = B1;
        best_B2Only = B2;
        best_B3Only = B3;
        best_B4Only = B4;
    
        best_chi2_BosonOnly = chi2_Boson;
        
    }else if (nBtags == 3 && nLooseBtags == 0 && nOtherJets > 0){
        for(int i=0; i<nBtags; i++) {
        //second loop
            for(int j=0; j<nBtags; j++) {
                if(i==j) continue;
                // third loop
                for(int k=0; k<nBtags; k++) {
                    if(k==j) continue;
                    if(k==i) continue;
                    //forth loop
                    for(int l=0; l<nOtherJets; l++) {
                        //add jets to obatin boson
                        current_boson1 = jets.at(BJet_order[i])+jets.at(BJet_order[j]);
                        current_boson2 = jets.at(BJet_order[k])+jets.at(OtherJet_order[l]);
                       
                        //calculate corresponding chi2
                        current_chi2_Boson = pow((current_boson1.M()-templateMass_Z)/templateSigma_Z, 2) + pow((current_boson2.M()-templateMass_Boson)/templateSigma_Boson, 2);
            
                        // save best solution
                        if(current_chi2_Boson < chi2_Boson)
                        {
                            chi2_Boson = current_chi2_Boson;
                            boson1 = current_boson1;
                            boson2 = current_boson2;
                            B1 = jets.at(BJet_order[i]);
                            B2 = jets.at(BJet_order[j]);
                            B3 = jets.at(BJet_order[k]);
                            B4 = jets.at(OtherJet_order[l]);
                        }
                    }
                }
            }
        }
        
        best_BosonOnly1 = boson1;
        best_BosonOnly2 = boson2;
        best_B1Only = B1;
        best_B2Only = B2;
        best_B3Only = B3;
        best_B4Only = B4;
    
        best_chi2_BosonOnly = chi2_Boson;
    }
}

void Chi2Reconstruction::ReconstructTTXSystem(const TLorentzVector& lepton, const std::vector<TLorentzVector>& jets, const std::vector<double>& csvs, TLorentzVector met, double CSVMwp)
{    
   
    // initializers
    //std::cout << "starting chi2 reco ========================================" << std::endl;
    InitReconstruction(jets, csvs, CSVMwp);
    ReconstructNeutrino(lepton, met);
    TLorentzVector metNew;

    int mistags = 0;
    int mistagged_topLepB = 0;
    int mistagged_topHadB = 0;
    int mistagged_topHadQ1 = 0;
    int mistagged_topHadQ2 = 0;
    int mistagged_BosonB1 = 0;
    int mistagged_BosonB2 = 0;

    double chi2_topLep = 10000;
    double chi2_topHad = 10000;
    double chi2_wHad   = 10000;
    double chi2_Boson  = 10000;

    // neurino solutions
    //std::cout << "looping over neutrino solutions" << std::endl;
    for( int ipznu=0; ipznu<2; ipznu++ )
    {
        // get met vector for neutrino solution
        metNew.SetXYZM(met.Px(),met.Py(),metPz[ipznu],0.0);

        // first make possible leptonic tops with b-tagged jets (includes making wLep)
        for (int i = 0; i < nJets; i++)
        {
            // leptonic W
            TLorentzVector wLep = metNew+lepton;
            // leptonic Top
            TLorentzVector topLep = metNew+lepton+jets.at(i);

            // chi2 for leptonic top
            chi2_topLep = pow((topLep.M()-templateMass_topLep)/templateSigma_topLep, 2);
            mistagged_topLepB = 0;
            if (csvs.at(i) < btagCut) mistagged_topLepB = 1;

            // next make possible wHads
            for (int j = 0; j < nJets; j++)
            {
                if ( j == i ) continue;
                for (int k = 0; k < nJets; k++)
                {
                    if( k == i || k == j ) continue;
                    
                    // hadronic W
                    TLorentzVector wHad = jets.at(j)+jets.at(k);
                    
                    // chi2 for hadronic W
                    chi2_wHad = pow((wHad.M()-templateMass_W)/templateSigma_W, 2);
                    mistagged_topHadQ1 = 0;                                   
                    mistagged_topHadQ2 = 0;
                    if (csvs.at(j) > btagCut) mistagged_topHadQ1 = 1;
                    if (csvs.at(k) > btagCut) mistagged_topHadQ2 = 1;

                    // next make possible hadronic top
                    for (int l = 0; l < nJets; l++)
                    {
                        if( l == i || l == j || l == k ) continue;
                        
                        // hadronic top
                        TLorentzVector topHad = wHad+jets.at(l);

                        // chi2 for hadronic top
                        chi2_topHad = pow((topHad.M()-templateMass_topHad)/templateSigma_topHad, 2);
                        mistagged_topHadB = 0;
                        if (csvs.at(l) < btagCut) mistagged_topHadB = 1;

                        // count mistags so far
                        mistags = mistagged_topHadB+mistagged_topLepB+mistagged_topHadQ1+mistagged_topHadQ2;
                        if (mistags > allowedMistags) continue;

                        bool bosonReconstructed = false;
                        TLorentzVector boson;
                        TLorentzVector B1;
                        TLorentzVector B2;
                        boson.SetPxPyPzE(0.,0.,0.,0.);
                        B1.SetPxPyPzE(0.,0.,0.,0.);
                        B2.SetPxPyPzE(0.,0.,0.,0.);

                        if (reconstructBoson)
                        {
                            chi2_Boson = 1e6;
                            double current_chi2_Boson = 1e5;
                            double current_mistags;
                            TLorentzVector current_boson;

                            // next make possible higgs bosons
                            for (int m = 0; m < nJets; m++) 
                            {
                                if( m == i || m == j || m == k || m == l ) continue;
                                for (int n = 0; n < nJets; n++)
                                {
                                    if( n == i || n == j || n == k || n == l || n == m ) continue;

                                    // Higgs
                                    current_boson = jets.at(m)+jets.at(n);

                                    // chi2 for higgs
                                    current_chi2_Boson = pow((current_boson.M()-templateMass_Boson)/templateSigma_Boson, 2);
                                    mistagged_BosonB1 = 0;
                                    mistagged_BosonB2 = 0;
                                    if (csvs.at(m) < btagCut) mistagged_BosonB1 = 1;
                                    if (csvs.at(n) < btagCut) mistagged_BosonB2 = 1;
                                    
                                    // count mistags
                                    current_mistags = mistags+mistagged_BosonB1+mistagged_BosonB2;
                                    if (current_mistags > allowedMistags) continue;

                                    // save best solution
                                    if(current_chi2_Boson < chi2_Boson)
                                    {
                                        bosonReconstructed = true;
                                        mistags = current_mistags;
                                        chi2_Boson = current_chi2_Boson;
                                        boson = current_boson;
                                        B1 = jets.at(m);
                                        B2 = jets.at(n);
                                    }
                                    
                                }
                            }
                        }
                        
                        // continue if no higgs could be reconstructed
                        if (reconstructBoson && !bosonReconstructed) continue;

                        // combine chi2
                        chi = chi2_topLep + chi2_wHad + chi2_topHad;
                        if (reconstructBoson) chi+=chi2_Boson;

                        // accept lowest chi
                        if( chi < minChi )
                        {
                            minChi = chi;
                            best_chi2_topHad = chi2_topHad;
                            best_chi2_topLep = chi2_topLep;
                            best_chi2_wHad = chi2_wHad;
                            best_chi2_Boson = chi2_Boson;

                            best_wLep = wLep;
                            best_wHad = wHad;                      
                            best_bLep = jets.at(i);
                            best_bHad = jets.at(l);
                            best_topLep = topLep;
                            best_topHad = topHad;
                            best_lepton = lepton;
                            best_mET = metNew;
                            best_Boson = boson;
                            best_B1 = B1;
                            best_B2 = B2;
                        }
                    }   
                }
            }
        }

        
    }
    //if( nJets >= 6 )
    //{
    //std::cout << "================================" << std::endl;
    //std::cout << "event printout for >=6jet events" << std::endl;
    //std::cout << "allowed mistags " << allowedMistags   << std::endl;
    //std::cout << "minChi          " << minChi           << std::endl;
    //std::cout << "bestChi_topHad  " << best_chi2_topHad << std::endl;
    //std::cout << "bestChi_topLep  " << best_chi2_topLep << std::endl;
    //std::cout << "bestChi_wHad    " << best_chi2_wHad   << std::endl;
    //std::cout << "bestChi_Boson   " << best_chi2_Boson  << std::endl;
    //std::cout << "hadtop pt       " << best_topHad.Pt() << std::endl;
    //std::cout << "higgs pt        " << best_Boson.Pt()  << std::endl;
    //std::cout << "================================" << std::endl;
    //}
    
}



///////////////////////// VARIABLE MAP INITIALIZATIONS ////////////////////////////////
std::map<std::string, double> Chi2Reconstruction::GetVariableMap_BosonSystemAngles() {
    std::map<std::string, double> varMap;
    
    // Deta
    varMap[tag+"Deta_"+bosonName+"_Lep"]            = -999.;
    varMap[tag+"Deta_"+bosonName+"_bLep"]           = -999.;
    varMap[tag+"Deta_"+bosonName+"_bHad"]           = -999.;
    varMap[tag+"Deta_"+bosonName+"_topHad"]         = -999.;
    varMap[tag+"Deta_"+bosonName+"_topLep"]         = -999.;
    
    // Dphi
    varMap[tag+"Dphi_"+bosonName+"_Lep"]            = -999.;
    varMap[tag+"Dphi_"+bosonName+"_bLep"]           = -999.;
    varMap[tag+"Dphi_"+bosonName+"_bHad"]           = -999.;
    varMap[tag+"Dphi_"+bosonName+"_topHad"]         = -999.;
    varMap[tag+"Dphi_"+bosonName+"_topLep"]         = -999.;
    
    // cosdtheta
    varMap[tag+"cosdTheta_"+bosonName+"_Lep"]      = -999.;
    varMap[tag+"cosdTheta_"+bosonName+"_bLep"]     = -999.;
    varMap[tag+"cosdTheta_"+bosonName+"_bHad"]     = -999.;
    varMap[tag+"cosdTheta_"+bosonName+"_topHad"]   = -999.;
    varMap[tag+"cosdTheta_"+bosonName+"_topLep"]   = -999.;

    return varMap;    
    }
    

std::map<std::string, double> Chi2Reconstruction::GetVariableMap_BosonOnlyAngles() {
    std::map<std::string, double> varMap;
    
    // Deta
    varMap[tag+"1_Deta"]            = -999.;
    varMap[tag+"2_Deta"]            = -999.;

    
    // Dphi
    varMap[tag+"1_Dphi"]            = -999.;
    varMap[tag+"2_Dphi"]            = -999.;
    varMap[tag+"1_cosdTheta"]     =-999.;
    varMap[tag+"2_cosdTheta"]     =-999.;

    
    varMap[tag+"1_Dr"]               = -999.;
    varMap[tag+"2_Dr"]               = -999.;

    
    return varMap;    
    }

std::map<std::string, double> Chi2Reconstruction::GetVariableMap_TopSystemAngles() {
    std::map<std::string, double> varMap;

    // Deta
    //varMap[tag+"Deta_Lep_bLep"]      = -999.;
    varMap[tag+"Deta_Lep_bHad"]      = -999.;
    varMap[tag+"Deta_Lep_topHad"]    = -999.;
    //varMap[tag+"Deta_Lep_wHad"]      = -999.;
    varMap[tag+"Deta_bLep_bHad"]     = -999.;
    varMap[tag+"Deta_bLep_topHad"]   = -999.;
    //varMap[tag+"Deta_bLep_wHad"]     = -999.;
    varMap[tag+"Deta_bLep_wLep"]     = -999.;
    varMap[tag+"Deta_topLep_bHad"]   = -999.;
    varMap[tag+"Deta_topLep_topHad"] = -999.;
    //varMap[tag+"Deta_topLep_wHad"]   = -999.;
    //varMap[tag+"Deta_bHad_wLep"]     = -999.;
    varMap[tag+"Deta_bHad_wHad"]     = -999.;
    varMap[tag+"Deta_topHad_wLep"]   = -999.;
    //varMap[tag+"Deta_wHad_wLep"]     = -999.;

    // Dphi
    //varMap[tag+"Dphi_Lep_bLep"]      = -999.;
    varMap[tag+"Dphi_Lep_bHad"]      = -999.;
    varMap[tag+"Dphi_Lep_topHad"]    = -999.;
    //varMap[tag+"Dphi_Lep_wHad"]      = -999.;
    varMap[tag+"Dphi_bLep_bHad"]     = -999.;
    varMap[tag+"Dphi_bLep_topHad"]   = -999.;
    //varMap[tag+"Dphi_bLep_wHad"]     = -999.;
    varMap[tag+"Dphi_bLep_wLep"]     = -999.;
    varMap[tag+"Dphi_topLep_bHad"]   = -999.;
    varMap[tag+"Dphi_topLep_topHad"] = -999.;
    //varMap[tag+"Dphi_topLep_wHad"]   = -999.;
    //varMap[tag+"Dphi_bHad_wLep"]     = -999.;
    varMap[tag+"Dphi_bHad_wHad"]     = -999.;
    varMap[tag+"Dphi_topHad_wLep"]   = -999.;
    //varMap[tag+"Dphi_wHad_wLep"]     = -999.;

    // cos dtheta
    //varMap[tag+"cosdTheta_Lep_bLep"]      = -999.;
    varMap[tag+"cosdTheta_Lep_bHad"]      = -999.;
    varMap[tag+"cosdTheta_Lep_topHad"]    = -999.;
    //varMap[tag+"cosdTheta_Lep_wHad"]      = -999.;
    varMap[tag+"cosdTheta_bLep_bHad"]     = -999.;
    varMap[tag+"cosdTheta_bLep_topHad"]   = -999.;
    //varMap[tag+"cosdTheta_bLep_wHad"]     = -999.;
    varMap[tag+"cosdTheta_bLep_wLep"]     = -999.;
    varMap[tag+"cosdTheta_topLep_bHad"]   = -999.;
    varMap[tag+"cosdTheta_topLep_topHad"] = -999.;
    //varMap[tag+"cosdTheta_topLep_wHad"]   = -999.;
    //varMap[tag+"cosdTheta_bHad_wLep"]     = -999.;
    varMap[tag+"cosdTheta_bHad_wHad"]     = -999.;
    varMap[tag+"cosdTheta_topHad_wLep"]   = -999.;
    //varMap[tag+"cosdTheta_wHad_wLep"]     = -999.;

    return varMap;
    }

std::map<std::string, double> Chi2Reconstruction::GetVariableMap_BosonSystemObjects() {
    std::map<std::string, double> varMap;

    // higgs
    varMap[tag+bosonName+"_M"]   = -999.;
    varMap[tag+bosonName+"_E"]   = -999.;
    varMap[tag+bosonName+"_Eta"] = -999.;
    varMap[tag+bosonName+"_Phi"] = -999.;
    varMap[tag+bosonName+"_Pt"]  = -999.;

    // b jets from higgs decay
    varMap[tag+bosonName+"_BJet1_E"] =   -999.;
    varMap[tag+bosonName+"_BJet1_Eta"] = -999.;
    varMap[tag+bosonName+"_BJet1_Phi"] = -999.;
    varMap[tag+bosonName+"_BJet1_Pt"] =  -999.;

    varMap[tag+bosonName+"_BJet2_E"] =   -999.;
    varMap[tag+bosonName+"_BJet2_Eta"] = -999.;
    varMap[tag+bosonName+"_BJet2_Phi"] = -999.;
    varMap[tag+bosonName+"_BJet2_Pt"] =  -999.;

    // chi2
    varMap[tag+"Chi2"+bosonName] = -999.;
    
    return varMap;
}

std::map<std::string, double> Chi2Reconstruction::GetVariableMap_BosonOnlyObjects() {
    std::map<std::string, double> varMap;
    
    
    // higgs
    varMap[tag+"1_M"]   = -999.;
    varMap[tag+"1_E"]   = -999.;
    varMap[tag+"1_Eta"] = -999.;
    varMap[tag+"1_Phi"] = -999.;
    varMap[tag+"1_Pt"]  = -999.;
    
    varMap[tag+"2_M"]   = -999.;
    varMap[tag+"2_E"]   = -999.;
    varMap[tag+"2_Eta"] = -999.;
    varMap[tag+"2_Phi"] = -999.;
    varMap[tag+"2_Pt"]  = -999.;

    // b jets from higgs decay
    varMap[tag+"BJet1_E"] =   -999.;
    varMap[tag+"BJet1_M"] = -999.;
    varMap[tag+"BJet1_Eta"] = -999.;
    varMap[tag+"BJet1_Phi"] = -999.;
    varMap[tag+"BJet1_Pt"] =  -999.;

    varMap[tag+"BJet2_E"] =   -999.;
    varMap[tag+"BJet2_M"] = -999.;
    varMap[tag+"BJet2_Eta"] = -999.;
    varMap[tag+"BJet2_Phi"] = -999.;
    varMap[tag+"BJet2_Pt"] =  -999.;
    
    varMap[tag+"BJet3_E"] =   -999.;
    varMap[tag+"BJet3_M"] = -999.;
    varMap[tag+"BJet3_Eta"] = -999.;
    varMap[tag+"BJet3_Phi"] = -999.;
    varMap[tag+"BJet3_Pt"] =  -999.;

    varMap[tag+"BJet4_E"] =   -999.;
    varMap[tag+"BJet4_M"] = -999.;
    varMap[tag+"BJet4_Eta"] = -999.;
    varMap[tag+"BJet4_Phi"] = -999.;
    varMap[tag+"BJet4_Pt"] =  -999.;

    // chi2
    varMap[tag+"Chi2"] = -999.;
    varMap[tag+"logChi2"] = -999;
    
    return varMap;
}

std::map<std::string, double> Chi2Reconstruction::GetVariableMap_TopSystemObjects() {
    std::map<std::string, double> varMap;

    // hadronic top system
    //varMap[tag+"TopHad_BJet_M"] = -999.;
    varMap[tag+"TopHad_BJet_E"] = -999.;
    varMap[tag+"TopHad_BJet_Eta"] = -999.;
    varMap[tag+"TopHad_BJet_Phi"] = -999.;
    varMap[tag+"TopHad_BJet_Pt"] = -999.;

    //varMap[tag+"TopHad_W_M"] = -999.;
    varMap[tag+"TopHad_W_E"] = -999.;
    varMap[tag+"TopHad_W_Eta"] = -999.;
    varMap[tag+"TopHad_W_Phi"] = -999.;
    varMap[tag+"TopHad_W_Pt"] = -999.;

    varMap[tag+"TopHad_M"] = -999.;
    varMap[tag+"TopHad_E"] = -999.;
    varMap[tag+"TopHad_Eta"] = -999.;
    varMap[tag+"TopHad_Phi"] = -999.;
    varMap[tag+"TopHad_Pt"] = -999.;


    // leptonic top system
    //varMap[tag+"TopLep_BJet_M"] = -999.;
    varMap[tag+"TopLep_BJet_E"] = -999.;
    varMap[tag+"TopLep_BJet_Eta"] = -999.;
    varMap[tag+"TopLep_BJet_Phi"] = -999.;
    varMap[tag+"TopLep_BJet_Pt"] = -999.;

    //varMap[tag+"TopLep_W_M"] = -999.;
    varMap[tag+"TopLep_W_E"] = -999.;
    varMap[tag+"TopLep_W_Eta"] = -999.;
    varMap[tag+"TopLep_W_Phi"] = -999.;
    varMap[tag+"TopLep_W_Pt"] = -999.;

    varMap[tag+"TopLep_M"] = -999.;
    varMap[tag+"TopLep_E"] = -999.;
    varMap[tag+"TopLep_Eta"] = -999.;
    varMap[tag+"TopLep_Phi"] = -999.;
    varMap[tag+"TopLep_Pt"] = -999.;

    // chi2 values
    varMap[tag+"Chi2Total"] = -999.;
    varMap[tag+"Chi2TopHad"] = -999.;
    varMap[tag+"Chi2TopLep"] = -999.;
    varMap[tag+"Chi2WHad"] = -999.;

    return varMap;
    }



double Chi2Reconstruction::GetDeltaPhi(double phi1, double phi2) {
    double dphi = TMath::Abs( phi1-phi2 );
    if(dphi >= TMath::Pi())
        dphi = 2.*TMath::Pi() - dphi;
    
    return dphi;
    }

///////////////////////// VARIABLE MAP FILLERS ////////////////////////////////
std::map<std::string,double> Chi2Reconstruction::FillVariableMap_BosonSystemAngles() {
    std::map<std::string, double> varMap;

    // Deta
    varMap[tag+"Deta_"+bosonName+"_Lep"]            = TMath::Abs( best_Boson.Eta() - best_lepton.Eta() );
    varMap[tag+"Deta_"+bosonName+"_bLep"]           = TMath::Abs( best_Boson.Eta() - best_bLep.Eta() );
    varMap[tag+"Deta_"+bosonName+"_bHad"]           = TMath::Abs( best_Boson.Eta() - best_bHad.Eta() );
    varMap[tag+"Deta_"+bosonName+"_topHad"]         = TMath::Abs( best_Boson.Eta() - best_topHad.Eta() );
    varMap[tag+"Deta_"+bosonName+"_topLep"]         = TMath::Abs( best_Boson.Eta() - best_topLep.Eta() );
    
    // Dphi
    varMap[tag+"Dphi_"+bosonName+"_Lep"]            = GetDeltaPhi( best_Boson.Phi(), best_lepton.Phi() );
    varMap[tag+"Dphi_"+bosonName+"_bLep"]           = GetDeltaPhi( best_Boson.Phi(), best_bLep.Phi() );
    varMap[tag+"Dphi_"+bosonName+"_bHad"]           = GetDeltaPhi( best_Boson.Phi(), best_bHad.Phi() );
    varMap[tag+"Dphi_"+bosonName+"_topHad"]         = GetDeltaPhi( best_Boson.Phi(), best_topHad.Phi() );
    varMap[tag+"Dphi_"+bosonName+"_topLep"]         = GetDeltaPhi( best_Boson.Phi(), best_topLep.Phi() );
    
    // cosdtheta
    varMap[tag+"cosdTheta_"+bosonName+"_Lep"]      = TMath::Cos( (best_Boson.Vect()).Angle(best_lepton.Vect()) );
    varMap[tag+"cosdTheta_"+bosonName+"_bLep"]     = TMath::Cos( (best_Boson.Vect()).Angle(best_bLep.Vect()) );
    varMap[tag+"cosdTheta_"+bosonName+"_bHad"]     = TMath::Cos( (best_Boson.Vect()).Angle(best_bHad.Vect()) );
    varMap[tag+"cosdTheta_"+bosonName+"_topHad"]   = TMath::Cos( (best_Boson.Vect()).Angle(best_topHad.Vect()) );
    varMap[tag+"cosdTheta_"+bosonName+"_topLep"]   = TMath::Cos( (best_Boson.Vect()).Angle(best_topLep.Vect()) );


    return varMap;
}

std::map<std::string,double> Chi2Reconstruction::FillVariableMap_TopSystemAngles() {
    std::map<std::string, double> varMap;

    // Deta
    //varMap[tag+"Deta_Lep_bLep"]      = TMath::Abs( best_lepton.Eta() - best_bLep.Eta() );
    varMap[tag+"Deta_Lep_bHad"]      = TMath::Abs( best_lepton.Eta() - best_bHad.Eta() );
    varMap[tag+"Deta_Lep_topHad"]    = TMath::Abs( best_lepton.Eta() - best_topHad.Eta() );
    //varMap[tag+"Deta_Lep_wHad"]      = TMath::Abs( best_lepton.Eta() - best_wHad.Eta() );
    varMap[tag+"Deta_bLep_bHad"]     = TMath::Abs( best_bLep.Eta() - best_bHad.Eta() );
    varMap[tag+"Deta_bLep_topHad"]   = TMath::Abs( best_bLep.Eta() - best_topHad.Eta() );
    //varMap[tag+"Deta_bLep_wHad"]     = TMath::Abs( best_bLep.Eta() - best_wHad.Eta() );
    varMap[tag+"Deta_bLep_wLep"]     = TMath::Abs( best_bLep.Eta() - best_wLep.Eta() );
    varMap[tag+"Deta_topLep_bHad"]   = TMath::Abs( best_topLep.Eta() - best_bHad.Eta() );
    varMap[tag+"Deta_topLep_topHad"] = TMath::Abs( best_topLep.Eta() - best_topHad.Eta() );
    //varMap[tag+"Deta_topLep_wHad"]   = TMath::Abs( best_topLep.Eta() - best_wHad.Eta() );
    //varMap[tag+"Deta_bHad_wLep"]     = TMath::Abs( best_bHad.Eta() - best_wLep.Eta() );
    varMap[tag+"Deta_bHad_wHad"]     = TMath::Abs( best_bHad.Eta() - best_wHad.Eta() );
    varMap[tag+"Deta_topHad_wLep"]   = TMath::Abs( best_topHad.Eta() - best_wLep.Eta() );
    //varMap[tag+"Deta_wHad_wLep"]     = TMath::Abs( best_wHad.Eta() - best_wLep.Eta() );

    // Dphi
    //varMap[tag+"Dphi_Lep_bLep"]      = GetDeltaPhi( best_lepton.Phi(), best_bLep.Phi() );
    varMap[tag+"Dphi_Lep_bHad"]      = GetDeltaPhi( best_lepton.Phi(), best_bHad.Phi() );
    varMap[tag+"Dphi_Lep_topHad"]    = GetDeltaPhi( best_lepton.Phi(), best_topHad.Phi() );
    //varMap[tag+"Dphi_Lep_wHad"]      = GetDeltaPhi( best_lepton.Phi(), best_wHad.Phi() );
    varMap[tag+"Dphi_bLep_bHad"]     = GetDeltaPhi( best_bLep.Phi(), best_bHad.Phi() );
    varMap[tag+"Dphi_bLep_topHad"]   = GetDeltaPhi( best_bLep.Phi(), best_topHad.Phi() );
    //varMap[tag+"Dphi_bLep_wHad"]     = GetDeltaPhi( best_bLep.Phi(), best_wHad.Phi() );
    varMap[tag+"Dphi_bLep_wLep"]     = GetDeltaPhi( best_bLep.Phi(), best_wLep.Phi() );
    varMap[tag+"Dphi_topLep_bHad"]   = GetDeltaPhi( best_topLep.Phi(), best_bHad.Phi() );
    varMap[tag+"Dphi_topLep_topHad"] = GetDeltaPhi( best_topLep.Phi(), best_topHad.Phi() );
    //varMap[tag+"Dphi_topLep_wHad"]   = GetDeltaPhi( best_topLep.Phi(), best_wHad.Phi() );
    //varMap[tag+"Dphi_bHad_wLep"]     = GetDeltaPhi( best_bHad.Phi(), best_wLep.Phi() );
    varMap[tag+"Dphi_bHad_wHad"]     = GetDeltaPhi( best_bHad.Phi(), best_wHad.Phi() );
    varMap[tag+"Dphi_topHad_wLep"]   = GetDeltaPhi( best_topHad.Phi(), best_wLep.Phi() );
    //varMap[tag+"Dphi_wHad_wLep"]     = GetDeltaPhi( best_wHad.Phi(), best_wLep.Phi() );

    // cos dtheta
    //varMap[tag+"cosdTheta_Lep_bLep"]      =  TMath::Cos( (best_lepton.Vect()).Angle(best_bLep.Vect()) );
    varMap[tag+"cosdTheta_Lep_bHad"]      =  TMath::Cos( (best_lepton.Vect()).Angle(best_bHad.Vect()) );
    varMap[tag+"cosdTheta_Lep_topHad"]    =  TMath::Cos( (best_lepton.Vect()).Angle(best_topHad.Vect()) );
    //varMap[tag+"cosdTheta_Lep_wHad"]      =  TMath::Cos( (best_lepton.Vect()).Angle(best_wHad.Vect()) );
    varMap[tag+"cosdTheta_bLep_bHad"]     =  TMath::Cos( (best_bLep.Vect()).Angle(best_bHad.Vect()) );
    varMap[tag+"cosdTheta_bLep_topHad"]   =  TMath::Cos( (best_bLep.Vect()).Angle(best_topHad.Vect()) );
    //varMap[tag+"cosdTheta_bLep_wHad"]     =  TMath::Cos( (best_bLep.Vect()).Angle(best_wHad.Vect()) );
    varMap[tag+"cosdTheta_bLep_wLep"]     =  TMath::Cos( (best_bLep.Vect()).Angle(best_wLep.Vect()) );
    varMap[tag+"cosdTheta_topLep_bHad"]   =  TMath::Cos( (best_topLep.Vect()).Angle(best_bHad.Vect()) );
    varMap[tag+"cosdTheta_topLep_topHad"] =  TMath::Cos( (best_topLep.Vect()).Angle(best_topHad.Vect()) );
    //varMap[tag+"cosdTheta_topLep_wHad"]   =  TMath::Cos( (best_topLep.Vect()).Angle(best_wHad.Vect()) );
    //varMap[tag+"cosdTheta_bHad_wLep"]     =  TMath::Cos( (best_bHad.Vect()).Angle(best_wLep.Vect()) );
    varMap[tag+"cosdTheta_bHad_wHad"]     =  TMath::Cos( (best_bHad.Vect()).Angle(best_wHad.Vect()) );
    varMap[tag+"cosdTheta_topHad_wLep"]   =  TMath::Cos( (best_topHad.Vect()).Angle(best_wLep.Vect()) );
    //varMap[tag+"cosdTheta_wHad_wLep"]     =  TMath::Cos( (best_wHad.Vect()).Angle(best_wLep.Vect()) );
    
    return varMap;
}

std::map<std::string,double> Chi2Reconstruction::FillVariableMap_BosonSystemObjects() {
    std::map<std::string, double> varMap;

    // higgs 
    varMap[tag+bosonName+"_M"]   = best_Boson.M();
    varMap[tag+bosonName+"_E"]   = best_Boson.E();
    varMap[tag+bosonName+"_Eta"] = best_Boson.Eta();
    varMap[tag+bosonName+"_Phi"] = best_Boson.Phi();
    varMap[tag+bosonName+"_Pt"]  = best_Boson.Pt();

    // b jets from higgs decay
    varMap[tag+bosonName+"_BJet1_E"] =   best_B1.E();
    varMap[tag+bosonName+"_BJet1_Eta"] = best_B1.Eta();
    varMap[tag+bosonName+"_BJet1_Phi"] = best_B1.Phi();
    varMap[tag+bosonName+"_BJet1_Pt"] =  best_B1.Pt();

    varMap[tag+bosonName+"_BJet2_E"] =   best_B2.E();
    varMap[tag+bosonName+"_BJet2_Eta"] = best_B2.Eta();
    varMap[tag+bosonName+"_BJet2_Phi"] = best_B2.Phi();
    varMap[tag+bosonName+"_BJet2_Pt"] =  best_B2.Pt();

    // chi2
    varMap[tag+"Chi2"+bosonName] = best_chi2_Boson;
    
    return varMap;
}




std::map<std::string,double> Chi2Reconstruction::FillVariableMap_TopSystemObjects() {
    std::map<std::string, double> varMap;

    // hadronic top system
    //varMap[tag+"TopHad_BJet_M"]   = best_bHad.M();
    varMap[tag+"TopHad_BJet_E"]   = best_bHad.E();
    varMap[tag+"TopHad_BJet_Eta"] = best_bHad.Eta();
    varMap[tag+"TopHad_BJet_Phi"] = best_bHad.Phi();
    varMap[tag+"TopHad_BJet_Pt"]  = best_bHad.Pt();

    //varMap[tag+"TopHad_W_M"]   = best_wHad.M();
    varMap[tag+"TopHad_W_E"]   = best_wHad.E();
    varMap[tag+"TopHad_W_Eta"] = best_wHad.Eta();
    varMap[tag+"TopHad_W_Phi"] = best_wHad.Phi();
    varMap[tag+"TopHad_W_Pt"]  = best_wHad.Pt();

    varMap[tag+"TopHad_M"]   = best_topHad.M();
    varMap[tag+"TopHad_E"]   = best_topHad.E();
    varMap[tag+"TopHad_Eta"] = best_topHad.Eta();
    varMap[tag+"TopHad_Phi"] = best_topHad.Phi();
    varMap[tag+"TopHad_Pt"]  = best_topHad.Pt();


    // leptonic top system
    //varMap[tag+"TopLep_BJet_M"]   = best_bLep.M();
    varMap[tag+"TopLep_BJet_E"]   = best_bLep.E();
    varMap[tag+"TopLep_BJet_Eta"] = best_bLep.Eta();
    varMap[tag+"TopLep_BJet_Phi"] = best_bLep.Phi();
    varMap[tag+"TopLep_BJet_Pt"]  = best_bLep.Pt();

    //varMap[tag+"TopLep_W_M"]   = best_wLep.M();
    varMap[tag+"TopLep_W_E"]   = best_wLep.E();
    varMap[tag+"TopLep_W_Eta"] = best_wLep.Eta();
    varMap[tag+"TopLep_W_Phi"] = best_wLep.Phi();
    varMap[tag+"TopLep_W_Pt"]  = best_wLep.Pt();

    varMap[tag+"TopLep_M"]   = best_topLep.M();
    varMap[tag+"TopLep_E"]   = best_topLep.E();
    varMap[tag+"TopLep_Eta"] = best_topLep.Eta();
    varMap[tag+"TopLep_Phi"] = best_topLep.Phi();
    varMap[tag+"TopLep_Pt"]  = best_topLep.Pt();

    // chi2 values
    varMap[tag+"Chi2Total"] = minChi;
    varMap[tag+"Chi2TopHad"] = best_chi2_topHad;
    varMap[tag+"Chi2TopLep"] = best_chi2_topLep;
    varMap[tag+"Chi2WHad"] = best_chi2_wHad;

    return varMap;
}

/*
std::map<std::string,double> Chi2Reconstruction::FillVariableMap_BosonOnlyAngles() {
    std::map<std::string, double> varMap;
    
    return varMap;
}*/

std::map<std::string,double> Chi2Reconstruction::FillVariableMap_BosonOnlyObjects() {
    std::map<std::string, double> varMap;
    
    // higgs 
    varMap[tag+"1_M"]   = best_BosonOnly1.M();
    varMap[tag+"1_E"]   = best_BosonOnly1.E();
    varMap[tag+"1_Eta"] = best_BosonOnly1.Eta();
    varMap[tag+"1_Phi"] = best_BosonOnly1.Phi();
    varMap[tag+"1_Pt"]  = best_BosonOnly1.Pt();
    
    varMap[tag+"2_M"]   = best_BosonOnly2.M();
    varMap[tag+"2_E"]   = best_BosonOnly2.E();
    varMap[tag+"2_Eta"] = best_BosonOnly2.Eta();
    varMap[tag+"2_Phi"] = best_BosonOnly2.Phi();
    varMap[tag+"2_Pt"]  = best_BosonOnly2.Pt();

    // b jets from higgs decay
    varMap[tag+"BJet1_E"] =   best_B1Only.E();
    varMap[tag+"BJet1_M"] =   best_B1Only.M();
    varMap[tag+"BJet1_Eta"] = best_B1Only.Eta();
    varMap[tag+"BJet1_Phi"] = best_B1Only.Phi();
    varMap[tag+"BJet1_Pt"] =  best_B1Only.Pt();

    varMap[tag+"BJet2_E"] =   best_B2Only.E();
    varMap[tag+"BJet2_M"] =   best_B2Only.M();
    varMap[tag+"BJet2_Eta"] = best_B2Only.Eta();
    varMap[tag+"BJet2_Phi"] = best_B2Only.Phi();
    varMap[tag+"BJet2_Pt"] =  best_B2Only.Pt();
    
    varMap[tag+"BJet3_E"] =   best_B3Only.E();
    varMap[tag+"BJet3_M"] =   best_B3Only.M();
    varMap[tag+"BJet3_Eta"] = best_B3Only.Eta();
    varMap[tag+"BJet3_Phi"] = best_B3Only.Phi();
    varMap[tag+"BJet3_Pt"] =  best_B3Only.Pt();
    
    varMap[tag+"BJet4_E"] =   best_B4Only.E();
    varMap[tag+"BJet4_M"] =   best_B4Only.M();
    varMap[tag+"BJet4_Eta"] = best_B4Only.Eta();
    varMap[tag+"BJet4_Phi"] = best_B4Only.Phi();
    varMap[tag+"BJet4_Pt"] =  best_B4Only.Pt();

    // chi2
    varMap[tag+"Chi2"] = best_chi2_BosonOnly;
    varMap[tag+"logChi2"] = TMath::Log(best_chi2_BosonOnly);
    
    return varMap;
}

std::map<std::string, double> Chi2Reconstruction::FillVariableMap_BosonOnlyAngles() {
    std::map<std::string, double> varMap;
    

    // Deta
    varMap[tag+"1_Deta"]            = fabs(best_B1Only.Eta() - best_B2Only.Eta());
    
    
    // Dphi
    varMap[tag+"1_Dphi"]            = GetDeltaPhi(best_B1Only.Phi(), best_B2Only.Phi());
    
    varMap[tag+"1_cosdTheta"]       = TMath::Cos( best_B1Only.Vect().Angle(best_B2Only.Vect()) );
    
    varMap[tag+"1_Dr"]              = best_B1Only.DeltaR(best_B2Only);
    

        varMap[tag+"2_Deta"]            = fabs(best_B3Only.Eta() - best_B4Only.Eta());
        varMap[tag+"2_Dphi"]            = GetDeltaPhi(best_B3Only.Phi(), best_B4Only.Phi());
        varMap[tag+"2_cosdTheta"]       = TMath::Cos( best_B3Only.Vect().Angle(best_B4Only.Vect()) );
        varMap[tag+"2_Dr"]              = best_B3Only.DeltaR(best_B4Only);

    
    return varMap;    
    }

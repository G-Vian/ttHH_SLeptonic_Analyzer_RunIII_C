#include "TTH/CommonClassifier/interface/MVAvars.h"

using namespace std;

MVAvars::MVAvars(const char* era)
{

    // mem = MEMClassifier(0, "btagDeepFlavB_", "2018");
    // mem = MEMClassifier();
    MEMClassifier mem(0, "btagDeepFlavB_", era);
    // ==================================================
    //init all variables potentially used
    
    // higgs specific variables
    variableMap["Reco_best_higgs_mass"] = -999.;
    variableMap["Reco_dEta_fn"] = -999.;


    // CSV variables
    variableMap["Evt_CSV_min"] = -999.;
    variableMap["Evt_CSV_min_tagged"] = -999.;
    variableMap["Evt_CSV_avg"] = -999.;
    variableMap["Evt_CSV_avg_tagged"] = -999.;
    variableMap["Evt_CSV_dev"] = -999.;
    variableMap["Evt_CSV_dev_tagged"] = -999.;

    // fox wolfram moments
    variableMap["Evt_h0"] = -999.;
    variableMap["Evt_h1"] = -999.;
    variableMap["Evt_h2"] = -999.;
    variableMap["Evt_h3"] = -999.;

    // event masses 
    variableMap["Evt_M_Total"] = -999.;
    variableMap["Evt_M3"] = -999.;
    variableMap["Evt_M3_oneTagged"] = -999.;

    // HT variables
    variableMap["Evt_HT"] = -999.;
    variableMap["Evt_HT_wo_MET"] = -999.;
    variableMap["Evt_HT_jets"] = -999.;
    variableMap["Evt_HT_tags"] = -999.;
    
    // missing transveral variables
    variableMap["Evt_MHT"] = -999.;
    variableMap["Evt_MET"] = -999.;
    variableMap["Evt_MTW"] = -999.;

    // btag likelihood ratios
    variableMap["Evt_blr"] = -999.;
    variableMap["Evt_blr_transformed"] = -999.;

    // Jet variables
    variableMap["Evt_M_JetsAverage"] = -999.;
    variableMap["Evt_Eta_JetsAverage"] = -999.;
    variableMap["Evt_Pt_JetsAverage"] = -999.;
    variableMap["Evt_E_JetsAverage"] = -999.;
    variableMap["Evt_M_TaggedJetsAverage"] = -999.;
    variableMap["Evt_Eta_TaggedJetsAverage"] = -999.;
    variableMap["Evt_Pt_TaggedJetsAverage"] = -999.;
    variableMap["Evt_E_TaggedJetsAverage"] = -999.;
    variableMap["Evt_M_UntaggedJetsAverage"] = -999.;
    variableMap["Evt_Eta_UntaggedJetsAverage"] = -999.;
    variableMap["Evt_Pt_UntaggedJetsAverage"] = -999.;
    variableMap["Evt_E_UntaggedJetsAverage"] = -999.;

    // dijet variables (tagged jets)
    variableMap["Evt_M2_closestTo125TaggedJets"] = -999.;
    variableMap["Evt_M2_closestTo91TaggedJets"] = -999.;
    variableMap["Evt_M2_TaggedJetsAverage"] = -999.;
    variableMap["Evt_Dr_TaggedJetsAverage"] = -999.;
    variableMap["Evt_Deta_TaggedJetsAverage"] = -999.;
    variableMap["Evt_M2_minDrTaggedJets"] = -999.;
    variableMap["Evt_Dr_minDrTaggedJets"] = -999.;
    variableMap["Evt_Pt_minDrTaggedJets"] = -999.;
    variableMap["Evt_Dr_maxDrTaggedJets"] = -999.;
    variableMap["Evt_Dr_closestTo91TaggedJets"] = -999.;

    // dijet variables (all jets)
    variableMap["Evt_M2_JetsAverage"] = -999.;
    variableMap["Evt_Dr_JetsAverage"] = -999.;
    variableMap["Evt_Deta_JetsAverage"] = -999.;
    variableMap["Evt_M2_minDrJets"] = -999.;
    variableMap["Evt_Dr_minDrJets"] = -999.;
    variableMap["Evt_Pt_minDrJets"] = -999.;
    variableMap["Evt_Dr_maxDrJets"] = -999.;

    // dijet variables (untagged jets)
    variableMap["Evt_M2_UntaggedJetsAverage"] = -999.;
    variableMap["Evt_Dr_UntaggedJetsAverage"] = -999.;
    variableMap["Evt_Deta_UntaggedJetsAverage"] = -999.;
    variableMap["Evt_M2_minDrUntaggedJets"] = -999.;
    variableMap["Evt_Dr_minDrUntaggedJets"] = -999.;
    variableMap["Evt_Pt_minDrUntaggedJets"] = -999.;
    variableMap["Evt_Dr_maxDrUntaggedJets"] = -999.;

    // max detas 
    variableMap["Evt_Deta_maxDetaJetJet"] = -999.;
    variableMap["Evt_Deta_maxDetaTagTag"] = -999.;
    variableMap["Evt_Deta_maxDetaJetTag"] = -999.;

    // lepton + jet variables
    variableMap["Evt_Dr_minDrLepTag"] = -999.;
    variableMap["Evt_Dr_minDrLepJet"] = -999.;
    variableMap["Evt_M_minDrLepTag"] = -999.;
    variableMap["Evt_M_minDrLepJet"] = -999.;


    // event shape variables
    variableMap["Evt_aplanarity"] = -999.;
    variableMap["Evt_aplanarity_jets"] = -999.;
    variableMap["Evt_aplanarity_tags"] = -999.;
    variableMap["Evt_sphericity"] = -999.;
    variableMap["Evt_sphericity_jets"] = -999.;
    variableMap["Evt_sphericity_tags"] = -999.;
    variableMap["Evt_transverse_sphericity"] = -999.;
    variableMap["Evt_transverse_sphericity_jets"] = -999.;
    variableMap["Evt_transverse_sphericity_tags"] = -999.;
    variableMap["Evt_JetPt_over_JetE"] = -999.;
    variableMap["Evt_TaggedJetPt_over_TaggedJetE"] = -999.;

    // reconstructed WLep variables
    variableMap["Reco_WLep_E"] = -999.;
    variableMap["Reco_WLep_Eta"] = -999.;
    variableMap["Reco_WLep_Phi"] = -999.;
    variableMap["Reco_WLep_Pt"] = -999.;
    variableMap["Reco_WLep_Mass"] = -999.;

    
}

MVAvars::~MVAvars()
{
}

void MVAvars::FillMVAvarMap(const std::vector<TLorentzVector> &selectedLeptonP4,
                      const std::vector<TLorentzVector> &selectedJetP4,
                      const std::vector<double> &selectedJetCSV,
                      const TLorentzVector &metP4)
{
    //check if CSVWP is set properly
    if(btagMcut<0){
        std::cerr << "Please set CSV Wp properly\n" << std::endl;
        throw std::exception();
    }
    // Reset all map entries to globalDefault value so that noting is left over from the last event
    ResetVariableMap();


    // ==================================================
    // construct object vectors etc

    int njets = selectedJetP4.size();
    int ntags = 0;
    for (uint i = 0; i < selectedJetCSV.size(); i++)
    {
        if (selectedJetCSV[i] > btagMcut)
            ntags++;
    }

    std::vector<double> selectedJetCSV_fixed;
    std::vector<double> looseSelectedJetCSV;

    for (uint i = 0; i < selectedJetCSV.size(); i++)
    {
        double tag = selectedJetCSV[i];
        if (std::isnan(tag))
        {
            tag = -.1;
        }
        else if (tag < 0)
        {
            tag = -.1;
        }
        else if (tag > 1)
        {
            tag = 1.;
        }
        selectedJetCSV_fixed.push_back(tag);
        //         looseSelectedJetCSV.push_back(tag);
    }

    std::vector<TLorentzVector> looseSelectedJetP4;
    std::vector<TLorentzVector> selectedTaggedJetP4;
    std::vector<TLorentzVector> selectedUntaggedJetP4;
    for (uint i = 0; i < selectedJetP4.size(); i++)
    {
        //         looseSelectedJetP4.push_back(selectedJetP4[i]);
        if (selectedJetCSV_fixed[i] > btagMcut)
        {
            selectedTaggedJetP4.push_back(selectedJetP4[i]);
        }
        else 
        {
            selectedUntaggedJetP4.push_back(selectedJetP4[i]);
        }
    }

    vector<vector<double>> jets_vvdouble;
    for (auto jet = selectedJetP4.begin(); jet != selectedJetP4.end(); jet++)
    {
        vector<double> pxpypzE;
        pxpypzE.push_back(jet->Px());
        pxpypzE.push_back(jet->Py());
        pxpypzE.push_back(jet->Pz());
        pxpypzE.push_back(jet->E());
        jets_vvdouble.push_back(pxpypzE);
    }

    vector<vector<double>> jets_vvdouble_tagged;
    for (auto jet = selectedTaggedJetP4.begin(); jet != selectedTaggedJetP4.end(); jet++)
    {
        vector<double> pxpypzE;
        pxpypzE.push_back(jet->Px());
        pxpypzE.push_back(jet->Py());
        pxpypzE.push_back(jet->Pz());
        pxpypzE.push_back(jet->E());
        jets_vvdouble_tagged.push_back(pxpypzE);
    }

    /// create jets_vvdouble_tagged
    vector<double> sortedCSV = selectedJetCSV_fixed;
    std::sort(sortedCSV.begin(), sortedCSV.end(), std::greater<double>());

    // ==================================================
    // calculate variables
    // aplanarity and sphericity

    double aplanarity, sphericity, transverse_sphericity;
    bdtvar.getSp(selectedLeptonP4[0], metP4, selectedJetP4, aplanarity, sphericity, transverse_sphericity);

    TLorentzVector dummy_Lepton;
    TLorentzVector dummy_met;

    double aplanarity_jets, sphericity_jets, transverse_sphericity_jets;
    bdtvar.getSp(dummy_Lepton, dummy_met, selectedJetP4, aplanarity_jets, sphericity_jets, transverse_sphericity_jets);

    double aplanarity_tags, sphericity_tags, transverse_sphericity_tags;
    bdtvar.getSp(dummy_Lepton, dummy_met, selectedTaggedJetP4, aplanarity_tags, sphericity_tags, transverse_sphericity_tags);

    // Fox Wolfram
    double h0, h1, h2, h3, h4;
    bdtvar.getFox(selectedJetP4, h0, h1, h2, h3, h4);

    // best higgs mass 1
    double minChi, dRbb;
    TLorentzVector bjet1, bjet2;
    double bestHiggsMass = -9.9;
    // if (category == "6j4t")
    if (njets >= 6 && ntags >= 4)
    {
        bestHiggsMass = bdtvar.getBestHiggsMass(selectedLeptonP4[0], metP4, selectedJetP4, selectedJetCSV_fixed, minChi, dRbb, bjet1, bjet2, looseSelectedJetP4, looseSelectedJetCSV, btagMcut);
    }

    // study top bb system
    TLorentzVector dummy_metv;
    double minChiStudy, chi2lepW, chi2leptop, chi2hadW, chi2hadtop, mass_lepW, mass_leptop, mass_hadW, mass_hadtop, dRbbStudy;
    double testquant1, testquant2, testquant3, testquant4, testquant5, testquant6, testquant7;

    TLorentzVector b1, b2;
    bdtvar.study_tops_bb_syst(
        metP4.Pt(), metP4.Phi(), dummy_metv, selectedLeptonP4[0], jets_vvdouble, selectedJetCSV_fixed, 
        minChiStudy, chi2lepW, chi2leptop, chi2hadW, chi2hadtop, 
        mass_lepW, mass_leptop, mass_hadW, mass_hadtop, 
        dRbbStudy, 
        testquant1, testquant2, testquant3, testquant4, testquant5, testquant6, testquant7, 
        b1, b2, btagMcut);

    // fn
    double dEta_fn = testquant6;
    // ptE ratios
    double pt_E_ratio = bdtvar.pt_E_ratio_jets(jets_vvdouble);
    double pt_E_ratio_tag = bdtvar.pt_E_ratio_jets(jets_vvdouble_tagged);

    // etamax
    double jet_jet_etamax = bdtvar.get_jet_jet_etamax(jets_vvdouble);
    double jet_tag_etamax = bdtvar.get_jet_tag_etamax(jets_vvdouble, selectedJetCSV_fixed, btagMcut);
    double tag_tag_etamax = bdtvar.get_tag_tag_etamax(jets_vvdouble, selectedJetCSV_fixed, btagMcut);

    // jet variables
    double ht_jets = 0;
    double ht_taggedjets = 0;
    double Mlj = 0;
    double dr_between_lep_and_closest_jet = 99;
    double mht_px = 0;
    double mht_py = 0;
    TLorentzVector sum_pt_vec = selectedLeptonP4[0];
    sum_pt_vec += metP4;
    for (auto jetvec = selectedJetP4.begin(); jetvec != selectedJetP4.end(); ++jetvec)
    {
        ht_jets += jetvec->Pt();
        mht_px += jetvec->Px();
        mht_py += jetvec->Py();
        sum_pt_vec += *jetvec;

        double drLep = selectedLeptonP4[0].DeltaR(*jetvec);
        if (drLep < dr_between_lep_and_closest_jet)
        {
            dr_between_lep_and_closest_jet = drLep;
            Mlj = (selectedLeptonP4[0] + *jetvec).M();
        }

    }

    mht_px += selectedLeptonP4[0].Px();
    mht_py += selectedLeptonP4[0].Py();
    double mass_of_everything = sum_pt_vec.M();
    double ht_wo_met = ht_jets + selectedLeptonP4[0].Pt();
    double ht = metP4.Pt() + ht_wo_met;
    double MHT = sqrt(mht_px * mht_px + mht_py * mht_py);

    double MTW = sqrt(2*(selectedLeptonP4[0].Pt()*metP4.Pt() - selectedLeptonP4[0].Px()*metP4.Px() - selectedLeptonP4[0].Py()*metP4.Py()));

    // mass of lepton and closest bt-tagged jet
    double Mlb = 0;
    double dr_between_lep_and_closest_tagged_jet = 999.;
    for (auto tagged_jet = selectedTaggedJetP4.begin(); tagged_jet != selectedTaggedJetP4.end(); tagged_jet++)
    {
        ht_taggedjets += tagged_jet->Pt();

        double drLep = selectedLeptonP4[0].DeltaR(*tagged_jet);
        if (drLep < dr_between_lep_and_closest_tagged_jet)
        {
            dr_between_lep_and_closest_tagged_jet = drLep;
            Mlb = (selectedLeptonP4[0] + *tagged_jet).M();
        }

    }

    // JET AVERAGES
    // jet variables
    double sumJetM = 0;
    double sumJetPt = 0;
    double sumJetE = 0;
    double sumJetEta = 0;
    int njets_all = 0;
    for (auto itjet = selectedJetP4.begin(); itjet != selectedJetP4.end(); itjet++)
    {
        njets_all ++;
        sumJetM += itjet->M();
        sumJetPt += itjet->Pt();
        sumJetE += itjet->E();
        sumJetEta += itjet->Eta();
    }
    double avgJetM = -99.;
    double avgJetPt = -99.;
    double avgJetE = -99.;
    double avgJetEta = -99.;
    if (njets_all>0)
    {
        avgJetM = sumJetM/njets_all;
        avgJetPt = sumJetPt/njets_all;
        avgJetE = sumJetE/njets_all;
        avgJetEta = sumJetEta/njets_all;
    }

    // tagged jet variables
    double sumTaggedJetM = 0;
    double sumTaggedJetPt = 0;
    double sumTaggedJetE = 0;
    double sumTaggedJetEta = 0;
    int njets_tagged = 0;
    for (auto itjet = selectedTaggedJetP4.begin(); itjet != selectedTaggedJetP4.end(); itjet++)
    {
        njets_tagged ++;
        sumTaggedJetM += itjet->M();
        sumTaggedJetPt += itjet->Pt();
        sumTaggedJetE += itjet->E();
        sumTaggedJetEta += itjet->Eta();
    }
    double avgTaggedJetM = -99.;
    double avgTaggedJetPt = -99.;
    double avgTaggedJetE = -99.;
    double avgTaggedJetEta = -99.;
    if (njets_tagged>0)
    {
        avgTaggedJetM = sumTaggedJetM/njets_tagged;
        avgTaggedJetPt = sumTaggedJetPt/njets_tagged;
        avgTaggedJetE = sumTaggedJetE/njets_tagged;
        avgTaggedJetEta = sumTaggedJetEta/njets_tagged;
    }

    // untagged jet variables
    double sumUntaggedJetM = 0;
    double sumUntaggedJetPt = 0;
    double sumUntaggedJetE = 0;
    double sumUntaggedJetEta = 0;
    int njets_untagged = 0;
    for (auto itjet = selectedUntaggedJetP4.begin(); itjet != selectedUntaggedJetP4.end(); itjet++)
    {
        njets_untagged ++;
        sumUntaggedJetM += itjet->M();
        sumUntaggedJetPt += itjet->Pt();
        sumUntaggedJetE += itjet->E();
        sumUntaggedJetEta += itjet->Eta();
    }
    double avgUntaggedJetM = -99.;
    double avgUntaggedJetPt = -99.;
    double avgUntaggedJetE = -99.;
    double avgUntaggedJetEta = -99.;
    if (njets_untagged>0)
    {
        avgUntaggedJetM = sumUntaggedJetM/njets_untagged;
        avgUntaggedJetPt = sumUntaggedJetPt/njets_untagged;
        avgUntaggedJetE = sumUntaggedJetE/njets_untagged;
        avgUntaggedJetEta = sumUntaggedJetEta/njets_untagged;
    }
    

    // DIJET SYSTEMS
    // tagged jets
    double closest_tagged_dijet_mass = -99;
    double minDrTagged = 99;
    double minPtTagged = -99;
    double maxDrTagged = -99;
    double sumDrTagged = 0;
    double sumDetaTagged = 0;
    double sumM2Tagged = 0;
    double DrTagged91 = 99.;


    int npairs_tagged = 0;
    double tagged_dijet_mass_closest_to_125 = -999;
    double tagged_dijet_mass_closest_to_91 = -999;
    for (auto tagged_jet1 = selectedTaggedJetP4.begin(); tagged_jet1 != selectedTaggedJetP4.end(); tagged_jet1++)
    {
        for (auto tagged_jet2 = tagged_jet1 + 1; tagged_jet2 != selectedTaggedJetP4.end(); tagged_jet2++)
        {
            double dr = tagged_jet1->DeltaR(*tagged_jet2);
            double m = (*tagged_jet1 + *tagged_jet2).M();
            double pt = (*tagged_jet1 + *tagged_jet2).Pt();
            sumDrTagged += dr;
            sumM2Tagged += m;
            double deta = fabs(tagged_jet1->Eta() - tagged_jet2->Eta());
            sumDetaTagged += deta;
            npairs_tagged++;
            if (dr < minDrTagged)
            {
                minDrTagged = dr;
                minPtTagged = pt;
                closest_tagged_dijet_mass = m;
            }
            if (dr > maxDrTagged)
            {
                maxDrTagged = dr;
            }
            if (fabs(tagged_dijet_mass_closest_to_125 - 125) > fabs(m - 125))
            {
                tagged_dijet_mass_closest_to_125 = m;
            }
            if (fabs(tagged_dijet_mass_closest_to_91 - 91) > fabs(m - 91))
            {
                tagged_dijet_mass_closest_to_91 = m;
                DrTagged91 = dr;
            }
        }
    }

    double avgDrTagged = -1;
    double avgDetaTagged = 0;
    double avgM2Tagged = 0;
    if (npairs_tagged != 0)
    {
        avgDrTagged = sumDrTagged / (double)npairs_tagged;
        avgDetaTagged = sumDetaTagged / (double)npairs_tagged;
        avgM2Tagged = sumM2Tagged / (double)npairs_tagged;
    }

    // all jets
    double closest_dijet_mass = -99;
    double minDrJets = 99;
    double minPtJets = -99;
    double maxDrJets = -99;
    double sumDrJets = 0;
    double sumDetaJets = 0;
    double sumM2Jets = 0;

    int npairs_jets = 0;
    for (auto jet1 = selectedJetP4.begin(); jet1 != selectedJetP4.end(); jet1++)
    {
        for (auto jet2 = jet1 + 1; jet2 != selectedJetP4.end(); jet2++)
        {
            double dr = jet1->DeltaR(*jet2);
            double m = (*jet1 + *jet2).M();
            double pt = (*jet1 + *jet2).Pt();
            sumDrJets += dr;
            sumM2Jets += m;
            double deta = fabs(jet1->Eta() - jet2->Eta());
            sumDetaJets += deta;
            npairs_jets++;
            if (dr < minDrJets)
            {
                minDrJets = dr;
                minPtJets = pt;
                closest_dijet_mass = m;
            }
            if (dr > maxDrJets)
            {
                maxDrJets = dr;
            }
        }
    }

    double avgDrJets = -1;
    double avgDetaJets = 0;
    double avgM2Jets = 0;
    if (npairs_jets != 0)
    {
        avgDrJets = sumDrJets / (double)npairs_jets;
        avgDetaJets = sumDetaJets / (double)npairs_jets;
        avgM2Jets = sumM2Jets / (double)npairs_jets;
    }

    // untagged jets
    double closest_untagged_dijet_mass = -99.;
    double minDrUntagged = 99;
    double minPtUntagged = -99;
    double maxDrUntagged = -99;
    double sumDrUntagged = 0;
    double sumDetaUntagged = 0;
    double sumM2Untagged = 0;

    int npairs_untagged = 0;
    for (auto jet1 = selectedUntaggedJetP4.begin(); jet1 != selectedUntaggedJetP4.end(); jet1++)
    {
        for (auto jet2 = jet1 + 1; jet2 != selectedUntaggedJetP4.end(); jet2++)
        {
            double dr = jet1->DeltaR(*jet2);
            double m = (*jet1 + *jet2).M();
            double pt = (*jet1 + *jet2).Pt();
            sumDrUntagged += dr;
            sumM2Untagged += m;
            double deta = fabs(jet1->Eta() - jet2->Eta());
            sumDetaUntagged += deta;
            npairs_untagged++;
            if (dr < minDrUntagged)
            {
                minDrUntagged = dr;
                minPtUntagged = pt;
                closest_untagged_dijet_mass = m;
            }
            if (dr > maxDrUntagged)
            {
                maxDrUntagged = dr;
            }
        }
    }

    double avgDrUntagged = -1;
    double avgDetaUntagged = 0;
    double avgM2Untagged = 0;
    if (npairs_untagged != 0)
    {
        avgDrUntagged = sumDrUntagged / (double)npairs_untagged;
        avgDetaUntagged = sumDetaUntagged / (double)npairs_untagged;
        avgM2Untagged = sumM2Untagged / (double)npairs_untagged;
    }


    // M3
    double m3 = -1.;
    double maxpt_for_m3 = -1;
    for (auto itJetVec1 = selectedJetP4.begin(); itJetVec1 != selectedJetP4.end(); ++itJetVec1)
    {
        for (auto itJetVec2 = itJetVec1 + 1; itJetVec2 != selectedJetP4.end(); ++itJetVec2)
        {
            for (auto itJetVec3 = itJetVec2 + 1; itJetVec3 != selectedJetP4.end(); ++itJetVec3)
            {

                TLorentzVector m3vec = *itJetVec1 + *itJetVec2 + *itJetVec3;

                if (m3vec.Pt() > maxpt_for_m3)
                {
                    maxpt_for_m3 = m3vec.Pt();
                    m3 = m3vec.M();
                }
            }
        }
    }

    // M3 one tagged
    double m3_onetagged = -1.;
    double maxpt_for_m3_onetagged = -1.;
    for (auto itTaggedJetVec1 = selectedTaggedJetP4.begin(); itTaggedJetVec1 != selectedTaggedJetP4.end(); ++itTaggedJetVec1)
    {
        for( auto itUntaggedJetVec1 = selectedUntaggedJetP4.begin(); itUntaggedJetVec1 != selectedUntaggedJetP4.end(); ++itUntaggedJetVec1)
        {
            for( auto itUntaggedJetVec2 = itUntaggedJetVec1+1; itUntaggedJetVec2 != selectedUntaggedJetP4.end(); ++itUntaggedJetVec2)
            {
                TLorentzVector m3vec_onetagged = *itTaggedJetVec1 + *itUntaggedJetVec1 + *itUntaggedJetVec2;
                if (m3vec_onetagged.Pt() > maxpt_for_m3_onetagged)
                {
                    maxpt_for_m3_onetagged = m3vec_onetagged.Pt();
                    m3_onetagged = m3vec_onetagged.M();
                }
            }
        }
    }

    double detaJetsAverage = 0;
    int nPairsJets = 0;
    for (auto itJetVec1 = selectedJetP4.begin(); itJetVec1 != selectedJetP4.end(); ++itJetVec1)
    {
        for (auto itJetVec2 = itJetVec1 + 1; itJetVec2 != selectedJetP4.end(); ++itJetVec2)
        {
            detaJetsAverage += fabs(itJetVec1->Eta() - itJetVec2->Eta());
            nPairsJets++;
        }
    }
    if (nPairsJets > 0)
    {
        detaJetsAverage /= (double)nPairsJets;
    }

    // btag variables


    // CSV averages
    double averageCSV_tagged = 0;
    double averageCSV_all = 0;
    double lowest_btag_all = 99;
    double lowest_btag_tagged = 99;
    // int njets = selectedJetP4.size();
    // int ntags = 0;
    for (auto itCSV = selectedJetCSV_fixed.begin(); itCSV != selectedJetCSV_fixed.end(); ++itCSV)
    {
        averageCSV_all += fmax(*itCSV, 0);
        lowest_btag_all = fmin(*itCSV, lowest_btag_all);
        if (*itCSV < btagMcut)
            continue;
        lowest_btag_tagged = fmin(*itCSV, lowest_btag_tagged);
        averageCSV_tagged += fmax(*itCSV, 0);
        // ntags++;
    }
    if (ntags > 0)
        averageCSV_tagged /= (double)ntags;
    else
        averageCSV_tagged = 0;
    if (selectedJetCSV_fixed.size() > 0)
        averageCSV_all /= selectedJetCSV_fixed.size();
    else
        averageCSV_all = 0;

    if (lowest_btag_all > 90)
        lowest_btag_all = -1;
    if (lowest_btag_tagged > 90)
        lowest_btag_tagged = -1;


    // squared CSV deviations from average
    double csvDev_all = 0;
    double csvDev_tagged = 0;
    for (auto itCSV = selectedJetCSV_fixed.begin(); itCSV != selectedJetCSV_fixed.end(); ++itCSV)
    {
        csvDev_all += pow(*itCSV - averageCSV_all, 2);
        if (*itCSV < btagMcut)
            continue;
        csvDev_tagged += pow(*itCSV - averageCSV_tagged, 2);
    }
    if (ntags > 0)
        csvDev_tagged /= (double)ntags;
    else
        csvDev_tagged = -1.;
    if (selectedJetCSV_fixed.size() > 0)
        csvDev_all /= selectedJetCSV_fixed.size();
    else
        csvDev_all = -1;
    


    //calculate blr_transformed
    double blr_transformed = -999.0;
    std::vector<unsigned int> out_best_perm;
    double out_P_4b = -1;
    double out_P_2b = -1;
    double eth_blr = -1;
    if (selectedJetP4.size() > 5)
    {
        eth_blr = mem.GetBTagLikelihoodRatio(selectedJetP4,
                                             selectedJetCSV_fixed,
                                             out_best_perm,
                                             out_P_4b,
                                             out_P_2b);
        blr_transformed = log(eth_blr / (1.0 - eth_blr));
    }
//    blr_transformed = log(eth_blr / (1.0 - eth_blr));


    // reconstruct Wlep
    TLorentzVector Wlep = GetLeptonicW(selectedLeptonP4.at(0), metP4);
    
    // ==================================================
    // Fill variable map
    
    // higgs specific variables
    variableMap["Reco_best_higgs_mass"] = bestHiggsMass;
    variableMap["Reco_dEta_fn"] = dEta_fn;

    // CSV variables
    variableMap["Evt_CSV_min"] = lowest_btag_all;
    variableMap["Evt_CSV_min_tagged"] = lowest_btag_tagged;
    variableMap["Evt_CSV_avg"] = averageCSV_all;
    variableMap["Evt_CSV_avg_tagged"] = averageCSV_tagged;
    variableMap["Evt_CSV_dev"] = csvDev_all;
    variableMap["Evt_CSV_dev_tagged"] = csvDev_tagged;

    // fox wolfram moments
    variableMap["Evt_h0"] = h0;
    variableMap["Evt_h1"] = h1;
    variableMap["Evt_h2"] = h2;
    variableMap["Evt_h3"] = h3;

    // event masses 
    variableMap["Evt_M_Total"] = mass_of_everything;
    variableMap["Evt_M3"] = m3;
    variableMap["Evt_M3_oneTagged"] = m3_onetagged;

    // HT variables
    variableMap["Evt_HT"] = ht;
    variableMap["Evt_HT_wo_MET"] = ht_wo_met;
    variableMap["Evt_HT_jets"] = ht_jets;
    variableMap["Evt_HT_tags"] = ht_taggedjets;
    
    // missing transveral variables
    variableMap["Evt_MHT"] = MHT;
    variableMap["Evt_MET"] = metP4.Pt();
    variableMap["Evt_MTW"] = MTW;

    // btag likelihood ratios
    variableMap["Evt_blr"] = eth_blr;
    variableMap["Evt_blr_transformed"] = blr_transformed;

    // Jet variables
    variableMap["Evt_M_JetsAverage"] = avgJetM;
    variableMap["Evt_Eta_JetsAverage"] = avgJetEta;
    variableMap["Evt_Pt_JetsAverage"] = avgJetPt;
    variableMap["Evt_E_JetsAverage"] = avgJetE;
    variableMap["Evt_M_TaggedJetsAverage"] = avgTaggedJetM;
    variableMap["Evt_Eta_TaggedJetsAverage"] = avgTaggedJetEta;
    variableMap["Evt_Pt_TaggedJetsAverage"] = avgTaggedJetPt;
    variableMap["Evt_E_TaggedJetsAverage"] = avgTaggedJetE;
    variableMap["Evt_M_UntaggedJetsAverage"] = avgUntaggedJetM;
    variableMap["Evt_Eta_UntaggedJetsAverage"] = avgUntaggedJetEta;
    variableMap["Evt_Pt_UntaggedJetsAverage"] = avgUntaggedJetPt;
    variableMap["Evt_E_UntaggedJetsAverage"] = avgUntaggedJetE;

    // dijet variables (tagged jets)
    variableMap["Evt_M2_closestTo125TaggedJets"] = tagged_dijet_mass_closest_to_125;
    variableMap["Evt_M2_closestTo91TaggedJets"] = tagged_dijet_mass_closest_to_91;
    variableMap["Evt_M2_TaggedJetsAverage"] = avgM2Tagged;
    variableMap["Evt_Dr_TaggedJetsAverage"] = avgDrTagged;
    variableMap["Evt_Deta_TaggedJetsAverage"] = avgDetaTagged;
    variableMap["Evt_M2_minDrTaggedJets"] = closest_tagged_dijet_mass;
    variableMap["Evt_Dr_minDrTaggedJets"] = minDrTagged;
    variableMap["Evt_Pt_minDrTaggedJets"] = minPtTagged;
    variableMap["Evt_Dr_maxDrTaggedJets"] = maxDrTagged;
    variableMap["Evt_Dr_closestTo91TaggedJets"] = DrTagged91;


    // dijet variables (all jets)
    variableMap["Evt_M2_JetsAverage"] = avgM2Jets;
    variableMap["Evt_Dr_JetsAverage"] = avgDrJets;
    variableMap["Evt_Deta_JetsAverage"] = avgDetaJets;
    variableMap["Evt_M2_minDrJets"] = closest_dijet_mass;
    variableMap["Evt_Dr_minDrJets"] = minDrJets;
    variableMap["Evt_Pt_minDrJets"] = minPtJets;
    variableMap["Evt_Dr_maxDrJets"] = maxDrJets;

    // dijet variables (untagged jets)
    variableMap["Evt_M2_UntaggedJetsAverage"] = avgM2Untagged;
    variableMap["Evt_Dr_UntaggedJetsAverage"] = avgDrUntagged;
    variableMap["Evt_Deta_UntaggedJetsAverage"] = avgDetaUntagged;
    variableMap["Evt_M2_minDrUntaggedJets"] = closest_untagged_dijet_mass;
    variableMap["Evt_Dr_minDrUntaggedJets"] = minDrUntagged;
    variableMap["Evt_Pt_minDrUntaggedJets"] = minPtUntagged;
    variableMap["Evt_Dr_maxDrUntaggedJets"] = maxDrUntagged;

    // max detas 
    variableMap["Evt_Deta_maxDetaJetJet"] = jet_jet_etamax;
    variableMap["Evt_Deta_maxDetaTagTag"] = tag_tag_etamax;
    variableMap["Evt_Deta_maxDetaJetTag"] = jet_tag_etamax;

    // lepton + jet variables
    variableMap["Evt_Dr_minDrLepTag"] = dr_between_lep_and_closest_tagged_jet;
    variableMap["Evt_Dr_minDrLepJet"] = dr_between_lep_and_closest_jet;
    variableMap["Evt_M_minDrLepTag"] = Mlb;
    variableMap["Evt_M_minDrLepJet"] = Mlj;

    // event shape variables
    variableMap["Evt_aplanarity"] = aplanarity;
    variableMap["Evt_aplanarity_jets"] = aplanarity_jets;
    variableMap["Evt_aplanarity_tags"] = aplanarity_tags;
    variableMap["Evt_sphericity"] = sphericity;
    variableMap["Evt_sphericity_jets"] = sphericity_jets;
    variableMap["Evt_sphericity_tags"] = sphericity_tags;
    variableMap["Evt_transverse_sphericity"] = transverse_sphericity;
    variableMap["Evt_transverse_sphericity_jets"] = transverse_sphericity_jets;
    variableMap["Evt_transverse_sphericity_tags"] = transverse_sphericity_tags;
    variableMap["Evt_JetPt_over_JetE"] = pt_E_ratio;
    variableMap["Evt_TaggedJetPt_over_TaggedJetE"] = pt_E_ratio_tag;

    // reconstructed WLep variables
    variableMap["Reco_WLep_E"] = Wlep.E();
    variableMap["Reco_WLep_Eta"] = Wlep.Eta();
    variableMap["Reco_WLep_Phi"] = Wlep.Phi();
    variableMap["Reco_WLep_Pt"] = Wlep.Pt();
    variableMap["Reco_WLep_Mass"] = Wlep.M();

}

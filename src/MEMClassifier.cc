#include "TTH/CommonClassifier/interface/MEMClassifier.h"

static const double mem_weight = 0.10;
/*This code is the interface for the MEM integrator. It receives a collection of objects which is processed 
in order to be passed on to the integrator. 
For the leptons and the MET, only the kinematics and the charge (for the lepton) needs to be provided.
For the jets, the situation is more complicated because the jet kinematics and therefore the MEM output
changes with different factorized JES sources, and might also cause the event category to change. When the
category doesn't change, the MEM is integrated over a vector of kinematics therefore one needs to provide 
the list of all jet variations. When the event category changes, the MEM integration is run again. In order
to consider the correct jets, a vector is needed for all jet properties, where 
each entry corresponds to a vector of all jets considered for the particular JES source. A vector of bools 
indicating if the event category changed for each of the sources is also needed.
*/


void MEMClassifier::setup_mem(
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
) {

    integrand->set_cfg(cfg);
    integrand->next_event();

    switch (hypo) {

    case HYPO_DEBUG: {
        break;
    }

    case SL_0W2H2T:
    case SL_1W2H2T:
    case SL_2W2H2T:
    case DL_0W2H2T: {
        setup_mem_impl(selectedLeptonP4,
                       selectedLeptonCharge,
                       selectedJetP4,
                       selectedJetCSV,
                       selectedJetType,
                       selectedJetVariation,
                       looseSelectedJetP4,
                       looseSelectedJetCSV,
                       metP4,
                       objs,
                       res,
                       changes_jet_category);
        break;
                    }
    case SL_2W2H2T_SJ: {
        setup_mem_sl_2w2h2t_sj(selectedLeptonP4,
                               selectedLeptonCharge,
                               selectedJetP4,
                               selectedJetCSV,
                               selectedJetType,
                               selectedJetVariation,
                               looseSelectedJetP4,
                               looseSelectedJetCSV,
                               metP4,
                               objs,
                               res,
                               changes_jet_category);
        break;
    }

    // Fallback
    default: {
        std::cerr << "Warning! Fallback case reached. Invalid MEM Hypo!!!" << std::endl;
        throw std::exception();
    }

    }

}


void MEMClassifier::setup_mem_impl(
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
) {
    
    if (!(selectedLeptonP4.size() == 1) && !(selectedLeptonP4.size() == 2)) {
        throw std::runtime_error("expected a single-lepton or dilepton event");
    }

    std::vector<unsigned int> best_perm;
    double blr_4b = 0.0;
    double blr_2b = 0.0;

    GetBTagLikelihoodRatio(
        selectedJetP4, selectedJetCSV, best_perm, blr_4b, blr_2b
    );

    assert(best_perm.size() >= 6);

    cout << "best_perm = [";
    for (auto ip : best_perm) {
        cout << ip << " ";
    }
    cout << "]" << endl;

    res.blr_4b = blr_4b;
    res.blr_2b = blr_2b;

    std::vector<MEM::Object*> tagged;
    std::vector<MEM::Object*> untagged;
    //use up to numMaxJets jets
    for (unsigned int ij=0; ij<std::min(selectedJetP4.size(), numMaxJets); ij++) {
    //for (auto ij : best_perm) {
        TLorentzVector p4 = selectedJetP4.at(ij);
        assert(p4.Pt() > 0);
        //Check if this jet was in the best 4b permutation, i.e. the first 4 indices of the permutation
        auto last = best_perm.begin() + 6;

        bool is_btagged = std::find(best_perm.begin(), last, ij) != last;

        MEM::Object* jet = make_jet(
                               p4.Pt(), p4.Eta(), p4.Phi(), p4.M(), is_btagged ? 1.0 : 0.0,
                               selectedJetCSV.at(ij),
                               false
                           );
        if (selectedJetVariation.size() > 0) {
            for (unsigned int ivar=0; ivar < selectedJetVariation.at(ij).size(); ivar++) {
                jet->p4_variations.push_back(selectedJetVariation.at(ij).at(ivar));
            }
        }

        if (is_btagged) {
            tagged.push_back(jet);
        } else {
            untagged.push_back(jet);
        }
        objs.push_back(jet);
        //integrand->push_back_object(jet);
    }
    for (unsigned int ij=0; ij<tagged.size(); ij++) {
        integrand->push_back_object(tagged[ij]);
    }
    for (unsigned int ij=0; ij<untagged.size(); ij++) {
        integrand->push_back_object(untagged[ij]);
    }

    assert(tagged.size() == 6);

    for (unsigned int il=0; il < selectedLeptonP4.size(); il++) {
        TLorentzVector lep_p4 = selectedLeptonP4.at(il);
        assert(lep_p4.Pt() > 0);
        MEM::Object* lep = make_lepton(lep_p4.Pt(), lep_p4.Eta(), lep_p4.Phi(), lep_p4.M(), selectedLeptonCharge[il]);
        objs.push_back(lep);
        integrand->push_back_object(lep);
    }

    assert(metP4.Pt() > 0);
    MEM::Object* met = new MEM::Object(metP4, MEM::ObjectType::MET );
    objs.push_back(met);
    integrand->push_back_object(met);
}


void MEMClassifier::setup_mem_sl_2w2h2t_sj(
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
) {

    if (selectedLeptonP4.size() != 1) {
        throw std::runtime_error("Expected a single-lepton event");
    }

    // Make sure we received a suitable set of input jets
    // Require:
    //  * exactly 6 jets
    //  * 2 light boosted
    //  * 1 b boosted
    //  * 3 resolved

    int n_light_boosted(0), n_b_boosted(0), n_resolved(0);

    assert(selectedJetP4.size() == 6);

    std::vector<MEM::Object*> tagged;
    std::vector<MEM::Object*> untagged;

    for (unsigned int ij=0; ij<selectedJetP4.size(); ij++) {

        TLorentzVector p4 = selectedJetP4.at(ij);
        assert(p4.Pt() > 0);

        switch (selectedJetType[ij]) {

        case RESOLVED: {
            MEM::Object* jet = make_jet(
                                   p4.Pt(), p4.Eta(), p4.Phi(), p4.M(), 1.0,
                                   selectedJetCSV.at(ij),
                                   false
                               );

            tagged.push_back(jet);
            n_resolved++;
            break;
        }

        case BOOSTED_LIGHT: {
            MEM::Object* jet = make_jet(
                                   p4.Pt(), p4.Eta(), p4.Phi(), p4.M(), 0.0,
                                   selectedJetCSV.at(ij),
                                   true
                               );
            untagged.push_back(jet);
            n_light_boosted++;
            break;
        }

        case BOOSTED_B: {
            MEM::Object* jet = make_jet(
                                   p4.Pt(), p4.Eta(), p4.Phi(), p4.M(), 1.0,
                                   selectedJetCSV.at(ij),
                                   true
                               );
            tagged.push_back(jet);
            n_b_boosted++;
            break;
        }
        }
    }

    assert(n_light_boosted==2);
    assert(n_b_boosted==1);
    assert(n_resolved==3);

    for (auto* jet : tagged) {
        objs.push_back(jet);
        integrand->push_back_object(jet);
    }

    for (auto* jet : untagged) {
        objs.push_back(jet);
        integrand->push_back_object(jet);
    }

    for (unsigned int il=0; il < selectedLeptonP4.size(); il++) {
        TLorentzVector lep_p4 = selectedLeptonP4.at(il);
        assert(lep_p4.Pt() > 0);
        MEM::Object* lep = make_lepton(lep_p4.Pt(), lep_p4.Eta(), lep_p4.Phi(), lep_p4.M(), selectedLeptonCharge[il]);
        objs.push_back(lep);
        integrand->push_back_object(lep);
    }

    assert(metP4.Pt() > 0);
    MEM::Object* met = new MEM::Object(metP4, MEM::ObjectType::MET );
    integrand->push_back_object(met);
}

std::vector<MEMResult> MEMClassifier::GetOutput(
    const std::vector<TLorentzVector>& selectedLeptonP4,
    const std::vector<double>& selectedLeptonCharge,
    const std::vector<vector<TLorentzVector>>& selectedJetP4,
    const std::vector<vector<double>>& selectedJetCSV,
    const std::vector<vector<JetType>>& selectedJetType,
    const std::vector<vector<vector<double>>>& selectedJetVariation,
    TLorentzVector& metP4,
    const std::vector<bool>& changes_jet_category
) {
    const int nleps = selectedLeptonP4.size();
    const int njets = selectedJetP4.at(0).size();

    auto hypo = SL_0W2H2T;
    if (nleps == 1) {
        if (njets >= 6) {
            hypo = SL_2W2H2T;
            cout << "njets=" << njets << " nleps=" << nleps << " chose hypo SL_2w2h2t" << endl;
        }
        else if (njets == 5) {
            hypo = SL_1W2H2T;
            cout << "njets=" << njets << " nleps=" << nleps << " chose hypo SL_1w2h2t" << endl;
        }
        else if (njets == 4) {
            hypo = SL_0W2H2T;
            cout << "njets=" << njets << " nleps=" << nleps << " chose hypo SL_0w2h2t" << endl;
        }
        else  {
            cout << "njets=" << njets << " nleps=" << nleps << " chose default hypo SL_0w2h2t" << endl;
        }
    } else if (nleps == 2) {
        hypo = DL_0W2H2T;
            cout << "njets=" << njets << " nleps=" << nleps << " chose hypo DL_0w2h2t" << endl;
    } else {
        throw std::runtime_error("Expected a single-lepton event or dilepton event");
    }

    return GetOutput(
        hypo,
        selectedLeptonP4,
        selectedLeptonCharge,
        selectedJetP4,
        selectedJetCSV,
        selectedJetType,
        selectedJetVariation,
        {},
        {},
        metP4,
        changes_jet_category
    );
}


std::vector<MEMResult> MEMClassifier::GetOutput(
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
    int ncalls
) {

    // Make sure vector sizes match up
    assert(selectedLeptonP4.size() == selectedLeptonCharge.size());
    assert(selectedJetP4.size() == selectedJetCSV.size());
    assert(selectedJetP4.size() == selectedJetType.size());
    assert(looseSelectedJetP4.size() == looseSelectedJetCSV.size());
    //This part needs to be looper over for all systematics

    vector<MEMResult> res;

    for (unsigned int i=0; i < selectedJetVariation.size();i++) {
        MEMResult res2;
        res.push_back(res2);
        //Change hypothesis for the JES sources if needed
        if (i==0 || changes_jet_category.at(i) == true) {
            auto newhypo = hypo;
            if (i>0) {
                const int nleps = selectedLeptonP4.size();
                const int njets = selectedJetP4.at(i).size();
                if (nleps == 1) {
                    if (njets >= 6) {
                        newhypo = SL_2W2H2T;
                        cout << "Changed hypo to SL_2w2h2t" << endl;
                    }
                    else if (njets == 5) {
                        newhypo = SL_1W2H2T;
                        cout << "Changed hypo to SL_1w2h2t" << endl;
                    }
                    else if (njets == 4) {
                        newhypo = SL_0W2H2T;
                        cout << "Changed hypo to SL_0w2h2t" << endl;
                    }
                } else if (nleps == 2) {
                    newhypo = DL_0W2H2T;
                        cout << "Changed hypo to DL_0w2h2t" << endl;
                }

            }
            std::vector<MEM::Object*> objs;


            setup_mem(
                newhypo,
                selectedLeptonP4,
                selectedLeptonCharge,
                selectedJetP4.at(i),
                selectedJetCSV.at(i),
                selectedJetType.at(i),
                selectedJetVariation.at(i),
                looseSelectedJetP4,
                looseSelectedJetCSV,
                metP4, objs, res.at(i),
                changes_jet_category
            );

            MEM::MEMOutput res_sig;
            MEM::MEMOutput res_bkg;

            vector<MEM::PSVar::PSVar> vars_to_integ;
            vector<MEM::PSVar::PSVar> vars_to_marginalize;
            auto fstate = MEM::FinalState::LH;
            switch (newhypo) {
            case SL_1W2H2T: {
                vars_to_marginalize = {MEM::PSVar::cos_qbar1, MEM::PSVar::phi_qbar1};
                break;
            }
            case SL_0W2H2T: {
                vars_to_marginalize = {MEM::PSVar::cos_q1, MEM::PSVar::phi_q1, MEM::PSVar::cos_qbar1, MEM::PSVar::phi_qbar1};
                break;
            }
            case DL_0W2H2T: {
                fstate = MEM::FinalState::LL;
                break;
            }
            default: {
                break;
            }
            }


            // float btagWP = 0;
            // if (string(this->btag_prefix) == "btagDeepCSV_") btagWP = 0.4941;
            // else if (string(this->btag_prefix) == "btagCSV_") btagWP = 0.8484;
            // else if (string(this->btag_prefix) == "btagDeepFlav_") btagWP = 0.2770;

            //Last check if event has enough jets
            int njets = selectedJetP4.at(i).size();
            int nleps = selectedLeptonP4.size();
            int nbjets = 0;
            for (unsigned int ij=0; ij < selectedJetP4.at(i).size(); ij++) {
                if (selectedJetCSV.at(i).at(ij)>btagWP) {
                    nbjets++;
                }
            }

            bool recalculated = 0;
            //Calculate MEM if event passes selection
            if(((nleps==2 && njets>=4 && nbjets>=3) || (nleps==1 && ((njets>=4 && nbjets>=4) || (njets>=6 && nbjets>=3))))) {
                recalculated = 1;
                res_sig = integrand->run(fstate, MEM::Hypothesis::TTH, vars_to_marginalize, vars_to_integ, ncalls);
                res_bkg = integrand->run(fstate, MEM::Hypothesis::TTBB, vars_to_marginalize, vars_to_integ, ncalls);
            }
            else {
                res_sig = MEM::MEMOutput();
                res_bkg = MEM::MEMOutput();
            }

            for(auto* o : objs) {
                delete o;
            }
            objs.clear();

            integrand->next_event();
            
            res.at(i).p_sig = res_sig.p;
            res.at(i).p_bkg = res_bkg.p;
            res.at(i).p_err_sig = res_sig.p_err;
            res.at(i).p_err_bkg = res_bkg.p_err;
            res.at(i).n_perm_sig = res_sig.num_perm;
            res.at(i).n_perm_bkg = res_bkg.num_perm;
            if (res.at(i).p_sig > 0) res.at(i).p = res.at(i).p_sig / (res.at(i).p_sig + mem_weight*res.at(i).p_bkg);
            else res.at(i).p = 0;

            //Get MEM output for the factorized JES sources
            if (selectedJetVariation.size() > 0) {
                const auto num_var = selectedJetVariation.size() - 1;

                for (unsigned int nvar = 0; nvar < num_var; nvar++) {
                    if (recalculated == 1) {
                        if (res_sig.variated.at(nvar)>0) res.at(i).p_variated.push_back(res_sig.variated.at(nvar) / (res_sig.variated.at(nvar) + mem_weight * res_bkg.variated.at(nvar)));
                        else res.at(i).p_variated.push_back(0);
                    }
                    else res.at(i).p_variated.push_back(0);
                }
            }
        }
    }
    return res;
}

MEM::Object* MEMClassifier::make_jet(double pt, double eta, double phi, double mass, double istagged, double csv, bool is_subjet) const {
    TLorentzVector lv;
    lv.SetPtEtaPhiM(pt, eta, phi, mass);
    MEM::Object* obj = new MEM::Object(lv, MEM::ObjectType::Jet);
    obj->addObs( MEM::Observable::BTAG, istagged); // 0 - jet is assumed to be from a light quark, 1 - a b quark
    obj->addObs( MEM::Observable::CSV, csv); //b-tagger
    //obj->addObs( MEM::Observable::PDGID, 0);  // NB: must NOT be specified!
    // attach the transfer functions corresponding to the jet

    if (is_subjet) {
        obj->addTransferFunction(MEM::TFType::bReco, getTransferFunction("sjb", lv.Eta()));
        obj->addTransferFunction(MEM::TFType::qReco, getTransferFunction("sjl", lv.Eta()));
    } else {
        obj->addTransferFunction(MEM::TFType::bReco, getTransferFunction("b", lv.Eta()));
        obj->addTransferFunction(MEM::TFType::qReco, getTransferFunction("l", lv.Eta()));
    }
    return obj;
}

MEM::Object* MEMClassifier::make_lepton(double pt, double eta, double phi, double mass, double charge) const {
    TLorentzVector lv;
    lv.SetPtEtaPhiM(pt, eta, phi, mass);
    MEM::Object* obj = new MEM::Object(lv, MEM::ObjectType::Lepton );
    obj->addObs( MEM::Observable::CHARGE, charge);
    return obj;
}

// Returns the transfer function corresponding to a jet flavour and eta
TF1* MEMClassifier::getTransferFunction(const char* flavour, double eta) const {
    int etabin = 0;
    if (std::abs(eta) > 1.0) {
        etabin = 1;
    }
    stringstream ss;
    ss << "tf_" << flavour << "_etabin" << etabin;
    const char* fname = ss.str().c_str();
    TF1* tf = (TF1*)(transfers->Get(fname));
    if (tf == nullptr) {
        cerr << "could not get transfer function " << fname << endl;
        cerr << flush;
        throw std::exception();
    }
    tf->SetNpx(10000);
    tf->SetRange(0, 500);
    return tf;
}

MEMClassifier::MEMClassifier(int verbosity, const char* _btag_prefix, const char* _era) : btag_prefix(_btag_prefix), era(_era) {

    TString era_ = "2018";
    if (TString(era).Contains("2018")){
        era_ = "2018";
        if (string(this->btag_prefix) == "btagDeepCSV_") btagWP = 0.4184;
        else if (string(this->btag_prefix) == "btagCSV_") btagWP = 0.8484; // careful no WP recommended here
        else if (string(this->btag_prefix) == "btagDeepFlav_") btagWP = 0.2770;
    }
    else if (TString(era).Contains("2017")){
        era_ = "2017";
        if (string(this->btag_prefix) == "btagDeepCSV_") btagWP = 0.4941;
        else if (string(this->btag_prefix) == "btagCSV_") btagWP = 0.8838;
        else if (string(this->btag_prefix) == "btagDeepFlav_") btagWP = 0.3033;
    }
    else if (TString(era).Contains("2016")){
        era_ = "2016";
        if (string(this->btag_prefix) == "btagDeepCSV_") btagWP = 0.6321;
        else if (string(this->btag_prefix) == "btagCSV_") btagWP = 0.8838;
        else if (string(this->btag_prefix) == "btagDeepFlav_") btagWP = 0.3093;
    }

    const string cmssw_path(std::getenv("CMSSW_BASE"));

    const string transfers_path = (
                                    string("file://") +
                                    cmssw_path +
                                    string("/src/TTH/CommonClassifier/data/transfer_") + 
                                    string(era_) + 
                                    string("_deepFlavour.root")
                                   ).c_str();

    const string btagfile_path = (
                                    string("file://") +
                                    cmssw_path +
                                    string("/src/TTH/CommonClassifier/data/btag_pdfs_")
                                    + string(era_) +
                                    string("_deepFlavour.root") 
                                 ).c_str();
    cout << "MEMClassifier: Loaded " << btagfile_path << endl;

    transfers = new TFile(transfers_path.c_str());
    assert(transfers != nullptr);

    btagfile = new TFile(btagfile_path.c_str());
    assert(btagfile != nullptr);

    cfg.defaultCfg();
    cfg.transfer_function_method = MEM::TFMethod::External;
    //Transfer functions for jet reconstruction efficiency
    cfg.set_tf_global(MEM::TFType::bLost, 0, getTransferFunction("beff", 0.0));
    cfg.set_tf_global(MEM::TFType::bLost, 1, getTransferFunction("beff", 2.0));
    cfg.set_tf_global(MEM::TFType::qLost, 0, getTransferFunction("leff", 0.0));
    cfg.set_tf_global(MEM::TFType::qLost, 1, getTransferFunction("leff", 2.0));
    cfg.add_distribution_global(MEM::DistributionType::DistributionType::csv_b, GetBTagPDF(this->btag_prefix, "b"));
    cfg.add_distribution_global(MEM::DistributionType::DistributionType::csv_c, GetBTagPDF(this->btag_prefix, "c"));
    cfg.add_distribution_global(MEM::DistributionType::DistributionType::csv_l, GetBTagPDF(this->btag_prefix, "l"));
    cfg.perm_pruning.push_back(MEM::Permutations::BTagged);
    cfg.perm_pruning.push_back(MEM::Permutations::QUntagged);
    cfg.perm_pruning.push_back(MEM::Permutations::QQbarBBbarSymmetry);
    cfg.num_jet_variations = 110;

    integrand = new MEM::Integrand(
        verbosity,
        cfg
    );

    blr = new MEM::JetLikelihood(); // TTH/MEIntegratorStandalone/interface/JetLikelihood.h
}

MEMClassifier::MEMClassifier( ) : MEMClassifier(0, "btagDeepFlavB_", "2018") {}

TH3D* MEMClassifier::GetBTagPDF(const char* prefix, const char* flavour) {
    assert(btagfile != nullptr);
    TH3D* ret = nullptr;
    ret = (TH3D*)(btagfile->Get((string(prefix)+string(flavour)+string("_pt_eta")).c_str()));
    assert(ret != nullptr);
    return ret;
}

double MEMClassifier::GetJetBProbability(const char* prefix, const char* flavour, double pt, double eta, double bdisc) {
    TH3D* h = GetBTagPDF(prefix, flavour);
    const int i = h->FindBin(pt, std::abs(eta), bdisc);
    const double c = h->GetBinContent(i);
    delete h;
    return c;
}

MEM::JetProbability MEMClassifier::GetJetBProbabilities(
    const TLorentzVector& p4, double bdisc
) {
    MEM::JetProbability jp;
    jp.setProbability(MEM::JetInterpretation::b, GetJetBProbability(this->btag_prefix, "b", p4.Pt(), p4.Eta(), bdisc));
    jp.setProbability(MEM::JetInterpretation::c, GetJetBProbability(this->btag_prefix, "c", p4.Pt(), p4.Eta(), bdisc));
    jp.setProbability(MEM::JetInterpretation::l, GetJetBProbability(this->btag_prefix, "l", p4.Pt(), p4.Eta(), bdisc));
    return jp; // setProbability in JetLikellihood.h
}

double MEMClassifier::GetBTagLikelihoodRatio(
    const std::vector<TLorentzVector>& selectedJetP4,
    const std::vector<double>& selectedJetCSV,
    std::vector<unsigned int>& out_best_perm,
    double& out_P_4b,
    double& out_P_2b
) {
//    assert(selectedJetP4.size() >= 4);
    assert(selectedJetP4.size() >= 6);

    for (unsigned int ij=0; ij < min(selectedJetP4.size(), numMaxJetsBLR); ij++) {
        blr->push_back_object(GetJetBProbabilities(selectedJetP4[ij], selectedJetCSV[ij]));
    }


    std::vector<unsigned int> best_perm_4b;
    std::vector<unsigned int> best_perm_2b;
//    double P_4b = blr->calcProbability(MEM::JetInterpretation::b, MEM::JetInterpretation::l, 4, best_perm_4b);
//    double P_2b = blr->calcProbability(MEM::JetInterpretation::b, MEM::JetInterpretation::l, 2, best_perm_2b);
    
    double P_4b = blr->calcProbability(MEM::JetInterpretation::b, MEM::JetInterpretation::l, 6, best_perm_4b);
    double P_2b = blr->calcProbability(MEM::JetInterpretation::b, MEM::JetInterpretation::l, 4, best_perm_2b);
    
    out_best_perm = best_perm_4b;
    blr->next_event();
    out_P_4b = P_4b;
    out_P_2b = P_2b;
    return P_4b / (P_4b + P_2b);
}

MEMClassifier::~MEMClassifier() {
    delete integrand;
    transfers->Close();
    delete blr;
    btagfile->Close();
}

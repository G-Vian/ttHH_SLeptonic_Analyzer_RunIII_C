// CorrectionsManager.cc
#include "CorrectionsManager.h"
#include <fstream>
#include <iostream>
#include <TRandom3.h>

// static flag for verbose logging:
static constexpr bool kVerbose = false;

CorrectionsManager::CorrectionsManager(const std::string& runYear,
                                       const std::string& dataEra,
                                       bool isData)
  : runYear_(runYear)
  , dataEra_(dataEra)
  , isData_(isData)
{
    loadJME_();          // always load MC JEC/JER; Data only if isData_
    loadPU_();           // always load PU (both MC & Data)
    loadBTag_();         // MC-only b-tag SF (or Data if you wish)
    loadGoldenJSON_();   // only Data
}

//--------------------------------------------------------------------------------------------------
// 1) Jet / MET corrections (JEC + JER)
//--------------------------------------------------------------------------------------------------
void CorrectionsManager::loadJME_() {
    const std::string base = "jsonpog-integration/POG/JME/" + runYear_ + "/";
    const std::string file = base + "jet_jerc.json.gz";
    if (kVerbose) std::cout << "[loadJME] Loading JME from " << file << "\n";

    auto cset = correction::CorrectionSet::from_file(file);
    if (kVerbose) {
        // compound()가 반환하는 map의 key들(= correction 이름)을 찍어봅니다.
        std::cout << "[loadJME] available corrections:\n";
        for (const auto &kv : cset->compound()) {
            std::cout << "  - " << kv.first << "\n";
        }
    }

    // 1) MC key
    std::string key_mc, key_res, key_sf;
    if (runYear_ == "2016PreVFP_UL" || runYear_ == "2016PostVFP_UL") {
        key_mc  = "Summer19UL16_V7_MC_L1L2L3Res_AK4PFchs";
        key_res = "Summer19UL16_JRV2_MC_PtResolution_AK4PFchs";
        key_sf  = "Summer19UL16_JRV2_MC_ScaleFactor_AK4PFchs";
    }
    else if (runYear_ == "2017_UL") {
        key_mc  = "Summer19UL17_V5_MC_L1L2L3Res_AK4PFchs";
        key_res = "Summer19UL17_JRV2_MC_PtResolution_AK4PFchs";
        key_sf  = "Summer19UL17_JRV2_MC_ScaleFactor_AK4PFchs";
    }
    else if (runYear_ == "2018_UL") {
        key_mc  = "Summer19UL18_V5_MC_L1L2L3Res_AK4PFchs";
        key_res = "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs";
        key_sf  = "Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs";
    }
    else {
        throw std::runtime_error("Unsupported runYear in loadJME_: " + runYear_);
    }


    // 2) 항상 MC용 JEC/JER 로드
    jec_MC_  = cset->compound().at(key_mc);
    jerRes_  = cset->at(key_res);
    jerSF_   = cset->at(key_sf);

    if (kVerbose) {
        std::cout
            << "  -> JEC_MC key  : " << key_mc  << "\n"
            << "  -> JER_Res key : " << key_res << "\n"
            << "  -> JER_SF key  : " << key_sf  << "\n";
    }

    // 3) 데이터면 Data용 JEC도 로드
    if (isData_) {
        std::string key_data;
        if      (runYear_ == "2016PreVFP_UL" || runYear_ == "2016PostVFP_UL")
            key_data = "Summer19UL16_Run" + dataEra_ + "_V7_DATA_L1L2L3Res_AK4PFchs";
        else if (runYear_ == "2017_UL")
            key_data = "Summer19UL17_Run" + dataEra_ + "_V5_DATA_L1L2L3Res_AK4PFchs";
        else /* 2018_UL */
            key_data = "Summer19UL18_Run" + dataEra_ + "_V5_DATA_L1L2L3Res_AK4PFchs";

        jec_Data_ = cset->compound().at(key_data);
        if (kVerbose) std::cout << "  -> JEC_Data key: " << key_data << "\n";
    }
}

//--------------------------------------------------------------------------------------------------
// 2) Pileup weights
//--------------------------------------------------------------------------------------------------
void CorrectionsManager::loadPU_() {
    const std::string file = "jsonpog-integration/POG/LUM/" + runYear_ + "/puWeights.json.gz";
    if (kVerbose) std::cout << "[loadPU] Loading PU from " << file << "\n";
    auto set = correction::CorrectionSet::from_file(file);

    // determine two‐digit year code
    std::string yearCode =
          (runYear_ == "2016PreVFP_UL"  || runYear_ == "2016PostVFP_UL") ? "16"
        : (runYear_ == "2017_UL") ? "17"
        : (runYear_ == "2018_UL") ? "18"
        : throw std::runtime_error("Unsupported runYear in loadPU_: " + runYear_);

    std::string key = "Collisions" + yearCode + "_UltraLegacy_goldenJSON";
    if (kVerbose) std::cout << "  -> PU key    : " << key << "\n";
    puCorr_ = set->at(key);
}

//--------------------------------------------------------------------------------------------------
// 3) b-tag SF  (typically only for MC)
//--------------------------------------------------------------------------------------------------
void CorrectionsManager::loadBTag_() {
    if (!isData_) {
        const std::string file = "jsonpog-integration/POG/BTV/" + runYear_ + "/btagging.json.gz";
        if (kVerbose) std::cout << "[loadBTag] Loading BTag from " << file << "\n";
        auto set = correction::CorrectionSet::from_file(file);

        const char* key = "deepJet_shape";
        btagCorr_ = set->at(key);
        if (kVerbose) std::cout << "  -> BTag key  : " << key << "\n";
    }
}

//--------------------------------------------------------------------------------------------------
// 4) Golden JSON  (only for Data)
//--------------------------------------------------------------------------------------------------
void CorrectionsManager::loadGoldenJSON_() {
    if (!isData_) return;

    // build path
    std::string suffix =
          (runYear_ == "2017_UL")          ? "UL2017_Collisions17"
        : (runYear_ == "2018_UL")          ? "UL2018_Collisions18"
        : (runYear_ == "2016PreVFP_UL")    ? "UL2016preVFP_Collisions16"
        : /*runYear_=="2016PostVFP_UL"*/     "UL2016postVFP_Collisions16";
    std::string txt = "jsonpog-integration/USER/GoldenJson/"
                    + runYear_ + "/Cert_294927-306462_13TeV_"
                    + suffix + "_GoldenJSON.txt";

    if (kVerbose) std::cout << "[loadGoldenJSON] Loading from " << txt << "\n";
    std::ifstream in(txt);
    if (!in) {
        std::cerr << "[loadGoldenJSON] ERROR opening " << txt << "\n";
        return;
    }
    nlohmann::json j; in >> j;
    for (auto& it : j.items()) {
        int run = std::stoi(it.key());
        for (auto& rng : it.value()) {
            goldenMask_[run].emplace_back(rng[0], rng[1]);
            if (kVerbose)
                std::cout << "  -> JSON run " << run
                          << " [" << rng[0] << "," << rng[1] << "]\n";
        }
    }
}

//--------------------------------------------------------------------------------------------------
//  API implementations
//--------------------------------------------------------------------------------------------------
double CorrectionsManager::getPUWeight(double nTrueInt,
                                       const std::string& var) const
{
    if (kVerbose) std::cout << "[getPUWeight] nTrueInt=" << nTrueInt
                            << " var=" << var << "\n";
    return puCorr_->evaluate({nTrueInt, var});
}

//double CorrectionsManager::getJEC(int run,
double CorrectionsManager::getJEC(
                                  double eta,
                                  double raw_pt,
                                  double area,
				  double rho) const
{
//    if (kVerbose) std::cout << "[getJEC] run="<<run
    if (kVerbose) std::cout << "[getJEC] -----  "
                            << " isData="<<isData_
                            << " eta="<<eta
                            << " raw_pt="<<raw_pt
                            << " area="<<area
			    << " rho="<<rho<<"\n";
    auto &corr = isData_ ? jec_Data_ : jec_MC_;
    return corr->evaluate({ area, eta, raw_pt, rho});
}

double CorrectionsManager::smearJER(double corr_pt,
                                    double gen_pt,
                                    double eta,
                                    double rho,
                                    const std::string& syst) const
{

    if(isData_) return corr_pt;

    double sf = 0, res = 0;
    // 1) JER scale factor
    try {
      // 1) JER scale factor: now pass both inputs (eta, systematic)
      if (kVerbose) std::cout << "[smearJER] calling jerSF_->evaluate({eta, syst}) -> {"
                         << eta << ", " << syst << "}\n";
      sf = jerSF_->evaluate({eta, syst});
      if (kVerbose) std::cout << "[smearJER] jerSF returned " << sf << "\n";
    }
    catch (const std::exception &e) {
      std::cerr << "[smearJER] ERROR in jerSF_->evaluate: " << e.what() << "\n";
      throw;
    }

    // 2) JER resolution
    try {
      if (kVerbose) std::cout << "[smearJER] calling jerRes_->evaluate({eta, corr_pt, rho})\n";
      res = jerRes_->evaluate({eta, corr_pt, rho});
      if (kVerbose) std::cout << "[smearJER] jerRes returned " << res << "\n";
    }
    catch (const std::exception &e) {
      std::cerr << "[smearJER] ERROR in jerRes_->evaluate: " << e.what() << "\n";
      throw;
    }

    if (kVerbose) std::cout << "[smearJER] corr_pt="<<corr_pt
                            << " gen_pt="<<gen_pt
                            << " eta="<<eta
                            << " rho="<<rho
                            << " sf="<<sf
                            << " res="<<res<<"\n";
    if (gen_pt >= 0) {
        return std::max(0.0, gen_pt + sf*(corr_pt - gen_pt));
    } else {
        double gaus   = TRandom3().Gaus(0,1);
        double factor = 1.0 + std::sqrt(sf*sf - 1.0)*res*gaus;
        return std::max(0.0, corr_pt * factor);
    }
}

double CorrectionsManager::getBTagSF(int hf,
                                     double eta,
                                     double pt,
                                     double disc,
                                     const std::string& sys) const
{
    if (!isData_) return 1.0;  // no SF applied for Data
    if (hf!=5 && hf!=4) hf=0;  // light
    if (kVerbose) std::cout << "[getBTagSF] hf="<<hf
                            << " eta="<<eta
                            << " pt="<<pt
                            << " disc="<<disc
                            << " sys="<<sys<<"\n";
    return btagCorr_->evaluate({sys, hf, std::fabs(eta), pt, disc});
}

bool CorrectionsManager::passGoldenJSON(int run, int lumi) const {
    if (!isData_) return true; // MC always “passes”
    if (kVerbose) std::cout << "[passGoldenJSON] run="<<run
                            << " lumi="<<lumi<<"\n";
    auto it = goldenMask_.find(run);
    if (it==goldenMask_.end()) return false;
    for (auto &p : it->second) {
        if (lumi>=p.first && lumi<=p.second) return true;
    }
    return false;
}

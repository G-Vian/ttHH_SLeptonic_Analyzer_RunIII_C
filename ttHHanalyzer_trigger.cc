#include "tnm.h"
#include <cmath> 
#include <algorithm>
#include <vector>
#include <map>
#include "TVector3.h"
#include "ttHHanalyzer_trigger.h"
#include <iostream>
#include "correction.h"// 🔧 


using namespace std;

#include "Logger.h"  // 🔧 
using namespace Logger; // 🔧 

void ttHHanalyzer::performAnalysis(){
    loop(noSys, false);
    /*    getbJetEffMap();
    initHistograms(kJES, false);
    initTree(kJES, false);
    loop(kJES, false);
    initHistograms(kJES, true);
    initTree(kJES, true);
    loop(kJES, true);

    initHistograms(kJER, false);
    initTree(kJER, false);
    loop(kJER, false);
    initHistograms(kJER, true);
    initTree(kJER, true);
    loop(kJER, true);

    initHistograms(kbTag, false);
    initTree(kbTag, false);
    loop(kbTag, false);
    initHistograms(kbTag, true);
    initTree(kbTag, true);
    loop(kbTag, true); */

}
void ttHHanalyzer::loop(sysName sysType, bool up) {
    int nevents = _ev->size();
    ////int nevents = 1000;

    std::cout << std::endl;
    std::cout << "[WARNING] This analyzer commented out [ \"WTF\" log ] in the header, Please check if you want!!!" << std::endl;

    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << "Before start, Let's check the analysis information" << std::endl;
    std::cout << "Run Year    ----> [  " << _year << "  ]" << std::endl;
    std::cout << "Data or MC  ----> [  " << _DataOrMC << "  ]" << std::endl;
    std::cout << "Sample Name ----> [  " << _sampleName << "  ]" << std::endl;

    std::string checklist = "[ tnm.cc ] & [ analyzer header ] & [ main ] & [ analyzer constructor ]";
    bool exitFlag = false;

    if (_year == "nothing") {
        std::cout << "[ERROR] year is not defined, Please check the " << checklist << std::endl;
        exitFlag = true;
    }
    if (_DataOrMC == "nothing") {
        std::cout << "[ERROR] Whether Data or MC is not defined, Please check the " << checklist << std::endl;
        exitFlag = true;
    }
    if (_sampleName == "nothing") {
        std::cout << "[ERROR] SampleName is not defined, Please check the " << checklist << std::endl;
        exitFlag = true;
    }

    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    if (exitFlag) std::exit(EXIT_FAILURE);

    std::string analysisInfo = _year + ", " + _DataOrMC + ", " + _sampleName;

    for (int entry = 0; entry < nevents; entry++) {
        event *currentEvent = new event;
        ////std::cout << "Processed events: " << entry << std::endl;
        _ev->read(entry);  // read an event into event buffer
        process(currentEvent, sysType, up);

        if (entry % 1000 == 0) {
            std::cout << "[INFO] Processed events of " << analysisInfo << ": " << entry << std::endl;
            currentEvent->summarize();
        }

        events.push_back(currentEvent);
    }
    // events.back()->summarize();
// 🔧 
   if(debugCorrections) std::cout<<"debug : Before [ writeHistos() ]"<<std::endl;
    writeHistos();
    if(debugCorrections) std::cout<<"debug : After writeHistos() & Before writeTree()"<<std::endl;
    writeTree();
    if(debugCorrections) std::cout<<"debug : After writeTree() & Before hcutFlow()"<<std::endl;
// 🔧 
	
 //   writeHistos(); // 🔧 
 //   writeTree(); // 🔧 


    for (const auto &x : cutflow) {
        std::cout << x.first  // string (key)
                  << ": " 
                  << x.second // string's value 
                  << std::endl;
    } 


    hCutFlow->Write(); // 🔧 
    if(debugCorrections) std::cout<<"debug : After hCutFlow() & Before hCutFlow_w()"<<std::endl; // 🔧 
    hCutFlow_w->Write();// 🔧 
    if(debugCorrections) std::cout<<"Problem : After hCutFlow_w()"<<std::endl; // 🔧 
}

void ttHHanalyzer::createObjects(event * thisEvent, sysName sysType, bool up){
cutflow["noCut"]+=1;
hCutFlow->Fill("noCut",1);
hCutFlow_w->Fill("noCut",_weight);
	
    _ev->fillObjects();

// This trigger paths are for the SL channel!
    if(_year == "2022" or _year == "2022EE" or _year == "2023" or _year == "2023B" or _year == "2024" ){

	if(_DataOrMC == "MC" or _DataOrMC == "Data"){thisEvent->setFilter(_ev->Flag_goodVertices ||				 
		  					  _ev->Flag_globalSuperTightHalo2016Filter ||		
		                                          _ev->Flag_EcalDeadCellTriggerPrimitiveFilter ||	
		                                          _ev->Flag_BadPFMuonFilter ||
                                                        _ev->Flag_BadPFMuonDzFilter ||	
                                                      _ev->Flag_hfNoisyHitsFilter ||
		                        		  _ev->Flag_eeBadScFilter ||				 
		                                          _ev->Flag_ecalBadCalibFilter);

							 
	    thisEvent->setTrigger(
				  _ev->HLT_Ele30_WPTight_Gsf ||
				  _ev->HLT_IsoMu24);
	}
    }
    
    
    thisEvent->setPV(_ev->PV_npvsGood);
    std::vector<eventBuffer::GenPart_s> genPart = _ev->GenPart;      
    std::vector<eventBuffer::Jet_s> jet = _ev->Jet;
    std::vector<eventBuffer::Muon_s> muonT = _ev->Muon;
    std::vector<eventBuffer::Electron_s> ele = _ev->Electron;
    std::vector<eventBuffer::FatJet_s> boostedJet = _ev->FatJet;
    objectGenPart * currentGenPart; 
    objectBoostedJet * currentBoostedJet;
    objectJet * currentJet;
    objectLep * currentMuon;
    objectLep * currentEle;
    int nVetoMuons = 0, nVetoEle = 0;
    objectMET * MET = new objectMET(_ev->PuppiMET_pt, 0, _ev->PuppiMET_phi, 0);
    float e = 1., es  = 1., pe = 1., pes = 1.;
    float me = 1., mes = 1., pme = 1.,  pmes = 1.;   
    thisEvent->setMET(MET);


    for(int i=0; i < boostedJet.size(); i++){
       	currentBoostedJet = new objectBoostedJet(boostedJet[i].pt, boostedJet[i].eta, boostedJet[i].phi, boostedJet[i].mass);
	currentBoostedJet->softDropMass = boostedJet[i].msoftdrop;
	if(currentBoostedJet->getp4()->Pt() > cut["boostedJetPt"] && fabs(currentBoostedJet->getp4()->Eta()) < fabs(cut["boostedJetEta"])){
	    //	    if((boostedJet[i].jetId & 4) == true){  	     
	    thisEvent->selectBoostedJet(currentBoostedJet);	
	    if(currentBoostedJet->getp4()->Pt() > cut["hadHiggsPt"]){
		if(boostedJet[i].particleNet_HbbvsQCD > cut["bTagDisc"]){
		    thisEvent->selectHadronicHiggs(currentBoostedJet);
		}
		//	}
	    }
	}
    }
    

    bool thereIsALeadLepton = false;
	
    for(int i = 0; i < muonT.size(); i++){
	if(fabs(muonT[i].eta) < cut["muonEta"] && muonT[i].tightId == true && muonT[i].pfRelIso04_all  < cut["muonIso"]){
	//if(fabs(muonT[i].eta) < cut["muonEta"] && muonT[i].mvaTTH > 0.15 && muonT[i].pfRelIso04_all  < cut["muonIso"]){
	    if(muonT[i].pt > cut["leadMuonPt"]){
		thereIsALeadLepton = true;
		break;
	    }
	}
    }
    if(!thereIsALeadLepton){
	for(int i = 0; i < ele.size(); i++){
	    if(fabs(ele[i].deltaEtaSC + ele[i].eta) < 1.4442 || fabs(ele[i].deltaEtaSC + ele[i].eta) > 1.5660){  //Electrons tracked neither in the barrel nor in the endcap are discarded.
		if(fabs(ele[i].eta) < cut["eleEta"] && ele[i].mvaIso_WP90 == true) { // && ele[i].pfRelIso03_all  < cut["eleIso"]){ 
		    if(ele[i].pt > cut["leadElePt"]){
			thereIsALeadLepton = true;
			break;
		    }
		}
	    }
	}
    }
    
    if(thereIsALeadLepton){ //we can add all leptons passing to the sublead selection to our containers
	for(int i = 0; i < muonT.size(); i++){
	    if(fabs(muonT[i].eta) < cut["muonEta"] && muonT[i].tightId == true && muonT[i].pfRelIso04_all < cut["muonIso"]){
	    //	    if(fabs(muonT[i].eta) < cut["muonEta"] && muonT[i].mvaTTH > 0.15 && muonT[i].pfRelIso04_all  < cut["muonIso"]){	
		if(muonT[i].pt > cut["subLeadMuonPt"]){
		    currentMuon = new objectLep(muonT[i].pt, muonT[i].eta, muonT[i].phi, 0.);
		    currentMuon->charge = muonT[i].charge;
		    currentMuon->miniPFRelIso = muonT[i].miniPFRelIso_all;
		    currentMuon->pfRelIso04 = muonT[i].pfRelIso04_all;
		    thisEvent->selectMuon(currentMuon);
		} // else if (muonT[i].pt > cut["vetoLepPt"]){	 
		// nVetoMuons++;
		//}
	    }
	}
	for(int i = 0; i < ele.size(); i++){
	    if(fabs(ele[i].deltaEtaSC + ele[i].eta) < 1.4442 || fabs(ele[i].deltaEtaSC + ele[i].eta) > 1.5660){  //Electrons tracked neither in the barrel nor in the endcap are discarded.
		if(fabs(ele[i].eta) < cut["eleEta"] && ele[i].mvaIso_WP90 == true) { // && ele[i].pfRelIso03_all  < cut["eleIso"]){ 
		    if(ele[i].pt > cut["subLeadElePt"]){
			currentEle = new objectLep(ele[i].pt, ele[i].eta, ele[i].phi, 0.);	 
			currentEle->charge = ele[i].charge;
			currentEle->miniPFRelIso = ele[i].miniPFRelIso_all;
			currentEle->pfRelIso03 = ele[i].pfRelIso03_all;
			thisEvent->selectEle(currentEle);
		    }// else if(ele[i].pt > cut["vetoLepPt"]){	 
		    //	nVetoEle++;
		    //}
		}
	    }
	}
    }
    thisEvent->orderLeptons();
// 🔧
    float dR = 0., deltaEta = 0., deltaPhi = 0.;
    for (int i = 0; i < (int)jet.size(); ++i) {
        const auto& jetRaw = jet[i];
    
        // 1) build the ORIGINAL jet (nanoAOD JEC already applied)
        objectJet* origJet = new objectJet(
            jetRaw.pt,
            jetRaw.eta,
            jetRaw.phi,
            jetRaw.mass
        );
        origJet->bTagCSV    = jetRaw.btagDeepFlavB;
        origJet->jetID      = jetRaw.jetId;
        origJet->jetPUid    = jetRaw.puId;
        origJet->hadFlav    = jetRaw.hadronFlavour;
        origJet->partonFlav = jetRaw.partonFlavour;
    
        // === A) Recover RAW pT / mass ===
        float ntuplePt  = origJet->getp4()->Pt();      // pt_nom
        float rawFactor = jetRaw.rawFactor;            // nanoAOD‐provideds
        float rawPt     = ntuplePt * (1.0 - rawFactor);

	////float oldraw = jetRaw.pt_nom * (1.0 - rawFactor);
	////std::cout<<"!!!!! diff = "<<rawPt - oldraw<<std::endl; //이런 로직이면 차이 0
        ////float rawMass   = origJet->getp4()->M() * rawFactor;
        float rawMass   = origJet->getp4()->M() * (1.0 - rawFactor);

        // store for validation (make sure objectJet has these members!)
        origJet->origPt   = ntuplePt;
        origJet->rawPt    = rawPt;
        origJet->rawMass  = rawMass;
    
        if (debugCorrections) {
            std::cout << "[JEC/JER] ntuplePt=" << ntuplePt
                      << " rawFactor=" << rawFactor
                      << " -> rawPt=" << rawPt << "\n";
        }
    
        // === B) Get gen‐matched pT for JER ===
        float genPt = -1.f;
        if (jetRaw.genJetIdx >= 0 && jetRaw.genJetIdx < (int)genJet.size()) {
            genPt = genJet[jetRaw.genJetIdx].pt;
        }
        origJet->genMatchedPt = genPt;
    
        // === C) Apply JEC ===
        float rho = _ev->fixedGridRhoFastjetAll;
        double jecSF = corrMgr->getJEC(
//            _ev->run,                        // UInt_t run number
            origJet->getp4()->Eta(),         // η
            rawPt,                           // raw pT
            jetRaw.area,                     // jet area
            rho                              // ρ (for L1FastJet)
        );
        float ptJEC       = rawPt * jecSF;
        origJet->jecSF    = jecSF;
        origJet->ptJEC    = ptJEC;
    
        if (debugCorrections) {
            std::cout << "[JEC] run=" << _ev->run
	              << " origPt(ntuple's corrJet)=" << ntuplePt
                      << " rawPt=" << rawPt
                      << " eta=" << origJet->getp4()->Eta()
                      << " jecSF=" << jecSF
                      << " -> ptJEC=" << ptJEC << "\n";
        }

        // === D) Apply JER ===
        double smearedPt = corrMgr->smearJER(
            ptJEC,                           // JEC‐corrected pT
            genPt,                           // gen‐matched pT (-1 if none)
            origJet->getp4()->Eta(),         // η
            rho,                             // ρ
	    //"nominal"
            "nom"
        );
        double jerSF      = (ptJEC > 0. ? smearedPt/ptJEC : 1.0);
        origJet->jerSF    = jerSF;
        origJet->ptJER    = smearedPt;
    
        if (debugCorrections) {
            std::cout << "[JER] ptJEC=" << ptJEC
                      << " genPt=" << genPt
                      << " rho="   << rho
                      << " jerSF=" << jerSF
                      << " -> ptJER=" << smearedPt << "\n";
        }
        
        objectJet* newJet = new objectJet(
            smearedPt,
            jetRaw.eta,
            jetRaw.phi,
            jetRaw.mass
        );

        newJet->origPt     = ntuplePt;
        newJet->ptJEC      = ptJEC;
        newJet->rawPt      = rawPt;
        newJet->rawMass    = rawMass;
        newJet->jerSF      = jerSF;
	newJet->ptJER      = smearedPt;
        newJet->jetPUid    = jetRaw.puId;
        newJet->hadFlav    = jetRaw.hadronFlavour;
        newJet->partonFlav = jetRaw.partonFlavour;

	// === E) hand off to event ===
        thisEvent->selectJet(newJet);
	delete origJet;

    }
    thisEvent->orderJets();
// 🔧	
    // Selecting b jets from genParticle info
    for (int i = 0; i < genPart.size(); i++){
       	currentGenPart = new objectGenPart(genPart[i].pt, genPart[i].eta, genPart[i].phi, genPart[i].mass);
	currentGenPart->hasHiggsMother = false;
	currentGenPart->hasTopMother = false;

      	if (bool((abs(genPart[i].pdgId) == 5) && (genPart[i].statusFlags & 256)) == true){
	    thisEvent->selectGenPart(currentGenPart);
	    int motherInd = genPart[i].genPartIdxMother;
	    if(abs(genPart[motherInd].pdgId) == 25 && (genPart[motherInd].statusFlags & 256)) {
		currentGenPart->hasHiggsMother = true;
	    } else if(abs(genPart[motherInd].pdgId) == 6 && (genPart[motherInd].statusFlags & 256)) {
		currentGenPart->hasTopMother = true;
	    } else {
		while((abs(genPart[motherInd].pdgId) != 25 || abs(genPart[motherInd].pdgId) != 6 || !(genPart[motherInd].statusFlags & 256)) && motherInd > 1){
		    motherInd = genPart[motherInd].genPartIdxMother;
		    if(abs(genPart[motherInd].pdgId) == 25 && (genPart[motherInd].statusFlags & 256)) {
			currentGenPart->hasHiggsMother = true;
		    } else if(abs(genPart[motherInd].pdgId) == 6 && (genPart[motherInd].statusFlags & 256)) {
			currentGenPart->hasTopMother = true;
		    }
		}
	    }
	    thisEvent->selectGenPart(currentGenPart);
	}
    }
}



bool ttHHanalyzer::selectObjects(event *thisEvent){
    //    std::cout << "bjet CSV: " << thisEvent->getSelbJets()->at(0)->bTagCSV << std::endl;
	
//    cutflow["noCut"]+=1;
//    hCutFlow->Fill("noCut",1);
//    hCutFlow_w->Fill("noCut",_weight);
//Checks if the event was accepted by the triggers
    if(cut["trigger"] > 0 && thisEvent->getTriggerAccept() == false){
	return false;
    }

    cutflow["nTrigger"]+=1;
    hCutFlow->Fill("nTrigger",1);
    hCutFlow_w->Fill("nTrigger",_weight);
	
//Checks if the event passes the filters
    if(cut["filter"] > 0 && thisEvent->getMETFilter() == false){
	return false;
    }

    cutflow["nFilter"]+=1;
    hCutFlow->Fill("nFilter",1);
    hCutFlow_w->Fill("nFilter",_weight);
//Checks if the event has a valid primary vertex 
    if(cut["pv"] < 0 && thisEvent->getPVvalue() == false){
	return false;
    } 

    cutflow["nPV"]+=1;
    hCutFlow->Fill("nPV",1);
    hCutFlow_w->Fill("nPV",_weight);


	
    if(!(thisEvent->getnSelJet()  >=  cut["nJets"] )){
	return false;
    }
    
    cutflow["njets>5"]+=1;
    hCutFlow->Fill("njets>5",1);
    hCutFlow_w->Fill("njets>5",_weight);
    
	
    if(!(thisEvent->getnbJet() >=  cut["nbJets"])){
	return false;
    }
    
    cutflow["nbjets>4"]+=1;
    hCutFlow->Fill("nbjets>4",1);
    hCutFlow_w->Fill("nbjets>4",_weight);
   
    if(!(thisEvent->getnSelLepton() == cut["nLeptons"])){
	return false;
    }
    
    cutflow["nlepton==1"]+=1;
    hCutFlow->Fill("nlepton==1",1);
    hCutFlow_w->Fill("nlepton==1",_weight);
	

	


    thisEvent->getStatsComb(thisEvent->getSelJets(), thisEvent->getSelLeptons(), ljetStat);
    thisEvent->getStatsComb(thisEvent->getSelbJets(), thisEvent->getSelLeptons(), lbjetStat);
 

    
/*    if (thisEvent->getSelLeptons()->size() == 2) {
        if (thisEvent->getSelLeptons()->at(0)->charge == thisEvent->getSelLeptons()->at(1)->charge) {
            return false;
         }
     }*/


//    cutflow["nOpositeChargedLep"]+=1;
//    hCutFlow->Fill("nOpositeChargedLep",1);
//    hCutFlow_w->Fill("nOpositeChargedLep",_weight);


    //    if(!(thisEvent->getnVetoLepton()  == cut["nVetoLeptons"])){
    //	return false;
    // }

    
   /* if(thisEvent->getnSelMuon()  == cut["nLeptons"]){
    	if(!((thisEvent->getSelMuonsMass() > 20) && (thisEvent->getSelMuonsMass() < 76 || thisEvent->getSelMuonsMass() > 106))){
	    return false;
	}
    }
    
    if(thisEvent->getnSelElectron()  == cut["nLeptons"]){
    	if(!((thisEvent->getSelElectronsMass() > 20) && (thisEvent->getSelElectronsMass() < 76 || thisEvent->getSelElectronsMass() > 106))){
    	    return false;
    	}
    }
    
    cutflow["nMassCut"]+=1;
    hCutFlow->Fill("nMassCut",1);
    hCutFlow_w->Fill("nMassCut",_weight);
    /*
    
        
/*    if(thisEvent->getnSelMuon()  == cut["nLeptons"] || thisEvent->getnSelElectron()  == cut["nLeptons"]){	
    	if(!(thisEvent->getMET()->getp4()->Pt() > cut["MET"] )){
    	    return false;
    	}
    }*/


	if (thisEvent->getSelElectrons()->size() == 1 && thisEvent->getSelMuons()->size() == 0) {
	    // Case 1: one electron, no muons
	    if (thisEvent->getSelElectrons()->at(0)->getp4()->Pt() >= cut["leadElePt"] &&
	        fabs(thisEvent->getSelElectrons()->at(0)->getp4()->Eta()) <= cut["eleEta"] ) {
	        
	//	cutflow["count_elec"]+=1;	        
	    } else {
	        return false;
	    }
	} 
	// Check if there is one muon and zero electrons
	else if (thisEvent->getSelElectrons()->size() == 0 && thisEvent->getSelMuons()->size() == 1) {
	    if (thisEvent->getSelMuons()->at(0)->getp4()->Pt() >= cut["leadMuonPt"] &&
	        fabs(thisEvent->getSelMuons()->at(0)->getp4()->Eta()) <= cut["muonEta"] ) {
	            
	///	cutflow["count_muon"]+=1;	        
	    } else {
	        return false;
	    }
	} 

/*	
// Check if there is one electron and one muon
else if (thisEvent->getSelElectrons()->size() == 1 && thisEvent->getSelMuons()->size() == 1) {
    // Case 3: One electron and one muon
    auto ele = thisEvent->getSelElectrons()->at(0);
    auto mu = thisEvent->getSelMuons()->at(0);
    
    // Determine which lepton has the highest transverse momentum (pT)
    if (ele->getp4()->Pt() > mu->getp4()->Pt()) {
        // Electron is leading
        if (ele->getp4()->Pt() < cut["leadElePt"] || fabs(ele->getp4()->Eta()) > cut["eleEta"] ||
            mu->getp4()->Pt() < cut["subLeadMuonPt"] || fabs(mu->getp4()->Eta()) > cut["muonEta"]) {
            return false;
        }
    } else {
        // Muon is leading
        if (mu->getp4()->Pt() < cut["leadMuonPt"] || fabs(mu->getp4()->Eta()) > cut["muonEta"] ||
            ele->getp4()->Pt() < cut["subLeadElePt"] || fabs(ele->getp4()->Eta()) > cut["eleEta"]) {
            return false;
        }
    }
}
*/

//////////////////////////////	
    if(!(thisEvent->getMET()->getp4()->Pt() > cut["MET"] )){
        return false;
    	}
/////////////////////////////

    cutflow["MET>20"]+=1;
    hCutFlow->Fill("MET>20",1);
    hCutFlow_w->Fill("MET>20",_weight);

    /*	std::cout << x.first  // string (key)
		  << ':' 
		  << x.second // string's value 
		  << std::endl;
		  } */
	

    return true;
}


void ttHHanalyzer::motherReco(const TLorentzVector & dPar1p4,const TLorentzVector & dPar2p4, const float mother1mass, float & _minChi2,float & _bbMassMin1){
    float bbMass1, chi2;
    bbMass1 = (dPar1p4+dPar2p4).M();
    chi2 = pow((bbMass1 - mother1mass),2)/pow((dPar1p4.Pt()+dPar2p4.Pt())/2.*0.03,0.5);
    if(_minChi2 > chi2){
	_minChi2      = chi2;
	_bbMassMin1   = bbMass1;
	_bpTHiggs1    = (dPar1p4+dPar2p4).Pt();
    }
} 


void ttHHanalyzer::diMotherReco(const TLorentzVector & dPar1p4,const TLorentzVector & dPar2p4,const TLorentzVector & dPar3p4,const TLorentzVector & dPar4p4, const float mother1mass, const float  mother2mass, float & _minChi2,float & _bbMassMin1, float & _bbMassMin2){
    float bbMass1, bbMass2, chi2;
    bbMass1 = (dPar1p4+dPar2p4).M();
    bbMass2 = (dPar3p4+dPar4p4).M();
    chi2 = pow((bbMass1 - mother1mass),2)/pow((dPar1p4.Pt()+dPar2p4.Pt())/2.*0.03,0.5) + pow((bbMass2 - mother2mass),2)/pow((dPar3p4.Pt()+dPar4p4.Pt())/2.*0.03,0.5);
    if(_minChi2 > chi2){
	_minChi2      = chi2;
	_bbMassMin1   = bbMass1;
	_bbMassMin2   = bbMass2;
	_bpTHiggs1    = (dPar1p4+dPar2p4).Pt();
	_bpTHiggs2    = (dPar3p4+dPar4p4).Pt();
    }
} 

void ttHHanalyzer::analyze(event *thisEvent){

    std::vector<objectJet*>* bJetsInv = thisEvent->getSelbJets(); 
    std::vector<objectJet*>* lbJetsInv = thisEvent->getLoosebJets(); 
    std::vector<objectJet*>* jetsInv = thisEvent->getSelJets(); 

    std::vector<TVector3> vectorsJet, vectorsBjet;
    // Event Shape Calculation & genbjet matching for mother particle
    for(int k = 0; k < jetsInv->size(); k++){
	vectorsJet.push_back(jetsInv->at(k)->getp4()->Vect());
    }
    for(int m = 0; m < bJetsInv->size(); m++){
	vectorsBjet.push_back(bJetsInv->at(m)->getp4()->Vect());
	if(thisEvent->getnGenPart() < 1) continue;
	bJetsInv->at(m)->matchedtoHiggs = false;	    	
	for(auto genParticle: *thisEvent->getGenParts()){
	    float dR = bJetsInv->at(m)->getp4()->DeltaR( *genParticle->getp4());
	    if(genParticle->hasHiggsMother == true && dR < 0.8){
	      if(dR > genParticle->dRmatched) {
		  //std::cout << "this was matched to a closer particle before" << std::endl;
	      }else if( genParticle->matched){
		  //std::cout << "this was matched before" << std::endl;
	      }
	      bJetsInv->at(m)->matchedtoHiggs   = true;	
	      bJetsInv->at(m)->matchedtoHiggsdR = dR;	
	      genParticle->matched = true;
	      genParticle->dRmatched = dR;
	      break;
	    }
	} 
	//	std::cout << bJetsInv->at(m)->matchedtoHiggsb << std::endl;  
    }

    _minChi2Higgs  = cLargeValue;
    _minChi2Z      = cLargeValue;
    _minChi2HiggsZ = cLargeValue;
    _minChi2SHiggsNotMatched = cLargeValue;
    _minChi2SHiggsMatched = cLargeValue;
    _minChi2HHNotMatched = cLargeValue;
    _minChi2HHMatched = cLargeValue;
    _bbMassMinSHiggsMatched = -1;
    _bbMassMinSHiggsNotMatched = -1;
    _bbMassMinHH1Matched = -1;
    _bbMassMinHH1NotMatched = -1;
    _bbMassMinHH2Matched = -1;
    _bbMassMinHH2NotMatched = -1;

    float tempminChi2 = cLargeValue, tmpMassMin1HiggsZ = 0., tmpMassMin2HiggsZ = 0.;
    float tempMinChi2SHiggs = cLargeValue, tempMinChi2SHiggs_r = cLargeValue, tmpMassMinSHiggs = 0.;
    float tempMinChi2SHiggsMatched = cLargeValue, tempMinChi2SHiggsMatched_r = cLargeValue, tmpMassMinSHiggsMatched = 0.;
    float tempMinChi2SHiggsNotMatched = cLargeValue, tempMinChi2SHiggsNotMatched_r = cLargeValue, tmpMassMinSHiggsNotMatched = 0.;
    //extract H
    for( int ibjet1 = 0; ibjet1 < bJetsInv->size(); ibjet1++){
	tempMinChi2SHiggs = cLargeValue;
	bJetsInv->at(ibjet1)->minChiHiggsIndex = -1;
	for( int ibjet2 = 1; ibjet2 < bJetsInv->size(); ibjet2++){
	    if( ibjet1 == ibjet2) continue;	   
	    tempMinChi2SHiggs_r = tempMinChi2SHiggs;
	    motherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4(),cHiggsMass, tempMinChi2SHiggs, tmpMassMinSHiggs);
	    if(tempMinChi2SHiggs_r > tempMinChi2SHiggs){
		bJetsInv->at(ibjet1)->minChiHiggsIndex = ibjet2;
	    }
	    if(bJetsInv->at(ibjet1)->matchedtoHiggs == true && bJetsInv->at(ibjet2)->matchedtoHiggs == true){
		motherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4(),cHiggsMass, _minChi2SHiggsMatched, _bbMassMinSHiggsMatched);		
	    } else {
		motherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4(),cHiggsMass, _minChi2SHiggsNotMatched, _bbMassMinSHiggsNotMatched);		
	    }
	}
	bJetsInv->at(ibjet1)->minChiHiggs = tempMinChi2SHiggs;
    }


    // HH & ZZ reco : 4 medium b jet case

    if(thisEvent->getnbJet() >  3){
    	for( int ibjet1 = 0; ibjet1 < bJetsInv->size(); ibjet1++){
   	    for( int ibjet2 = ibjet1+1; ibjet2 < bJetsInv->size(); ibjet2++){
    		if( ibjet1 == ibjet2) continue;
    		for( int ibjet3 = 1; ibjet3 < bJetsInv->size(); ibjet3++){
    		    if(ibjet1 == ibjet3 || ibjet2 == ibjet3) continue;
    		    for( int ibjet4 = ibjet3+1; ibjet4 < bJetsInv->size(); ibjet4++){
    			if(ibjet1 == ibjet4 || ibjet2 == ibjet4 || ibjet3 == ibjet4 ) continue;		
			if(bJetsInv->at(ibjet1)->matchedtoHiggs == true && bJetsInv->at(ibjet2)->matchedtoHiggs == true && bJetsInv->at(ibjet3)->matchedtoHiggs == true && bJetsInv->at(ibjet4)->matchedtoHiggs == true){
			    diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHMatched, _bbMassMinHH1Matched, _bbMassMinHH2Matched);
			} else {
			    diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHNotMatched, _bbMassMinHH1NotMatched, _bbMassMinHH2NotMatched);
			}
			diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
				     , cHiggsMass, cHiggsMass, _minChi2Higgs, _bbMassMin1Higgs, _bbMassMin2Higgs);
			diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cZMass, _minChi2Z, _bbMassMin1Z, _bbMassMin2Z);  
			diMotherReco(*bJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() //ZH 
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cZMass, _minChi2HiggsZ, _bbMassMin1HiggsZ, _bbMassMin2HiggsZ); 
		    }
		}
	    }
	}
	// HH & ZZ reco : 3 medium + 1 loose b jet case
    } else if(thisEvent->getnbJet() == 3 && thisEvent->getnbLooseJet() > 3){
	for( int ibjet1 = 0; ibjet1 < lbJetsInv->size(); ibjet1++){
	    for( int ibjet2 = 0; ibjet2 < bJetsInv->size(); ibjet2++){
		if( lbJetsInv->at(ibjet1) == bJetsInv->at(ibjet2)) continue;
		for( int ibjet3 = 0; ibjet3 < bJetsInv->size(); ibjet3++){
		    if(lbJetsInv->at(ibjet1) == bJetsInv->at(ibjet3) || ibjet2 == ibjet3) continue;
		    for( int ibjet4 = ibjet3+1; ibjet4 < bJetsInv->size(); ibjet4++){
			if(lbJetsInv->at(ibjet1) == bJetsInv->at(ibjet4) || ibjet2 == ibjet4 || ibjet3 == ibjet4 ) continue;		
			if( bJetsInv->at(ibjet2)->matchedtoHiggs == true && bJetsInv->at(ibjet3)->matchedtoHiggs == true && bJetsInv->at(ibjet4)->matchedtoHiggs == true){
			    diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHMatched, _bbMassMinHH1Matched, _bbMassMinHH2Matched);
			} else {
			    diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
					 , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4() 
					 , cHiggsMass, cHiggsMass, _minChi2HHNotMatched, _bbMassMinHH1NotMatched, _bbMassMinHH2NotMatched);
			}

			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cHiggsMass, _minChi2Higgs, _bbMassMin1Higgs, _bbMassMin2Higgs); 
			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cZMass, _minChi2Z, _bbMassMin1Z, _bbMassMin2Z); 
			// ZH combinatorics 
			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH H(lbb)Z(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cZMass, _minChi2HiggsZ, _bbMassMin1HiggsZ, _bbMassMin2HiggsZ);  
			tempminChi2 = _minChi2HiggsZ; tmpMassMin1HiggsZ = _bbMassMin1HiggsZ; tmpMassMin2HiggsZ = _bbMassMin2HiggsZ;
			diMotherReco(*lbJetsInv->at(ibjet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH Z(lbb)H(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cHiggsMass, _minChi2HiggsZ, _bbMassMin2HiggsZ, _bbMassMin1HiggsZ); 
			// pick the lowest minChi2 for the two combinatorics
			if(tempminChi2 < _minChi2HiggsZ){
			    _minChi2HiggsZ = tempminChi2; _bbMassMin1HiggsZ = tmpMassMin1HiggsZ; _bbMassMin2HiggsZ = tmpMassMin2HiggsZ;
			}
		    }
		}
	    }
	}
    } else if(thisEvent->getnbJet() == 3){   	// HH & ZZ reco : 3 medium b jet + 1 jet case 
	for( int ijet1 = 0; ijet1 < jetsInv->size(); ijet1++){
	    for( int ibjet2 = 0; ibjet2 < bJetsInv->size(); ibjet2++){
		if( jetsInv->at(ijet1) == bJetsInv->at(ibjet2)) continue;
		for( int ibjet3 = 0; ibjet3 < bJetsInv->size(); ibjet3++){
		    if(jetsInv->at(ijet1) == bJetsInv->at(ibjet3) || ibjet2 == ibjet3) continue;
		    for( int ibjet4 = ibjet3+1; ibjet4 < bJetsInv->size(); ibjet4++){
			if(jetsInv->at(ijet1) == bJetsInv->at(ibjet4) || ibjet2 == ibjet4 || ibjet3 == ibjet4 ) continue;		
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cHiggsMass, _minChi2Higgs, _bbMassMin1Higgs, _bbMassMin2Higgs);
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4()
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cZMass, _minChi2Z, _bbMassMin1Z, _bbMassMin2Z);
			// ZH combinatorics 
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH H(jb)Z(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cHiggsMass, cZMass, _minChi2HiggsZ, _bbMassMin1HiggsZ, _bbMassMin2HiggsZ);
			tempminChi2 = _minChi2HiggsZ; tmpMassMin1HiggsZ = _bbMassMin1HiggsZ; tmpMassMin2HiggsZ = _bbMassMin2HiggsZ;
			diMotherReco(*jetsInv->at(ijet1)->getp4(), *bJetsInv->at(ibjet2)->getp4() // ZH Z(jb)H(bb)
				     , *bJetsInv->at(ibjet3)->getp4(), *bJetsInv->at(ibjet4)->getp4()
				     , cZMass, cHiggsMass, _minChi2HiggsZ, _bbMassMin2HiggsZ, _bbMassMin1HiggsZ);
			// pick the lowest minChi2 for the two combinatorics
			if(tempminChi2 < _minChi2HiggsZ){
			    _minChi2HiggsZ = tempminChi2; _bbMassMin1HiggsZ = tmpMassMin1HiggsZ; _bbMassMin2HiggsZ = tmpMassMin2HiggsZ;
			}
		    }
		}
	    }
	}
    } 


    //    std::map<std::string, float> testVars = HypoComb->GetBestPermutation(getLepP4(thisEvent),getJetP4(thisEvent),getJetCSV(thisEvent),*(thisEvent->getMET()->getp4()));
    //    HypoComb.GetBestPermutation(getLepP4(thisEvent),getJetP4(thisEvent),getJetCSV(thisEvent),*(thisEvent->getMET()->getp4()));
    //    std::cout<< "BLR: " << testVars["Evt_blr"] << std::endl;

    thisEvent->eventShapeJet  = new EventShape(vectorsJet);
    thisEvent->eventShapeBjet  = new EventShape(vectorsBjet);
}


void ttHHanalyzer::process(event* thisEvent, sysName sysType, bool up){
// 🔧
    // 1) Golden JSON filter for data
    if (_DataOrMC == "Data") {
        int run  = _ev->run;
        int lumi = _ev->luminosityBlock;
        if (!corrMgr->passGoldenJSON(run, lumi)) {
            if (debugCorrections) {
                std::cout << "[GoldenJSON] Skipping event: run=" 
                          << run << " lumi=" << lumi << std::endl;
            }
            return;
        }
        if (debugCorrections) {
            std::cout << "[GoldenJSON] Accepted event: run=" 
                      << run << " lumi=" << lumi << std::endl;
        }
    }

    _PUWeight = 1.0;
    // 2) Pileup reweighting for MC
    if (_DataOrMC != "Data") {
        float nTrue = _ev->Pileup_nTrueInt;
        double puW  = corrMgr->getPUWeight(nTrue, "nominal");
        _PUWeight   *= puW;
        if (debugCorrections) {
            std::cout << "[PU] nTrue=" << nTrue 
                      << " weight=" << puW 
                      << " -> evtWeight=" << _weight << std::endl;
        }
    }

    // 3) L1 pre-firing weight
    _L1PrefiringWeight = 1.0;
    if (_DataOrMC != "Data") {
	float preFireWgt = _ev->L1PreFiringWeight_Nom;
	_L1PrefiringWeight  *= preFireWgt;
        if (debugCorrections) {
	    std::cout << "[L1PreFiring] weight=" << preFireWgt <<std::endl;
	}
    }

    // 4) Build physics objects in thisEvent from the raw buffer
    createObjects(thisEvent, sysType, up);

    // 5) Final analysis steps: kinematics, histograms, tree
    analyze(thisEvent);
    fillHistos(thisEvent);
    fillTree(thisEvent);

}// 🔧



void ttHHanalyzer::fillHistos(event * thisEvent){

    //    for(int i=0; i < cutflow.size(); i++){
    //	std::cout << "cutflow size: " << cutflow[i] << std::endl; 
	//	hCutFlow->SetBinContent(1, cutflow[0]);
    //    }
    
    int i = 1; 
    for (const auto& x : cutflow){
	hCutFlow->SetBinContent(i, x.second);
	i++;
    }



    thisEvent->getCentrality(thisEvent->getSelJets(), thisEvent->getSelbJets(), jbjetCent);
    thisEvent->getCentrality(thisEvent->getSelJets(), thisEvent->getSelLeptons(), jlepCent);
    thisEvent->getStats(thisEvent->getSelJets(), jetStat);
    thisEvent->getStats(thisEvent->getSelbJets(), bjetStat);
    thisEvent->getStatsComb(thisEvent->getSelJets(), thisEvent->getSelbJets(), bjStat);
    thisEvent->getMaxPTComb(thisEvent->getSelJets(), thisEvent->getSelbJets(), jbbMaxs);
    thisEvent->getMaxPTSame(thisEvent->getSelJets(), jjjMaxs);

    thisEvent->getStatsComb(thisEvent->getSelJets(), thisEvent->getSelLeptons(), ljetStat);
    thisEvent->getStatsComb(thisEvent->getSelbJets(), thisEvent->getSelLeptons(), lbjetStat);
    
    thisEvent->getFoxWolfram(thisEvent->getSelJets(), jetFoxWolfMom);
    thisEvent->getFoxWolfram(thisEvent->getSelbJets(), bjetFoxWolfMom);
 

    //    std::cout << "Number of Hadronic Higgs: " << thisEvent->getnHadronicHiggs() << std::endl;

    hjetNumber->Fill(thisEvent->getnSelJet(),_weight*thisEvent->getbTagSys());
    hHadronicHiggsNumber->Fill(thisEvent->getnHadronicHiggs(),_weight*thisEvent->getbTagSys());
    hBjetNumber->Fill(thisEvent->getnbJet(),_weight*thisEvent->getbTagSys());
    hLightJetNumber->Fill(thisEvent->getnLightJet(),_weight*thisEvent->getbTagSys());

    hjetAverageMass->Fill(thisEvent->getSumSelJetMass()/(float)thisEvent->getnSelJet(),_weight*thisEvent->getbTagSys());
    hHadronicHiggsAverageMass->Fill(thisEvent->getSumSelHadronicHiggsMass()/(float)thisEvent->getnHadronicHiggs(),_weight*thisEvent->getbTagSys());
    hBjetAverageMass->Fill(thisEvent->getSumSelbJetMass()/(float)thisEvent->getnbJet(),_weight*thisEvent->getbTagSys());
    hLightJetAverageMass->Fill(thisEvent->getSumSelLightJetMass()/(float)thisEvent->getnLightJet(),_weight*thisEvent->getbTagSys());
    hBjetAverageMassSqr->Fill((thisEvent->getSumSelbJetMass()*thisEvent->getSumSelbJetMass())/(float)thisEvent->getnbJet(), _weight*thisEvent->getbTagSys());

    if(thisEvent->getnHadronicHiggs() > 0){ 
	hHadronicHiggsSoftDropMass1->Fill(thisEvent->getSelHadronicHiggses()->at(0)->softDropMass,_weight*thisEvent->getbTagSys());
    }

    if(thisEvent->getnHadronicHiggs() > 1){ 
	hHadronicHiggsSoftDropMass2->Fill(thisEvent->getSelHadronicHiggses()->at(1)->softDropMass,_weight*thisEvent->getbTagSys());
    }


    hjetHT->Fill(thisEvent->getSumSelJetScalarpT(),_weight*thisEvent->getbTagSys());
    hBjetHT->Fill(thisEvent->getSumSelbJetScalarpT(),_weight*thisEvent->getbTagSys());
    hHadronicHiggsHT->Fill(thisEvent->getSumSelHadronicHiggsScalarpT(),_weight*thisEvent->getbTagSys());
    hLightJetHT->Fill(thisEvent->getSumSelLightJetScalarpT(),_weight*thisEvent->getbTagSys());
    hAvgDeltaRjj->Fill(jetStat.meandR,_weight*thisEvent->getbTagSys());
    hAvgDeltaEtajj->Fill(jetStat.meandEta,_weight*thisEvent->getbTagSys());
    hminDeltaRjj->Fill(jetStat.mindR,_weight*thisEvent->getbTagSys());
    hminDeltaRMassjj->Fill(jetStat.mindRMass,_weight*thisEvent->getbTagSys());
    hminDeltaRpTjj->Fill(jetStat.mindRpT,_weight*thisEvent->getbTagSys());
    hAvgDeltaRbb->Fill(bjetStat.meandR,_weight*thisEvent->getbTagSys());
    hAvgDeltaRbj->Fill(bjStat.meandR,_weight*thisEvent->getbTagSys());
    hAvgDeltaEtabj->Fill(bjStat.meandEta,_weight*thisEvent->getbTagSys());
    hminDeltaRbj->Fill(bjStat.mindR,_weight*thisEvent->getbTagSys());
    hminDeltaRMassbj->Fill(bjStat.mindRMass,_weight*thisEvent->getbTagSys());
    hminDeltaRpTbj->Fill(bjStat.mindRpT,_weight*thisEvent->getbTagSys());    
    hAvgDeltaEtabb->Fill(bjetStat.meandEta,_weight*thisEvent->getbTagSys());
    hminDeltaRbb->Fill(bjetStat.mindR,_weight*thisEvent->getbTagSys());
    hminDeltaRMassbb->Fill(bjetStat.mindRMass,_weight*thisEvent->getbTagSys());
    hminDeltaRpTbb->Fill(bjetStat.mindRpT,_weight*thisEvent->getbTagSys());
    hmaxDeltaEtabb->Fill(bjetStat.maxdEta,_weight*thisEvent->getbTagSys());
    hmaxDeltaEtajj->Fill(jetStat.maxdEta,_weight*thisEvent->getbTagSys());
    hAvgDeltaEtabj->Fill(bjStat.meandEta,_weight*thisEvent->getbTagSys());
    hmaxDeltaEtabj->Fill(bjStat.maxdEta,_weight*thisEvent->getbTagSys());
    hmaxPTmassjbb->Fill(jbbMaxs.maxPTmass, _weight*thisEvent->getbTagSys());
    hmaxPTmassjjj->Fill(jjjMaxs.maxPTmass, _weight*thisEvent->getbTagSys());

    hInvMassHSingleMatched->Fill(_bbMassMinSHiggsMatched,_weight*thisEvent->getbTagSys());
    hInvMassHSingleNotMatched->Fill(_bbMassMinSHiggsNotMatched,_weight*thisEvent->getbTagSys());
    hChi2HiggsSingleMatched->Fill(_minChi2SHiggsMatched,_weight*thisEvent->getbTagSys());
    hChi2HiggsSingleNotMatched->Fill(_minChi2SHiggsNotMatched,_weight*thisEvent->getbTagSys());

    hInvMassHH1Matched->Fill(_bbMassMinHH1Matched,_weight*thisEvent->getbTagSys());
    hInvMassHH1NotMatched->Fill(_bbMassMinHH1NotMatched,_weight*thisEvent->getbTagSys());
    hInvMassHH2Matched->Fill(_bbMassMinHH2Matched,_weight*thisEvent->getbTagSys());
    hInvMassHH2NotMatched->Fill(_bbMassMinHH2NotMatched,_weight*thisEvent->getbTagSys());
    hChi2HHMatched->Fill(_minChi2HHMatched,_weight*thisEvent->getbTagSys());
    hChi2HHNotMatched->Fill(_minChi2HHNotMatched,_weight*thisEvent->getbTagSys());


    hInvMassH1->Fill(_bbMassMin1Higgs,_weight*thisEvent->getbTagSys());
    hInvMassH2->Fill(_bbMassMin2Higgs,_weight*thisEvent->getbTagSys());
    hPTH1->Fill(_bpTHiggs1,_weight*thisEvent->getbTagSys());
    hPTH2->Fill(_bpTHiggs2,_weight*thisEvent->getbTagSys());
    hInvMassZ1->Fill(_bbMassMin1Z,_weight*thisEvent->getbTagSys());
    hInvMassZ2->Fill(_bbMassMin2Z,_weight*thisEvent->getbTagSys());
    hInvMassHZ1->Fill(_bbMassMin1HiggsZ,_weight*thisEvent->getbTagSys());
    hInvMassHZ2->Fill(_bbMassMin2HiggsZ,_weight*thisEvent->getbTagSys());
    if(fabs(_bbMassMin1Higgs-cHiggsMass) < fabs(_bbMassMin2Higgs-cHiggsMass)){
	hInvMassH1mChi->Fill(_bbMassMin1Higgs,_weight*thisEvent->getbTagSys());
	hInvMassH2mChi->Fill(_bbMassMin2Higgs,_weight*thisEvent->getbTagSys());
    } else  {
	hInvMassH2mChi->Fill(_bbMassMin1Higgs,_weight*thisEvent->getbTagSys());
	hInvMassH1mChi->Fill(_bbMassMin2Higgs,_weight*thisEvent->getbTagSys());
    }
    hChi2Higgs->Fill(_minChi2Higgs,_weight*thisEvent->getbTagSys());
    hChi2Z->Fill(_minChi2Z,_weight*thisEvent->getbTagSys());
    hChi2HiggsZ->Fill(_minChi2HiggsZ,_weight*thisEvent->getbTagSys());

    
    hmet->Fill(thisEvent->getMET()->getp4()->Pt(),_weight*thisEvent->getbTagSys());
    //    hmetPhi->Fill(thisEvent->getMET()->getp4()->Phi(),_weight*thisEvent->getbTagSys());
    // hmetEta->Fill(thisEvent->getMET()->getp4()->Eta(),_weight*thisEvent->getbTagSys());
    
    hAplanarity->Fill(thisEvent->eventShapeJet->getAplanarity(), _weight*thisEvent->getbTagSys());
    hSphericity->Fill(thisEvent->eventShapeJet->getSphericity(), _weight*thisEvent->getbTagSys());
    hTransSphericity->Fill(thisEvent->eventShapeJet->getTransSphericity(), _weight*thisEvent->getbTagSys());
    hCvalue->Fill(thisEvent->eventShapeJet->getC(), _weight*thisEvent->getbTagSys());
    hDvalue->Fill(thisEvent->eventShapeJet->getD(), _weight*thisEvent->getbTagSys());
    hCentralityjb->Fill(jbjetCent.centrality, _weight*thisEvent->getbTagSys());    
    hCentralityjl->Fill(jlepCent.centrality, _weight*thisEvent->getbTagSys());    

    hH0->Fill(jetFoxWolfMom.h0, _weight*thisEvent->getbTagSys());
    hH1->Fill(jetFoxWolfMom.h1, _weight*thisEvent->getbTagSys());
    hH2->Fill(jetFoxWolfMom.h2, _weight*thisEvent->getbTagSys());
    hH3->Fill(jetFoxWolfMom.h3, _weight*thisEvent->getbTagSys());
    hH4->Fill(jetFoxWolfMom.h4, _weight*thisEvent->getbTagSys());
    hR1->Fill(jetFoxWolfMom.r1, _weight*thisEvent->getbTagSys());
    hR2->Fill(jetFoxWolfMom.r2, _weight*thisEvent->getbTagSys());
    hR3->Fill(jetFoxWolfMom.r3, _weight*thisEvent->getbTagSys());
    hR4->Fill(jetFoxWolfMom.r4, _weight*thisEvent->getbTagSys()); 


    hBjetH0->Fill(bjetFoxWolfMom.h0, _weight*thisEvent->getbTagSys());
    hBjetH1->Fill(bjetFoxWolfMom.h1, _weight*thisEvent->getbTagSys());
    hBjetH2->Fill(bjetFoxWolfMom.h2, _weight*thisEvent->getbTagSys());
    hBjetH3->Fill(bjetFoxWolfMom.h3, _weight*thisEvent->getbTagSys());
    hBjetH4->Fill(bjetFoxWolfMom.h4, _weight*thisEvent->getbTagSys());
    hBjetR1->Fill(bjetFoxWolfMom.r1, _weight*thisEvent->getbTagSys());
    hBjetR2->Fill(bjetFoxWolfMom.r2, _weight*thisEvent->getbTagSys());
    hBjetR3->Fill(bjetFoxWolfMom.r3, _weight*thisEvent->getbTagSys());
    hBjetR4->Fill(bjetFoxWolfMom.r4, _weight*thisEvent->getbTagSys()); 


    hBjetAplanarity->Fill(thisEvent->eventShapeBjet->getAplanarity(), _weight*thisEvent->getbTagSys());
    hBjetSphericity->Fill(thisEvent->eventShapeBjet->getSphericity(), _weight*thisEvent->getbTagSys());
    hBjetTransSphericity->Fill(thisEvent->eventShapeBjet->getTransSphericity(), _weight*thisEvent->getbTagSys());
    hBjetCvalue->Fill(thisEvent->eventShapeBjet->getC(), _weight*thisEvent->getbTagSys());
    hBjetDvalue->Fill(thisEvent->eventShapeBjet->getD(), _weight*thisEvent->getbTagSys());


    for(int ih=0; ih < thisEvent->getnSelJet() && ih < nHistsJets; ih++){
	hjetsPTs.at(ih)->Fill(thisEvent->getSelJets()->at(ih)->getp4()->Pt(),_weight*thisEvent->getbTagSys());
	hjetsEtas.at(ih)->Fill(thisEvent->getSelJets()->at(ih)->getp4()->Eta(),_weight*thisEvent->getbTagSys());
	hjetsBTagDisc.at(ih)->Fill(getJetCSV(thisEvent).at(ih),_weight*thisEvent->getbTagSys());
    }

    for(int ih=0; ih < thisEvent->getnbJet() && ih < nHistsbJets; ih++){
	hbjetsPTs.at(ih)->Fill(thisEvent->getSelbJets()->at(ih)->getp4()->Pt(),_weight*thisEvent->getbTagSys());
	hbjetsEtas.at(ih)->Fill(thisEvent->getSelbJets()->at(ih)->getp4()->Eta(),_weight*thisEvent->getbTagSys());
	hbjetsBTagDisc.at(ih)->Fill(getbJetCSV(thisEvent).at(ih),_weight*thisEvent->getbTagSys());
    }


    for(int ih=0; ih < thisEvent->getnLightJet() && ih < nHistsLightJets; ih++){
	hLightJetsPTs.at(ih)->Fill(thisEvent->getSelLightJets()->at(ih)->getp4()->Pt(),_weight*thisEvent->getbTagSys());
	hLightJetsEtas.at(ih)->Fill(thisEvent->getSelLightJets()->at(ih)->getp4()->Eta(),_weight*thisEvent->getbTagSys());
	hLightJetsBTagDisc.at(ih)->Fill(getlightJetCSV(thisEvent).at(ih),_weight*thisEvent->getbTagSys());
    }

    hleptonNumber->Fill(thisEvent->getnSelLepton(),_weight*thisEvent->getbTagSys());
    hElecNumber->Fill(thisEvent->getnSelElectron(),_weight*thisEvent->getbTagSys());
    hMuonNumber->Fill(thisEvent->getnSelMuon(),_weight*thisEvent->getbTagSys());

//    if(thisEvent->getnSelMuon() == 2){
//	hDiMuonMass->Fill(thisEvent->getSelMuonsMass(),_weight*thisEvent->getbTagSys());
//	hDiMuonPT->Fill(thisEvent->getSelMuonsPT(),_weight*thisEvent->getbTagSys());
//	hDiMuonEta->Fill(thisEvent->getSelMuonsEta(),_weight*thisEvent->getbTagSys());
//    }

//    if(thisEvent->getnSelElectron() == 2){
//	hDiElectronMass->Fill(thisEvent->getSelElectronsMass(),_weight*thisEvent->getbTagSys());
//	hDiElectronPT->Fill(thisEvent->getSelElectronsPT(),_weight*thisEvent->getbTagSys());
//	hDiElectronEta->Fill(thisEvent->getSelElectronsEta(),_weight*thisEvent->getbTagSys());
//    }

    hleptonHT->Fill(thisEvent->getSelLeptonHT(),_weight*thisEvent->getbTagSys());
    hST->Fill(thisEvent->getSelLeptonST(),_weight*thisEvent->getbTagSys());
    hLeptonPT1->Fill(thisEvent->getSelLeptons()->at(0)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
    hLeptonEta1->Fill(thisEvent->getSelLeptons()->at(0)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
   // hLeptonPT2->Fill(thisEvent->getSelLeptons()->at(1)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
  //  hLeptonEta2->Fill(thisEvent->getSelLeptons()->at(1)->getp4()->Eta(), _weight*thisEvent->getbTagSys());


    if(thisEvent->getnSelMuon() > 0){
	hMuonPT1->Fill(thisEvent->getSelMuons()->at(0)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hMuonEta1->Fill(thisEvent->getSelMuons()->at(0)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    }
    
    if(thisEvent->getnSelElectron() > 0){
	hElePT1->Fill(thisEvent->getSelElectrons()->at(0)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hEleEta1->Fill(thisEvent->getSelElectrons()->at(0)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    }
  /*  
    if(thisEvent->getnSelMuon() > 1){
	hMuonPT2->Fill(thisEvent->getSelMuons()->at(1)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hMuonEta2->Fill(thisEvent->getSelMuons()->at(1)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    }
    
    if(thisEvent->getnSelElectron() > 1){
	hElePT2->Fill(thisEvent->getSelElectrons()->at(1)->getp4()->Pt(), _weight*thisEvent->getbTagSys());
	hEleEta2->Fill(thisEvent->getSelElectrons()->at(1)->getp4()->Eta(), _weight*thisEvent->getbTagSys());
    } 
*/
    hLepCharge1->Fill(thisEvent->getSelLeptons()->at(0)->charge, _weight*thisEvent->getbTagSys());
//    hLepCharge2->Fill(thisEvent->getSelLeptons()->at(1)->charge, _weight*thisEvent->getbTagSys());

   int nSel = thisEvent->getnSelJet(); // 🔧
   for(int ih = 0; ih < nSel && ih < nHistsJets; ++ih){
     float origPt    = thisEvent->getSelJets()->at(ih)->origPt;   // nanoAOD JEC 적용 후 pT
     float smearedPt = thisEvent->getSelJets()->at(ih)->getp4()->Pt();           // JEC+JER 적용된 최종 pT
     float ptJEC     = thisEvent->getSelJets()->at(ih)->ptJEC;

     // A) 원본 pT
     //h_origPt.at(ih)->Fill(origPt, _weight);
     h_origPt.at(ih)->Fill(origPt);

     h_JECPt.at(ih)->Fill(ptJEC);

     // B) 스메어된 pT
     //h_smearedPt.at(ih)->Fill(smearedPt, _weight);
     h_smearedPt.at(ih)->Fill(smearedPt);
   }// 🔧

	
}



void ttHHanalyzer::writeHistos(){
    _of->file->cd();
    _histoDirs.at(0)->cd();
    for(int ih=0; ih<nHistsJets; ih++){
	hjetsPTs.at(ih)->Write();
	hjetsEtas.at(ih)->Write();
	hjetsBTagDisc.at(ih)->Write();
    }
    for(int ih=0; ih<nHistsbJets; ih++){
	hbjetsPTs.at(ih)->Write();
	hbjetsEtas.at(ih)->Write();
	hbjetsBTagDisc.at(ih)->Write();
    }

    for(int ih=0; ih<nHistsLightJets; ih++){
	hLightJetsPTs.at(ih)->Write();
	hLightJetsEtas.at(ih)->Write();
	hLightJetsBTagDisc.at(ih)->Write();
    }

	   for(int ih=0; ih<nHistsJets; ih++){ // 🔧
//        hjetsPTs.at(ih)->Write();
//        hjetsEtas.at(ih)->Write();
//        hjetsBTagDisc.at(ih)->Write();
  
        // JEC-only 원본 pT
        h_origPt.at(ih)->Write();
	h_JECPt.at(ih)->Write();

        // JEC+JER 스메어된 pT
        h_smearedPt.at(ih)->Write();

    } // 🔧
    hInvMassHSingleMatched->Write();
    hInvMassHSingleNotMatched->Write();
    hChi2HiggsSingleNotMatched->Write();
    hChi2HiggsSingleMatched->Write();
    hInvMassHH1Matched->Write();
    hInvMassHH1NotMatched->Write();
    hInvMassHH2Matched->Write();
    hInvMassHH2NotMatched->Write();
    hChi2HHNotMatched->Write();
    hChi2HHMatched->Write();

    hjetNumber->Write();
    hBjetNumber->Write();
    hHadronicHiggsNumber->Write();
    hLightJetNumber->Write();
    hjetAverageMass->Write();
    hBjetAverageMass->Write();
    hHadronicHiggsAverageMass->Write();
    hLightJetAverageMass->Write();
    hBjetAverageMassSqr->Write();
    hHadronicHiggsSoftDropMass1->Write();
    hHadronicHiggsSoftDropMass2->Write();
    hjetHT->Write();
    hBjetHT->Write();
    hHadronicHiggsHT->Write();
    hLightJetHT->Write();
    hAvgDeltaRjj->Write();
    hminDeltaRjj->Write();
    hminDeltaRMassjj->Write();
    hminDeltaRpTjj->Write();
    hAvgDeltaRbb->Write();
    hAvgDeltaEtajj->Write();
    hAvgDeltaEtabb->Write();
    hminDeltaRbb->Write();
    hminDeltaRMassbb->Write();
    hminDeltaRpTbb->Write();
    hmaxDeltaEtabb->Write();
    hmaxDeltaEtajj->Write();
    hAvgDeltaRbj->Write();
    hAvgDeltaEtabj->Write();
    hminDeltaRbj->Write();
    hminDeltaRMassbj->Write();
    hminDeltaRpTbj->Write();
    hmaxDeltaEtabj->Write();
    hmaxPTmassjbb->Write();
    hmaxPTmassjjj->Write();

    hPTH1->Write();
    hPTH2->Write();
    hInvMassH1->Write();
    hInvMassH2->Write();
    hInvMassH1mChi->Write();
    hInvMassH2mChi->Write();
    hInvMassHZ1->Write();
    hInvMassHZ2->Write();
    hInvMassZ1->Write();
    hInvMassZ2->Write();
    hChi2Higgs->Write();
    hChi2HiggsZ->Write();
    hChi2Z->Write();


    hmet->Write();
    // hmetPhi->Write();
    // hmetEta->Write();

    hAplanarity->Write();
    hSphericity->Write();
    hTransSphericity->Write();
    hCvalue->Write();
    hDvalue->Write();
    hCentralityjb->Write();
    hCentralityjl->Write();

    hH0->Write();
    hH1->Write();
    hH2->Write();
    hH3->Write();
    hH4->Write();
    hR1->Write();
    hR2->Write();
    hR3->Write();
    hR4->Write(); 

    hBjetH0->Write();
    hBjetH1->Write();
    hBjetH2->Write();
    hBjetH3->Write();
    hBjetH4->Write();
    hBjetR1->Write();
    hBjetR2->Write();
    hBjetR3->Write();
    hBjetR4->Write(); 

    hBjetAplanarity->Write();
    hBjetSphericity->Write();
    hBjetTransSphericity->Write();
    hBjetCvalue->Write();
    hBjetDvalue->Write();

    _histoDirs.at(1)->cd();
    
    hLepCharge1->Write();
  //  hLepCharge2->Write();

    hleptonNumber->Write();
    hElecNumber->Write();
    hMuonNumber->Write();
   // hDiMuonMass->Write();
   // hDiElectronMass->Write();
    hDiMuonPT->Write();
    hDiElectronPT->Write();
    hDiMuonEta->Write();
    hDiElectronEta->Write();
    hleptonHT->Write();
    hST->Write();
    hLeptonEta1->Write();
    hLeptonPT1->Write();
 //   hLeptonEta2->Write();
 //   hLeptonPT2->Write();
    
    hMuonEta1->Write();
    hMuonPT1->Write();
    hEleEta1->Write();
    hElePT1->Write();
//    hMuonEta2->Write();
//    hMuonPT2->Write();
//    hEleEta2->Write();
//    hElePT2->Write(); 
}
void ttHHanalyzer::fillTree(event * thisEvent){
// 🔧
////////////////////////////////////////////////////////////////////////////////////////
    // For Trigger Path
    passTrigger_HLT_IsoMu24 = _ev->HLT_IsoMu24; // Reference Muon Trigger
    passTrigger_HLT_Ele30_WPTight_Gsf = _ev->HLT_Ele30_WPTight_Gsf;
    //passTrigger_HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59 = _ev->HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59;
    //passTrigger_HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94 = _ev->HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94;
    //passTrigger_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 = _ev->HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5;
    passTrigger_6J1T_B    = _ev->HLT_PFHT430_SixJet40_BTagCSV_p080;
    passTrigger_6J1T_CDEF = _ev->HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5;
    passTrigger_6J2T_B    = _ev->HLT_PFHT380_SixJet32_DoubleBTagCSV_p075;
    passTrigger_6J2T_CDEF = _ev->HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2;
    passTrigger_4J3T_B    = _ev->HLT_HT300PT30_QuadJet_75_60_45_40_TripeCSV_p07;
    passTrigger_4J3T_CDEF = _ev->HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0;

    nMuons = thisEvent->getnSelMuon();
    nElecs = thisEvent->getnSelElectron();
    nJets = thisEvent->getnSelJet();
    nbJets = thisEvent->getnbJet();
    HT = thisEvent->getSumSelJetScalarpT();

    // Fill the jet information [ It will fill the nJets && maximum 30th jets ]
    for (int i = 0; i < nJets && i < 30; ++i) {
        jetPt[i] = thisEvent->getSelJets()->at(i)->getp4()->Pt();
        jetEta[i] = thisEvent->getSelJets()->at(i)->getp4()->Eta();
        bTagScore[i] = thisEvent->getSelJets()->at(i)->bTagCSV;
        hadFlavs[i] = thisEvent->getSelJets()->at(i)->hadFlav;
        partonFlavs[i] = thisEvent->getSelJets()->at(i)->partonFlav;
        // In the Event Buffer, you can find the struct jet_s which accesses to defined branch variables in NanoAODv9
        // And then in the .hh files, you also can check variables like  bTagCSV or hadFlav as objectJet class's  member variable
        // Then you can get the info at NanoAODv9 from jet_s structure and put them into the .hh's objectJet class member variable in the iteration of jets
        // Finally here, you can store variables in the objectJet class's member variables into tree

    } 
    for (int i = nJets; i < 30; ++i) {
        jetPt[i] = -999;
        jetEta[i] = -999;
        bTagScore[i] = -999;
        hadFlavs[i] = -999;
        partonFlavs[i] = -999;
    }

    eventNumber = _ev->event;
    runNumber = _ev->run;
 
    SampleWeight = _SampleWeight;
    PUWeight = _PUWeight;
    L1PrefiringWeight = _L1PrefiringWeight;

////////////////////////////////////////////////////////////////////////////////////////
   // 🔧
	
   
    jetPT1 = thisEvent->getSelJets()->at(0)->getp4()->Pt();
    jetPT2 = thisEvent->getSelJets()->at(1)->getp4()->Pt();
    jetPT3 = thisEvent->getSelJets()->at(2)->getp4()->Pt();
    jetPT4 = thisEvent->getSelJets()->at(3)->getp4()->Pt();
    jetEta1 = thisEvent->getSelJets()->at(0)->getp4()->Eta();
    jetEta2 = thisEvent->getSelJets()->at(1)->getp4()->Eta();
    jetEta3 = thisEvent->getSelJets()->at(2)->getp4()->Eta();
    jetEta4 = thisEvent->getSelJets()->at(3)->getp4()->Eta();
    jetBTagDisc1 = getJetCSV(thisEvent).at(0);
    jetBTagDisc2 = getJetCSV(thisEvent).at(1);
    jetBTagDisc3 = getJetCSV(thisEvent).at(2);
    jetBTagDisc4 = getJetCSV(thisEvent).at(3);


    bjetPT1 = thisEvent->getSelbJets()->at(0)->getp4()->Pt();
    bjetPT2 = thisEvent->getSelbJets()->at(1)->getp4()->Pt();
    bjetPT3 = thisEvent->getSelbJets()->at(2)->getp4()->Pt();
    bjetEta1 = thisEvent->getSelbJets()->at(0)->getp4()->Eta();
    bjetEta2 = thisEvent->getSelbJets()->at(1)->getp4()->Eta();
    bjetEta3 = thisEvent->getSelbJets()->at(2)->getp4()->Eta();
    bjetBTagDisc1 = getbJetCSV(thisEvent).at(0);
    bjetBTagDisc2 = getbJetCSV(thisEvent).at(1);
    bjetBTagDisc3 = getbJetCSV(thisEvent).at(2);
    bbjetHiggsMatched1 = thisEvent->getSelbJets()->at(0)->matchedtoHiggs;
    bbjetHiggsMatched2 = thisEvent->getSelbJets()->at(1)->matchedtoHiggs;
    bbjetHiggsMatched3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggs;
    bbjetHiggsMatcheddR1 = thisEvent->getSelbJets()->at(0)->matchedtoHiggsdR;
    bbjetHiggsMatcheddR2 = thisEvent->getSelbJets()->at(1)->matchedtoHiggsdR;
    bbjetHiggsMatcheddR3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggsdR;
    bbjetMinChiHiggsIndex1 = thisEvent->getSelbJets()->at(0)->minChiHiggsIndex;
    bbjetMinChiHiggsIndex2 = thisEvent->getSelbJets()->at(1)->minChiHiggsIndex;
    bbjetMinChiHiggsIndex3 = thisEvent->getSelbJets()->at(2)->minChiHiggsIndex;


    if(thisEvent->getnSelJet() > 4){
    	jetPT5 = thisEvent->getSelJets()->at(4)->getp4()->Pt();
	jetEta5 = thisEvent->getSelJets()->at(4)->getp4()->Eta();
	jetBTagDisc5 = getJetCSV(thisEvent).at(4);
    } else{
	jetPT5 = -6;
	jetEta5 = -6;
	jetBTagDisc5 = -6;
    }
    
    if(thisEvent->getnSelJet() > 5){
	jetPT6 = thisEvent->getSelJets()->at(5)->getp4()->Pt();
	jetEta6 = thisEvent->getSelJets()->at(5)->getp4()->Eta();
	jetBTagDisc6 = getJetCSV(thisEvent).at(5);
    } else{
	jetPT6 = -6;
	jetEta6 = -6;
	jetBTagDisc6 = -6;
    }


    if(thisEvent->getnSelJet() > 6){
	jetPT7 = thisEvent->getSelJets()->at(6)->getp4()->Pt();
	jetEta7 = thisEvent->getSelJets()->at(6)->getp4()->Eta();
	jetBTagDisc7 = getJetCSV(thisEvent).at(6);
    } else{
	jetPT7 = -6;
	jetEta7 = -6;
	jetBTagDisc7 = -6;
    }

    if(thisEvent->getnSelJet() > 7){
	jetPT8 = thisEvent->getSelJets()->at(7)->getp4()->Pt();
	jetEta8 = thisEvent->getSelJets()->at(7)->getp4()->Eta();
	jetBTagDisc8 = getJetCSV(thisEvent).at(7);
    } else{
	jetPT8 = -6;
	jetEta8 = -6;
	jetBTagDisc8 = -6;
    }

    if(thisEvent->getnbJet() > 2){
	bjetPT3 = thisEvent->getSelbJets()->at(2)->getp4()->Pt();
	bjetEta3 = thisEvent->getSelbJets()->at(2)->getp4()->Eta();
	bjetBTagDisc3 = getbJetCSV(thisEvent).at(2);
	bbjetHiggsMatched3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggs;
	bbjetHiggsMatcheddR3 = thisEvent->getSelbJets()->at(2)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex3 = thisEvent->getSelbJets()->at(2)->minChiHiggsIndex;
    } else {
	bjetPT3 = -6;
	bjetEta3 = -6;
	bjetBTagDisc3 = -6;
	bbjetHiggsMatched3 = 0;
	bbjetHiggsMatcheddR3 = -6;
	bbjetMinChiHiggsIndex3 = -6;
	} 


    if(thisEvent->getnbJet() > 3){
	bjetPT4 = thisEvent->getSelbJets()->at(3)->getp4()->Pt();
	bjetEta4 = thisEvent->getSelbJets()->at(3)->getp4()->Eta();
	bjetBTagDisc4 = getbJetCSV(thisEvent).at(3);
	bbjetHiggsMatched4 = thisEvent->getSelbJets()->at(3)->matchedtoHiggs;
	bbjetHiggsMatcheddR4 = thisEvent->getSelbJets()->at(3)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex4 = thisEvent->getSelbJets()->at(3)->minChiHiggsIndex;
    } else {
	bjetPT4 = -6;
	bjetEta4 = -6;
	bjetBTagDisc4 = -6;
	bbjetHiggsMatched4 = 0;
	bbjetHiggsMatcheddR4 = -6;
	bbjetMinChiHiggsIndex4 = -6;
    }
    
    if(thisEvent->getnbJet() > 4){
	bjetPT5 = thisEvent->getSelbJets()->at(4)->getp4()->Pt();
	bjetEta5 = thisEvent->getSelbJets()->at(4)->getp4()->Eta();
	bjetBTagDisc5 = getbJetCSV(thisEvent).at(4);
	bbjetHiggsMatched5 = thisEvent->getSelbJets()->at(4)->matchedtoHiggs;
	bbjetHiggsMatcheddR5 = thisEvent->getSelbJets()->at(4)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex5 = thisEvent->getSelbJets()->at(4)->minChiHiggsIndex;
    } else {
	bjetPT5 = -6;
	bjetEta5 = -6;
	bjetBTagDisc5 = -6;
	bbjetHiggsMatched5 = 0;
	bbjetHiggsMatcheddR5 = -6;
	bbjetMinChiHiggsIndex5 = -6;
    }

    if(thisEvent->getnbJet() > 5){
	bjetPT6 = thisEvent->getSelbJets()->at(5)->getp4()->Pt();
	bjetEta6 = thisEvent->getSelbJets()->at(5)->getp4()->Eta();
	bjetBTagDisc6 = getbJetCSV(thisEvent).at(5);
	bbjetHiggsMatched6 = thisEvent->getSelbJets()->at(5)->matchedtoHiggs;
	bbjetHiggsMatcheddR6 = thisEvent->getSelbJets()->at(5)->matchedtoHiggsdR;
	bbjetMinChiHiggsIndex6 = thisEvent->getSelbJets()->at(5)->minChiHiggsIndex;
    } else {
	bjetPT6 = -6;
	bjetEta6 = -6;
	bjetBTagDisc6 = -6;
	bbjetHiggsMatched6 = 0;
	bbjetHiggsMatcheddR6 = -6;
	bbjetMinChiHiggsIndex6 = -6;
    }

    if(thisEvent->getnbJet() > 6){
	bjetPT7 = thisEvent->getSelbJets()->at(6)->getp4()->Pt();
	bjetEta7 = thisEvent->getSelbJets()->at(6)->getp4()->Eta();
	bjetBTagDisc7 = getbJetCSV(thisEvent).at(6);
	//bbjetHiggsMatched7 = thisEvent->getSelbJets()->at(6)->matchedtoHiggs;
	//bbjetHiggsMatcheddR7 = thisEvent->getSelbJets()->at(6)->matchedtoHiggsdR;
    } else {
	bjetPT7 = -6;
	bjetEta7 = -6;
	bjetBTagDisc7 = -6;
	//bbjetHiggsMatched7 = 0;
	//bbjetHiggsMatcheddR7 = -6;
    }


    if(thisEvent->getnbJet() > 7){
	bjetPT8 = thisEvent->getSelbJets()->at(7)->getp4()->Pt();
	bjetEta8 = thisEvent->getSelbJets()->at(7)->getp4()->Eta();
	bjetBTagDisc8 = getbJetCSV(thisEvent).at(7);
	//bbjetHiggsMatched8 = thisEvent->getSelbJets()->at(7)->matchedtoHiggs;
	//bbjetHiggsMatcheddR8 = thisEvent->getSelbJets()->at(7)->matchedtoHiggsdR;
    } else {
	bjetPT8 = -6;
	bjetEta8 = -6;
	bjetBTagDisc8 = -6;
	//bbjetHiggsMatched8 = 0;
	//bbjetHiggsMatcheddR8 = 0;
	} 


    if(thisEvent->getnLightJet() > 0){
    	lightjetPT1 = thisEvent->getSelLightJets()->at(0)->getp4()->Pt();
	lightjetEta1 = thisEvent->getSelLightJets()->at(0)->getp4()->Eta();
	lightjetBTagDisc1 = getlightJetCSV(thisEvent).at(0);
    } else{
	lightjetPT1 = -6;
	lightjetEta1 = -6;
	lightjetBTagDisc1 = -6;
    }

    if(thisEvent->getnLightJet() > 1){
    	lightjetPT2 = thisEvent->getSelLightJets()->at(1)->getp4()->Pt();
	lightjetEta2 = thisEvent->getSelLightJets()->at(1)->getp4()->Eta();
	lightjetBTagDisc2 = getlightJetCSV(thisEvent).at(1);
    } else{
	lightjetPT2 = -6;
	lightjetEta2 = -6;
	lightjetBTagDisc2 = -6;
    }

    if(thisEvent->getnLightJet() > 2){
    	lightjetPT3 = thisEvent->getSelLightJets()->at(2)->getp4()->Pt();
	lightjetEta3 = thisEvent->getSelLightJets()->at(2)->getp4()->Eta();
	lightjetBTagDisc3 = getlightJetCSV(thisEvent).at(2);
    } else{
	lightjetPT3 = -6;
	lightjetEta3 = -6;
	lightjetBTagDisc3 = -6;
    }


    weight= _weight;
    jetAverageMass = thisEvent->getSumSelJetMass()/thisEvent->getnSelJet();
    bjetAverageMass = thisEvent->getSumSelbJetMass()/thisEvent->getnbJet();
    lightJetAverageMass = thisEvent->getSumSelLightJetMass()/thisEvent->getnLightJet();
    bjetAverageMassSqr = (thisEvent->getSumSelbJetMass()*thisEvent->getSumSelbJetMass())/thisEvent->getnbJet();
    met = thisEvent->getMET()->getp4()->Pt();
    averageDeltaRjj = jetStat.meandR;
    averageDeltaRbb = bjetStat.meandR;
    averageDeltaRbj = bjStat.meandR;
    averageDeltaEtajj = jetStat.meandEta;
    averageDeltaEtabb = bjetStat.meandEta;
    averageDeltaEtabj = bjStat.meandEta;
    minDeltaRjj = jetStat.mindR;
    minDeltaRbb = bjetStat.mindR;
    minDeltaRbj = bjStat.mindR;
    maxDeltaEtabb = bjetStat.maxdEta;
    maxDeltaEtajj = jetStat.maxdEta;
    maxDeltaEtabj = bjStat.maxdEta;
    minDeltaRMassjj = jetStat.mindRMass;
    minDeltaRMassbb = bjetStat.mindRMass;
    minDeltaRMassbj = bjStat.mindRMass;
    minDeltaRpTjj = jetStat.mindRpT;
    minDeltaRpTbb = bjetStat.mindRpT;
    minDeltaRpTbj = bjStat.mindRpT;
    maxPTmassjjj = jjjMaxs.maxPTmass;
    maxPTmassjbb = jbbMaxs.maxPTmass;
    H0 = jetFoxWolfMom.h0;
    H1 = jetFoxWolfMom.h1;
    H2 = jetFoxWolfMom.h2;
    H3 = jetFoxWolfMom.h3;
    H4 = jetFoxWolfMom.h4;
    bH0 = bjetFoxWolfMom.h0;
    bH1 = bjetFoxWolfMom.h1;
    bH2 = bjetFoxWolfMom.h2;
    bH3 = bjetFoxWolfMom.h3;
    bH4 = bjetFoxWolfMom.h4;
    R1 = jetFoxWolfMom.r1;
    R2 = jetFoxWolfMom.r2;
    R3 = jetFoxWolfMom.r3;
    R4 = jetFoxWolfMom.r4;
    bR1 = bjetFoxWolfMom.r1;
    bR2 = bjetFoxWolfMom.r2;
    bR3 = bjetFoxWolfMom.r3;
    bR4 = bjetFoxWolfMom.r4;

    jetHT = thisEvent->getSumSelJetScalarpT();  
    bjetHT = thisEvent->getSumSelbJetScalarpT();
    lightjetHT = thisEvent->getSumSelLightJetScalarpT();
    jetNumber = thisEvent->getnSelJet();
    bjetNumber = thisEvent->getnbJet();
    lightjetNumber = thisEvent->getnLightJet();
    invMassZ1 = _bbMassMin1Z; //ZZ
    invMassZ2 = _bbMassMin2Z;
    chi2Z = _minChi2Z;
    invMassH1 = _bbMassMin1Higgs; //HH
    invMassH2 = _bbMassMin2Higgs;
    chi2Higgs = _minChi2Higgs;
    chi2HiggsZ = _minChi2HiggsZ; //ZH
    invMassHiggsZ1 = _bbMassMin1HiggsZ;
    invMassHiggsZ2 = _bbMassMin2HiggsZ;
    PTH1 = _bpTHiggs1;
    PTH2 = _bpTHiggs2;


    centralityjb = jbjetCent.centrality; 
    centralityjl = jlepCent.centrality; 
    aplanarity = thisEvent->eventShapeJet->getAplanarity();
    sphericity = thisEvent->eventShapeJet->getSphericity();
    transSphericity = thisEvent->eventShapeJet->getTransSphericity();
    cValue = thisEvent->eventShapeJet->getC();
    dValue = thisEvent->eventShapeJet->getD();
    baplanarity = thisEvent->eventShapeBjet->getAplanarity();
    bsphericity = thisEvent->eventShapeBjet->getSphericity();
    btransSphericity = thisEvent->eventShapeBjet->getTransSphericity();
    bcValue = thisEvent->eventShapeBjet->getC();
    bdValue = thisEvent->eventShapeBjet->getD();

    leptonPT1 = thisEvent->getSelLeptons()->at(0)->getp4()->Pt();
  //  leptonPT2 = thisEvent->getSelLeptons()->at(1)->getp4()->Pt();
    leptonEta1 = thisEvent->getSelLeptons()->at(0)->getp4()->Eta();
    //leptonEta2 = thisEvent->getSelLeptons()->at(1)->getp4()->Eta();
    leptonCharge1 = thisEvent->getSelLeptons()->at(0)->charge;
  //  leptonCharge2 = thisEvent->getSelLeptons()->at(1)->charge;
    leptonHT = thisEvent->getSelLeptonHT();
    ST = thisEvent->getSelLeptonST();


    if(thisEvent->getnSelMuon() > 0){
	muonPT1 = thisEvent->getSelMuons()->at(0)->getp4()->Pt();
	muonEta1 = thisEvent->getSelMuons()->at(0)->getp4()->Eta();
    } else {
	muonPT1 = -6;
	muonEta1 = -6;
    }
/*
    if(thisEvent->getnSelMuon() > 1){
	muonPT2 = thisEvent->getSelMuons()->at(1)->getp4()->Pt();
	muonEta2 = thisEvent->getSelMuons()->at(1)->getp4()->Eta();
	diMuonMass = thisEvent->getSelMuonsMass();
    } else {
	muonPT2 = -6;
	muonEta2 = -6;
	diMuonMass = -6;
    } 
*/
    if(thisEvent->getnSelElectron() > 0){
	elePT1 = thisEvent->getSelElectrons()->at(0)->getp4()->Pt();
	eleEta1 = thisEvent->getSelElectrons()->at(0)->getp4()->Eta();
    } else {
	elePT1 = -6;
	eleEta1 = -6;
    }
/*
    if(thisEvent->getnSelElectron() > 1){
	elePT2 = thisEvent->getSelElectrons()->at(1)->getp4()->Pt();
	eleEta2 = thisEvent->getSelElectrons()->at(1)->getp4()->Eta();
	diElectronMass = thisEvent->getSelElectronsMass();
    } else {
	elePT2 = -6;
	eleEta2 = -6;
	diElectronMass = -6;
    } 

  */  
    _inputTree->Fill();
}

void ttHHanalyzer::writeTree(){
    _of->file->cd();
    _treeDirs.at(0)->cd();
    _inputTree->Write();
    //    _inputTree->Delete();
    
}
//changed from here
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



int main(int argc, char** argv){
    commandLine cl(argc, argv);
    vector<string> filenames = fileNames(cl.filelist);
    double weight = cl.externalweight;   // Get global weight 
 
    // Create tree reader
    itreestream stream(filenames, "Events");
    if ( !stream.good() ) error("can't read root input files");

    eventBuffer ev(stream);
    std::cout << " Output filename: " << cl.outputfilename << std::endl;
    ////ttHHanalyzer analysis(cl.outputfilename, &ev, weight, true)
  
    // If you want to check or modify arguments,
    // Please check the [ src/tnm.cc ]
    // Arguments structure --> filelist, outputDirName, weight, Year, Data or MC, sampleName
    ttHHanalyzer analysis(cl.outputfilename, &ev, weight, true, cl.year, cl.isData, cl.sampleName);
    analysis.performAnalysis();

    ev.close();
    //    of.close();
    return 0;
}

#ifndef _CommonBDTvars_h
#define _CommonBDTvars_h

#include <vector>
#include <map>
#include "TVector.h"
#include "TLorentzVector.h"

typedef std::map<std::string, std::string> mparams;
typedef std::vector< TLorentzVector > vecTLorentzVector;
typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<std::vector<int> > vvint;
typedef std::vector<std::string> vstring;
typedef std::vector<double> vdouble;
typedef std::vector<int> vint;


using namespace std;

class CommonBDTvars{

	// === Functions === //
	public: 
		// Constructor(s) and destructor
		CommonBDTvars();
		virtual ~CommonBDTvars();
		


	//Algorithms 
   
		void getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, double &aplanarity, double &sphericity, double &tsphericity);
		void getFox(vecTLorentzVector jets, double &h0, double &h1, double &h2, double &h3, double &h4);
		double getBestHiggsMass(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, vdouble btag, double &minChi, double &dRbb, TLorentzVector &bjet1, TLorentzVector &bjet2, vecTLorentzVector loose_jets, vdouble loose_btag, double CSVMwp);
   
    
	
	
		void convert_jets_to_TLVs(vvdouble jets, vecTLorentzVector &vect_of_jet_TLVs);
		void vect_of_tagged_TLVs(vvdouble jets, vdouble jetCSV, vecTLorentzVector &vect_of_btag_TLVs, double CSVMwp);
		double get_jet_jet_etamax (vvdouble jets);
		double get_jet_tag_etamax (vvdouble jets, vdouble jetCSV, double CSVMwp);
		double get_tag_tag_etamax (vvdouble jets, vdouble jetCSV, double CSVMwp);
		
		
		
		
		double study_tops_bb_syst (double MET, double METphi, TLorentzVector &metv, TLorentzVector lepton, vvdouble jets, vdouble csv, double &minChi, double &chi2lepW, double &chi2leptop, double &chi2hadW, double &chi2hadtop, double &mass_lepW, double &mass_leptop, double &mass_hadW, double &mass_hadtop, double &dRbb, double &testquant1, double &testquant2, double &testquant3, double &testquant4, double &testquant5, double &testquant6, double &testquant7, TLorentzVector &b1, TLorentzVector &b2, double CSVMwp);
		double getBestHiggsMass2(TLorentzVector lepton, TLorentzVector &met, vecTLorentzVector jets, vdouble btag, double &minChi, double &dRbb, TLorentzVector &bjet1, TLorentzVector &bjet2, double &chi2lepW, double &chi2leptop, double &chi2hadW, double &chi2hadtop, double &mass_lepW, double &mass_leptop, double &mass_hadW, double &mass_hadtop, TLorentzVector &toplep, TLorentzVector &tophad, double CSVMwp);
		double get_median_bb_mass(vvdouble jets, vdouble jetCSV, double CSVMwp);
		double pt_E_ratio_jets(vvdouble jets);
		
		double JetDelta_EtaAvgEta(vvdouble jet_vect_TLV, vdouble jet_CSV, std::string JetorTag, std::string JetorTag_Avg, double CSVMwp );

	
	private:

		// Parameter management
	private:
  
		// Old functions
	public:

	protected:
	
	//double CSVLwp, CSVMwp, CSVTwp;

	private:


	// === Variables === //
	public:

	protected:

	private:

}; // End of class prototype

#endif // _CommonBDTvars_h






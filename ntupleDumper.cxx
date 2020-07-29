#include "stdlib.h"
#include "ECALELF_ntuple_processor.cc"
#include <string>
#include "TString.h"

using namespace std;

int ntupleDumper(string m_basedir, TString jobname, Bool_t m_isData, string s_runlist="") {

	string time_bins = m_basedir + "inputs/Run2_UL_timebins_all.root";
	string lumi_file = m_basedir + "inputs/lumi_3years_csv.dat";

	if(jobname.Contains("2016")) {
	  lumi_file = m_basedir + "inputs/lumi_2016_csv.dat";
	  time_bins = m_basedir + "inputs/UL2016_timebins_all.root";
	}
	else if(jobname.Contains("2017")) {
	  lumi_file = m_basedir + "inputs/lumi_2017_csv.dat";
	  time_bins = m_basedir + "inputs/UL2017_timebins_all.root";
        }
        else if(jobname.Contains("2018")) {
 	  lumi_file = m_basedir + "inputs/lumi_2018_csv.dat";
	  time_bins = m_basedir + "inputs/UL2018_timebins_all.root";
        }

	int isFindBins = 0;
	string writedir, filein, fileout;
        TString m_runlist = s_runlist;

        filein = m_basedir +"inputs/"+ jobname + ".txt";
        writedir = m_basedir+"outputs/"+jobname + "/";
	if(m_runlist!="") fileout = writedir + jobname + "_"+m_runlist+".root";
	else fileout = writedir + jobname + "_all.root";

        mkdir(writedir);

	vector<UInt_t> m_runFilter;
	if(m_runlist!=""){
		if(m_runlist.Contains("_")){
			TObjArray *tx = m_runlist.Tokenize("_");
			tx->Print();
			for (Int_t i = 0; i < tx->GetEntries(); i++) {
				string sstring = ((TObjString *)(tx->At(i)))->String().Data();
				cout<<"Filter with run "<<sstring<<endl;
				UInt_t istring = stoi(sstring);
				m_runFilter.push_back(istring);
			}
		}
		else {
			cout<<"Filter with run "<<s_runlist<<endl;
			m_runFilter.push_back(stoi(s_runlist));
		}
	}

        ECALELF_ntuple_processor ecal(filein, fileout, m_isData, 10000, time_bins, lumi_file, 0, m_runFilter);
	ecal.analyseNtuples();
        return 1;
}

#ifndef ECALELF_NTUPLE_PROCESSOR
#define ECALELF_NTUPLE_PROCESSOR

#include "sstream"
#include "iostream"
#include "fstream"
#include "string.h"
#include "map"
#include "functional"
#include "algorithm"
#include "vector"
#include "math.h"
#include "bitset"
#include "TROOT.h"
#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"
#include "TClass.h"
#include "TTree.h"
#include "TChain.h"
#include "TColor.h"
#include "TDatime.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TBrowser.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "boost/function.hpp"
#include <boost/algorithm/string/replace.hpp>
#include "boost/filesystem.hpp"
#include "statistic.cc"
#include "ECALELF_ntuple_processor.h"

#include <sys/stat.h> // for mkdir

#include "extra_tools.cc"

Bool_t is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

ECALELF_ntuple_processor::ECALELF_ntuple_processor(std::string _ntuple_list_file, std::string _file_out, Bool_t _isData, UInt_t _eventsPerBin,
	std::string _time_bins_file, std::string _lumi_file, Bool_t _noTrain = 0, std::vector<UInt_t> _runFilter = {}):
ntuple_list_file(_ntuple_list_file), file_out(_file_out), isData(_isData),
eventsPerBin(_eventsPerBin), time_bins_file(_time_bins_file), lumi_file(_lumi_file), noTrain(_noTrain), runFilter(_runFilter) {

	std::cout << "[" << __LINE__ << "]" << "\t"  << " Ntuples list: " << ntuple_list_file << std::endl << "\t" <<
	" Output file: " << file_out << std::endl << "\t"  <<
	" Events PerBin: " << eventsPerBin << std::endl << "\t"  <<
	" Time bins file: " << time_bins_file << std::endl << "\t"  <<
	" Lumi file: " << lumi_file << std::endl << std::endl << std::endl;

        std::vector<std::string> ntuples = getNonemptyLinesWithFilenameKeyword(_ntuple_list_file, "EGamma-*.root");
        std::cout << "[" << __LINE__ << "]" << "Ntuples to analyse :" << endl;
        for (auto & it : ntuples) {
                ntuple_list.push_back(it);
                std::cout << "\t" << it << std::endl;
        }

        std::sort(ntuple_list.begin(), ntuple_list.end());

        std::vector<std::string> eleIDtrees = getNonemptyLinesWithFilenameKeyword(_ntuple_list_file, "recal_*.root");
        std::cout << "[" << __LINE__ << "]" << "recal Tree to analyse :" << endl;
        for (auto & it : eleIDtrees) {
                eleIDtree_list.push_back(it);
                std::cout << "\t" << it << std::endl;
        }

        std::sort(eleIDtree_list.begin(), eleIDtree_list.end());

        if(eleIDtree_list.size() != ntuples.size()) {
                cout<<"ERRROR! Unequal eleIDtrees and ntuples"<<endl;
                exit (EXIT_FAILURE);
        }

        fill_hists = (isData) ? (&ECALELF_ntuple_processor::FillAll) : (&ECALELF_ntuple_processor::Fill1D);
};

void ECALELF_ntuple_processor::analyseNtuples() {
	setHistograms();
	analyse();
	saveHistograms();
};

void ECALELF_ntuple_processor::getTimedGraphs() {
	std::cout << "[" << __LINE__ << "]" << " Making time evolution graphs from file " << file_out << std::endl;
	TFile file(file_out.c_str(), "READ");
	TIter next(file.GetListOfKeys());
	TKey *key;
	std::vector<TGraphAsymmErrors*> graphs;
	while ((key = (TKey*) next())) {
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom("TH2D")) continue;
		TH2D *hist = (TH2D*) key->ReadObj();
		std::string variable_title = hist->GetYaxis()->GetTitle();
		histogramTemplate *variable_info = NULL;
		for (auto & it : plot_variables) {
			if (variable_title.find(it.xtitle + it.xunit) != std::string::npos) {
				variable_info = &it;
				break;
			}
		}
		if (variable_info == NULL) continue;
		std::cout << "[" << __LINE__ << "]" << " Found TH2D " << hist->GetName() << std::endl;
		if (variable_info->makeMedianGraph) {
			std::cout << "[" << __LINE__ << "] \t" << " Making median vs time graph..." << std::endl;
			TGraphAsymmErrors* scale = stat.D2toXstatistic(hist, "median", 100.0, graphResamplings, 95.0);
			scale->GetXaxis()->SetTimeDisplay(1);
			scale->GetXaxis()->SetTimeFormat("#splitline{%d-%b}{%H:%M}%F1970-01-01 00:00:00s0");
			graphs.push_back(scale);
		}
		if (variable_info->makeSigma60Graph) {
			std::cout << "[" << __LINE__ << "] \t" << " Making StdDev60\% vs time graph..." << std::endl;
			TGraphAsymmErrors* res60 = stat.D2toXstatistic(hist, "stdev", 60.0, graphResamplings, 95.0);
			res60->GetXaxis()->SetTimeDisplay(1);
			res60->GetXaxis()->SetTimeFormat("#splitline{%d-%b}{%H:%M}%F1970-01-01 00:00:00s0");
			graphs.push_back(res60);
		}
	}
	file.Close();
	TFile file2(file_out.c_str(), "UPDATE");
	file2.cd();
	for (auto & it : graphs) {
		it->Write();
	}
	file2.Close();
	std::cout << "[" << __LINE__ << "] \t" << " Graphs written to file " << file_out << std::endl;
};


void ECALELF_ntuple_processor::getSummary() {
	TFile file(file_out.c_str(), "READ");
	TIter next(file.GetListOfKeys());
	TKey *key;

	std::string summary_file_name = file_out;
	boost::replace_all(summary_file_name, ".root", "_summary.txt");
	std::ofstream statsfile;
	statsfile.open(summary_file_name);

	errors *scale = NULL, *res60 = NULL, *res = NULL;
	Double_t x_val = -999, x_eL = -999, x_eH = -999, x_r60 = -999, x_r60_eL = -999, x_r60_eH = -999;

	statsfile << "category,scale,scale_errLow,scale_errHigh,resolution60, resolution60_errLow,resolution60_errHigh,stat" << std::endl;

	while ((key = (TKey*) next())) {
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom("TH1D")) continue;
		TH1D *hist = (TH1D*) key->ReadObj();

		std::string variable_title = hist->GetXaxis()->GetTitle();
		histogramTemplate *variable_info = NULL;
		for (auto & it : plot_variables) {
			if (variable_title.find(it.xtitle + it.xunit) != std::string::npos) {
				variable_info = &it;
				break;
			}
		}
		if (variable_info == NULL) continue;

		std::cout << "[" << __LINE__ << "] " << "Analysing " << hist->GetName() << std::endl;

		x_val = -999, x_eL = -999, x_eH = -999, x_r60 = -999, x_r60_eL = -999, x_r60_eH = -999;

		if (variable_info->makeMedianGraph) {
			scale = stat.getHistStatistic(hist, "median", 100.0, histResamplings, 95.0);
			if (scale != NULL) {
				x_val = scale->x;
				x_eL = scale->errLow;
				x_eH = scale->errHigh;
			}
		}
		if (variable_info->makeSigma60Graph) {
			res60 = stat.getHistStatistic(hist, "stdev", 60.0, histResamplings, 95.0);
			if (res60 != NULL) {
				x_r60 = res60->x;
				x_r60_eL = res60->errLow;
				x_r60_eH = res60->errHigh;
			}
		}

		statsfile << hist->GetName() << "," << x_val << "," <<
		x_eL << "," << x_eH << "," <<
		x_r60 << "," << x_r60_eL << "," <<
		x_r60_eH << "," <<
		std::fixed << hist->GetEntries() << std::endl;
		delete scale;
		scale = NULL;
		delete res60;
		res60 = NULL;
		hist->Delete();
	}

	statsfile.close();
	file.Close();
	std::cout << "[" << __LINE__ << "]" << " Summary written to file: " << summary_file_name << std::endl;
};

void ECALELF_ntuple_processor::readLumi() {
	std::cout << "[" << __LINE__ << "]" << " Reading lumi file..." << std::endl;
	CSVReader reader(lumi_file);
	std::vector<std::vector<std::string> > data = reader.getData();
	UInt_t low_edge, fillNum;

	UInt_t first_row;
	for (UInt_t it = 2; it < data.size() - 3; it++) {
		try {
			low_edge = std::stoi(data[it][2]);
		} catch (...) {
			std::cout << "[" << __LINE__ << "]" << " Error reading " << lumi_file << " in line" << it << std::endl;
			break;
		}
		if (low_edge >= begin_time - max_bin_width) {
			first_row = it;
			break;
		}
	}

	std::vector<Double_t> bins, bin_contents;
	bins.reserve(data.size() + 1);
	bin_contents.reserve(data.size() + 1);
	Double_t content, fillLum = 0;

	low_edge = std::stoi(data[first_row][2]);
	fillNum = std::stoi(data[first_row][0].erase(0, 7));
	fillMap[low_edge] = fillNum;

	std::cout << "[" << __LINE__ << "]" << " Fill info:" << std::endl;
	std::cout << "\t" << "Time" << "\t\t\t" << "Fill" << "\t" << "Luminosity[/pb]" << std::endl;
	for (UInt_t it = first_row + 1; it < data.size() - 3; it++) {
		try {
			low_edge = std::stoi(data[it][2]);
			fillNum = std::stoi(data[it][0].erase(0, 7));
			content = std::stod(data[it][6]);
			fillLum += content;

		} catch (...) {
			std::cout << "[" << __LINE__ << "]" << " Error reading " << lumi_file << " in line" << it << std::endl;
			break;
		}
		bin_contents.push_back(content);
		bins.push_back(low_edge);
		if (fillMap.rbegin()->second < fillNum) {
			fillMap[low_edge] = fillNum;
			auto rit = fillMap.rbegin();
			rit++;
			std::cout << "\t" << TDatime(rit->first).AsSQLString()  << "\t" << rit->second << "\t" << fillLum / 1000000 << std::endl;
			fillLum = 0;
		}
        // if (end_time + < bins.back()) break;
	}
	lumi_hist = new TH1D("Luminosity", "Luminosity", bins.size() - 1, bins.data());
	for (UInt_t i = 0; i < lumi_hist->GetNbinsX(); i++) {
		lumi_hist->SetBinContent(i + 1, bin_contents[i] / 1000000000);
	}
	lumi_hist->Sumw2();
	lumi_hist->GetXaxis()->SetTimeDisplay(1);
	lumi_hist->GetXaxis()->SetTimeFormat("#splitline{%d-%b}{%H:%M}%F1970-01-01 00:00:00s0");
	lumi_hist->GetXaxis()->SetNdivisions(510, kFALSE);
	lumi_hist->Scale(1000000, "width");
	lumi_hist->GetYaxis()->SetTitle("Luminosity [10^{33} cm^{-2} s^{-1}]");
	lumi_hist->GetYaxis()->CenterTitle();
	UInt_t color = TColor::GetColor(lumi_color.c_str());
	lumi_hist->SetFillColor(color);
	lumi_hist->SetLineColor(color);
	lumi_hist->SetMarkerSize(0);
	lumi_hist->SetStats(0);
	lumi_hist->SetDrawOption("bar hist");
	lumi_hist->SetMinimum(0.000001);
	lumi_hist->SetMaximum(30);
	lumi_hist->GetXaxis()->SetRangeUser(begin_time, end_time);
	lumi_hist->GetXaxis()->SetLimits(begin_time, end_time);
	std::cout << "[" << __LINE__ << "]" << " Lumisection times read from file " << lumi_file << std::endl;
};

UInt_t ECALELF_ntuple_processor::FindFill(UInt_t t) {
	for (auto & it : fillMap) {
		if (it.first >= t) {
			return it.second;
		}
	}
	cout << "[" << __LINE__ << "]"  << " Fill NOT found for time " << t << endl;
	return 0;
};

Bool_t ECALELF_ntuple_processor::InSameFill(Double_t t1, Double_t t2) {
	UInt_t _t1 = std::round(t1), _t2 = std::round(t2);
	UInt_t fill1 = FindFill(_t1), fill2 = FindFill(_t2);
	Bool_t IsSameFill = (fill1 == fill2) || (fill1 == 0) || (fill2 == 0);
	return IsSameFill;
};

void ECALELF_ntuple_processor::optiIseTimeBinning(TH1D *&tHist) {
	std::cout << "[" << __LINE__ << "]" << " Optimising time bins for " << tHist->GetName() << std::endl;
	std::vector<Double_t> bins, contents;
	UInt_t firstBin = tHist->FindFirstBinAbove(0, 1), lastBin = tHist->FindLastBinAbove(0, 1);
	bins.push_back(tHist->GetXaxis()->GetBinLowEdge(firstBin));
	Double_t binSum = 0;
	Bool_t metThreshold, tooWideBin, differentFill, insertBin;
	for (UInt_t j = firstBin; j < lastBin + 1; j++) {
		binSum += tHist->GetBinContent(j);
		UInt_t lowEdge = round(bins.back());
		UInt_t upEdgeCand = round(tHist->GetXaxis()->GetBinUpEdge(j));

		metThreshold = (binSum >= eventsPerBin);
		tooWideBin = (upEdgeCand - lowEdge > max_bin_width);

		insertBin = metThreshold || tooWideBin;

		if (insertBin) {
			bins.push_back(upEdgeCand);
			contents.push_back(binSum);
			binSum = 0;
		}
	}
	if (bins.back() < tHist->GetXaxis()->GetBinUpEdge(lastBin)) {
		bins.push_back(tHist->GetXaxis()->GetBinUpEdge(lastBin));
	};
	TH1D* outHist = new TH1D("", "", bins.size() - 1, bins.data());
	outHist->GetXaxis()->CenterTitle(true);
	outHist->GetXaxis()->SetTimeDisplay(1);
	outHist->GetXaxis()->SetLabelOffset(0.03);;
	outHist->GetXaxis()->SetTimeFormat("#splitline{%d-%b}{%H:%M}%F1970-01-01 00:00:00s0");
	outHist->GetYaxis()->CenterTitle(true);
	outHist->GetYaxis()->SetTitle("# of events");
	for (UInt_t i = 0; i < contents.size(); i++) {
		outHist->SetBinContent(i + 1, contents[i]);
	};
	std::string name = tHist->GetName();
	std::string title = tHist->GetTitle();
	tHist->Delete();
	tHist = NULL;
	outHist->SetName(name.c_str());
	outHist->SetTitle(title.c_str());
	tHist = outHist;
	std::cout << "[" << __LINE__ << "]" << name << ": Number of bins = " << bins.size() - 1 << std::endl;
};

void ECALELF_ntuple_processor::constructTimeHistograms() {
	std::cout << "[" << __LINE__ << "]" << "Making fine-binned time histograms... " << std::endl;
	for (auto & itC : categories) {
		TH1D *hist = new TH1D(itC.name.c_str(), itC.title.c_str(),  end_time - begin_time, begin_time, end_time);
		itC.stats = hist;
	}
};

void ECALELF_ntuple_processor::setHistograms() {
	TFile *bins_file = NULL;
	if (isData) {
		bins_file = new TFile(time_bins_file.c_str(), "READ");
		if (bins_file->IsZombie()) {
			std::cout << "[" << __LINE__ << "]" << " Zombie time bins file " << time_bins_file << std::endl;
			exit (EXIT_FAILURE);
		}
	}
	TH1D *timeBinHist = NULL;

	outFile = new TFile(file_out.c_str(), "RECREATE");
	for (auto & _cat : categories) {
		std::vector<Double_t> t_bins;
		if (isData) {
			std::cout<<"cat "<<_cat.name<<std::endl;
			timeBinHist = (TH1D*) bins_file->Get(_cat.name.c_str());
			t_bins = stat.histBins(timeBinHist);
			std::cout << "[" << __LINE__ << "]" << " # of bins  in " << _cat.name << ": " << t_bins.size() <<std::endl<<"\t\tevents= "<< timeBinHist->Integral()<< std::endl;
			timeBinHist->Delete();
			timeBinHist = NULL;
		}
		for (auto & _var : plot_variables) {
			TH1D* hist1D = setVarHistograms(_cat, _var);
			hist1D->SetDirectory(outFile->GetDirectory(""));
			TH2D* hist2D = NULL;
			if (isData){
				hist2D = setTimedHistogram(_cat, _var, t_bins);
				hist2D->SetDirectory(outFile->GetDirectory(""));
			}
			histogram_set *hist_set = new histogram_set(_var.var, hist1D, hist2D);
			_cat.histograms.push_back(hist_set);
		}
	}
	if (isData) bins_file->Close();
};

TH2D *ECALELF_ntuple_processor::setTimedHistogram(category &cat, const histogramTemplate &_template, const std::vector<Double_t> &t_bins) {
	std::string name = stat.removeNonAlpha(_template.xtitle) + "_vs_time_in_" + cat.name;
	TH2D *hist2D = new TH2D(name.c_str(), cat.title.c_str(), t_bins.size() - 1, t_bins.data(),  _template.nbins, _template.xmin, _template.xmax);
	hist2D->GetYaxis()->SetTitle((_template.xtitle + _template.xunit).c_str());
	hist2D->GetXaxis()->SetTitle("Time");
	hist2D->GetXaxis()->SetTimeDisplay(1);
	hist2D->GetXaxis()->SetTimeFormat("#splitline{%d-%b}{%H:%M}%F1970-01-01 00:00:00s0");
	std::cout << "[" << __LINE__ << "]" << "TH2D created: " << name << std::endl;
	return hist2D;
};

TH1D *ECALELF_ntuple_processor::setVarHistograms(category &cat, const histogramTemplate &_template) {
	std::string name = stat.removeNonAlpha(_template.xtitle) + "_in_" + cat.name;
	TH1D *hist1D = new TH1D(name.c_str(), cat.title.c_str(), _template.nbins, _template.xmin, _template.xmax);
	hist1D->GetXaxis()->SetTitle((_template.xtitle + _template.xunit).c_str());
	hist1D->GetYaxis()->SetTitle(("# of events/" + _template.xdiv).c_str());
	std::cout << "[" << __LINE__ << "]" << "TH1D created: " << name << std::endl;
	return hist1D;
};

void ECALELF_ntuple_processor::Fill1D() {
	for (auto & _cat : categories) {
		if (!(*_cat.Pass)) continue;
		for (auto & _var : _cat.histograms) {
			_var->hist->Fill(*(_var->var), mcGenWeight);
		}
	}
};

void ECALELF_ntuple_processor::FillAll() {
	for (auto & _cat : categories) {
		if (!(*_cat.Pass)) continue;
		for (auto & _var : _cat.histograms) {
			_var->hist->Fill(*(_var->var), mcGenWeight);
			_var->histTimed->Fill(eventTime, *(_var->var));
		}
	}
};

void ECALELF_ntuple_processor::analyse() {

        ntupleChain = new TChain("selected");
        for (auto &it : ntuple_list) {
                ntupleChain->Add(it.c_str());
                std::cout<<"\tAdded ntuple "<<it<<std::endl;
        }

	ntupleChain->SetBranchStatus("runNumber", 1);
	ntupleChain->SetBranchAddress("runNumber", &runNumber);
        ntupleChain->SetBranchStatus("eventNumber", 1);
	ntupleChain->SetBranchAddress("eventNumber", &eventNumber);
	ntupleChain->SetBranchStatus("eventTime", 1);
	ntupleChain->SetBranchAddress("eventTime", &eventTime);

	if (!isData){
		ntupleChain->SetBranchStatus("mcGenWeight", 1);
		ntupleChain->SetBranchAddress("mcGenWeight", &mcGenWeight);
	}

	ntupleChain->SetBranchStatus("eleID", 1);
	ntupleChain->SetBranchAddress("eleID", &eleID);

	ntupleChain->SetBranchStatus("energy_ECAL_ele", 1);
	ntupleChain->SetBranchAddress("energy_ECAL_ele", &energy_ECAL_ele);

	ntupleChain->SetBranchStatus("pAtVtxGsfEle", 1);
	ntupleChain->SetBranchAddress("pAtVtxGsfEle", &pAtVtxGsfEle);

	ntupleChain->SetBranchStatus("R9Ele", 1);
	ntupleChain->SetBranchAddress("R9Ele", &R9Ele);

	ntupleChain->SetBranchStatus("etaEle", 1);
	ntupleChain->SetBranchAddress("etaEle", &etaEle);

	ntupleChain->SetBranchStatus("etaSCEle", 1);
	ntupleChain->SetBranchAddress("etaSCEle", &etaSCEle);

	ntupleChain->SetBranchStatus("phiSCEle", 1);
	ntupleChain->SetBranchAddress("phiSCEle", &phiSCEle);

	ntupleChain->SetBranchStatus("phiEle", 1);
	ntupleChain->SetBranchAddress("phiEle", &phiEle);

	ntupleChain->SetBranchStatus("invMass_ECAL_ele", 1);
	ntupleChain->SetBranchAddress("invMass_ECAL_ele", &invMass_ECAL_ele);

	ntupleChain->SetBranchStatus("xSeedSC", 1);
	ntupleChain->SetBranchAddress("xSeedSC", &xSeedSC);

	ntupleChain->SetBranchStatus("ySeedSC", 1);
	ntupleChain->SetBranchAddress("ySeedSC", &ySeedSC);

	ntupleChain->SetBranchStatus("invMass", 1);
	ntupleChain->SetBranchAddress("invMass", &invMass);

        eleIDchain = new TChain("recal_tree");
        for (auto &it : eleIDtree_list) {
                eleIDchain->Add(it.c_str());
                std::cout<<"\tAdded recal_tree "<<it<<std::endl;
        }

        eleIDchain->SetBranchStatus("runNumber", 1);
        eleIDchain->SetBranchAddress("runNumber", &runNumber2);
        eleIDchain->SetBranchStatus("eventNumber", 1);
        eleIDchain->SetBranchAddress("eventNumber", &eventNumber2);

        eleIDchain->SetBranchStatus("invMass_ECAL_ele_recal", 1);
        eleIDchain->SetBranchAddress("invMass_ECAL_ele_recal", &invMass_ECAL_ele_recal);

//	eleIDchain->Show(10);
//	eleIDchain->BuildIndex("runNumber","eventNumber");
        ntupleChain->AddFriend(eleIDchain);

//        eleIDchain->Show(10);
        std::cout << "[" << __LINE__ << "]" << " Analysing ntuples " <<std::endl;
        analyseTree(ntupleChain);

};

void ECALELF_ntuple_processor::analyseTree(TChain *tree) {
	for (UInt_t i = 0; i < tree->GetEntries(); i++) {
		tree->GetEntry(i);
		if (!ApplyCuts()) continue;
		(this->*fill_hists)();
	}

};

void ECALELF_ntuple_processor::setNtupleTree(TChain * _tree) {
//	_tree->SetBranchStatus("*", 0);

	_tree->SetBranchStatus("runNumber", 1);
	_tree->SetBranchAddress("runNumber", &runNumber);
    //_tree->SetBranchAddress("lumiBlock", &lumiBlock);
	_tree->SetBranchAddress("eventNumber", &eventNumber);
	_tree->SetBranchStatus("eventTime", 1);
	_tree->SetBranchAddress("eventTime", &eventTime);

	if (!isData){
		_tree->SetBranchStatus("mcGenWeight", 1);
		_tree->SetBranchAddress("mcGenWeight", &mcGenWeight);
	}

//	_tree->SetBranchStatus("isTrain", 1);
//	_tree->SetBranchAddress("isTrain", &isTrain);


    //_tree->SetBranchAddress("nPU", &nPU);
	_tree->SetBranchStatus("eleID", 1);
	_tree->SetBranchAddress("eleID", &eleID);

	_tree->SetBranchStatus("energy_ECAL_ele", 1);
	_tree->SetBranchAddress("energy_ECAL_ele", &energy_ECAL_ele);

	_tree->SetBranchStatus("pAtVtxGsfEle", 1);
	_tree->SetBranchAddress("pAtVtxGsfEle", &pAtVtxGsfEle);

	_tree->SetBranchStatus("R9Ele", 1);
	_tree->SetBranchAddress("R9Ele", &R9Ele);

	_tree->SetBranchStatus("etaEle", 1);
	_tree->SetBranchAddress("etaEle", &etaEle);

	_tree->SetBranchStatus("etaSCEle", 1);
	_tree->SetBranchAddress("etaSCEle", &etaSCEle);

	_tree->SetBranchStatus("phiSCEle", 1);
	_tree->SetBranchAddress("phiSCEle", &phiSCEle);

	_tree->SetBranchStatus("phiEle", 1);
	_tree->SetBranchAddress("phiEle", &phiEle);

	_tree->SetBranchStatus("invMass_ECAL_ele", 1);
	_tree->SetBranchAddress("invMass_ECAL_ele", &invMass_ECAL_ele);

	_tree->SetBranchStatus("xSeedSC", 1);
	_tree->SetBranchAddress("xSeedSC", &xSeedSC);

	_tree->SetBranchStatus("ySeedSC", 1);
	_tree->SetBranchAddress("ySeedSC", &ySeedSC);

	_tree->SetBranchStatus("invMass", 1);
	_tree->SetBranchAddress("invMass", &invMass);
};


void ECALELF_ntuple_processor::setEleIDTree(TChain * _tree) {
        _tree->SetBranchStatus("*", 0);

        _tree->SetBranchStatus("invMass_ECAL_ele_recal", 1);
        _tree->SetBranchAddress("invMass_ECAL_ele_recal", &invMass_ECAL_ele_recal);

};


void ECALELF_ntuple_processor::saveHistograms() {
	outFile->Write();
	outFile->Close();
	std::cout << "[" << __LINE__ << "]" << " Histograms saved to file: " << file_out << std::endl;
};

Bool_t ECALELF_ntuple_processor::ApplyCuts() {
	mPassBasicCuts = 0;

	if(!runFilter.empty()){
		if(std::find(runFilter.begin(), runFilter.end(), runNumber) == runFilter.end()) return mPassBasicCuts;
	}

	if (!((eleID[0] & eleIDmask) && (eleID[1] & eleIDmask))) return mPassBasicCuts;

	Et[0] = pAtVtxGsfEle[0] / std::cosh(etaEle[0]);
	Et[1] = pAtVtxGsfEle[1] / std::cosh(etaEle[1]);

	if ((Et[0] < EtLowLimit) || (Et[1] < EtLowLimit)) return mPassBasicCuts;

	iL = (Et[0] > Et[1]) ? 0 : 1;
	iS = 1 - iL;

	EtL = Et[iL];
	EtS = Et[iS];

	etaEle_eL = etaSCEle[iL];
	phiEle_eL = phiSCEle[iL];
	absEtaL = std::abs(etaSCEle[iL]);
	absEtaS = std::abs(etaSCEle[iS]);

	if(absEtaL > 2.5 || absEtaS > 2.5) return mPassBasicCuts;

	if (invMass < 70. || invMass > 110.) return mPassBasicCuts;
        
	R9L = R9Ele[iL];
	R9S = R9Ele[iS];
	ebL = (absEtaL < BETRetaMin);
	ebS = (absEtaS < BETRetaMin);
	hR9L = (R9L > highR9);
	hR9S = (R9S > highR9);
	EoPL = energy_ECAL_ele[iL] / pAtVtxGsfEle[iL];
	EoPS = energy_ECAL_ele[iS] / pAtVtxGsfEle[iS];

	etaSC_L = etaSCEle[iL];
	etaSC_S = etaSCEle[iS];

	phiSC_L = phiSCEle[iL];
	phiSC_S = phiSCEle[iS];

	zSide[0] = (etaSCEle[0] > 0) ? 1 : -1;
	zSide[1] = (etaSCEle[1] > 0) ? 1 : -1;


	mPassBasicCuts = 1;

	IsHR9 = (R9L > 0.96) && (R9S > 0.96);
	IsLR9 = (R9L < 0.96) && (R9S < 0.96);

	IsEb = (absEtaL < BETRetaMin) && (absEtaS < BETRetaMin);
	IsEbHR9 = IsEb && IsHR9;
	IsEbLR9 = IsEb && IsLR9;

	IsEe = (absEtaL > BETRetaMax) && (absEtaS > BETRetaMax) && (absEtaL < InnerEta) && (absEtaS < InnerEta);
	IsEeHR9 = IsEe && IsHR9;
	IsEeLR9 = IsEe && IsLR9;

	eL_IsInEe = (absEtaL > BETRetaMax) && (absEtaL < InnerEta) ;
	eS_IsInEb= (absEtaS < BETRetaMin);

	return mPassBasicCuts;
};


ECALELF_ntuple_processor::~ECALELF_ntuple_processor() {
};

#endif

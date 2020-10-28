#include "extra_tools.cc"
#include "statistic.cc"
#include "TStyle.h"
#include "TROOT.h"
#include "THStack.h"
#include <TLegend.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TF1.h>
#include "TText.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TBox.h"
#include "TGaxis.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string histOptionsFile="histOptions.txt";
parseOptions histOpts;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

string m_writeDir="outputs/plots/";


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plotHist(std::string _histName, TString sample1, TString sample2);
void plotAll();
statistic st;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plotHist(std::string _histName, TString sample1, TString sample2){
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	sample ul2018eop2dIOVN(sample1, "eop_2d_median_IOVN", 21, "#238b45");
	sample ul2018eop2dIOVNm1(sample2, "eop_2d_median_IOVNm1", 22, "#023858");

	std::vector<sample*> data = {&ul2018eop2dIOVNm1};
	std::vector<sample*> mc = {&ul2018eop2dIOVN};

	std::vector<TObject*> _deleteThese;

	std::cout<<"Plotting "<<_histName<<std::endl;

	TCanvas canvas((_histName + "_canvas").c_str(), "", 2400, 1800);
	canvas.Draw();

	TPad pad0("pad0", "", histOpts.getFloat("pad0x1"), histOpts.getFloat("pad0y1"), histOpts.getFloat("pad0x2"), histOpts.getFloat("pad0y2"));
	TPad pad1("pad1", "", histOpts.getFloat("pad1x1"), histOpts.getFloat("pad1y1"), histOpts.getFloat("pad1x2"), histOpts.getFloat("pad1y2"));
	TPad pad2("pad2", "", histOpts.getFloat("pad2x1"), histOpts.getFloat("pad2y1"), histOpts.getFloat("pad2x2"), histOpts.getFloat("pad2y2"));
	TPad pad3("pad3", "", histOpts.getFloat("pad3x1"), histOpts.getFloat("pad3y1"), histOpts.getFloat("pad3x2"), histOpts.getFloat("pad3y2"));

	pad0.SetMargin(histOpts.getFloat("pad0marginL"), histOpts.getFloat("pad0marginR"), histOpts.getFloat("pad0marginB"), histOpts.getFloat("pad0marginT"));
	pad1.SetMargin(histOpts.getFloat("pad1marginL"), histOpts.getFloat("pad1marginR"), histOpts.getFloat("pad1marginB"), histOpts.getFloat("pad1marginT"));
	pad2.SetMargin(histOpts.getFloat("pad2marginL"), histOpts.getFloat("pad2marginR"), histOpts.getFloat("pad2marginB"), histOpts.getFloat("pad2marginT"));
	pad3.SetMargin(histOpts.getFloat("pad3marginL"), histOpts.getFloat("pad3marginR"), histOpts.getFloat("pad3marginB"), histOpts.getFloat("pad3marginT"));

	pad0.SetFillStyle(0);
	pad1.SetFillStyle(0);
	pad2.SetFillStyle(0);
	pad3.SetFillStyle(0);

	pad0.SetFillColor(0);
	pad1.SetFillColor(0);
	pad2.SetFillColor(0);
	pad3.SetFillColor(0);

	pad0.SetFrameFillStyle(0);
	pad1.SetFrameFillStyle(0);
	pad2.SetFrameFillStyle(0);
	pad3.SetFrameFillStyle(0);

	canvas.cd();


	pad1.Draw();
	pad1.cd();
	TLegend legend(histOpts.getFloat("pad2marginL")-0.04, 0.0, 0.99, 1.0);
	legend.SetTextSize(histOpts.getFloat("legendTextSize"));
	legend.SetNColumns(histOpts.getInt("legendNcols"));
	legend.SetFillStyle(1001);
	legend.SetBorderSize(0);
	legend.Draw();


	std::string _title, _xTitle, _yTitle;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	THStack mainStack("mainStack","");

	//Add mc
	for(Int_t index = 0; index < mc.size(); index++){
		sample * _mc = mc[index];

		TH1D * hist_mc = (TH1D*) getHistFromFile(_histName, _mc->ntuple, 0);

		Float_t meanHist = hist_mc->GetMean();
		Float_t sigma60 = st.getHistStatistic(hist_mc, "stdev", 60.);
		cout<<"data sigma60 is "<<sigma60<<endl;


		hist_mc->Rebin(5);
//		if(_histName.find("mee") != std::string::npos) hist_mc->Rebin(5);
//		if(_histName.find("R9eL") != std::string::npos) hist_mc->Rebin(10);
//		if(_histName.find("etaSC") != std::string::npos) hist_mc->Rebin(10);



		_mc->assignAtt(hist_mc, histOpts.getFloat("markerSize"), histOpts.getFloat("lineWidth"));

		hist_mc->Scale(histOpts.getFloat("scale")/hist_mc->Integral());

		mainStack.Add(hist_mc, "HIST");

		_deleteThese.push_back(hist_mc);

//                std::string legend_text = _mc->legend + " (#mu=" + to_string_with_precision(meanHist, 2) + ", #sigma^{60}=" + to_string_with_precision(sigma60, 4) + ")";
		std::string legend_text = _mc->legend;
		TLegendEntry* legendEntry =  legend.AddEntry(hist_mc, legend_text.c_str(), "L");
		legendEntry->SetTextColor(TColor::GetColor(_mc->color.c_str()));

	}

		//Add data
	for(Int_t index = 0; index < data.size(); index++){
		sample * _data = data[index];

		std::string tmp_histName = _histName;
		// if(_data->ntuple.find("March") != std::string::npos) boost::replace_all(_histName, "R92", "R9");

		TH1D * hist_data = (TH1D*) getHistFromFile(_histName, _data->ntuple, 0);

		// hist_data->Scale(-1.);

		Float_t meanHist = hist_data->GetMean();
		Float_t sigma60 = st.getHistStatistic(hist_data, "stdev", 60.);
                cout<<"mc sigma60 is "<<sigma60<<endl;

		hist_data->Rebin(5);
//		if(_histName.find("mee") != std::string::npos) hist_data->Rebin(5);
//		if(_histName.find("R9eL") != std::string::npos) hist_data->Rebin(10);
//		if(_histName.find("etaSC") != std::string::npos) hist_data->Rebin(10);

		_data->assignAtt(hist_data, histOpts.getFloat("markerSize"), histOpts.getFloat("lineWidth"));
		hist_data->SetFillStyle(histOpts.getInt("errorFillStyle"));
		hist_data->SetFillColor(TColor::GetColor(_data->color.c_str()));
		hist_data->Scale(histOpts.getFloat("scale")/hist_data->Integral());

		mainStack.Add(hist_data, "PE");
		mainStack.Add(hist_data, "E2");

		_deleteThese.push_back(hist_data);

		if(index==0){
			_title = hist_data->GetTitle();
			_xTitle = hist_data->GetXaxis()->GetTitle();
			_yTitle = hist_data->GetYaxis()->GetTitle();
		}

//                std::string legend_text = _data->legend + " (#mu=" + to_string_with_precision(meanHist, 2) + ", #sigma^{60}=" + to_string_with_precision(sigma60, 4) + ")";
		std::string legend_text = _data->legend;
		TLegendEntry* legendEntry =  legend.AddEntry(hist_data, legend_text.c_str(), "LFPE");
		legendEntry->SetTextColor(TColor::GetColor(_data->color.c_str()));
		_histName = tmp_histName;
	}

	canvas.cd();

	std::vector<Float_t> xLimits = getXlimits({(TH1D*) mainStack.GetStack()->Last()}, histOpts.getFloat("binThreshold"));

	pad2.Draw();
	pad2.cd();
	pad2.SetGrid(1,1);
	mainStack.Draw("NOSTACK");


	mainStack.GetXaxis()->SetRangeUser(xLimits[1], xLimits[0]);
	mainStack.GetXaxis()->SetLabelSize(0);
	mainStack.GetXaxis()->SetTitleSize(0);

	mainStack.GetYaxis()->SetTitle(_yTitle.c_str());
	((TGaxis*)mainStack.GetYaxis())->SetMaxDigits(histOpts.getInt("pad2yMaxDigits"));
	mainStack.GetYaxis()->SetLabelSize(histOpts.getFloat("pad2axislabelsize"));
	mainStack.GetYaxis()->SetTitleSize(histOpts.getFloat("pad2axisTitleSize"));
	mainStack.GetYaxis()->SetTitleOffset(histOpts.getFloat("pad2ytitleoffset"));
	mainStack.GetYaxis()->SetNdivisions(histOpts.getInt("nDivsPad2Y"), kTRUE);
	mainStack.GetXaxis()->SetNdivisions(histOpts.getInt("nDivsPad2X"), kTRUE);

	if(_histName.find("R9eL") != std::string::npos) pad2.SetLogy();
	if(_histName.find("R9eS") != std::string::npos) pad2.SetLogy();
	if(_histName.find("HEL") != std::string::npos) pad2.SetLogy();
	if(_histName.find("sigmaietaietaL") != std::string::npos) pad2.SetLogy();
	if(_histName.find("SigmaETECALRecHit") != std::string::npos) pad2.SetLogy();
        if(_histName.find("hOverE") != std::string::npos) pad2.SetLogy();


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	THStack ratioStack("ratioStack","");


	Double_t _ratioSpread = 0.2;
	Double_t _ratioMean = 1.;

	for(Int_t index = 0; index < data.size(); index++){

		TH1D * hist_data = (TH1D*) getHistFromFile(_histName, data[index]->ntuple, 0);

		std::string _ratioHistName =(std::string) "data_mc_ratio_"  + (std::string) hist_data->GetName();
		hist_data->SetName(_ratioHistName.c_str());

		hist_data->Rebin(5);
//		if(_histName.find("mee") != std::string::npos) hist_data->Rebin(5);
//		if(_histName.find("R9eL") != std::string::npos) hist_data->Rebin(10);
//		if(_histName.find("etaSC") != std::string::npos) hist_data->Rebin(10);

		hist_data->Scale(histOpts.getFloat("scale")/hist_data->Integral());
		data[index]->assignAtt(hist_data, histOpts.getFloat("markerSize"), histOpts.getFloat("lineWidth"));


		TH1D * hist_mc = (TH1D*) getHistFromFile(_histName, mc.at(0)->ntuple, 0);

		hist_mc->Rebin(5);
//		if(_histName.find("mee") != std::string::npos) hist_mc->Rebin(5);
//		if(_histName.find("R9eL") != std::string::npos) hist_mc->Rebin(10);
//		if(_histName.find("etaSC") != std::string::npos) hist_mc->Rebin(10);

		hist_mc->Scale(histOpts.getFloat("scale")/hist_mc->Integral());
		hist_data->Divide(hist_mc);
		hist_mc->SetName("to_delete");
		hist_mc->Delete();

		hist_data->SetLineColor(TColor::GetColor(data.at(index)->color.c_str()));
		hist_data->SetMarkerColor(TColor::GetColor(data.at(index)->color.c_str()));
		hist_data->SetFillStyle(histOpts.getInt("errorFillStyle"));
		hist_data->SetFillColor(TColor::GetColor(data.at(index)->color.c_str()));

		_deleteThese.push_back(hist_data);

		ratioStack.Add(hist_data, "E2");
		ratioStack.Add(hist_data, "PE");

		Float_t tmpSpread =weightedYspread(hist_data);

		Float_t tmpMean =weightedYmean(hist_data);

		_ratioSpread = (_ratioSpread < tmpSpread) ? tmpSpread : _ratioSpread;
		_ratioMean = (_ratioMean < tmpMean) ? tmpMean : _ratioMean;
	}

	Double_t ratioMin = _ratioMean - histOpts.getFloat("ratioNspread") * _ratioSpread;
	Double_t ratioMax = _ratioMean + histOpts.getFloat("ratioNspread") * _ratioSpread;

	if(ratioMin < histOpts.getFloat("ratioMin")) ratioMin = histOpts.getFloat("ratioMin");
	if(ratioMin > histOpts.getFloat("ratioAtLeast")) ratioMin = histOpts.getFloat("ratioAtLeast");

	if(ratioMax > histOpts.getFloat("ratioMax")) ratioMax = histOpts.getFloat("ratioMax");
	if(ratioMax < histOpts.getFloat("ratioAtMost")) ratioMax = histOpts.getFloat("ratioAtMost");

	if(ratioMin != ratioMin) ratioMin = histOpts.getFloat("ratioMin");
	if(ratioMax != ratioMax) ratioMax = histOpts.getFloat("ratioMax");

	cout<<"\t\t Rmin = "<<ratioMin<<"\tRmax"<<ratioMax<<endl;

	canvas.cd();
	pad3.Draw();
	pad3.cd();
	pad3.SetGrid(1,1);

	TF1 unityLine("unityLine","1.",xLimits[1],xLimits[0]);
	unityLine.SetLineWidth(histOpts.getFloat("lineWidth"));
	unityLine.SetLineColor(TColor::GetColor(histOpts.get("unitylinecolor").c_str()));
	unityLine.Draw();

	ratioStack.Draw("SAME NOSTACK A");

	unityLine.GetXaxis()->SetTitle(_xTitle.c_str());
	unityLine.GetXaxis()->SetTitleSize(histOpts.getFloat("pad3axisTitleSize"));
	unityLine.GetXaxis()->SetTitleOffset(histOpts.getFloat("pad3xtitleoffset"));
	unityLine.GetXaxis()->SetLabelSize(histOpts.getFloat("pad3axislabelsize"));
	unityLine.GetXaxis()->SetRangeUser(xLimits[1], xLimits[0]);
	unityLine.GetYaxis()->SetTitle(histOpts.get("ratioPlotYtitle").c_str());
	unityLine.GetYaxis()->SetTitleSize(histOpts.getFloat("pad3axisTitleSize"));
	unityLine.GetYaxis()->SetTitleOffset(histOpts.getFloat("pad3ytitleoffset"));
	unityLine.GetYaxis()->SetLabelSize(histOpts.getFloat("pad3axislabelsize"));
	unityLine.GetYaxis()->SetRangeUser(ratioMin, ratioMax);
	unityLine.GetYaxis()->SetLimits(ratioMin, ratioMax);
	unityLine.GetYaxis()->SetNdivisions(histOpts.getInt("nDivsPad3Y"), kTRUE);
	unityLine.GetXaxis()->SetNdivisions(histOpts.getInt("nDivsPad3X"), kTRUE);



	canvas.cd();
	pad0.Draw();
	pad0.cd();
	TPaveText pt(0.0, .0, 1., 1.);
	pt.SetTextSize(histOpts.getFloat("titleSize"));
	pt.AddText(_title.c_str());
	pt.Draw();

	gPad->RedrawAxis();
	gPad->RedrawAxis("g");
	gPad->Update();

	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();


//        std::string _writePath = histOpts.get("writeDir") + "/" + _histName + ".png";
	std::string _writePath = m_writeDir + "/" + _histName + ".png";


	TBox canvasFrame(0.,0.,1.,1.);
	canvasFrame.SetFillStyle(0);
	// canvasFrame.SetFillColor(0);
	canvasFrame.SetLineColor(kBlack);
	canvasFrame.SetLineWidth(histOpts.getFloat("lineWidth"));
	canvas.cd();
	canvasFrame.Draw();

	canvas.SaveAs(_writePath.c_str());

	for(auto & toDelete: _deleteThese){
		toDelete->Delete();
	}
};


void plotAll(TString sample1, TString sample2){

	gROOT->SetBatch();

	histOpts.parseIt(histOptionsFile);

//	mkdir(histOpts.get("writeDir"));
        mkdir(m_writeDir);

	std::vector<std::string> histList = getObjectList(data[0]->ntuple, "TH1D");
	// plotHist(histList[0]);
	for(auto & hist : histList){
		if(hist.find("S4") != std::string::npos || hist.find("s4") != std::string::npos) continue;
		plotHist(hist, sample1, sample2);
	}
};

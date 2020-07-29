#define BETRetaMin 1.4442
#define BETRetaMax 1.566
#define MidEta 2.0
#define InnerEta 2.5
#define EtLowLimit 20.0
#define max_bin_width 7.*24.*3600.
#define highR9 0.94
#define eleIDmask 0x0004 //loose94XV1Run2017
#define graphResamplings 200
#define histResamplings 200

std::string lumi_color = "#F4D03F";

struct histogramTemplate {
	Float_t *var = nullptr;
	Float_t xmin;
	Float_t xmax;
	UInt_t nbins;
	Bool_t makeMedianGraph;
	Bool_t makeSigma60Graph;
	std::string xtitle;
	std::string xunit;
	std::string xdiv;
	std::string ytitle;

	histogramTemplate(Float_t *_var, Float_t _xmin, Float_t _xmax, UInt_t _nbins, Bool_t _makeMedianGraph = 1, Bool_t _makeSigma60Graph = 0, std::string _xtitle = "", std::string _xunit = "",
		std::string _xdiv = "", std::string _ytitle = "") :
	var(_var), xmin(_xmin), xmax(_xmax), nbins(_nbins), makeMedianGraph(_makeMedianGraph), makeSigma60Graph(_makeSigma60Graph), xtitle(_xtitle), xunit(_xunit), xdiv(_xdiv),
	ytitle(_ytitle) {
	};

	histogramTemplate() {
	};

	~histogramTemplate() {
	};

	std::vector<Float_t> getBinVector() {
		Float_t binsize = (xmax - xmin) / (Float_t) nbins;
		std::vector<Float_t> binEdges;
		if (binsize < 0) return {};
		for (UInt_t i = 0; i < nbins + 1; i++) {
			binEdges.push_back(xmin + i * binsize);
		}
		return binEdges;
	};
};

struct histogram_set {
	Float_t *var = nullptr;
	TH1D* hist = nullptr;
	TH2D* histTimed = nullptr;

	histogram_set(Float_t *_var, TH1D* _hist = nullptr, TH2D* _histTimed = nullptr) :
	var(_var), hist(_hist), histTimed(_histTimed) {
	};

	histogram_set() {
	};

	~histogram_set() {
	};

};

struct category {
	std::string name;
	std::string title;
	//PassCategory selector;
	Bool_t *Pass;
	TH1D *stats;
	category(std::string _name, std::string _title, Bool_t *_Pass, TH1D *_stats = NULL): name(_name), title(_title), Pass(_Pass), stats(_stats) {
	};
	category() {};
	std::vector<histogram_set*> histograms;
};


class ECALELF_ntuple_processor {
public:
	ECALELF_ntuple_processor(std::string _ntuple_list_file, std::string _file_out, Bool_t _isData, UInt_t _eventsPerBin, std::string _time_bins_file, std::string _lumi_file, Bool_t _noTrain, std::vector<UInt_t> _runFilter);
	std::vector<std::string> ntuple_list;
        std::vector<std::string> eleIDtree_list;
        Bool_t isData, noTrain;
	~ECALELF_ntuple_processor();

	void findTimeBins();
	void analyseNtuples();
	void getTimedGraphs();
	void getSummary();

	//typedef Bool_t(ECALELF_ntuple_processor::*PassCategory)();

private:
	/************Inputs*************/
	std::string ntuple_list_file;
        std::string eleIDtree_list_file;
	std::string time_bins_file;
	std::string lumi_file;
	std::string file_out;
	UInt_t eventsPerBin;
	std::vector<UInt_t> runFilter;
	/******************************/

	TFile *outFile = nullptr;

	void getTimeLimits();
	UInt_t begin_time, end_time;

	void readLumi();
	TH1D* get_luminosity_hist();
	TH1D* lumi_hist = NULL;
	std::map<UInt_t, UInt_t> fillMap;

	void setNtupleTree(TChain* _tree);
        void setEleIDTree(TChain* _tree);
	UInt_t eventTime;
	UInt_t runNumber;
        UInt_t runNumber2;
        ULong64_t eventNumber;
        ULong64_t eventNumber2;


        Char_t isTrain;

	UShort_t lumiBlock;
	Float_t mcGenWeight = 1.0;
	UChar_t nPU;
	UInt_t eleID[3];
	Float_t energy_ECAL_ele[3], R9Ele[3], etaEle[3], phiEle[3], pAtVtxGsfEle[3];
	Float_t esEnergySCEle[3];
	Float_t invMass, invMass_ECAL_ele, invMass_5x5SC, invMass_rawSC, etaEle_eL, phiEle_eL;
	Short_t xSeedSC[3], ySeedSC[3];
	Float_t etaSCEle[3];
	Float_t phiSCEle[3];
	Float_t ele1E, ele2E;
        Int_t zSide[2];

        Float_t S4[3], etaWidth[3], phiWidth[3], hOverE[3], deltaEtaSuperClusterTrackAtVtx[3], deltaPhiSuperClusterTrackAtVtx[3];
        Float_t invMass_ECAL_ele_recal;

	Float_t invMassTkr;
	Float_t invMassTkrTnP;
	Float_t invMassEleTnP;


	Float_t etaSC_L, etaSC_S;
	Float_t phiSC_L, phiSC_S;

	Int_t seedRing[3];


	void constructTimeHistograms();
	void skimTree(TTree* _tree);
	UInt_t FindFill(UInt_t t);
	Bool_t InSameFill(Double_t t1, Double_t t2);
	void optiIseTimeBinning(TH1D *&tHist);

	statistic stat;
	void setHistograms();
	TH2D *setTimedHistogram(category &cat, const histogramTemplate &_template, const std::vector<Double_t> &t_bins);
	TH1D *setVarHistograms(category &cat, const histogramTemplate &_template);

	void Fill1D();
	void FillAll();

	typedef void(ECALELF_ntuple_processor::*fill_func)();
	fill_func fill_hists;

	void analyse();
        void analyseTree(TChain *tree);

	Bool_t ApplyCuts();
	Bool_t mPassBasicCuts;
	Bool_t IsHR9;
	Bool_t IsLR9;
	Bool_t IsEb;
	Bool_t IsEbHR9;
	Bool_t IsEbLR9;
	Bool_t IsEe;
	Bool_t IsEeHR9;
	Bool_t IsEeLR9;

	Bool_t IsOuterEe;
	Bool_t IsOuterEeHR9;
	Bool_t IsOuterEeLR9;
	Bool_t IsMidEe;
	Bool_t IsMidEeHR9;
	Bool_t IsMidEeLR9;
	Bool_t IsInnerEe;
	Bool_t IsInnerEeLR9;
	Bool_t IsInnerEeHR9;

	Bool_t eL_IsInEe;
	Bool_t eS_IsInEb;

	Float_t energyECALeleS;

	void saveHistograms();


	std::vector<category> categories = {
		{"EB", "EB", &IsEb},
		{"EE", "EE", &IsEe},
	};

	std::vector<histogramTemplate> plot_variables = {
//		{&invMass_ECAL_ele, 40.0, 140.0, 1000, 1, 1, "m_{ee}_default", "[GeV]", "0.01 GeV"},
                {&invMass_ECAL_ele_recal, 40.0, 140.0, 1000, 1, 1, "m_{ee}", "[GeV]", "0.01 GeV"},
	};

	/*other event variables*/
	UInt_t ID;
	UInt_t era;
	UInt_t CatKey;
	UInt_t EBCatKey;
	UInt_t R9CatKey;
	Float_t Et[2];
	Float_t m12;
	Float_t R9L, R9S;
	Int_t iL, iS;
	Float_t absEtaL, absEtaS;
	Float_t EtL, EtS;
	Float_t EoPL, EoPS;
	Bool_t ebL, ebS, hR9L, hR9S;

        TChain *ntupleChain = nullptr;
        TChain *eleIDchain = nullptr;
};

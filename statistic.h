struct errors {

    errors(Double_t y, Double_t e1, Double_t e2) : x(y), errLow(e1), errHigh(e2) {
    }

    errors() {
    }
    Double_t x;
    Double_t errLow;
    Double_t errHigh;
};

class statistic {
public:
    /*Function to get N percent minimum interval*/
    TH1D *getNpercentMinHist(TH1 *hist, Float_t N);

    /*Umbrella functions for calling TH1 statistic functions*/
    Double_t getHistStatistic(TH1 *hist, std::string statistic, Float_t N);
    errors* getHistStatistic(TH1 *hist, std::string statistic, Float_t N, Int_t Nresamples, Double_t confidenceLevel);

    /*Functions to get from TH2* a TH1D* or TGraphAsymErrors* of y-statistic vs  x*/
    TH1D* D2toXstatistic(TH2 *h2, std::string statistic, Float_t N);
    TGraphAsymmErrors* D2toXstatistic(TH2 *h2, std::string statistic, Float_t N, Int_t Nresamples, Double_t confidenceLevel);

    /*Function to get from TH3 a TH2D* of */
    //    TH2D* D3toD2(TH3* threeD, std::string statistic, Int_t N = 100);

    /*Function for percentile of a vector of doubles*/
    Double_t getPercentile(vector<Double_t> list, Double_t P);

    /*Function to make basic checks on the histogram*/
    Bool_t checkHist(TH1 * hist);

    /*Functions to calculate TH1 statistic*/
    Double_t getHistStats(TH1 * hist);
    Double_t getHistMedian(TH1 * hist);
    Double_t getHistMean(TH1 * hist);
    Double_t getHistRMS(TH1 * hist);
    Double_t getHistStDev(TH1 * hist);

    /*Tools*/
    TGraphAsymmErrors* buildGraph(std::vector<errors*> xCollection, std::vector<errors*>yCollection);
    void copyErrorsToDouble(std::vector<errors*> errorList, std::vector<Double_t> & vector, Int_t pos = 0);
    //    void print_array(const Double_t* array);
    void clear_errors_vector(std::vector<errors*> & vec);
    std::vector<Double_t> array2vec(const Double_t *array);
    std::vector<Double_t> histBins(TH1* hist);
    std::string to_string_with_precision(Double_t a_value, const int n) ;
    std::string removeNonAlpha(std::string word);

private:
    typedef Double_t(statistic::*H_statistic)(TH1*);

    std::map<std::string, H_statistic > function_map = {
        {
            "median", &statistic::getHistMedian
        },
        {
            "stats", &statistic::getHistStats
        },
        {
            "mean", &statistic::getHistMean
        },
        {
            "rms", &statistic::getHistRMS
        },
        {
            "stdev", &statistic::getHistStDev
        }
    };
};
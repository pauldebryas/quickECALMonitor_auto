import os
import argparse
import json

parser = argparse.ArgumentParser(description=r'"submitNtuples.py" aims at proceeding the data, taking as an input a json file in which are stored the path of the ntuple list that you want to proceed.')

parser.add_argument("-f", nargs='+', help="Path to the json files where the data library is.", default=["/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/102X_dataRun2_Sep2018Rereco_Run2018D_SingleIOVrun320500_normalized_v1/EGamma-Run2018A-ZSkim-17Sep2018-v2/campaign_summary.json"])

parser.add_argument("-d", type=str, help="The data you wants to read.", default="102X_dataRun2_Sep2018Rereco_Run2018D_SingleIOVrun320500_normalized_v1")

parser.add_argument("-t", help="Type of the data file (merged or unmerged).", choices=['mergedFiles','unmergedFiles'], default="mergedFiles")

args = parser.parse_args()

nb_of_json_file = len(args.f) 

tot_path = []
tot_main = []
tot_eleID = []
tot_extraCalib = []
tot_extraStudy = []


if(args.t=="mergedFiles"):
    for i in range(nb_of_json_file):
	    with open(args.f[i]) as json_file:
		    files_path = json.load(json_file)
		    path = files_path[args.d]['path']
		    main = files_path[args.d]['mergedFiles']['main']
		    eleID= files_path[args.d]['mergedFiles']['eleID']
		    extraCalib = files_path[args.d]['mergedFiles']['extraCalib']
		    extraStudy = files_path[args.d]['mergedFiles']['extraStudy']
            tot_path.append(path)
            tot_main.append(main)
            tot_eleID.append(eleID)
            tot_extraCalib.append(extraCalib)
            tot_extraStudy.append(extraStudy)            

elif(args.t=="unmergedFiles"):
    for i in range(nb_of_json_file):
	    main = []
	    eleID = []
	    extraCalib = []
	    extraStudy = []
	    with open(args.f[i]) as json_file:
		    files_path = json.load(json_file)
		    path = files_path[args.d]['path']
		    for i in files_path[args.d]['unmergedFiles']:
			main.append(i['main'])
			eleID.append(i['eleID'])
			extraCalib.append(i['extraCalib'])
			extraStudy.append(i['extraStudy'])
            tot_path.append(path)
            tot_main.append(main)
            tot_eleID.append(eleID)
            tot_extraCalib.append(extraCalib)
            tot_extraStudy.append(extraStudy)

print(main)
'''
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
string = r'root -b -q -l "ntupleDumper.cxx("' + dname + r'", "ul2018_eop_2d_median_IOVN",1, "inputs/Run2_UL_timebins_all.root","inputs/lumi_3years_csv.dat")" &> logs/ul2018_eop_2d_median_IOVN_all.log'
print(string)

#os.system('cmsenv')
#os.system('source /eos/project-c/cms-ecal-calibration/ecal-venv/bin/activate')
#print("Current working directory " + os.getcwd())
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname) # python equivalent to the line workdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"/  
#print("Current working directory " + os.getcwd())

#os.system('mkdir -vp logs')

#os.system('root -b -q -l ntupleDumper.cxx() &> logs/ul2018_eop_2d_median_IOVN_all.log &')
#os.system('root -b -q -l ')

root -b -q -l "ntupleDumper.cxx(\"$workdir\", \"ul2018_eop_2d_median_IOVN\",1, \"inputs/Run2_UL_timebins_all.root\",\"inputs/lumi_3years_csv.dat\")" &> logs/ul2018_eop_2d_median_IOVN_all.log &
root -b -q -l "ntupleDumper.cxx(\"$workdir\", \"ul2018_eop_2d_median_IOVNm1\",1,\"inputs/Run2_UL_timebins_all.root\",\"inputs/lumi_3years_csv.dat\")" &> logs/ul2018_eop_2d_median_IOVNm1_all.log &
os.system('')
os.system('')
'''

cmsenv

workdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"/  # pass the directory of this script to workdir

mkdir -vp logs

root -b -q -l "ntupleDumper.cxx(\"$workdir\", \"ul2018_eop_2d_median_IOVN\",1, \"inputs/Run2_UL_timebins_all.root\",\"inputs/lumi_3years_csv.dat\")" &> logs/ul2018_eop_2d_median_IOVN_all.log &
root -b -q -l "ntupleDumper.cxx(\"$workdir\", \"ul2018_eop_2d_median_IOVNm1\",1,\"inputs/Run2_UL_timebins_all.root\",\"inputs/lumi_3years_csv.dat\")" &> logs/ul2018_eop_2d_median_IOVNm1_all.log &


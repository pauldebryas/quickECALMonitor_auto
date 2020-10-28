# quickECALMonitor

This is a package of codes for producing simple ECAL monitoring plots like mee, R9 etc. distributions.

# Authors

Muhammad Aamir Shazad, Rajdeep Mohan Chatterjee, Shilpi Jain, Jin Wang, and maybe some others

# Usage 

1. Initialize the CMSSW enviroment:

```
cmsrel CMSSW_9_4_9
cd CMSSW_9_4_9/src
cmsenv

```

2. Clone example codes:

```
git clone https://gitlab.cern.ch/wangjin/quickECALMonitor.git

```

3. Retrieve the ntuple list for the samples that you want to check. In "findNtuples.sh" you should give the path of the samples.

```
. findNtuples.sh

```

4. Dump the necessary variables you want to check, e.g. mee. You can modify "ECALELF_ntuple_processor" class to change selections and variables.

```
. dumpVariables.sh

```

5. Compare performance and extract resolution, change the file names and path in "plotAll.cxx", plot parameters can be modified with "histOption.txt".

```
root -b -q -l "plotAll.cxx(\"outputs/ul2018_eop_2d_median_IOVN/ul2018_eop_2d_median_IOVN_all.root\",\"outputs/ul2018_eop_2d_median_IOVNm1/ul2018_eop_2d_median_IOVNm1_all.root\")"

```


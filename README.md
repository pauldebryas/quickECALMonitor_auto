# quickECALMonitor

This is a package of codes for producing simple ECAL monitoring plots like mee, R9 etc. distributions.

# Authors

Muhammad Aamir Shazad, Rajdeep Mohan Chatterjee, Shilpi Jain, Jin Wang, and maybe some others

# Usage 

1. Initialize the CMSSW enviroment

```
cmsrel CMSSW_9_4_9
cd CMSSW_9_4_9/src
cmsenv

```

2. copy example codes

```
git clone https://gitlab.cern.ch/wangjin/quickECALMonitor.git

```

3. get the ntuple list for the samples that you want to check. in findNtuples.sh you should give the path of the samples 

```
. findNtuples.sh

```

4. dump the necessary variables you want to check, e.g. mee. You can modify "ECALELF_ntuple_processor" class to change selections and variables

```
. dumpVariables.sh

```

5. compare performance and extract resolution, change the file names and path in "plotHists.C", plot parameters can be modified with "histOption.txt"

```
root -b -q -l plotAll.cxx 

```


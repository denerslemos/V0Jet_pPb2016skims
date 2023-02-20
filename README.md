# Jet V0 skim code using HTCondor

Code to produce jets from the CMS HiForest and V0 skims from Dener. 

## Intructions

Setup CMSSW (just for root versioning)
```
export SCRAM_ARCH=slc7_amd64_gcc10
cmsrel CMSSW_12_5_0
cd CMSSW_12_5_0/src
cmsenv
```
Inside of the src folder, download the code using
```
git clone git@github.com:denerslemos/V0Jet_pPb2016skims.git
cd V0Jet_pPb2016skims
```
Before compile the code you must check the [sub_skim.sh](https://github.com/denerslemos/V0Jet_pPb2016skims/blob/main/sub_skim.sh) lines 4 (CMSSW/src) and 6 (.../V0Jet_pPb2016skims) and replace by your own folders.

Once this steps are done you can compile the code with
```
g++ -O2 V0Jet_pPbSkim.C `root-config --libs` `root-config --cflags` -o V0Jet_pPbSkim
```
This will create the executable: ```V0Jet_pPbSkim``` 

After that you will submit jobs by using HTCondor_submit.py, just text me to see how to submit jobs :).

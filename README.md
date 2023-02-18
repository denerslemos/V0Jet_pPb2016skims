# XeXe skim code using HTCondor

Code to produce jets and tracks from the CMS HiForest. Setup here: https://github.com/denerslemos/XeXe_repo

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
git clone git@github.com:denerslemos/XeXe2017skims.git
cd XeXe2017skims
```
Before compile the code you must check the [sub_skim.sh](https://github.com/denerslemos/XeXe2017skims/blob/main/sub_skim.sh) lines 4 (CMSSW/src) and 6 (.../XeXe2017skims) and replace by your own folders.

Once this steps are done you can compile the code with
```
g++ -O2 XeXeSkim.C `root-config --libs` `root-config --cflags` -o XeXeSkim
```
This will create the executable: ```XeXeSkim``` 

After that you will need your VOMS certificate, do it using
```
voms-proxy-init -rfc -voms cms --out voms_proxy.txt --hours 200
```
that creates a certificate file valid for 200 hours: ```voms_proxy.txt```

Now you can submit the condor jobs using the python script, [```HTCondor_submit.py```](https://github.com/denerslemos/XeXe2017skims/blob/main/HTCondor_submit.py):

```
python HTCondor_submit.py -i input_text_file -o output_name_file -m X -s Y
```

- input_text_file: is the text file (use it without the .txt extension) with inputs and can be found in the folders [MC_SAMPLES](https://github.com/denerslemos/XeXe2017skims/tree/main/MC_SAMPLES) or [DATA_SAMPLES](https://github.com/denerslemos/XeXe2017skims/tree/main/DATA_SAMPLES) each .root input will be a job

- output_name_file: output file name (use it without the .root extension), it will automatically include a counter for each input. You can use paths to save on EOS.

- X: 0 for data and 1 for MC

- Y: name for the submission files, I have used HTcondor_sub_ + some information from the sample, pthat, MB, ...

It will automatically include a counter for each input

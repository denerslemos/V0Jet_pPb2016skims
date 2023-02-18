#!/bin/bash

echo "Setup CMSSW (ROOT version)"
cd /afs/cern.ch/user/d/ddesouza/CMSSW_12_5_0/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/UIC/V0Jet/V0Jet_pPb2016skims
mkdir -p cond
echo "Submit skim jobs at "
echo PWD: $PWD

./V0Jet_pPbSkim $1 $2 $3 $4

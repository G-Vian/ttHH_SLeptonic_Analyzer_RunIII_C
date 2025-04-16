#!/bin/bash

dasgoclient --query="file dataset=/TTHHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM " | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_23.txt
#dasgoclient --query="file dataset=/TTZZTo4b_5f_LO_TuneCP5_13TeV_madgraph_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> filelist_TTZZ17.txt
#dasgoclient --query="file dataset=/TTTo2L2Nu-noSC_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> filelist_TT_DL17.txt

#how to use :
#  chmod +x makeFileList.sh
# voms-proxy-init --voms cms
# ./makeFileList.sh



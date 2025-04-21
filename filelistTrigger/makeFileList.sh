#!/bin/bash

dasgoclient --query="file dataset=/TTHHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM " | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_23.txt
dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_23.txt
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_23.txt
dasgoclient --query="file dataset=/TTZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_23.txt
dasgoclient --query="file dataset=/TTZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZZ_23.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v3/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_23.txt
#how to use :
#  chmod +x makeFileList.sh
# voms-proxy-init --voms cms
# ./makeFileList.sh



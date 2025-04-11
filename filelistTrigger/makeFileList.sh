#!/bin/bash

dasgoclient --query="file dataset=/TTHHTo4b_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> filelist_TTHH17.txt
dasgoclient --query="file dataset=/TT4b_TuneCP5_13TeV_madgraph_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> filelist_TT4b17.txt
dasgoclient --query="file dataset=/TTZHTo4b_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9_ext1-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> filelist_TTZH17.txt
dasgoclient --query="file dataset=/TTZZTo4b_5f_LO_TuneCP5_13TeV_madgraph_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> filelist_TTZZ17.txt
dasgoclient --query="file dataset=/TTTo2L2Nu-noSC_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> filelist_TT_DL17.txt



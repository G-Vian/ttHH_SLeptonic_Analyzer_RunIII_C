#!/bin/bash


dasgoclient --query="file dataset=/TTHHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_22.txt
dasgoclient --query="file dataset=/TTHHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_22EE.txt
dasgoclient --query="file dataset=/TTZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZZ_22.txt
dasgoclient --query="file dataset=/TTZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZZ_22EE.txt
dasgoclient --query="file dataset=/TTZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_22.txt
dasgoclient --query="file dataset=/TTZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_22EE.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22_1.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22_2.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22EE_1.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22EE_2.txt

dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_22EE.txt
dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_22.txt

dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE_1.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE_2.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE_3.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-JMENano12p5_132X_mcRun3_2022_realistic_postEE_v4-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE_4.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-JMENano12p5_132X_mcRun3_2022_realistic_postEE_v4_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE_5.txt

dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22_1.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22_2.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22_3.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-JMENano12p5_132X_mcRun3_2022_realistic_v3_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22_4.txt



dasgoclient --query="file dataset=//ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22_1.txt
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22_2.txt
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22EE_1.txt
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22EE_2.txt


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



#!/bin/sh

#### set tagNTuple value from crab submit jobs and then ./create_fileList_NTuple.sh tagNTuple

pathCrab_Signal=$1
pathCrab_Background=$2

#tagNTuple='20190505'
#pathCrab_Signal='/cms/store/user/sbhowmik/GluGluHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack/PhaseIIMTDTDRAutumn18MiniAOD_'${tagNTuple}'/*/*'
#pathCrab_Background='/cms/store/user/sbhowmik/NeutrinoGun_E_10GeV/PhaseIIMTDTDRAutumn18MiniAOD_'${tagNTuple}'/*/*'

unset JAVA_HOME

hdfs dfs -ls $pathCrab_Signal/*root | awk '{print $8}' | grep root | sed 's/^\//\/hdfs\//g' > fileList_NTuple_Signal.list

hdfs dfs -ls $pathCrab_Background/*root | awk '{print $8}' | grep root | sed 's/^\//\/hdfs\//g' > fileList_NTuple_Background.list


#ls -ltr $pathCrab_Signal/*root | awk '{print $9}' | grep root  > fileList_NTuple_Signal.list

#ls -ltr $pathCrab_Background/*root | awk '{print $9}' | grep root  > fileList_NTuple_Background.list

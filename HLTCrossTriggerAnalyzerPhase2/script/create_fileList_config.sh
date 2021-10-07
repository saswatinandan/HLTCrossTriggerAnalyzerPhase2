#!/bin/sh

filePath=$1
#filePath='/home/sbhowmik/L1TauTrigger/L1TauAnalyzerPhase2/CMSSW_10_6_1_patch2/src/L1TauAnalyzerPhase2/L1TauAnalyzerPhase2'

unset JAVA_HOME

echo $filePath
ls -ltr $filePath/test*py | awk '{print $9}' | grep py  > fileList_config.list

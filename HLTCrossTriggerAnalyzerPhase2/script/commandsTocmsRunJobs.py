import os, subprocess, sys


# ----------- *** Start Modification *** -------------------------------------

tagRootTree = 'VBFHToTauTau_DeepTau_20210611'
#tagRootTree = 'MinBias_DeepTau_20210611'

pathCrab_Signal = '/cms/store/user/sbhowmik/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8/VBFHToTauTau_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1117_'+tagNTuple+'/*/*'
#pathCrab_Background = '/cms/store/user/sbhowmik/MinBias_TuneCP5_14TeV-pythia8/MinBias_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1117_'+tagNTuple+'/*/*'

pathRootTree = '/home/sbhowmik/RootTree/HLTTau/Phase2/'

workingDir = os.getcwd()

#sampleType=["Signal", "Background"]
sampleType=["Signal"]
#sampleType=["Background"] 
#algoType=["HLTTau"]
algoType=["L1andHLTTau"]
#fileType=["rootTree", "bdt", "hist"] 
fileType=["rootTree"]



# ------------ *** End Modification *** --------------------------------------



# ------------ Define command to execute -------------------------------------
def run_cmd(command):
  print "executing command = '%s'" % command
  p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = p.communicate()
  return stdout




# -----------Create file list of crab output root files ------------

run_cmd('rm %s' % "*list")
scriptFile = os.path.join(workingDir, "script", "create_fileList_NTuple.sh")
run_cmd('bash %s %s %s' % (scriptFile, pathCrab_Signal, pathCrab_Background))

# -----------Create cfg file to run analyzer ------------

run_cmd('rm %s' % "test_*.py")
for i in range (0, len(sampleType)):
  for j in range (0, len(algoType)):
    scriptFile = os.path.join(workingDir, "script", "create_test_TauAnalyzer.sh")
    run_cmd('bash %s %s %s %s' % (scriptFile, tagRootTree, sampleType[i], algoType[j]))

# -----------cmsRun analyzer files -----------.   

run_cmd('rm %s' % "*.root")
for i in range (0, len(sampleType)):
  for j in range (0, len(algoType)):
    scriptFile = os.path.join(workingDir, "submit_jobs_cmsRun_"+algoType[j]+"_"+sampleType[i]+".sh")
    run_cmd('bash %s' % scriptFile)


# -----------manage output files -----------.
'''
for i in range (0, len(sampleType)):
  for j in range (0, len(algoType)):
    for k in range (0, len(fileType)):
      rootFiles=os.path.join(fileType[k]+"_test_"+algoType[j]+"Analyzer_"+sampleType[i]+"_"+tagRootTree+"_part_*.root")
      haddFile=os.path.join(fileType[k]+"_test_"+algoType[j]+"Analyzer_"+sampleType[i]+"_"+tagRootTree+".root")
      run_cmd('rm %s' % haddFile)
      run_cmd('hadd %s %s' % (haddFile, rootFiles))
      run_cmd('mv %s %s' % (haddFile, pathRootTree))

run_cmd('rm %s' % "*.root")
testDir = os.path.join(workingDir, "test")
run_cmd('mv %s %s' % ("*.list", testDir))
run_cmd('mv %s %s' % ("test_*.py", testDir))
run_cmd('mv %s %s' % ("submit_jobs_cmsRun*.sh", testDir))
run_cmd('rm %s' % "*.log")
'''




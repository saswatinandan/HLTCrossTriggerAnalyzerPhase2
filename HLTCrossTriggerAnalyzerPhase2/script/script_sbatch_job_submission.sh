#!/bin/bash
#SBATCH --job-name=HLT_Tau_Analyzer   # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=test_L1HPSPFTau_Background_part_1_%j.log   # Standard output and error log
#SBATCH --error=test_L1HPSPFTau_Background_part_1_%j.log   # Standard output and error log
#jobdir=$1
#outputdir=$2
inputFile=$1

pwd; hostname; date

export OMP_NUM_THREADS=8

echo "Running HLT Tau Analyzer script on a single CPU with 4 cores "

# unset JAVA_HOME, because hadoop commands might not work 
# this is especially true if one has sourced necessary files for the GRID proxy
echo 'Unsetting JAVA_HOME=$JAVA_HOME'
unset JAVA_HOME


#set cms enviroment
source /cvmfs/cms.cern.ch/cmsset_default.sh
#cd ${jobdir}
cd /home/sbhowmik/HLTTau/HLTTauAnalyzerPhase2/new/CMSSW_11_1_7/src/HLTTauAnalyzerPhase2/HLTTauAnalyzerPhase2/

eval $(scramv1 runtime -sh) # same as cmsenv

#cd ${outputdir}
#/usr/bin/time --verbose cmsRun /home/sbhowmik/L1TauTrigger/L1TauAnalyzerPhase2/CMSSW_10_6_1_patch2/src/L1TauAnalyzerPhase2/L1TauAnalyzerPhase2/test/test_L1HPSPFTau_Background_part_1.py 
/usr/bin/time --verbose cmsRun ${inputFile}


date

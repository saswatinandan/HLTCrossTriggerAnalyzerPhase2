# HLTCrossTriggerAnalyzerPhase2


```
cmsrel CMSSW_11_1_8
cd CMSSW_11_1_8/src
cmsenv
git cms-init
git cms-merge-topic -u  fwyzard:HLT_Phase2_menu

git clone https://github.com/sandeepbhowmik1/HLTProducerPhase2 $CMSSW_BASE/src/HLTrigger/Phase2HLTPFTaus
git clone https://github.com/sandeepbhowmik1/HLTDataformatsPhase2 $CMSSW_BASE/src/DataFormats/Phase2HLTPFTaus

scram b -j 10

```

# Download the HLTCrossTriggerAnalyzerPhase2/HLTCrossTriggerAnalyzerPhase2 code

Clone this repository into your `$CMSSW_BASE/src/` folder.

```
git clone https://github.com/sandeepbhowmik1/HLTCrossTriggerAnalyzerPhase2

```


import FWCore.ParameterSet.Config as cms
process = cms.Process('Analyze')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "file:/hdfs/local/snandan/hlttautimestudy/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/produce_particleflowwithtimecut_QCD_wgenmatch.root"
#            "file:/home/snandan/hlt_trigger/cmssw11/CMSSW_11_1_8/src/HLTrigger/Phase2HLTPFTaus/test/test.root"
#        "file:/hdfs/local/snandan/hlttautimestudy/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/produce_particleflowwithtimecut_VBF_withgenmatch.root"
    )
)
process.analysisSequence = cms.Sequence()

process.load("HLTCrossTriggerAnalyzerPhase2.HLTCrossTriggerAnalyzerPhase2.HLTTauIsolationStudy_cff")
process.analysisSequence += process.AnalyzerSeq
process.p = cms.Path(
    process.analysisSequence
    )
process.schedule = cms.Schedule(process.p)
suffix = 'VBF'
suffix += ''
process.TFileService=cms.Service('TFileService',
                                 #fileName=cms.string("HLTTauIsolationStudy_%s_isolation.root" %suffix))
                                 fileName=cms.string("test_QCD_015.root"))

process.options.numberOfThreads = cms.untracked.uint32(10)

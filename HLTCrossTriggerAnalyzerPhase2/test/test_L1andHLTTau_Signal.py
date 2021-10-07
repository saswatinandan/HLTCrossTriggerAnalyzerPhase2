import FWCore.ParameterSet.Config as cms
process = cms.Process('Analyze')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
	)
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
            'file:/home/sbhowmik/HLTTau/HLTTauProducerPhase2/CrossTrigger/CMSSW_11_1_8/src/HLTrigger/Phase2HLTPFTaus/test/NTuple_produce_HLT_CrossTrigger_10events.root'
	)
)
process.analysisSequence = cms.Sequence()

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.tauGenJets.GenParticles = cms.InputTag('genParticles')
process.analysisSequence += process.tauGenJets
process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.analysisSequence += process.tauGenJetsSelectorAllHadrons

process.load("HLTCrossTriggerAnalyzerPhase2.HLTCrossTriggerAnalyzerPhase2.L1andHLTTauAnalyzer_cff")
process.L1andHLTTauAnalyzer.histRootFileName = cms.string("hist_test_L1andHLTTauAnalyzer_Signal_VBFHToTauTau_DeepTau_20211007.root")
process.L1andHLTTauAnalyzer.bdtRootFileName = cms.string("bdt_test_L1andHLTTauAnalyzer_Signal_VBFHToTauTau_DeepTau_20211007.root")
process.analysisSequence += process.AnalyzerSeq
process.p = cms.Path(
	process.analysisSequence
)
process.schedule = cms.Schedule(process.p)
process.TFileService=cms.Service('TFileService',fileName=cms.string("rootTree_test_L1andHLTTauAnalyzer_Signal_VBFHToTauTau_DeepTau_20211007.root"))

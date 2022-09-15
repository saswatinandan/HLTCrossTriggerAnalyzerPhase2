import FWCore.ParameterSet.Config as cms

HLTTauIsolationStudy = cms.EDAnalyzer("HLTTauIsolationStudy",
                      debug    = cms.untracked.bool(False),
                      hltTauToken = cms.InputTag("hltUpdatedPatHpsPFTaus8HitsMaxDeltaZWithOfflineVerticesWithTimeCut015"),
                      hltTauSumChargedIsoToken = cms.string("chargedIsoPtSumdR03"),
                      genTauToken = cms.InputTag("tauGenJetsSelectorAllHadrons"),
)


AnalyzerSeq = cms.Sequence(
    HLTTauIsolationStudy
)

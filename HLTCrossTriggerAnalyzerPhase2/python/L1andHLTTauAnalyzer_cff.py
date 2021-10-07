import FWCore.ParameterSet.Config as cms

genMatchedTaus = cms.EDFilter("genMatchTauFilter",
        taus = cms.InputTag("slimmedTaus")
    )

goodTaus = cms.EDFilter("PATTauRefSelector",
        #src = cms.InputTag("slimmedTaus"),
        src = cms.InputTag("genMatchedTaus"),
        cut = cms.string(
        'pt > 20 && abs(eta) < 2.4 '
        '&& abs(charge) > 0 && abs(charge) < 2 '
        '&& tauID("decayModeFinding") > 0.5 '
        #'&& tauID("chargedIsoPtSum") < 2.5'
        '&& tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5'
        #'&& tauID("byMediumIsolationMVArun2v1DBoldDMwLT") > 0.5 '
        #'&& tauID("againstMuonTight3") > 0.5 '
        #'&& tauID("againstElectronVLooseMVA6") > 0.5 '
        ),
        filter = cms.bool(False)
)

genVertexProducer = cms.EDProducer("GenVertexProducer",
  #src = cms.InputTag('prunedGenParticles'),
  src = cms.InputTag('genParticles'),
  pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
) 

L1andHLTTauAnalyzer = cms.EDAnalyzer("L1andHLTTauAnalyzer",
                                     debug              = cms.untracked.bool(False),
                                isGenTau           = cms.untracked.bool(True),
                                isRecoTau          = cms.untracked.bool(False),
                                isHLTTau          = cms.untracked.bool(True),
                                isL1HPSTau          = cms.untracked.bool(True),
                                isL1NNTau          = cms.untracked.bool(True),
                                min_pt             = cms.untracked.double(0),
                                max_eta            = cms.untracked.double(3.0),
                                src_evtWeight      = cms.InputTag("stitchingWeight"),
                                genTagToken        = cms.InputTag("generator"),
                                genTauToken        = cms.InputTag("tauGenJetsSelectorAllHadrons"),
                                recoVertexToken    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                recoTauToken       = cms.InputTag("slimmedTaus"),
                                recoGMTauToken     = cms.InputTag("goodTaus"),
                                hltTauToken        = cms.InputTag("hltUpdatedPatHpsPFTaus8HitsMaxDeltaZWithOfflineVertices"),
                                hltTauSumChargedIsoToken = cms.string("byDeepTau2017v2VSjetraw"),
                                #hltTauSumChargedIsoToken = cms.string("chargedIsoPtSumHGCalFix"), 
                                l1hpsTauToken      = cms.InputTag("HPSPFTauProducerPF"),
                                #l1hpsTauToken      = cms.InputTag("L1HPSPFTauProducerWithStripsWithoutPreselectionPF"),
                                #l1nnTauToken       = cms.InputTag("L1NNTauProducer","L1PFTausNN"),
                                l1nnTauToken       = cms.InputTag("L1NNTauProducerPuppi","L1PFTausNN"),
                                treeName           = cms.string("L1andHLTTauAnalyzer"),
                                fillBDT            = cms.untracked.bool(True),
                                bdtRootFileName    = cms.string("bdt_test_L1andHLTTauAnalyzer.root"),
                                treeBDTName        = cms.string("L1andHLTTauAnalyzer"),
                                createHistRoorFile = cms.untracked.bool(False),
                                histRootFileName   = cms.string("hist_test_L1andHLTTauAnalyzer.root"),
                                applyBDT           = cms.untracked.bool(True),
                                bdtInputFileName   = cms.string("L1andHLTTauAnalyzerPhase2/L1andHLTTauAnalyzerPhase2/data/L1HPSPFTau_XGB_testVars_default_6Var_20200304_4.xml"),  # _4
)


AnalyzerSeq = cms.Sequence(
    #genVertexProducer +
    #genMatchedTaus +
    #goodTaus       +
    L1andHLTTauAnalyzer
)

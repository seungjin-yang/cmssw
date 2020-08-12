import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy


gemOfflineDQMTightGlbMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string(
        '(pt > 20)'
        '&& isGlobalMuon'
        '&& globalTrack.isNonnull'
        '&& passed(\'CutBasedIdTight\')'
    ),
    filter = cms.bool(False)
)


gemOfflineDQMStaMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string(
        '(pt > 20)'
        '&& isStandAloneMuon'
        '&& outerTrack.isNonnull'
    ),
    filter = cms.bool(False)
)


gemOfflineDQMCosmicMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muonsFromCosmics'),
    cut = cms.string(
        'outerTrack.isNonnull'
        '&& (pt > 20)'
        '&& (outerTrack.hitPattern.cscStationsWithValidHits > 0)'
    ),
    filter = cms.bool(False)
)

gemOfflineDQMCosmicMuons1Leg = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muonsFromCosmics1Leg'),
    cut = cms.string(
        'outerTrack.isNonnull'
#        '&& (outerTrack.hitPattern.cscStationsWithValidHits > 0)'
    ),
    filter = cms.bool(False)
)

# propagatorName = cms.untracked.string("SteppingHelixPropagatorAny"),
gemEfficiencyAnalyzerTightGlb = DQMEDAnalyzer('GEMEfficiencyAnalyzer',
    MuonServiceProxy,
    muonTag = cms.InputTag('gemOfflineDQMTightGlbMuons'),
    recHitTag = cms.InputTag('gemRecHits'),
    muonType = cms.untracked.string("Global"),
    residualXCut = cms.double(5.0),
    ptBinning = cms.untracked.vdouble(20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 200.),
    etaNbins = cms.untracked.int32(7),
    etaLow = cms.untracked.double(1.5),
    etaUp = cms.untracked.double(2.2),
    folder = cms.untracked.string('GEM/GEMEfficiency/TightGlobalMuon'),
    logCategory = cms.untracked.string('GEMEfficiencyAnalyzerTightGlb'),
)

gemEfficiencyAnalyzerSta = gemEfficiencyAnalyzerTightGlb.clone()
gemEfficiencyAnalyzerSta.muonTag = cms.InputTag("gemOfflineDQMStaMuons")
gemEfficiencyAnalyzerSta.muonType = cms.untracked.string("Standalone")
gemEfficiencyAnalyzerSta.folder = cms.untracked.string('GEM/GEMEfficiency/StandaloneMuon')
gemEfficiencyAnalyzerSta.logCategory = cms.untracked.string('GEMEfficiencyAnalyzerSta')

gemEfficiencyAnalyzerCosmic = gemEfficiencyAnalyzerSta.clone()
gemEfficiencyAnalyzerCosmic.muonTag = cms.InputTag("gemOfflineDQMCosmicMuons")
gemEfficiencyAnalyzerCosmic.muonType = cms.untracked.string("Cosmic")
# gemEfficiencyAnalyzerCosmic.ptBinning = cms.untracked.vdouble(20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 200.)
gemEfficiencyAnalyzerCosmic.folder = cms.untracked.string('GEM/GEMEfficiency/CosmicMuon')
gemEfficiencyAnalyzerCosmic.logCategory = cms.untracked.string('GEMEfficiencyAnalyzerCosmic')
gemEfficiencyAnalyzerCosmic.residualXCut = cms.double(5.0)
gemEfficiencyAnalyzerCosmic.etaLow = cms.untracked.double(0.0)
gemEfficiencyAnalyzerCosmic.etaUp = cms.untracked.double(3.0)
gemEfficiencyAnalyzerCosmic.etaNbins = cms.untracked.int32(30)


from Configuration.Eras.Modifier_phase2_GEM_cff import phase2_GEM
for _each in [gemEfficiencyAnalyzerTightGlb, gemEfficiencyAnalyzerSta]:
    phase2_GEM.toModify(_each, etaNbins=cms.untracked.int32(15), etaHigh=cms.untracked.double(3.0))

gemEfficiencyAnalyzerTightGlbSeq = cms.Sequence(
    cms.ignore(gemOfflineDQMTightGlbMuons) *
    gemEfficiencyAnalyzerTightGlb)

gemEfficiencyAnalyzerStaSeq = cms.Sequence(
    cms.ignore(gemOfflineDQMStaMuons) *
    gemEfficiencyAnalyzerSta)

gemEfficiencyAnalyzerCosmicSeq = cms.Sequence(
    cms.ignore(gemOfflineDQMCosmicMuons) *
    gemEfficiencyAnalyzerCosmic)

gemEfficiencyAnalyzerSeq = cms.Sequence(
    gemEfficiencyAnalyzerTightGlbSeq *
    gemEfficiencyAnalyzerStaSeq *
    gemEfficiencyAnalyzerCosmicSeq)

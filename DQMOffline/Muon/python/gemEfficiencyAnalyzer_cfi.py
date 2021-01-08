import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy
from DQMOffline.Muon.gemEfficiencyAnalyzerDefault_cfi import *


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


gemEfficiencyAnalyzerTight = gemEfficiencyAnalyzerDefault.clone()
gemEfficiencyAnalyzerTight.ServiceParameters = MuonServiceProxy.ServiceParameters.clone()
gemEfficiencyAnalyzerTight.muonTag = cms.InputTag('gemOfflineDQMTightGlbMuons')
gemEfficiencyAnalyzerTight.folder = cms.untracked.string('GEM/GEMEfficiency/TightGlobalMuon')
gemEfficiencyAnalyzerTight.logCategory = cms.untracked.string('GEMEfficiencyAnalyzerTight')

gemEfficiencyAnalyzerSTA = gemEfficiencyAnalyzerDefault.clone()
gemEfficiencyAnalyzerSTA.ServiceParameters = MuonServiceProxy.ServiceParameters.clone()
gemEfficiencyAnalyzerSTA.muonTag = cms.InputTag("gemOfflineDQMStaMuons")
gemEfficiencyAnalyzerSTA.useGlobalMuon = cms.untracked.bool(False)
gemEfficiencyAnalyzerSTA.folder = cms.untracked.string('GEM/GEMEfficiency/StandaloneMuon')
gemEfficiencyAnalyzerSTA.logCategory = cms.untracked.string('GEMEfficiencyAnalyzerSTA')

from Configuration.Eras.Modifier_phase2_GEM_cff import phase2_GEM
phase2_GEM.toModify(gemEfficiencyAnalyzerTight, etaNbins=cms.untracked.int32(15), etaHigh=cms.untracked.double(3.0))
phase2_GEM.toModify(gemEfficiencyAnalyzerSTA, etaNbins=cms.untracked.int32(15), etaHigh=cms.untracked.double(3.0))

gemEfficiencyAnalyzerTightSeq = cms.Sequence(
    cms.ignore(gemOfflineDQMTightGlbMuons) *
    gemEfficiencyAnalyzerTight)

gemEfficiencyAnalyzerSTASeq = cms.Sequence(
    cms.ignore(gemOfflineDQMStaMuons) *
    gemEfficiencyAnalyzerSTA)

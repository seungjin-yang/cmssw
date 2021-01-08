import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy
from DQMOffline.Muon.gemEfficiencyAnalyzerDefault_cfi import *

gemEfficiencyAnalyzerCosmics = gemEfficiencyAnalyzerDefault.clone()
gemEfficiencyAnalyzerCosmics.muonTag = cms.InputTag("muons")
gemEfficiencyAnalyzerCosmics.ServiceParameters = MuonServiceProxy.ServiceParameters.clone()
gemEfficiencyAnalyzerCosmics.isCosmics = cms.untracked.bool(True)
gemEfficiencyAnalyzerCosmics.useGlobalMuon = cms.untracked.bool(False)
gemEfficiencyAnalyzerCosmics.folder = cms.untracked.string('GEM/GEMEfficiency/StandaloneMuon')
gemEfficiencyAnalyzerCosmics.logCategory = cms.untracked.string('GEMEfficiencyAnalyzerCosmics')

from Configuration.Eras.Modifier_phase2_GEM_cff import phase2_GEM
phase2_GEM.toModify(gemEfficiencyAnalyzerCosmics, etaNbins=cms.untracked.int32(15), etaHigh=cms.untracked.double(3.0))

gemEfficiencyAnalyzerCosmicsSeq = cms.Sequence(gemEfficiencyAnalyzerCosmics)

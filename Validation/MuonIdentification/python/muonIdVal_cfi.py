import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
muonIdVal = DQMEDAnalyzer('MuonIdVal',
    inputMuonCollection           = cms.InputTag("muons"),
    inputDTRecSegment4DCollection = cms.InputTag("dt4DSegments"),
    inputCSCSegmentCollection     = cms.InputTag("cscSegments"),
    inputGEMSegmentCollection     = cms.InputTag("gemSegments"),
    inputMuonTimeExtraValueMap    = cms.InputTag("muons"),
    inputMuonCosmicCompatibilityValueMap = cms.InputTag("muons","cosmicsVeto"),
    inputMuonShowerInformationValueMap = cms.InputTag("muons","muonShowerInformation"), 
    useTrackerMuons               = cms.untracked.bool(True),
    useGlobalMuons                = cms.untracked.bool(True),
    useTrackerMuonsNotGlobalMuons = cms.untracked.bool(True),
    useGlobalMuonsNotTrackerMuons = cms.untracked.bool(True),
    useGEM                        = cms.untracked.bool(False),
    makeEnergyPlots               = cms.untracked.bool(True),
    makeTimePlots                 = cms.untracked.bool(True),
    make2DPlots                   = cms.untracked.bool(True),
    makeAllChamberPlots           = cms.untracked.bool(True),
    makeCosmicCompatibilityPlots  = cms.untracked.bool(True),
    makeShowerInformationPlots    = cms.untracked.bool(True),
    baseFolder                    = cms.untracked.string("Muons/MuonIdentificationV")
)

# fastsim has no cosmic muon veto in place
from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(muonIdVal, makeCosmicCompatibilityPlots = False)

from Configuration.Eras.Modifier_run3_GEM_cff import run3_GEM
run3_GEM.toModify(muonIdVal, useGEM=cms.untracked.bool(True))

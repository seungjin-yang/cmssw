import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester


gemEfficiencyHarvesterTightGlb = DQMEDHarvester('GEMEfficiencyHarvester',
    folder = cms.untracked.string('GEM/GEMEfficiency/TightGlobalMuon'),
    logCategory = cms.untracked.string('GEMEfficiencyHarvesterTightGlb')
)


gemEfficiencyHarvesterSta = DQMEDHarvester('GEMEfficiencyHarvester',
    folder = cms.untracked.string('GEM/GEMEfficiency/StandaloneMuon'),
    logCategory = cms.untracked.string('GEMEfficiencyHarvesterSta')
)


gemEfficiencyHarvesterCosmic = DQMEDHarvester('GEMEfficiencyHarvester',
    folder = cms.untracked.string('GEM/GEMEfficiency/CosmicMuon'),
    logCategory = cms.untracked.string('GEMEfficiencyHarvesterCosmic')
)


gemEfficiencyHarvesterSeq = cms.Sequence(gemEfficiencyHarvesterTightGlb *
                                         gemEfficiencyHarvesterSta *
                                         gemEfficiencyHarvesterCosmic)

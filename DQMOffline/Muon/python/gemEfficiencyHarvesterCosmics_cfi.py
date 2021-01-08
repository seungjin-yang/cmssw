import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from DQMOffline.Muon.gemEfficiencyHarvesterDefault_cfi import gemEfficiencyHarvesterDefault
from DQMOffline.Muon.gemEfficiencyAnalyzerCosmics_cfi import gemEfficiencyAnalyzerCosmics

gemEfficiencyHarvesterCosmics = gemEfficiencyHarvesterDefault.clone()
gemEfficiencyHarvesterCosmics.folder = cms.untracked.string(gemEfficiencyAnalyzerCosmics.folder.value())
gemEfficiencyHarvesterCosmics.logCategory = cms.untracked.string('GEMEfficiencyHarvesterTight')

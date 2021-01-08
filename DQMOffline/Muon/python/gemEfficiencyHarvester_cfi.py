import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from DQMOffline.Muon.gemEfficiencyHarvesterDefault_cfi import gemEfficiencyHarvesterDefault
from DQMOffline.Muon.gemEfficiencyAnalyzer_cfi import gemEfficiencyAnalyzerTight
from DQMOffline.Muon.gemEfficiencyAnalyzer_cfi import gemEfficiencyAnalyzerSTA

gemEfficiencyHarvesterTight = gemEfficiencyHarvesterDefault.clone()
gemEfficiencyHarvesterTight.folder = cms.untracked.string(gemEfficiencyAnalyzerTight.folder.value())
gemEfficiencyHarvesterTight.logCategory = cms.untracked.string('GEMEfficiencyHarvesterTight')

gemEfficiencyHarvesterSTA = gemEfficiencyHarvesterDefault.clone()
gemEfficiencyHarvesterSTA.folder = cms.untracked.string(gemEfficiencyAnalyzerSTA.folder.value())
gemEfficiencyHarvesterSTA.logCategory = cms.untracked.string('GEMEfficiencyHarvesterSTA')

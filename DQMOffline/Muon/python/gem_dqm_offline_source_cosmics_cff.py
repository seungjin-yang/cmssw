import FWCore.ParameterSet.Config as cms

from DQMOffline.Muon.gemOfflineMonitorDefault_cfi import *
from DQMOffline.Muon.gemEfficiencyAnalyzerCosmics_cfi import *

gemOfflineMonitorCosmics = gemOfflineMonitorDefault.clone()

gemSourcesCosmics = cms.Sequence(
    gemOfflineMonitorCosmics *
    gemEfficiencyAnalyzerCosmicsSeq)

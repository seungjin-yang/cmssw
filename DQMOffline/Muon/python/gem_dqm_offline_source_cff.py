import FWCore.ParameterSet.Config as cms

from DQMOffline.Muon.gemOfflineMonitorDefault_cfi import *
from DQMOffline.Muon.gemEfficiencyAnalyzer_cfi import *

gemOfflineMonitor = gemOfflineMonitorDefault.clone()

gemSources = cms.Sequence(
    gemOfflineMonitor *
    gemEfficiencyAnalyzerTightSeq *
    gemEfficiencyAnalyzerSTASeq
)

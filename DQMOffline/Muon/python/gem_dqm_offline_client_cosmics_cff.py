import FWCore.ParameterSet.Config as cms

from DQMOffline.Muon.gemOfflineHarvesterDefault_cfi import *
from DQMOffline.Muon.gemEfficiencyHarvesterCosmics_cfi import *

gemOfflineHarvesterCosmics = gemOfflineHarvesterDefault.clone()

gemClientsCosmics = cms.Sequence(
    gemOfflineHarvesterCosmics *
    gemEfficiencyHarvesterCosmics)

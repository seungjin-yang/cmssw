import FWCore.ParameterSet.Config as cms

from DQMOffline.Muon.gemOfflineHarvesterDefault_cfi import *
from DQMOffline.Muon.gemEfficiencyHarvester_cfi import *

gemOfflineHarvester = gemOfflineHarvesterDefault.clone()

gemClients = cms.Sequence(
    gemOfflineHarvester *
    gemEfficiencyHarvesterTight *
    gemEfficiencyHarvesterSTA)

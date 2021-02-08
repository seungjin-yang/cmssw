import FWCore.ParameterSet.Config as cms
from EventFilter.GEMRawToDigi.muonGEMDigisDefault_cfi import muonGEMDigisDefault as _muonGEMDigisDefault
from L1Trigger.L1TGEM.simMuonGEMPadDigis_cfi import simMuonGEMPadDigis as _simMuonGEMPadDigis
from L1Trigger.L1TGEM.simMuonGEMPadDigiClusters_cfi import simMuonGEMPadDigiClusters as _simMuonGEMPadDigiClusters

muonGEMDigis = _muonGEMDigisDefault.clone()
muonGEMPadDigis = _simMuonGEMPadDigis.clone(InputCollection=cms.InputTag("muonGEMDigis"))
muonGEMPadDigiClusters = _simMuonGEMPadDigiClusters.clone(InputCollection=cms.InputTag("muonGEMPadDigis"))

from Configuration.Eras.Modifier_run2_GEM_2017_cff import run2_GEM_2017
from Configuration.Eras.Modifier_run3_GEM_cff import run3_GEM
from Configuration.Eras.Modifier_phase2_GEM_cff import phase2_GEM

run2_GEM_2017.toModify(muonGEMDigis, useDBEMap = True)
run3_GEM.toModify(muonGEMDigis, useDBEMap = True)
phase2_GEM.toModify(muonGEMDigis, useDBEMap = False)

muonGEMDigiTask = cms.Task(
    muonGEMDigis,
    muonGEMPadDigis,
    muonGEMPadDigiClusters)

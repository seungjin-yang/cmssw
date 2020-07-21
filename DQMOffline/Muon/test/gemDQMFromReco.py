# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step10 --conditions auto:phase1_2021_realistic -s DQM:@none --geometry DB:Extended --era Run3 --mc --eventcontent DQM --filein file:/hdfs/store/user/seyang/GEM/Validation-SW/GEMOfflineMonitor-CosmicMuon/CMSSW_11_2_0_pre2/11654.0_region+1/step3_firstRun-1/step3_000.root --fileout file:inDQM.root -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('DQM',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('DQMOffline.Muon.gemOfflineMonitor_cfi')
process.load('DQMOffline.Muon.gemEfficiencyAnalyzer_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)


from FWCore.ParameterSet import VarParsing
import glob
import os.path

options = VarParsing.VarParsing ('analysis')
options.register(
    "inPath",
    "/eos/cms/store/user/hyunyong/cosmics/",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "a path to a RECO file or a directory having such files")
options.parseArguments()

if os.path.isdir(options.inPath):
    pathname = os.path.join(options.inPath, 'step3*.root')
    entries = ["file:" + each for each in glob.glob(pathname) if 'in' not in os.path.basename(each)]
    fileNames = cms.untracked.vstring(*entries)
elif os.path.isfile(options.inPath):
    fileNames = cms.untracked.vstring("file:" + options.inPath)
else:
    raise IOError('No such file or directory: {}'.format(options.inPath))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = fileNames,
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step10 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


process.gemOfflineMonitor.doDigi = cms.untracked.bool(False)

# Additional output definition

# Other statements
#process.mix.playback = True
#process.mix.digitizers = cms.PSet()
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

# Path and EndPath definitions
process.gem_dqm_step = cms.EndPath(process.gemOfflineMonitor + process.gemEfficiencyAnalyzerSeq)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.gem_dqm_step,
    process.endjob_step,
    process.DQMoutput_step)

# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: runGEM --data --era Run3 --conditions auto:run3_data_promptlike --scenario cosmics --geometry DB:Extended --step DQM --eventcontent DQM --datatier DQMIO --filein file:/hdfs/store/user/seyang/store/data/Commissioning2020/Cosmics/RAW-RECO/CosmicTP-PromptReco-v1/000/337/973/00000/0377F3AC-6086-CB47-946C-641886BEE2FA.root -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('DQM',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineCosmics_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('DQMOffline.Muon.gemOfflineMonitor_cfi')
process.load('DQMOffline.Muon.gemEfficiencyAnalyzer_cfi')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

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
    "/eos/cms/store/data/Commissioning2020/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/337/973/00000/",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "")
options.parseArguments()

if os.path.isdir(options.inPath):
    pathname = os.path.join(options.inPath, '*.root')
    entries = ["file:" + each for each in glob.glob(pathname)]
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
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('runGEM nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('GEM_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_promptlike', '')
process.gemOfflineMonitor.doDigi = cms.untracked.bool(False)


# Path and EndPath definitions
process.gem_dqmoffline_step = cms.EndPath(process.gemOfflineMonitor + process.gemEfficiencyAnalyzerSeq)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.gem_dqmoffline_step, process.DQMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)



# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

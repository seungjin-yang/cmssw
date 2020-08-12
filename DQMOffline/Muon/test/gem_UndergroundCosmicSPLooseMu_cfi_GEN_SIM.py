# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: UndergroundCosmicSPLooseMu_cfi --conditions auto:phase1_2021_cosmics -n 10 --era Run3 --eventcontent RAWSIM --relval 9000,50 -s GEN,SIM --scenario cosmics --datatier GEN-SIM --beamspot Run3RoundOptics25ns13TeVLowSigmaZ --geometry DB:Extended --fileout file:step1.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('SIM',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('SimGeneral.MixingModule.mixCosmics_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimNOBEAM_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource")

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
    annotation = cms.untracked.string('UndergroundCosmicSPLooseMu_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_cosmics', '')


process.cosmicInGEM = cms.EDFilter("CosmicGenFilterHelix",
    charges = cms.vint32(1, -1),
    doMonitor = cms.untracked.bool(False),
    maxZ = cms.double(-550),
    minP = cms.double(0.0),
    minPt = cms.double(0.0),
    minZ = cms.double(-600),
    pdgIds = cms.vint32(-13, 13),
    propagator = cms.string('SteppingHelixPropagatorAlong'),
    radius = cms.double(260.0),
    src = cms.InputTag("generator","unsmeared")
)


process.generator = cms.EDProducer("CosMuoGenProducer",
    AcptAllMu = cms.bool(False),
    ClayWidth = cms.double(50000.0),
    ElossScaleFactor = cms.double(1.0),
    MTCCHalf = cms.bool(False),
    MaxEnu = cms.double(10000.0),
    MaxP = cms.double(3000.0),
    MaxPhi = cms.double(360.0),
    MaxT0 = cms.double(12.5),
    MaxTheta = cms.double(84.0),
    MinEnu = cms.double(10.0),
    MinP = cms.double(10.0),
    MinP_CMS = cms.double(-1.0),
    MinPhi = cms.double(0.0),
    MinT0 = cms.double(-12.5),
    MinTheta = cms.double(0.0),
    MultiMuon = cms.bool(False),
    MultiMuonFileFirstEvent = cms.int32(1),
    MultiMuonFileName = cms.string('CORSIKAmultiMuon.root'),
    MultiMuonNmin = cms.int32(2),
    NuProdAlt = cms.double(7500000.0),
    PlugVx = cms.double(0.0),
    PlugVz = cms.double(-14000.0),
    RadiusOfTarget = cms.double(8000.0),
    RhoAir = cms.double(0.001214),
    RhoClay = cms.double(2.3),
    RhoPlug = cms.double(2.5),
    RhoRock = cms.double(2.5),
    RhoWall = cms.double(2.5),
    TIFOnly_constant = cms.bool(False),
    TIFOnly_linear = cms.bool(False),
    TrackerOnly = cms.bool(False),
    Verbosity = cms.bool(False),
    ZCentrOfTarget = cms.double(0.0),
    ZDistOfTarget = cms.double(15000.0)
)


process.SteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.ProductionFilterSequence = cms.Sequence(process.generator+process.cosmicInGEM)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.ProductionFilterSequence)



# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

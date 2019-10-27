import os
import FWCore.ParameterSet.Config as cms

# Set parameters externally 
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'skipEvents',
    0,
    VarParsing.multiplicity.singleton,VarParsing.varType.int,
    'Path of local input files'
)

params.register(
    'reportFreq',
    1,
    VarParsing.multiplicity.singleton,VarParsing.varType.int,
    'Path of local input files'
)

params.register(
    'localDatasetPath',
    '',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Path of local input files'
)

params.register(
    'Dataset',
    '',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Get input files of a published dataset'
)

params.register(
    'DBSInstance',
    '',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'DBS instance, e.g. prod/phys03 for USER samples'
)

params.register(
    'JSONFile',
    '',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'JSON file'
)

params.register(
    'GlobalTag',
    '94X_mc2017_realistic_v11',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Global tag'
)

params.register(
    'addEventInfo',
    True,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to save the event coordinates in the tree'
)

# Define the process
process = cms.Process("AV")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = params.reportFreq

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if params.GlobalTag != '' : process.GlobalTag.globaltag = params.GlobalTag
else : process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True)
)

# How many events to process
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(params.maxEvents))

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("tree.root")
)

# Define the inputs
process.source = cms.Source("PoolSource",
    fileNames  = cms.untracked.vstring(params.inputFiles),
    skipEvents = cms.untracked.uint32(params.skipEvents)
)

if params.localDatasetPath != '':
    for fname in os.listdir(params.localDatasetPath):
        path = params.localDatasetPath + "/" + fname
        if os.path.getsize(path) == 0 :
            print path + " is invalid ... skipping"
        else :
            path = "file:" + path
            process.source.fileNames += [path]

if params.Dataset != '':
    dbsinstance = " instance=" + params.DBSInstance if params.DBSInstance != '' else ""
    query = "das_client -query=\"file dataset=" + params.Dataset + dbsinstance + "\""
    fnames = os.popen(query).readlines()
    for fname in fnames:
        process.source.fileNames += [fname.rstrip()]

# Apply JSON if specified and if running on data 
import FWCore.PythonUtilities.LumiList as LumiList
if params.JSONFile != '' :
    process.source.lumisToProcess = LumiList.LumiList(filename = params.JSONFile).getVLuminosityBlockRange()

# The treemaker
process.tree = cms.EDAnalyzer('MuonHitsMapper',
    addEventInfo  = cms.bool(params.addEventInfo),
    muons         = cms.InputTag("muons"),
    gens          = cms.InputTag("genParticles"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
)

# Defining the path
process.p = cms.Path(process.tree)

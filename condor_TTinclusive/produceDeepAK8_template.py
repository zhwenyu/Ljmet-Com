import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register('isMC', '', VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Is MC') 
options.inputFiles = [CONDOR_FILELIST]
options.isMC = CONDOR_ISMC
options.maxEvents = -1
options.parseArguments()

# ---------------------------------------------------------
process = cms.Process("test")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),  
   wantSummary=cms.untracked.bool(False)
)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source('PoolSource',
    fileNames=cms.untracked.vstring(options.inputFiles)
)

# ---------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v13', '')
if options.isMC == False: process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6')
print 'Using global tag', process.GlobalTag.globaltag
# ---------------------------------------------------------
# set up TransientTrackBuilder
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)
# ---------------------------------------------------------
# recluster Puppi jets
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags'
]
JETCorrLevels = ['L2Relative', 'L3Absolute']

# ---------------------------------------------------------
#process.deepntuplizer = cms.EDAnalyzer('FatJetNNTestAnalyzer',
process.deepntuplizer = cms.EDProducer('FatJetNNTestProducer',
                                jets=cms.untracked.InputTag('slimmedJetsAK8'),
                                subjets=cms.untracked.InputTag(''),
                                jetR=cms.untracked.double(0.8),
                                useReclusteredJets=cms.untracked.bool(False),
                                decorrMode=cms.untracked.int32(0),
                                datapath=cms.untracked.string('NNKit/data/ak8'),
                                )
process.p = cms.Path(process.deepntuplizer)

process.output = cms.OutputModule(
                "PoolOutputModule",
                fileName = cms.untracked.string('CONDOR_MEDIATOR'),
                #fileName = cms.untracked.string('Tprime1700test_deepak8.root'),
                outputCommands = cms.untracked.vstring(['keep *','keep *_deepntuplizer_*_test']),
                )

process.output_step = cms.EndPath(process.output)
process.scedule = cms.Schedule(
    process.p,
    process.output_step)


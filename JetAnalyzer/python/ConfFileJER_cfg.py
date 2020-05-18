import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("AK4jets")

# QG likelihood
process.load("JetNtupleProducerTool.JetAnalyzer.QGLikelihood_cfi")
process.load("RecoJets.JetProducers.QGTagger_cfi")
process.QGTagger.srcJets = cms.InputTag("patJetsReapplyJEC")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")

# File service
process.load("CommonTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName=cms.string("JetNtuple_RunIISummer16_13TeV_MCJER.root")

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring("/store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/0C645070-A197-2648-B6E3-9AA2D7545A4F.root"),
                            skipEvents = cms.untracked.uint32(0)
                        )

process.jec = cms.ESSource("PoolDBESSource",
         DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
         timetype = cms.string('runnumber'),
         toGet = cms.VPSet(cms.PSet(record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string('JetCorrectorParametersCollection_Autumn18_V19_MC_AK4PFchs'),
                                    label  = cms.untracked.string('AK4PFchs')
                                    )
                 ),
         connect = cms.string('sqlite:Autumn18_V19_MC.db')
    )
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

### JECs =====================================================================================================
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet', 
            'L2Relative', 
            'L3Absolute'],
  payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
process.patJetsReapplyJEC = updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
    )

### JER
process.load('Configuration.StandardSequences.Services_cff')
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *

process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        toGet = cms.VPSet(
            # Resolution
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Autumn18_V7b_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),
            # Scale factors
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Autumn18_V7b_MC_SF_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
                ),
            ),
        connect = cms.string('sqlite:Autumn18_V7b_MC.db')
        )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


process.AK4jets = cms.EDAnalyzer("JetAnalyzer",
	## jet, PF and generator level collections ##
        jets = cms.InputTag("patJetsReapplyJEC"),
	pfCands = cms.InputTag("packedPFCandidates"),
	genJets = cms.InputTag("slimmedGenJets"),
	genEventInfo = cms.InputTag("generator"),
	## good primary vertices ##
	vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        payload = cms.string('AK4PFchs'),
	confGoodVtxNdof = cms.double(4),
	confGoodVtxZ = cms.double(24),
	confGoodVtxRho = cms.double(2),
	## pileup and rhos ##
	pileupInfo = cms.InputTag("slimmedAddPileupInfo"),
	pfRhoAll = cms.InputTag("fixedGridRhoFastjetAll"),
	pfRhoCentral = cms.InputTag("fixedGridRhoFastjetCentral"),
	pfRhoCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
	pfRhoCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
)

# Choose how many events to process (-1 = all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Report execution progress
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.p = cms.Path(process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC + process.QGTagger + process.AK4jets)

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("AK4jets")

# QG likelihood
process.load("JetNtupleProducerTool.JetAnalyzer.QGLikelihood_cfi")
process.load("RecoJets.JetProducers.QGTagger_cfi")
process.QGTagger.srcJets = cms.InputTag("slimmedJets")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")

# File service
process.load("CommonTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName=cms.string("JetNtuple_RunIISummer16_13TeV_MCNoJEC.root")

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring("/store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/0C645070-A197-2648-B6E3-9AA2D7545A4F.root"),
                            skipEvents = cms.untracked.uint32(0)
                        )

process.AK4jets = cms.EDAnalyzer("JetAnalyzer",
	## jet, PF and generator level collections ##
        jets = cms.InputTag("slimmedJets"),
	pfCands = cms.InputTag("packedPFCandidates"),
	genJets = cms.InputTag("slimmedGenJets"),
	genEventInfo = cms.InputTag("generator"),
	## good primary vertices ##
	vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
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

process.p = cms.Path(process.QGTagger + process.AK4jets)

import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList.json').getVLuminosityBlockRange()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#QG likelihood
process.load("Jetter.MiniAnalyzer.QGLikelihood_cfi")
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets=cms.InputTag("slimmedJets")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")

#Choose how many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000000) )

#Select the MiniAOD file to process
# With pileup:
# fileList = FileUtils.loadListFromFile('RunIISummer16MiniAODv2_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_filelist.txt')
# Without pileup:
# fileList = FileUtils.loadListFromFile('RunIISummer16MiniAODv2_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_filelist.txt')
#/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer16MiniAODv2-NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM

#'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/08278E4E-E4EF-E611-8BD7-FA163E3ABA64.root'
#'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DCD55D97-61F0-E611-9C1B-FA163EAD13C1.root'

process.source = cms.Source("PoolSource",
	# fileNames = cms.untracked.vstring(*fileList)
	# fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/08278E4E-E4EF-E611-8BD7-FA163E3ABA64.root')
	fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/024535A1-A5EF-E611-9C1C-FA163EE26C02.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/027011D6-A5EF-E611-8317-02163E01765E.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/0C375534-97EF-E611-A05E-02163E012DD7.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DEAC1B5D-79EF-E611-A80C-FA163E62DB22.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DCD55D97-61F0-E611-9C1B-FA163EAD13C1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DAE1461A-BFEF-E611-A7B0-FA163E306E8B.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DA47F5BD-84EF-E611-9717-02163E00C5B8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/D4C1FD3F-49F5-E611-9C57-ECF4BBE15B60.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A86E2919-BFEF-E611-A980-FA163E5DEB61.root')
)

process.demo = cms.EDAnalyzer('MiniAnalyzer',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
    packed = cms.InputTag("packedGenParticles"),
    genJets = cms.InputTag("slimmedGenJets"),
	genEventInfo = cms.InputTag("generator")

)


process.p = cms.Path(process.QGTagger + process.demo)

process.MessageLogger.cerr.FwkReport.reportEvery = 10

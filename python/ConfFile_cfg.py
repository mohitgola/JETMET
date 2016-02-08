import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
 #       '/store/user/khurana/ZprimeToZhToZlephbb_narrow_M-600_13TeV-madgraph/crab_ZprimeToZhToZlephbb_narrow_M-600_13TeV-madgraph_MC25ns_ReMiniAOD_NewMiniIso_20151120/151120_160912/0000/NCUGlobalTuples_1.root'
#'root://cms-xrd-global.cern.ch//store/user/khurana/ZprimeToZhToZlephbb_narrow_M-800_13TeV-madgraph/crab_ZprimeToZhToZlephbb_narrow_M-800_13TeV-madgraph_MC25ns_ReMiniAOD_NewMiniIso_20151120/151120_161007/0000/NCUGlobalTuples_1.root'
#'root:/store/user/khurana/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_NewMiniIso_20151120/151120_154908/0000/NCUGlobalTuples_1.root'
 #'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/DoubleElectron_FlatPt-1To300/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/0068838E-F508-E511-85C1-0026189438F9.root'
'root://xrootd.unl.edu//store/mc/Spring14miniaod/DYJetsToEEMuMu_M-120To200_13TeV-madgraph/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/80A19415-A40E-E411-8B04-001D09FDE18B.root' 
)
)

process.demo = cms.EDAnalyzer('MiniAnalyzer', 
jets = cms.InputTag("slimmedJets"),
mets = cms.InputTag("slimmedMETs")
) 


process.p = cms.Path(process.demo)

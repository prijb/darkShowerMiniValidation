import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

do_sig = True

filelist_sig = FileUtils.loadListFromFile("CheckDS/data/sigList.txt")
filelist_bkg = FileUtils.loadListFromFile("CheckDS/data/bkgListQCD.txt")
#filelist_bkg = FileUtils.loadListFromFile("data/bkgListZero.txt")
readFiles_sig = cms.untracked.vstring(*filelist_sig)
readFiles_bkg = cms.untracked.vstring(*filelist_bkg)

if do_sig:
	anlzr_name = 'CheckDS_sig'
	fname = readFiles_sig 	
	outname = 'file:signal.root'
else:
	anlzr_name = 'CheckDS_QCD'
	fname = readFiles_bkg 
	outname = 'file:background.root'	

L1Info = ["L1_DoubleMu0", "L1_DoubleMu0_Mass_Min1", "L1_DoubleMu0_OQ",
 "L1_DoubleMu0_SQ", "L1_DoubleMu0_SQ_OS", "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
 "L1_DoubleMu0er1p5_SQ", "L1_DoubleMu0er1p5_SQ_OS",
 "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu0er1p5_SQ_dR_Max1p4",
 "L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", "L1_DoubleMu0er2p0_SQ_dR_Max1p4",
 "L1_DoubleMu10_SQ", "L1_DoubleMu18er2p1", "L1_DoubleMu4_SQ_OS",
 "L1_DoubleMu4_SQ_OS_dR_Max1p2", "L1_DoubleMu4p5_SQ_OS",
 "L1_DoubleMu4p5_SQ_OS_dR_Max1p2", "L1_DoubleMu4p5er2p0_SQ_OS",
 "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", "L1_DoubleMu9_SQ", "L1_DoubleMu_12_5",
 "L1_DoubleMu_15_5_SQ", "L1_DoubleMu_15_7", "L1_DoubleMu_15_7_Mass_Min1",
 "L1_DoubleMu_15_7_SQ",
 ]

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.demo = cms.EDAnalyzer(anlzr_name,
l1Seeds =  cms.vstring(L1Info),
doL1 = cms.bool(True),
 AlgInputTag = cms.InputTag("gtStage2Digis"),
        l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
ReadPrescalesFromFile=cms.bool(False)
)


process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load(
    "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff"
)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorsNoErrorPropagation_cff')
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = fname
    fileNames = cms.untracked.vstring("/store/user/mcitron/ProjectMetis/HiddenValley_vector_m_2_ctau_10_xiO_1_xiL_1_privateMC_11X_MINIAODSIM_v7_generationForBParkingAddFullGenAddDisp/output_1.root")
    #fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18MiniAODv2/MinBias_TuneCP2_13TeV-pythia8/MINIAODSIM/FSUL18_pilot_106X_upgrade2018_realistic_v16_L1v1_ext1-v2/50000/24AD2EE9-05EF-DA4B-8660-4996045D0D71.root")
    #fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18MiniAODv2/MinBias_TuneCP5_13TeV-pythia8/MINIAODSIM/FSUL18_pilot_106X_upgrade2018_realistic_v16_L1v1-v2/2430000/79986FDF-72AB-CD45-80ED-EB4AB63A90E9.root")
)

print("Processing:", outname)
print("File used:")
print(process.source.fileNames)

print("-----------------------------")

process.TFileService = cms.Service("TFileService",
                                               fileName = cms.string(outname)
                                                                                  )
'''
slimmedMuonsWithUserData = cms.EDProducer("PATMuonUserDataEmbedder",
     src = cms.InputTag("slimmedMuons"),
     userFloats = cms.PSet(
        miniIsoChg = cms.InputTag("isoForMu:miniIsoChg"),
        miniIsoAll = cms.InputTag("isoForMu:miniIsoAll"),
        ptRatio = cms.InputTag("ptRatioRelForMu:ptRatio"),
        ptRel = cms.InputTag("ptRatioRelForMu:ptRel"),
        jetNDauChargedMVASel = cms.InputTag("ptRatioRelForMu:jetNDauChargedMVASel"),
     ),
     userCands = cms.PSet(
        jetForLepJetVar = cms.InputTag("ptRatioRelForMu:jetForLepJetVar") # warning: Ptr is null if no match is found
     ),
)
finalMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("slimmedMuonsWithUserData"),
    cut = cms.string("pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt'))")
)
'''

try:
  process.p = cms.Path(process.demo)
except RuntimeError as e:
  print("Error occurred while processing file:", e)

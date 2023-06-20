// -*- C++ -*-
//
// Package:    CheckGen/CheckDS
// Class:      CheckDS
//
/**\class CheckDS CheckDS.cc CheckGen/CheckDS/plugins/CheckDS.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Matthew Citron <mcitron@ucsb.edu> 10/19/2017
//         Created:  Wed, 21 Sep 2022 03:32:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "TH2.h"
#include "TTree.h"
#include <Math/Vector4D.h>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TMath.h"
#include <Math/Vector4D.h>

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"

#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"




//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class CheckDS_QCD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit CheckDS_QCD(const edm::ParameterSet&);
        ~CheckDS_QCD();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // ----------member data ---------------------------
        // reco::GenParticle
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;  //used to select what tracks to read from configuration file
        edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;  //used to select what tracks to read from configuration file
        edm::EDGetTokenT<reco::BeamSpot> bsToken_;  //used to select what tracks to read from configuration file

        bool doL1_;

        edm::EDGetToken algToken_;
        std::shared_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
        std::vector<std::string> l1Seeds_;
        TTree * outputTree;
        //L1 and muon output vectors for tree
        std::vector<std::string> * l1_name = new std::vector<std::string>();
        std::vector<bool> * l1_result = new std::vector<bool>();
        std::vector<int> * l1_prescale = new std::vector<int>();
        std::vector<float> * muon_pt = new std::vector<float>();
        std::vector<float> * muon_dz = new std::vector<float>();
        std::vector<float> * muon_phi = new std::vector<float>();
        std::vector<float> *  muon_eta = new std::vector<float>();
        std::vector<bool> *  muon_isLoose = new std::vector<bool>();
        std::vector<float> *  muon_dxy = new std::vector<float>();
		std::vector<float> * muon_dxyErr = new std::vector<float>();

        TH1D * pdgIdHist;
        TH1D * ptHist;
        TH1D * massHist;
        TH2D * outHistMass2DOpen;
        TH2D * outHistMass2DClose;
        TH2D * outHistMass2D;
        TH2D * outHistMass2DTrans;
        TH1D * outHistPhi1D;
        TH2D * outHistPhiMass2D;
        TH2D * outHistDPhiMass2D;
        TH2D * outHistMass2DOpenDG;
        TH2D * outHistMass2DCloseDG;
        TH2D * outHistMass2DDG;
        TH2D * outHistMass2DDGTrans;
        TH1D * outHistPhi1DDG;
        TH2D * outHistPhiMass2DDG;
        TH2D * outHistDPhiMass2DDG;
        TH1D * ctauHist;
        TH1D * dispHist;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CheckDS_QCD::CheckDS_QCD(const edm::ParameterSet& iConfig)
    :
        genToken_(consumes<std::vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"))),
        muonToken_(consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"))),
        bsToken_(consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot")))
       // dgMuonToken_(consumes<std::vector<reco::Track>>(edm::InputTag("displacedGlobalMuons")))
{
    doL1_ = iConfig.getParameter<bool>("doL1");
    //now do what ever initialization is needed
    if (doL1_){
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("AlgInputTag"));
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1GtUtils_ = std::make_shared<l1t::L1TGlobalUtil>(iConfig, consumesCollector());
    }
    edm::Service<TFileService> fs;
    outputTree = fs->make<TTree>("outputTree","outputTree");
    outputTree->Branch("l1_name",&l1_name);
    outputTree->Branch("l1_result",&l1_result);
    outputTree->Branch("l1_prescale",&l1_prescale);
    outputTree->Branch("muon_pt",&muon_pt);
    outputTree->Branch("muon_phi",&muon_phi);
    outputTree->Branch("muon_eta",&muon_eta);
    outputTree->Branch("muon_isLoose",&muon_isLoose);
    outputTree->Branch("muon_dxy",&muon_dxy);
    outputTree->Branch("muon_dxyErr",&muon_dxyErr);
    outputTree->Branch("muon_dz",&muon_dz);

    pdgIdHist = fs->make<TH1D>("pdgId","pdgId",300,0,300);
    ptHist = fs->make<TH1D>("pt","pt",200,0,50);
    massHist = fs->make<TH1D>("mass","mass",200,1.5,2.5);
    outHistMass2D = fs->make<TH2D>("dxyVsMass",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DTrans = fs->make<TH2D>("dxyVsMassTrans",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DOpen = fs->make<TH2D>("dxyVsMassOpen",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DClose = fs->make<TH2D>("dxyVsMassClose",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistPhiMass2D = fs->make<TH2D>("phiVsMass",";phi;mass",50,-3.2,3.2,50,1.5,2.5);
    outHistDPhiMass2D = fs->make<TH2D>("dphiVsMass",";dphi;mass",50,-0.1,0.1,50,1.5,2.5);
    outHistPhi1D = fs->make<TH1D>("phi",";phi;mass",100,-0.05,0.05);

    outHistMass2DDG = fs->make<TH2D>("dxyVsMassDG",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DDGTrans = fs->make<TH2D>("dxyVsMassDGTrans",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DOpenDG = fs->make<TH2D>("dxyVsMassOpenDG",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DCloseDG = fs->make<TH2D>("dxyVsMassCloseDG",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistPhiMass2DDG = fs->make<TH2D>("phiVsMassDG",";phi;mass",50,-3.2,3.2,50,1.5,2.5);
    outHistDPhiMass2DDG = fs->make<TH2D>("dphiVsMassDG",";dphi;mass",50,-0.1,0.1,50,1.5,2.5);
    outHistPhi1DDG = fs->make<TH1D>("phiDG",";phi;mass",100,-0.05,0.05);

    ctauHist = fs->make<TH1D>("ctau","ctau",200,0,1000);
    dispHist = fs->make<TH1D>("disp","disp",200,0,1000);

}


CheckDS_QCD::~CheckDS_QCD()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
    void
CheckDS_QCD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::ESHandle<Propagator> propagator_h;
    l1_result->clear();
    l1_name->clear();
    l1_prescale->clear();
    muon_pt->clear();
    muon_phi->clear();
    muon_eta->clear();
    muon_isLoose->clear();
    muon_dxy->clear();
    muon_dxyErr->clear();
    muon_dz->clear();

    if (doL1_){
        l1GtUtils_->retrieveL1(iEvent, iSetup, algToken_);
        for (auto const& l1seed:l1Seeds_){
            bool l1htbit = 0;
            double prescale = -1;
            l1GtUtils_->getFinalDecisionByName(l1seed, l1htbit);
            l1GtUtils_->getPrescaleByName(l1seed, prescale);
            l1_result->push_back(l1htbit);
            l1_name->push_back(l1seed);
            l1_prescale->push_back(prescale);
        }
    }
    // iSetup.get<TrackingComponentsRecord>().get( "SteppingHelixPropagatorAlong",
    //         propagator_h );
    // propagator_ = propagator_h.product();

    // edm::ESHandle<MagneticField> field_h;
    // iSetup.get<IdealMagneticFieldRecord>().get( field_h );
    // bField_ = field_h.product();

    // edm::ESHandle<GlobalTrackingGeometry> geometry_h;
    // iSetup.get<GlobalTrackingGeometryRecord>().get( geometry_h );
    // geometry_ = geometry_h.product();
    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);



    Handle<std::vector<reco::GenParticle>> genParticles;
    Handle<std::vector<pat::Muon>> muons;
   // Handle<std::vector<reco::Track>> dgMuons;
    Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(genToken_, genParticles);
    iEvent.getByToken(bsToken_, beamSpot);
    iEvent.getByToken(muonToken_, muons);
   // iEvent.getByToken(dgMuonToken_, dgMuons);
   // std::vector<reco::TransientTrack> dgMuonsTrans = (*theB).build (dgMuons);
    std::vector<float> outputVecPx;
    std::vector<float> outputVecPy;
    std::vector<float> outputVecPz;
    std::vector<float> outputVecE;
    std::vector<float> outputVecIp;
    std::vector<float> outputVecIp2;
    std::vector<int> outputVecId;
    std::vector<int> outputVecStatus;
    std::vector<pat::Muon> outputVecMatchedMuons;
    //std::vector<reco::Track> outputVecMatchedDGMuons;
    //std::vector<reco::TransientTrack> outputVecMatchedDGMuonsTrans;
    std::vector<reco::GenParticle> outputVecGen;
    for(pat::Muon muon : *muons){
        muon_pt->push_back(muon.pt());
        muon_phi->push_back(muon.phi());
        muon_eta->push_back(muon.eta());
        muon_isLoose->push_back(muon.isLooseMuon());
        if (muon.innerTrack().isNonnull()){
            muon_dxy->push_back(muon.innerTrack()->dxy(*beamSpot));
            muon_dxyErr->push_back(muon.innerTrack()->dxyError());
            muon_dz->push_back(muon.innerTrack()->dz(beamSpot->position()));
        }
        else{
            muon_dxy->push_back(-1.);
            muon_dxyErr->push_back(-1.);
            muon_dz->push_back(-1.);
        }
    }
    outputTree->Fill();
}
        // ------------ method called once each job just before starting event loop  ------------
        void
            CheckDS_QCD::beginJob()
            {
            }

        // ------------ method called once each job just after ending the event loop  ------------
        void
            CheckDS_QCD::endJob()
            {
            }

        // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
        void
            CheckDS_QCD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
                //The following says we do not know what parameters are allowed so do no validation
                // Please change this to state exactly what you do use, even if it is no parameters
                edm::ParameterSetDescription desc;
                desc.setUnknown();
                descriptions.addDefault(desc);

                //Specify that only 'tracks' is allowed
                //To use, remove the default given above and uncomment below
                //ParameterSetDescription desc;
                //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
                //descriptions.addDefault(desc);
            }

        //define this as a plug-in
        DEFINE_FWK_MODULE(CheckDS_QCD);

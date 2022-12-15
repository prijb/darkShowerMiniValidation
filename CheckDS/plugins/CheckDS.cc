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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "TH2.h"
#include <Math/Vector4D.h>

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
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



class CheckDS : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit CheckDS(const edm::ParameterSet&);
      ~CheckDS();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      // reco::GenParticle
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<std::vector<reco::Track>> dgMuonToken_;  //used to select what tracks to read from configuration file
      const MagneticField*          bField_;
      const Propagator*             propagator_;
      const GlobalTrackingGeometry* geometry_;

      TH1D * pdgIdHist;
      TH1D * ptHist;
      TH1D * massHist;
      TH2D * outHistMass2DOpen;
      TH2D * outHistMass2DClose;
      TH2D * outHistMass2D;
      TH1D * outHistPhi1D;
      TH2D * outHistPhiMass2D;
      TH2D * outHistDPhiMass2D;
      TH2D * outHistMass2DOpenDG;
      TH2D * outHistMass2DCloseDG;
      TH2D * outHistMass2DDG;
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
CheckDS::CheckDS(const edm::ParameterSet& iConfig)
 :
    genToken_(consumes<std::vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"))),
    muonToken_(consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"))),
    dgMuonToken_(consumes<std::vector<reco::Track>>(edm::InputTag("displacedGlobalMuons")))

{
   //now do what ever initialization is needed
    edm::Service<TFileService> fs;
    pdgIdHist = fs->make<TH1D>("pdgId","pdgId",300,0,300);
    ptHist = fs->make<TH1D>("pt","pt",200,0,50);
    massHist = fs->make<TH1D>("mass","mass",200,1.5,2.5);
    outHistMass2D = fs->make<TH2D>("dxyVsMass",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DOpen = fs->make<TH2D>("dxyVsMassOpen",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DClose = fs->make<TH2D>("dxyVsMassClose",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistPhiMass2D = fs->make<TH2D>("phiVsMass",";phi;mass",50,-3.2,3.2,50,1.5,2.5);
    outHistDPhiMass2D = fs->make<TH2D>("dphiVsMass",";dphi;mass",50,-0.1,0.1,50,1.5,2.5);
    outHistPhi1D = fs->make<TH1D>("phi",";phi;mass",100,-0.05,0.05);

    outHistMass2DDG = fs->make<TH2D>("dxyVsMassDG",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DOpenDG = fs->make<TH2D>("dxyVsMassOpenDG",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistMass2DCloseDG = fs->make<TH2D>("dxyVsMassCloseDG",";dxy;mass",50,0,20,50,1.5,2.5);
    outHistPhiMass2DDG = fs->make<TH2D>("phiVsMassDG",";phi;mass",50,-3.2,3.2,50,1.5,2.5);
    outHistDPhiMass2DDG = fs->make<TH2D>("dphiVsMassDG",";dphi;mass",50,-0.1,0.1,50,1.5,2.5);
    outHistPhi1DDG = fs->make<TH1D>("phiDG",";phi;mass",100,-0.05,0.05);

    ctauHist = fs->make<TH1D>("ctau","ctau",200,0,1000);
    dispHist = fs->make<TH1D>("disp","disp",200,0,1000);

}


CheckDS::~CheckDS()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
    void
CheckDS::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::ESHandle<Propagator> propagator_h;

    iSetup.get<TrackingComponentsRecord>().get( "SteppingHelixPropagatorAlong",
            propagator_h );
    propagator_ = propagator_h.product();

    edm::ESHandle<MagneticField> field_h;
    iSetup.get<IdealMagneticFieldRecord>().get( field_h );
    bField_ = field_h.product();

    edm::ESHandle<GlobalTrackingGeometry> geometry_h;
    iSetup.get<GlobalTrackingGeometryRecord>().get( geometry_h );
    geometry_ = geometry_h.product();



    Handle<std::vector<reco::GenParticle>> genParticles;
    Handle<std::vector<pat::Muon>> muons;
    Handle<std::vector<reco::Track>> dgMuons;
    iEvent.getByToken(genToken_, genParticles);
    iEvent.getByToken(muonToken_, muons);
    iEvent.getByToken(dgMuonToken_, dgMuons);
    std::vector<float> outputVecPx;
    std::vector<float> outputVecPy;
    std::vector<float> outputVecPz;
    std::vector<float> outputVecE;
    std::vector<float> outputVecIp;
    std::vector<float> outputVecIp2;
    std::vector<int> outputVecId;
    std::vector<int> outputVecStatus;
    std::vector<pat::Muon> outputVecMatchedMuons;
    std::vector<reco::Track> outputVecMatchedDGMuons;
    std::vector<reco::GenParticle> outputVecGen;
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> beamSpotVertex;
    for(auto genParticle : *genParticles){
        const reco::GenParticle * genParticleMother = dynamic_cast<const reco::GenParticle *>(genParticle.mother());    
        while (genParticleMother && !(abs(genParticleMother->pdgId()) == 4900113)){
            // if (genParticleMother->mother() == 0) break;
            genParticleMother =  dynamic_cast<const reco::GenParticle *>(genParticleMother->mother());
        }
        if (genParticleMother){
            // std::cout << genParticle.pdgId() << " " << genParticle.status() << " " << genParticle.mother()->pdgId() << std::endl;
            pdgIdHist->Fill(abs(genParticle.pdgId()));
            double genParticleBeta = genParticleMother->p()/genParticleMother->energy();
            double genParticleGamma = 1./TMath::Sqrt(1.-genParticleBeta*genParticleBeta);
            ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> neutralinoPos = genParticle.vertex();
            ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> parentPos = genParticleMother->vertex();
            beamSpotVertex = parentPos;
            ptHist->Fill(genParticle.pt());
            double displacementGluino = TMath::Sqrt((neutralinoPos-parentPos).Mag2());
            double ctau = displacementGluino*10 / (genParticleBeta*genParticleGamma);
            ctauHist->Fill(ctau);
            dispHist->Fill(displacementGluino); 
            if (genParticle.status() == 1){
                outputVecId.push_back(genParticle.pdgId());
                outputVecStatus.push_back(genParticle.status());
                outputVecIp.push_back(genParticleMother->energy());
                outputVecIp2.push_back(genParticleMother->pz());
                outputVecPx.push_back(genParticle.px());
                outputVecPy.push_back(genParticle.py());
                outputVecPz.push_back(genParticle.pz());
                outputVecE.push_back(genParticle.energy());
                outputVecGen.push_back(genParticle);
                pat::Muon matchedMuon;
                reco::Track matchedDGMuon;
                float minDeltaR = 1000;
                float minDeltaRDG = 1000;
                if (abs(genParticle.pdgId()) == 13 && genParticle.pt() > 3){
                    for(pat::Muon muon : *muons){
                        if (deltaR(genParticle,muon) < 0.05 && deltaR(genParticle,muon) < minDeltaR){
                            if (muon.charge() == genParticle.charge()){
                                // std::cout << "HERE" << std::endl;
                                // std::cout << deltaR(genParticle,muon) << std::endl;
                                minDeltaR = deltaR(genParticle,muon);
                                matchedMuon = muon;
                            }
                        }
                    }
                    for(reco::Track dgMuon : *dgMuons){
                        if (deltaR(genParticle,dgMuon) < 0.05 && deltaR(genParticle,dgMuon) < minDeltaRDG){
                            if (dgMuon.charge() == genParticle.charge()){
                                // std::cout << "HERE" << std::endl;
                                // std::cout << deltaR(genParticle,muon) << std::endl;
                                minDeltaRDG = deltaR(genParticle,dgMuon);
                                matchedDGMuon = dgMuon;
                            }
                        }
                    }
                }
                if (minDeltaR >500){
                    outputVecMatchedMuons.push_back(pat::Muon());
                }
                else{
                    outputVecMatchedMuons.push_back(matchedMuon);
                }
                if (minDeltaRDG >500){
                    outputVecMatchedDGMuons.push_back(reco::Track());
                }
                else{
                    outputVecMatchedDGMuons.push_back(matchedDGMuon);
                }
            }
        }
    }
    std::vector<float> seenIp;
    seenIp.clear();
    std::vector<float> seenIp2;
    seenIp2.clear();
    for (uint iV = 0; iV < outputVecPz.size(); iV++){
        // std::cout << outputVecPx[iV] << " " << outputVecPy[iV] << " " << outputVecPz[iV] << " " << outputVecE[iV] << " " << outputVecId[iV] << " " << outputVecStatus[iV] << " " << outputVecIp[iV] << std::endl;
        bool skip = false;
        float overalX=0;
        float overalY=0;
        float overalZ=0;
        float overalE=0;
        for (uint iV2 = 0; iV2 < seenIp.size(); iV2++){
            if (fabs(outputVecIp[iV] - seenIp[iV2]) < 0.0001 && fabs(outputVecIp2[iV] - seenIp2[iV2]) < 0.0001) skip = true;
        }
        if (skip) continue;
        for (uint iV2 = 0; iV2 < outputVecPz.size(); iV2++){
            // if (iV == iV2) continue;
            // std::cout << outputVecIp[iV] << " " << outputVecIp[iV2] << std::endl;
            if (fabs(outputVecIp[iV] - outputVecIp[iV2])<0.0001 && fabs(outputVecIp2[iV] - outputVecIp2[iV2])<0.0001){
                overalX += outputVecPx[iV2];
                overalY += outputVecPy[iV2];
                overalZ += outputVecPz[iV2];
                overalE += outputVecE[iV2];
                if (iV != iV2){
                    reco::GenParticle genMu1 = outputVecGen[iV];
                    reco::GenParticle genMu2 = outputVecGen[iV2];
                    if (abs(genMu1.pdgId()) == 13 && abs(genMu2.pdgId()) == 13){
                        // std::cout << iV << " " <<abs(genMu1.pdgId()) << " " << outputVecMatchedMuons[iV].pt() << std::endl;
                        // std::cout << iV2 << " " <<abs(genMu2.pdgId()) << " " << outputVecMatchedMuons[iV2].pt() << std::endl;
                        if (outputVecMatchedMuons[iV].pt() > 0 && outputVecMatchedMuons[iV2].pt() > 0){
                            //muons matched
                            pat::Muon mu1 = outputVecMatchedMuons[iV];
                            pat::Muon mu2 = outputVecMatchedMuons[iV2];
                            float dxy = TMath::Sqrt((genMu1.vertex() - beamSpotVertex).Perp2());
                            // std::cout <<  "Gen DXY " << TMath::Sqrt((genMu1.vertex() - beamSpotVertex).Perp2()) << std::endl;
                            float dphi = deltaPhi(genMu1.vertex(),mu1.p4()+mu2.p4());
                            // std::cout <<  "Gen Mass " << TMath::Sqrt((genMu1.p4()+genMu2.p4()).mag2()) << std::endl;
                            // recoMomDir = mu1.p4()+mu2.p4();
                            float deltaPhiMu = mu1.phi() - mu2.phi();
                            if (deltaPhiMu > 3.141) deltaPhiMu -= 2*3.141;
                            if (deltaPhiMu < -3.141) deltaPhiMu += 2*3.141;
                            // if deltaPhiMu < -2*3.141: deltaPhiMu += 2*3.141
                            outHistMass2D->Fill(dxy,TMath::Sqrt((mu1.p4()+mu2.p4()).mag2()));
                            if ((mu1.charge() > 0 && deltaPhiMu > 0) || (mu1.charge() < 0 && deltaPhiMu < 0)){
                                outHistMass2DOpen->Fill(dxy,TMath::Sqrt((mu1.p4()+mu2.p4()).mag2()));
                            }
                            if ((mu1.charge() < 0 && deltaPhiMu > 0) || (mu1.charge() > 0 && deltaPhiMu < 0)){
                                outHistMass2DClose->Fill(dxy,TMath::Sqrt((mu1.p4()+mu2.p4()).mag2()));
                            }
                            if (dxy > 4){
                                // std::cout << "Gen Mu 1" << " " << genMu1.pt() << " " << genMu1.phi() << " " << genMu1.eta() << std::endl;
                                // std::cout << "Rec Mu 1" << " " << mu1.pt() << " " << mu1.phi() << " " << mu1.eta() << std::endl;
                                // std::cout << "Gen Mu 2" <<" "  << genMu2.pt() << " " << genMu2.phi() << " " << genMu2.eta() << std::endl;
                                // std::cout << "Rec Mu 2" << " " << mu2.pt() << " " << mu2.phi() << " " << mu2.eta() << std::endl;
                                // std::cout <<  "Reco Mass " << TMath::Sqrt((mu1.p4()+mu2.p4()).mag2()) << std::endl;
                                outHistPhiMass2D->Fill(genMu1.vertex().phi(),TMath::Sqrt((mu1.p4()+mu2.p4()).mag2()));
                                outHistPhi1D->Fill(deltaPhi(mu1,genMu1));
                                outHistPhi1D->Fill(deltaPhi(mu2,genMu2));
                                outHistDPhiMass2D->Fill(dphi,TMath::Sqrt((mu1.p4()+mu2.p4()).mag2()));
                            }
                        }
                        if (outputVecMatchedDGMuons[iV].pt() > 0 && outputVecMatchedDGMuons[iV2].pt() > 0){
                            //muons matched
                            reco::Track mu1Track = outputVecMatchedDGMuons[iV];
                            reco::Track mu2Track = outputVecMatchedDGMuons[iV2];
                            std::cout << mu1Track.innerPosition() << std::endl;
                            std::cout << mu2Track.innerPosition() << std::endl;
                            std::cout <<  genMu1.vertex() << std::endl;
                            std::cout << std::endl;
                            float dxy = TMath::Sqrt((genMu1.vertex() - beamSpotVertex).Perp2());
                            ROOT::Math::PtEtaPhiMVector mu1 = ROOT::Math::PtEtaPhiMVector(mu1Track.pt(),mu1Track.eta(),mu1Track.phi(),0.107);
                            ROOT::Math::PtEtaPhiMVector mu2 = ROOT::Math::PtEtaPhiMVector(mu2Track.pt(),mu2Track.eta(),mu2Track.phi(),0.107);

                            // reco::Muon mu1 = reco::Muon(mu1Track.charge(),mu1Mom,genMu1.vertex());
                            // reco::Muon mu2 = reco::Muon(mu2Track.charge(),mu2Mom,genMu1.vertex());
                            float absZ_ = genMu1.vertex().z();
                            float radius_ =  TMath::Sqrt(genMu1.vertex().Perp2());

                            GlobalPoint vertex1 = GlobalPoint(
                                    mu1Track.innerPosition().x(), mu1Track.innerPosition().y(), mu1Track.innerPosition().z() );
                            GlobalVector momentum1 = GlobalVector( 
                                    mu1Track.innerMomentum().x(), mu1Track.innerMomentum().y(), mu1Track.innerMomentum().z() );

                            GlobalTrajectoryParameters trackParams( vertex1, momentum1, mu1Track.charge(), bField_ );
                            FreeTrajectoryState fts( trackParams );

                            TrajectoryStateOnSurface propagatedInfo = propagator_->propagate(fts, *Cylinder::build( radius_, Surface::PositionType( 0, 0, 0 ), Surface::RotationType() ) );
                            if ( propagatedInfo.globalPosition().z() > absZ_ ) {
                                propagatedInfo = propagator_->propagate(
                                        fts, *Plane::build( Surface::PositionType( 0, 0, absZ_ ),
                                            Surface::RotationType() ) );
                            } else if ( propagatedInfo.globalPosition().z() < -absZ_ ) {
                                propagatedInfo = propagator_->propagate(
                                        fts, *Plane::build( Surface::PositionType( 0, 0, -absZ_ ),
                                            Surface::RotationType() ) );
                            }
                            std::cout<< propagatedInfo.globalPosition().x() << std::endl;



                            // std::cout <<  "Gen DXY " << TMath::Sqrt((genMu1.vertex() - beamSpotVertex).Perp2()) << std::endl;
                            float dphi = deltaPhi(genMu1.vertex(),mu1+mu2);
                            // std::cout <<  "Gen Mass " << TMath::Sqrt((genMu1.p4()+genMu2.p4()).mag2()) << std::endl;
                            // recoMomDir = mu1.p4()+mu2.p4();
                            float deltaPhiMu = mu1.phi() - mu2.phi();
                            if (deltaPhiMu > 3.141) deltaPhiMu -= 2*3.141;
                            if (deltaPhiMu < -3.141) deltaPhiMu += 2*3.141;
                            // if deltaPhiMu < -2*3.141: deltaPhiMu += 2*3.141
                            outHistMass2DDG->Fill(dxy,TMath::Sqrt((mu1+mu2).mag2()));
                            if ((mu1Track.charge() > 0 && deltaPhiMu > 0) || (mu1Track.charge() < 0 && deltaPhiMu < 0)){
                                outHistMass2DOpenDG->Fill(dxy,TMath::Sqrt((mu1+mu2).mag2()));
                            }
                            if ((mu1Track.charge() < 0 && deltaPhiMu > 0) || (mu1Track.charge() > 0 && deltaPhiMu < 0)){
                                outHistMass2DCloseDG->Fill(dxy,TMath::Sqrt((mu1+mu2).mag2()));
                            }
                            if (dxy > 4){
                                // std::cout << "Gen Mu 1" << " " << genMu1.pt() << " " << genMu1.phi() << " " << genMu1.eta() << std::endl;
                                // std::cout << "Rec Mu 1" << " " << mu1.pt() << " " << mu1.phi() << " " << mu1.eta() << std::endl;
                                // std::cout << "Gen Mu 2" <<" "  << genMu2.pt() << " " << genMu2.phi() << " " << genMu2.eta() << std::endl;
                                // std::cout << "Rec Mu 2" << " " << mu2.pt() << " " << mu2.phi() << " " << mu2.eta() << std::endl;
                                // std::cout <<  "Reco Mass " << TMath::Sqrt((mu1+mu2).mag2()) << std::endl;
                                outHistPhiMass2DDG->Fill(genMu1.vertex().phi(),TMath::Sqrt((mu1+mu2).mag2()));
                                outHistPhi1DDG->Fill(deltaPhi(mu1,genMu1));
                                outHistPhi1DDG->Fill(deltaPhi(mu2,genMu2));
                                outHistDPhiMass2DDG->Fill(dphi,TMath::Sqrt((mu1+mu2).mag2()));
                            }
                        }
                    }
                }
            }

        }
        seenIp.push_back(outputVecIp[iV]);
        seenIp2.push_back(outputVecIp2[iV]);
        float overalM = TMath::Sqrt(overalE*overalE-overalX*overalX-overalY*overalY-overalZ*overalZ);
        massHist->Fill(overalM);
    }

}


// ------------ method called once each job just before starting event loop  ------------
    void
CheckDS::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
    void
CheckDS::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CheckDS::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(CheckDS);

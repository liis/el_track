// -*- C++ -*-
//
// Package:    MakeTrackValTree
// Class:      MakeTrackValTree
// 
/**\class MakeTrackValTree MakeTrackValTree.cc MakeTree_2/MakeTrackValTree/src/MakeTrackValTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Liis Rebane
//         Created:  Fri Apr 11 14:28:05 EEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//-------------------------------------------------

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
//#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimGeneral/TrackingAnalysis/interface/SimHitTPAssociationProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
//
// class declaration
//

#define MAXPART 100

class MakeTrackValTree : public edm::EDAnalyzer {
   public:
      explicit MakeTrackValTree(const edm::ParameterSet&);
      ~MakeTrackValTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

  // ------------------ thresholds for TPselector ----------
  double ptMinTP, minRapidityTP, maxRapidityTP, tipTP, lipTP, signalOnlyTP, chargedOnlyTP, stableOnlyTP; 
  std::vector<int> pdgIdTP;
  int minHitTP;
  // -------------------- for tree -------------
  TTree *trackValTree_;
  int run_nr_, evt_nr_, lumi_nr_;
  int np_gen_, np_reco_, np_reco_toGen_,np_fake_;
  int is_reco_matched_[MAXPART], is_gen_matched_[MAXPART];
  int gen_pdgId_[MAXPART];
  double reco_pt_[MAXPART], reco_eta_[MAXPART], reco_phi_[MAXPART], fake_pt_[MAXPART], fake_eta_[MAXPART], fake_phi_[MAXPART];
  double gen_pt_[MAXPART],  gen_eta_[MAXPART], gen_phi_[MAXPART];

  bool is_gsf_;
  edm::InputTag track_label_gsf_, track_label_, el_seed_label_;

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
MakeTrackValTree::MakeTrackValTree(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  is_gsf_ = iConfig.getParameter<bool>("isGSF");
  track_label_gsf_ = iConfig.getParameter<edm::InputTag>("trackLabelGSF");
  track_label_ = iConfig.getParameter<edm::InputTag>("trackLabel");
  el_seed_label_ = iConfig.getParameter<edm::InputTag>("elSeedLabel");

  ptMinTP =   iConfig.getParameter<double>("ptMinTP");
  minRapidityTP = iConfig.getParameter<double>("minRapidityTP");
  maxRapidityTP = iConfig.getParameter<double>("maxRapidityTP");
  tipTP =  iConfig.getParameter<double>("tipTP");
  lipTP =  iConfig.getParameter<double>("lipTP");
  minHitTP =  iConfig.getParameter<int>("minHitTP");
  signalOnlyTP =  iConfig.getParameter<bool>("signalOnlyTP");
  chargedOnlyTP =  iConfig.getParameter<bool>("chargedOnlyTP");
  stableOnlyTP =   iConfig.getParameter<bool>("stableOnlyTP");
  pdgIdTP =  iConfig.getParameter<std::vector<int> >("pdgIdTP");
}


MakeTrackValTree::~MakeTrackValTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MakeTrackValTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   run_nr_ = iEvent.id().run();
   lumi_nr_ = iEvent.id().luminosityBlock();
   evt_nr_ = iEvent.id().event();

   for (unsigned int i = 0; i < MAXPART; i++){
     //----- TPs ---------
     is_reco_matched_[i] = -10;
     gen_eta_[i] = -99;
     gen_phi_[i] = -99;
     gen_pt_[i] = -99;
     gen_pdgId_[i] = -99;

     //------ reco tracks --
     is_gen_matched_[i] = -10;
     reco_pt_[i] = -99;
     reco_eta_[i] = -99;
     reco_phi_[i] = -99;
     fake_pt_[i] = -99;
     fake_eta_[i] = -99;
     fake_phi_[i] = -99;
     
   }


   edm::InputTag track_label; 
   if(is_gsf_)
     track_label = track_label_gsf_;
   else
     track_label = track_label_;

   edm::Handle<edm::View<reco::Track> > trackCollection; //reconstructed tracks
   iEvent.getByLabel(track_label, trackCollection);
   
   edm::Handle<edm::View<reco::ElectronSeed> > elSeedCollection; 
   iEvent.getByLabel(el_seed_label_, elSeedCollection);

   edm::Handle<TrackingParticleCollection>  TPCollection ; //simulated tracks
   iEvent.getByLabel("mix","MergedTrackTruth",TPCollection);

   const TrackingParticleCollection tPCeff = *(TPCollection.product());
   TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP);


   edm::ESHandle<TrackAssociatorBase> myAssociator;
   iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", myAssociator);
   reco::SimToRecoCollection simRecColl = myAssociator->associateSimToReco(trackCollection, TPCollection, &iEvent, &iSetup);
   reco::RecoToSimCollection recSimColl = myAssociator->associateRecoToSim(trackCollection, TPCollection, &iEvent, &iSetup);
   

   //   reco::SimToRecoCollectionSeed simRecCollSeed = myAssociator->associateSimToReco(elSeedCollection, TPCollection, &iEvent, &iSetup);
   // < -- FIXME, associator for seeds, otherwise match manually

   std::cout<<" nr. electron seeds: " << elSeedCollection->size()<<std::endl;
   std::cout<<" nr. tracks: " << trackCollection->size()<<std::endl;
   
   np_gen_ = 0;
   for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
     TrackingParticleRef tpr(TPCollection, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
     
     if( !tpSelector(*tp) ) continue;
     if( tp->pt() < 1 ) continue; // such tracks are irrelevant for electrons
     
     //     std::vector<PSimHit> simhits=tp->trackPSimHit(DetId::Tracker); //sim hits no longer saved for TP-s in 7_1_x
     int nr_simhits = 1;
     gen_eta_[np_gen_] = tp->eta();
     gen_phi_[np_gen_] = tp->phi();
     gen_pt_[np_gen_] = tp->pt();
     gen_pdgId_[np_gen_] = tp->pdgId();
     
     //--------check for matched reco track defined by AssociatorByHits (efficiency and fake-rate plots)-----------
     const reco::Track* matchedTrackPointer=0;
     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
     int nSharedHits = 0; int nRecoTrackHits=0; double matchedRecoTrackPt =0; double matchedRecoTrackEta=0;

     if(simRecColl.find(tpr) != simRecColl.end()){
       rt = simRecColl[tpr]; // find sim-to-reco association
       if ( rt.size()!=0 ) {
	 matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track                                                 
	 nSharedHits = rt.begin()->second;
	 nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
	 matchedRecoTrackPt = matchedTrackPointer->pt();
	 matchedRecoTrackEta = matchedTrackPointer->eta();
	 is_reco_matched_[np_gen_] = 1;
       }
     }
     else
       is_reco_matched_[np_gen_] = 0;

     std::cout<<"Found tracking particle with pt = "<<tp->pt()<<", eta = "<<tp->eta()<<", nr simhits = "<< nr_simhits<<std::endl;
     std::cout<<"nr shared hits = "<<nSharedHits<<", nr reco hits "<<nRecoTrackHits<<", matched reco track pt = "<<matchedRecoTrackPt<<", matched reco track eta = "<<matchedRecoTrackEta<<std::endl;

     np_gen_++; // count tracking particles that pass the standard quality cuts
   } // <-- end loop over tracking particles

   std::cout<<"nr. tracking particles: " <<np_gen_ << std::endl;

   //------------------------ loop over reco tracks --------------------------
   np_reco_ = 0;
   np_reco_toGen_ = 0;
   np_fake_ = 0;
   for( edm::View<reco::Track>::size_type i=0; i<trackCollection->size(); i++){
     edm::RefToBase<reco::Track> track(trackCollection, i);
     
     if(track->pt() < 1 ) continue; // reduce noise, irrelevant for electrons

     reco_pt_[np_reco_] = track->pt();
     reco_phi_[np_reco_] = track->phi();
     reco_eta_[np_reco_] = track->eta();
     //     reco_nrRecoHits_[np_reco_] = track->numberOfValidHits();

     //---------------------check association maps to sim-tracks---------------------
     std::vector<std::pair<TrackingParticleRef, double> > tp;

     if(recSimColl.find(track) != recSimColl.end()){ // if matched
       tp = recSimColl[track];
       if(tp.size() != 0) {
	 is_gen_matched_[np_reco_] = 1;
	 np_reco_toGen_++;
	 TrackingParticleRef matchedTrackingParticle = tp.begin()->first;
	 //	 reco_nrSharedHits_[np_reco_] = tp.begin()->second;
	 // reco_nrSimHits_[np_reco_] = (matchedTrackingParticle->trackPSimHit(DetId::Tracker) ).size();
       }
     }

     else{  // if fake
       is_gen_matched_[np_reco_] = 0;
       fake_eta_[np_fake_] = track->eta();
       fake_pt_[np_fake_] = track->pt();
       fake_phi_[np_fake_] = track->phi();

       np_fake_++;
     }
     np_reco_++;     
   } // <-- end loop over reco tracks

   trackValTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTrackValTree::beginJob()
{
  edm::Service<TFileService> fs;
  trackValTree_ = fs->make<TTree>("trackValTree", "trackValTree");

  trackValTree_->Branch("run_nr", &run_nr_, "run_nr/I");
  trackValTree_->Branch("evt_nr", &evt_nr_, "evt_nr/I");
  trackValTree_->Branch("lumi_nr", &lumi_nr_, "lumi_nr/I");


  //---------------- reco tracks -------------------
  trackValTree_->Branch("np_reco", &np_reco_, "np_reco/I");
  trackValTree_->Branch("np_reco_toGen", &np_reco_toGen_, "np_reco_toGen/I");
  trackValTree_->Branch("np_fake", &np_fake_, "np_fake/I");

  trackValTree_->Branch("reco_gen_matched", is_gen_matched_, "reco_gen_matched[np_reco]/I");
  trackValTree_->Branch("reco_pt", reco_pt_, "reco_pt[np_reco]/D");
  trackValTree_->Branch("reco_eta", reco_eta_, "reco_eta[np_reco]/D");
  trackValTree_->Branch("reco_phi", reco_phi_, "reco_phi[np_reco]/D");
  trackValTree_->Branch("fake_pt", fake_pt_, "fake_pt[np_fake]/D");
  trackValTree_->Branch("fake_eta", fake_eta_, "fake_eta[np_fake]/D");
  trackValTree_->Branch("fake_phi", fake_phi_, "fake_phi[np_fake]/D");

  // ----------------- TPs ----------------------------
  trackValTree_->Branch("np_gen", &np_gen_, "np_gen/I");
  trackValTree_->Branch("gen_reco_matched", is_reco_matched_, "gen_reco_matched[np_gen]/I");
  trackValTree_->Branch("gen_pdgId", gen_pdgId_, "gen_pdgId[np_gen]/I");
  trackValTree_->Branch("gen_pt", gen_pt_, "gen_pt[np_gen]/D");
  trackValTree_->Branch("gen_eta", gen_eta_, "gen_eta[np_gen]/D");
  trackValTree_->Branch("gen_phi", gen_phi_, "gen_phi[np_gen]/D");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MakeTrackValTree::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MakeTrackValTree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MakeTrackValTree::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MakeTrackValTree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MakeTrackValTree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeTrackValTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeTrackValTree);

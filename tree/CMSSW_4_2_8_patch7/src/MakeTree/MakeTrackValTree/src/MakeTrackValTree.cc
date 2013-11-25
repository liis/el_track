// -*- C++ -*-
//
// Package:    MakeTrackValTree
// Class:      MakeTrackValTree
// 
/**\class MakeTrackValTree MakeTrackValTree.cc MakeTree/MakeTrackValTree/src/MakeTrackValTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Liis Rebane (ETHZ) [liis]
//         Created:  Mon Nov 18 16:27:35 CET 2013
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


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
//#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"

#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
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
  double ptMinTP, minRapidityTP, maxRapidityTP, tipTP, lipTP, signalOnlyTP, chargedOnlyTP, stableOnlyTP, minAbsEtaTP, maxAbsEtaTP;
  std::vector<int> pdgIdTP;
  int minHitTP;

  TTree *trackValTree_;

  int np_gen_, np_gen_toReco_, np_reco_, np_reco_toGen_;
  //  bool is_reco_matched_[MAXPART], 
  //  bool is_gen_matched_[MAXPART];
  int is_reco_matched_[MAXPART], is_gen_matched_[MAXPART];
  double gen_pt_[MAXPART], gen_eta_[MAXPART], reco_pt_[MAXPART], reco_eta_[MAXPART];
  int gen_pdgId_[MAXPART], gen_nrSimHits_[MAXPART], gen_nrSharedHits_[MAXPART], gen_nrRecoHits_[MAXPART], reco_nrRecoHits_[MAXPART], reco_nrSimHits_[MAXPART], reco_nrSharedHits_[MAXPART];
  
  bool debug_;

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
  //--------------get cut thresholds for sim tracks----------------
  debug_ = iConfig.getParameter<bool>("debug");
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
  minAbsEtaTP =  iConfig.getParameter<double>("minAbsEtaTP");
  maxAbsEtaTP =  iConfig.getParameter<double>("maxAbsEtaTP");
  //-----------------------------------------------------------

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
  //reinitialize at each event
  for (unsigned int i = 0; i < MAXPART; i++){
    gen_eta_[i] = -999;
    gen_pt_[i] = -999;
    gen_pdgId_[i] = -999;
    gen_nrSimHits_[i] = -999;
    gen_nrRecoHits_[i] = -999;
    gen_nrSharedHits_[i] = -999;
    is_reco_matched_[i] = 0;
    is_gen_matched_[i] = 0;
    reco_eta_[i] = -999;
    reco_pt_[i] = -999;
    reco_nrSimHits_[i] = -10;
    reco_nrRecoHits_[i] = -10;
    reco_nrSharedHits_[i] = -10;
  }

  edm::Handle<edm::View<reco::Track> >   trackCollection; //reconstructed tracks
  iEvent.getByLabel("cutsRecoTracksHp", trackCollection);
  //  iEvent.getByLabel("generalTracks", trackCollection);

  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ; //simulated tracks
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());

  edm::ESHandle<TrackAssociatorBase> theAssociator; //create track associators for MC truth matching
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits" ,theAssociator);

  reco::SimToRecoCollection simRecColl = theAssociator->associateSimToReco(trackCollection, TPCollectionHeff, &iEvent);
  reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(trackCollection, TPCollectionHeff, &iEvent);

  //  edm::Handle<reco::SimToRecoCollection> simtorecoCollection;
  //iEvent.getByLabel("trackingParticleRecoTrackAssociation", simtorecoCollection);

  if(debug_)
    std::cout <<"Number of simulated tracks = "<<tPCeff.size() << std::endl;

  TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,               
								 signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP,minAbsEtaTP,maxAbsEtaTP);  
  //----------------------Loop over sim tracks--------------------------------
  np_gen_ = 0;
  np_gen_toReco_ = 0;
  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

    if( !tpSelector(*tp) ) continue;

    gen_eta_[np_gen_] = tp->eta();
    gen_pt_[np_gen_] = tp->pt();
    gen_pdgId_[np_gen_] = tp->pdgId();

    std::vector<PSimHit> simhits=tp->trackPSimHit(DetId::Tracker);        
    gen_nrSimHits_[np_gen_] =  simhits.size();

    //--------check for matched reco tracks-----------

    const reco::Track* matchedTrackPointer=0;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
    int nSharedHits(0), nRecoTrackHits(0);

    if(simRecColl.find(tpr) != simRecColl.end()){
      rt = simRecColl[tpr]; // find sim-to-reco association
      if ( rt.size()!=0 ) {
        np_gen_toReco_++; //count the number of simTracks that have a recoTrack associated
	
	
	matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track                                                 
	nSharedHits = rt.begin()->second;
	nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
	
	//	if( (nRecoTrackHits < nSharedHits) || ( (int)simhits.size() < nSharedHits) ){
	//  std::cout<<"Found a matched reco track with:"<<"nr reco track hits = "<<nRecoTrackHits<<" and nr shared hits = "<<nSharedHits<<", nr sim track hits = "<<simhits.size()<<std::endl;
	//}
	
	if(debug_){
	  std::cout<<"Found a matched reco track with:"<<"nr reco track hits = "<<nRecoTrackHits<<" and nr shared hits = "<<nSharedHits<<", nr sim track hits = "<<simhits.size()<<", eta = "<<tp->eta()<<", phi = "<< tp->phi()<<std::endl;
	}
	is_reco_matched_[np_gen_] = 1;
	gen_nrSharedHits_[np_gen_] = nSharedHits;
	gen_nrRecoHits_[np_gen_] = nRecoTrackHits; 
      }
    }
     np_gen_++; //count selected sim tracks passing the selection (important to keep in the end of the loop for tree items)      
  }
  //--------------------Loop over reco tracks-----------------------------------
  np_reco_ = 0;
  np_reco_toGen_ = 0;
  
  for( edm::View<reco::Track>::size_type i=0; i<trackCollection->size(); i++){
    edm::RefToBase<reco::Track> track(trackCollection, i);
    reco_pt_[np_reco_] = track->pt();
    reco_eta_[np_reco_] = track->eta();
    reco_nrRecoHits_[np_reco_] = track->numberOfValidHits();

    const TrackingParticle* matchedTrackingParticlePointer = 0;
    std::vector<std::pair<TrackingParticleRef, double> > tp;
    int nSharedHits(0), nSimTrackHits(0);
    if(recSimColl.find(track) != recSimColl.end()){
      tp = recSimColl[track];
      if(tp.size() != 0) {
	is_gen_matched_[np_reco_] = 1;
	np_reco_toGen_++;
	matchedTrackingParticlePointer = tp.begin()->first.get();
	reco_nrSharedHits_[np_reco_] = tp.begin()->second;
	
	reco_nrSimHits_[np_reco_] = (matchedTrackingParticlePointer->trackPSimHit(DetId::Tracker) ).size();

	std::cout<<"Found a matched sim track with:"<<"nr sim track hits = "<<reco_nrSimHits_[np_reco_]<<" and nr shared hits = "<<tp.begin()->second<<", nr reco track hits = "<<track->numberOfValidHits()<<", eta = "<<track->eta()<<", phi = "<<track->phi()<<std::endl;
      }
    }
    np_reco_++;
  }
  //---------------------------------------------------------------------------
  if(debug_){
    std::cout<<"-----------------------------------"<<std::endl;
    std::cout<<"Total simulated = "<<np_gen_<<std::endl;
    std::cout<<"Total simulated SimToReco assoc. = "<<np_gen_toReco_<<std::endl;
    std::cout<<"Total reconstructed = "<<np_reco_<<std::endl;
    std::cout<<"Total reconstructed RecoToSim assoc. = "<<np_reco_toGen_<<std::endl;
    std::cout<<"Total fakes = "<<np_reco_ - np_reco_toGen_<<std::endl;
    std::cout<<"-------------------------------------"<<std::endl;
  }
  trackValTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTrackValTree::beginJob()
{
  edm::Service<TFileService> fs;
  trackValTree_ = fs->make<TTree>("trackValTree","trackValTree");
 
  trackValTree_->Branch("np_reco", &np_reco_, "np_reco/I");
  trackValTree_->Branch("np_reco_toGen", &np_reco_toGen_, "np_reco_toGen/I");

  trackValTree_->Branch("reco_gen_matched", is_gen_matched_, "reco_gen_matched[np_reco]/I");
  trackValTree_->Branch("reco_pt", reco_pt_, "reco_pt[np_reco]/D");
  trackValTree_->Branch("reco_eta", reco_eta_, "reco_eta[np_reco]/D");
  trackValTree_->Branch("reco_nrSimHits", reco_nrSimHits_, "reco_nrSimHits[np_reco]/I");
  trackValTree_->Branch("reco_nrRecoHits", reco_nrRecoHits_, "reco_nrRecoHits[np_reco]/I");
  trackValTree_->Branch("reco_nrSharedHits", reco_nrSharedHits_, "reco_nrSharedHits[np_reco]/I");

  trackValTree_->Branch("np_gen",&np_gen_,"np_gen/I"); // needs to be filled in order to fill x[np]
  trackValTree_->Branch("np_gen_toReco", &np_gen_toReco_, "np_gen_ToReco/I");

  trackValTree_->Branch("gen_reco_matched", is_reco_matched_, "gen_reco_matched[np_gen]/I");
  trackValTree_->Branch("gen_pdgId", gen_pdgId_, "gen_pdgId[np_gen]/I");
  trackValTree_->Branch("gen_pt", gen_pt_, "gen_pt[np_gen]/D");
  trackValTree_->Branch("gen_eta", gen_eta_, "gen_eta[np_gen]/D");
  trackValTree_->Branch("gen_nrSimHits", gen_nrSimHits_, "gen_nrSimHits[np_gen]/I");
  trackValTree_->Branch("gen_nrRecoHits", gen_nrRecoHits_, "gen_nrRecoHits[np_gen]/I");
  trackValTree_->Branch("gen_nrSharedHits", gen_nrSharedHits_, "gen_nrSharedHits[np_gen]/I");  
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


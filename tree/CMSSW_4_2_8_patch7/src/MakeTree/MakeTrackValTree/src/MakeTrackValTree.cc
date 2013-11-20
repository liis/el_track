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

  int np_gen_, np_gen_toReco_, np_reco_;
  bool is_reco_matched_[MAXPART];
  double gen_pt_[MAXPART], gen_eta_[MAXPART];


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
  debug_ = true;

  //reinitialize at each event
  for (unsigned int i = 0; i < MAXPART; i++){
    gen_eta_[i] = -9999;
    gen_pt_[i] = -9999;
    is_reco_matched_[i] = false;
  }

  edm::Handle<edm::View<reco::Track> >   trackCollection; //reconstructed tracks
  iEvent.getByLabel("generalTracks", trackCollection);
  np_reco_ = 0;

  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ; //simulated tracks
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());

  edm::ESHandle<TrackAssociatorBase> theAssociator; //create track associators for MC truth matching
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHitsRecoDenom" ,theAssociator);
  reco::SimToRecoCollection simRecColl = theAssociator->associateSimToReco(trackCollection, TPCollectionHeff, &iEvent);
  reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(trackCollection, TPCollectionHeff, &iEvent);


  edm::Handle<reco::SimToRecoCollection> simtorecoCollection;
  iEvent.getByLabel("trackingParticleRecoTrackAssociation", simtorecoCollection);

  if(debug_)
    std::cout <<"Number of simulated tracks = "<<tPCeff.size() << std::endl;

  TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,               
								 signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP,minAbsEtaTP,maxAbsEtaTP);  
  np_gen_ = 0;
  np_gen_toReco_ = 0;
  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

    if( !tpSelector(*tp) ) continue;
      
    gen_eta_[np_gen_] = tp->eta();
    gen_pt_[np_gen_] = tp->pt();
    std::vector<PSimHit> simhits=tp->trackPSimHit(DetId::Tracker);        
    std::cout<<"number of hits on sim track "<<simhits.size()<<std::endl;
    
    //--------check for matched reco tracks-----------

    const reco::Track* matchedTrackPointer=0;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
    int nSharedHits(0), nRecoTrackHits(0);
    if(simRecColl.find(tpr) != simRecColl.end()){
      rt = simRecColl[tpr];
      matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track                                                 
      nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
      nSharedHits = rt.begin()->second;
      std::cout<<"nr reco track hits = "<<nRecoTrackHits<<std::endl;
      std::cout<<"nr shared hits = "<<nSharedHits<<std::endl;

      if ( rt.size()!=0 ) {
	np_gen_toReco_++; //This counter counts the number of simTracks that have a recoTrack associated
	is_reco_matched_[np_gen_] = true;
	//	matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track
	//	nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
	//	nSharedHits = rt.begin()->second;
	//std::cout<<"nr reco track hits = "<<nRecoTrackHits<<std::endl;
	//std::cout<<"nr shared hits = "<<nSharedHits<<std::endl;
      }
    }
    np_gen_++;
  }
  
  if(debug_){
    std::cout<<"nr of selected gen tracks = "<<np_gen_<<std::endl;
    std::cout<<"nr of reco matched gen tracks = "<<np_gen_toReco_<<std::endl;
  }
  trackValTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTrackValTree::beginJob()
{
  edm::Service<TFileService> fs;
  trackValTree_ = fs->make<TTree>("trackValTree","trackValTree");
 
  trackValTree_->Branch("np_reco_", &np_reco_, "np_reco_/I");

  trackValTree_->Branch("np_gen",&np_gen_,"np_gen/I"); // needs to be filled in order to fill x[np]
  trackValTree_->Branch("np_gen_toReco", &np_gen_toReco_, "np_gen_ToReco/I");

  trackValTree_->Branch("gen_reco_matched", &is_reco_matched_, "gen_reco_matched/B");
  trackValTree_->Branch("gen_pt", gen_pt_, "gen_pt[np_gen]/D");
  trackValTree_->Branch("gen_eta", gen_eta_, "gen_eta[np_gen]/D");
  
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



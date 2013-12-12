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
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h" 
#include "DataFormats/SiStripDetId/interface/TIBDetId.h" 
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"


#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"

#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TTree.h"


//
// class declaration
//
#define MAXPART 100
#define MAXHIT 50

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
      //------functions----------
      std::vector<int> getHitPosition(DetId, bool);
      GlobalVector getHitMomentum(std::vector<PSimHit>::const_iterator, edm::ESHandle<TrackerGeometry>, bool);
      std::vector<PSimHit> getGoodHits(std::vector<PSimHit>, bool);
    
  // ----------member data ---------------------------
  double ptMinTP, minRapidityTP, maxRapidityTP, tipTP, lipTP, signalOnlyTP, chargedOnlyTP, stableOnlyTP, minAbsEtaTP, maxAbsEtaTP;
  std::vector<int> pdgIdTP;
  int minHitTP;

  TTree *trackValTree_;

  int np_gen_, np_gen_toReco_, np_reco_, np_reco_toGen_, np_fake_, run_nr_, evt_nr_, lumi_nr_;
  int is_reco_matched_[MAXPART], is_gen_matched_[MAXPART];
  
  float gen_hit_pt_[MAXPART][MAXHIT], gen_hit_eta_[MAXPART][MAXHIT], gen_hit_phi_[MAXPART][MAXHIT];
  int gen_hit_subdetector_[MAXPART][MAXHIT], gen_hit_layer_[MAXPART][MAXHIT];

  double gen_pt_[MAXPART], gen_eta_[MAXPART], gen_phi_[MAXPART], gen_matched_pt_[MAXPART], gen_matched_eta_[MAXPART], gen_matched_phi_[MAXPART], gen_dxy_[MAXPART], gen_dz_[MAXPART];
  double reco_pt_[MAXPART], reco_eta_[MAXPART], reco_phi_[MAXPART], fake_pt_[MAXPART], fake_eta_[MAXPART], fake_phi_[MAXPART];
  int gen_pdgId_[MAXPART], gen_nrSimHits_[MAXPART], gen_nrRecoHits_[MAXPART], reco_nrRecoHits_[MAXPART], reco_nrSimHits_[MAXPART];
  double reco_nrSharedHits_[MAXPART], gen_nrSharedHits_[MAXPART];
  
  bool debug_, is_gsf_;

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
  is_gsf_ = iConfig.getParameter<bool>("isGSF");

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

GlobalVector MakeTrackValTree::getHitMomentum(std::vector<PSimHit>::const_iterator it_hit, edm::ESHandle<TrackerGeometry> tracker, bool hitdebug = false){
  DetId dId = DetId(it_hit->detUnitId() );
  LocalVector local_p = it_hit->momentumAtEntry();
  const GeomDetUnit* detunit = tracker->idToDetUnit(dId.rawId());
  GlobalVector global_p = detunit->toGlobal(local_p);
  if( hitdebug )
    std::cout<<"hit pt = "<<it_hit->momentumAtEntry().perp()<<", eta = "<<it_hit->momentumAtEntry().eta()<<", phi = "<<it_hit->momentumAtEntry().phi()<<std::endl;

  return global_p;
}

std::vector<int> MakeTrackValTree::getHitPosition(DetId dId, bool hitdebug = false) {
  //  DetId dId = DetId(it_hit->detUnitId() );

  int layerNumber = 0;
  int subdetId = static_cast<int>(dId.subdetId());

  if(subdetId == PixelSubdetector::PixelBarrel){
    PXBDetId pxbid(dId.rawId());
    layerNumber = pxbid.layer();
    if(hitdebug)
      std::cout<<"Hit at pixel barrel layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }

  else if(subdetId == PixelSubdetector::PixelEndcap){
    PXFDetId pxfid(dId.rawId());
    layerNumber = pxfid.disk();
    if(hitdebug)
   std::cout<<"Hit at pixel endcap layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
      
  else if( subdetId == StripSubdetector::TIB){
    TIBDetId tibid(dId.rawId());
    layerNumber = tibid.layer();
    if(hitdebug)
      std::cout<<"Hit at TIB layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
    
  else if( subdetId == StripSubdetector::TOB){
    TOBDetId tobid(dId.rawId());
    layerNumber = tobid.layer();
    if(hitdebug)
      std::cout<<"Hit at TOB layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
  else if( subdetId == StripSubdetector::TID){
    TIDDetId tidid(dId.rawId());
    layerNumber = tidid.wheel();
    if(hitdebug)
      std::cout<<"Hit at TID layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
  else if( subdetId == StripSubdetector::TEC ){
    TECDetId tecid(dId.rawId());
    layerNumber = tecid.wheel();
    if(hitdebug)
      std::cout<<"Hit at TEC layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }

  std::vector<int> hit_position;
  hit_position.push_back(subdetId);
  hit_position.push_back(layerNumber);

  return hit_position;
}

std::vector<PSimHit>  MakeTrackValTree::getGoodHits(std::vector<PSimHit> simhits, bool hitdebug = false) 
//loop over all hits in the event and omit subsequent hit in the same subdetector layer 
{
  
  std::vector<PSimHit> good_hits;
  for(std::vector<PSimHit>::const_iterator it_hit = simhits.begin(); it_hit != simhits.end(); it_hit++ ){
    DetId dId = DetId(it_hit->detUnitId() ); 
    std::vector<int> hitposition = getHitPosition(dId, hitdebug);

    int hit_subdetector = hitposition.at(0);
    int hit_layer = hitposition.at(1);
      
    if( !good_hits.size() )
      good_hits.push_back(*it_hit);
    else{ //compare to the previous hit position
      DetId dId_good = DetId(good_hits.back().detUnitId() );

      int good_hit_subdetector = (getHitPosition(dId_good)).at(0);
      int good_hit_layer = (getHitPosition(dId_good)).at(1);

      if( hit_subdetector == good_hit_subdetector && hit_layer == good_hit_layer){
	good_hits.pop_back(); //replace the previous hit with a later one
	good_hits.push_back(*it_hit);
      }
      else
	good_hits.push_back(*it_hit);
    }
  }
 
  if(hitdebug)
    std::cout<<"good hits size = "<<good_hits.size()<<std::endl;
  
  return good_hits;
}

// ------------ method called for each event  ------------
void
MakeTrackValTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  run_nr_ = iEvent.id().run();
  lumi_nr_ = iEvent.id().luminosityBlock();
  evt_nr_ = iEvent.id().event();

  //reinitialize at each event
  for (unsigned int i = 0; i < MAXPART; i++){
    is_reco_matched_[i] = 0;
    gen_eta_[i] = -999;
    gen_phi_[i] = -999;
    gen_pt_[i] = -999;
    gen_pdgId_[i] = -999;
    gen_dxy_[i] = -10;
    gen_dz_[i] = -10;
    gen_nrSimHits_[i] = -10;
    gen_nrRecoHits_[i] = -10;
    gen_nrSharedHits_[i] = -10;
    
    gen_matched_eta_[i] = -999;
    gen_matched_pt_[i] = -999;
    gen_matched_phi_[i] = -999;

    is_gen_matched_[i] = 0;
    reco_eta_[i] = -999;
    reco_phi_[i] = -999;
    reco_pt_[i] = -999;
    reco_nrSimHits_[i] = -10;
    reco_nrRecoHits_[i] = -10;
    reco_nrSharedHits_[i] = -10;
  
    fake_eta_[i] = -999;
    fake_phi_[i] = -999;
    fake_pt_[i] = -999;

    for(unsigned int h = 0; h < MAXHIT; h++){
      gen_hit_pt_[i][h] = -999;
      gen_hit_eta_[i][h] = -999;
      gen_hit_phi_[i][h] = -999;
      gen_hit_layer_[i][h] = -10;
      gen_hit_subdetector_[i][h] = -10;
    }
  }

  edm::Handle<edm::View<reco::Track> > trackCollection; //reconstructed tracks
  //  iEvent.getByLabel("electronGsfTracks", trackCollection);
  iEvent.getByLabel("elGsfTracksWithQuality", trackCollection);

  edm::Handle<edm::View<reco::Track> > trackCollectionUsual;
  //  iEvent.getByLabel("generalTracks", trackCollectionUsual);
  iEvent.getByLabel("cutsRecoTracksHp", trackCollectionUsual);

  edm::Handle<edm::View<reco::GsfTrack> > trackCollectionFullElectron;
  iEvent.getByLabel("electronGsfTracks", trackCollectionFullElectron);
                     
    for( edm::View<reco::Track>::size_type i=0; i<trackCollectionUsual->size(); i++){
    edm::RefToBase<reco::Track> track(trackCollectionUsual, i);
    std::cout<<"general track with presel pT = "<<track->pt()<<std::endl;
  }

  for( edm::View<reco::Track>::size_type i=0; i<trackCollection->size(); i++){
    edm::RefToBase<reco::Track> track(trackCollection, i);
    std::cout<<"track with el_label with presel pT = "<<track->pt()<<std::endl;
  }

  for( edm::View<reco::GsfTrack>::size_type i=0; i<trackCollectionFullElectron->size(); i++){
    edm::RefToBase<reco::GsfTrack> track(trackCollectionFullElectron, i);
    std::cout<<"GSF track with no presel el_label pT = "<<track->pt()<<std::endl;
  }
  
  
  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ; //simulated tracks
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());
  if(debug_)
    std::cout <<"Number of simulated tracks = "<<tPCeff.size() << std::endl;

  TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP,minAbsEtaTP,maxAbsEtaTP);
       
  //  Handle<reco::SimToRecoCollection> simtorecocollection;
  // event.getByLabel

  reco::SimToRecoCollection simRecColl;
  reco::RecoToSimCollection recSimColl;

  //  edm::ESHandle<TrackAssociatorBase> theAssociator;
  //iSetup.get<TrackAssociatorRecord>().get("trackAssociatorByHits", theAssociator);

  edm::Handle<reco::SimToRecoCollection > simtorecoCollection;
  iEvent.getByLabel("trackingParticleRecoTrackAsssociation", simtorecoCollection);
  simRecColl = *(simtorecoCollection.product());

  edm::Handle<reco::RecoToSimCollection > rectosimCollection;
  iEvent.getByLabel("trackingParticleRecoTrackAsssociation",  rectosimCollection);
  recSimColl = *(rectosimCollection.product());

  //edm::ESHandle<TrackAssociatorBase> theAssociator; //create track associators for MC truth matching
  //iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits" ,theAssociator);
  
  
  //reco::SimToRecoCollection simRecColl = theAssociator->associateSimToReco(trackCollection, TPCollectionHeff, &iEvent);
  //reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(trackCollection, TPCollectionHeff, &iEvent);


  edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
  iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP", parametersDefinerTP);   

  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);


  //----------------------Loop over sim tracks--------------------------------
  np_gen_ = 0;
  np_gen_toReco_ = 0;


  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

    if( !tpSelector(*tp) ) continue;
    if( tp->pt() < 2 ) continue; // reduce noise. such tracks are irrelevant for electrons - never makes it to ECAL

    ParticleBase::Point vertex = parametersDefinerTP->vertex(iEvent,iSetup,*tp);
    gen_dxy_[np_gen_] = -vertex.x()*sin(tp->phi())+vertex.y()*cos(tp->phi()); 
    gen_dz_[np_gen_] = vertex.z()-(vertex.x()*tp->px()+vertex.y()*tp->py())/tp->pt()*tp->pz()/tp->pt();

    gen_eta_[np_gen_] = tp->eta();
    gen_phi_[np_gen_] = tp->phi();
    gen_pt_[np_gen_] = tp->pt();
    gen_pdgId_[np_gen_] = tp->pdgId();

    bool hitdebug = false;


    std::vector<PSimHit> simhits=tp->trackPSimHit(DetId::Tracker);
    double nr_simhits = simhits.size();

    if(hitdebug){
      std::cout<<"TRACKING PARTICLE with pt = "<<tp->pt()<<std::endl;
      std::cout<<"Sim hits size = "<<simhits.size()<<std::endl;
    }

    gen_nrSimHits_[np_gen_] =  nr_simhits;
    std::vector<PSimHit> goodhits=getGoodHits(simhits, hitdebug); //omit subsequent hits in the same subdetector layer  

    unsigned int hit_count = 0;
    for(std::vector<PSimHit>::const_iterator it_hit = goodhits.begin(); it_hit != goodhits.end(); it_hit++){

      GlobalVector global_p = getHitMomentum(it_hit,tracker, hitdebug);      
      float hit_eta_global = global_p.eta();
      float hit_phi_global = global_p.phi();
      float pt_at_entry = global_p.perp();
      float hit_pt_eff = pt_at_entry/tp->pt();

      DetId dId = DetId(it_hit->detUnitId() ); 
      std::vector<int> hit_position= getHitPosition(dId, hitdebug);
      int subdetector = hit_position.at(0);
      int layerNumber = hit_position.at(1);

      //-----------fill tree entries for each hit----------
      gen_hit_pt_[np_gen_][hit_count] = pt_at_entry;
      gen_hit_eta_[np_gen_][hit_count] = hit_eta_global;
      gen_hit_phi_[np_gen_][hit_count] = hit_phi_global;

      gen_hit_layer_[np_gen_][hit_count] = layerNumber;
      gen_hit_subdetector_[np_gen_][hit_count] = subdetector;        
      //---------------------------------------------------
	   
      if(hitdebug)
	std::cout<<"hit efficiency at hit "<<hit_count<<" = "<<hit_pt_eff<<std::endl;
      
      hit_count++;
    }

    if(hitdebug)
      std::cout<<"Nr hits counted = "<<hit_count<<std::endl;

    //--------check for matched reco tracks-----------

    const reco::Track* matchedTrackPointer=0;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
    int nSharedHits(0), nRecoTrackHits(0);

    if(simRecColl.find(tpr) != simRecColl.end()){
      rt = simRecColl[tpr]; // find sim-to-reco association
      if ( rt.size()!=0 ) {	
	matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track                                                 
	nSharedHits = rt.begin()->second;
	nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
	
	/*	if( (nRecoTrackHits < nSharedHits) || ( (int)simhits.size() < nSharedHits) ){
	  std::cout<<"SIMTORECO NR_RECO < NR_SHARED"<<"nr reco track hits = "<<nRecoTrackHits<<" and nr shared hits = "<<nSharedHits<<", nr sim track hits = "<<simhits.size()<<std::endl;
	}
	*/
	if(debug_){
	  std::cout<<"Found a simtrack with a matched reco track with:"<<"nr reco track hits = "<<nRecoTrackHits<<" quality = "<<nSharedHits<<", nr sim track hits = "<<simhits.size()<<", eta = "<<tp->eta()<<", phi = "<< tp->phi()<<std::endl;
	}
	is_reco_matched_[np_gen_] = 1;
	gen_nrSharedHits_[np_gen_] = nSharedHits;
	gen_nrRecoHits_[np_gen_] = nRecoTrackHits; 
	
	gen_matched_pt_[np_gen_] = tp->pt();
	gen_matched_eta_[np_gen_] = tp->eta();
	gen_matched_phi_[np_gen_] = tp->phi();

        np_gen_toReco_++; //count the number of simTracks that have a recoTrack associated
      }
    }

    

    
    np_gen_++; //count selected sim tracks passing the selection (important to keep in the end of the loop for tree items)      
  }
  //--------------------Loop over reco tracks-----------------------------------
  np_reco_ = 0;
  np_reco_toGen_ = 0;
  np_fake_ = 0;

  for( edm::View<reco::Track>::size_type i=0; i<trackCollection->size(); i++){
    edm::RefToBase<reco::Track> track(trackCollection, i);

    reco_pt_[np_reco_] = track->pt();
    reco_phi_[np_reco_] = track->phi();
    reco_eta_[np_reco_] = track->eta();
    reco_nrRecoHits_[np_reco_] = track->numberOfValidHits();

    std::vector<std::pair<TrackingParticleRef, double> > tp;

    if(recSimColl.find(track) != recSimColl.end()){ // if matched
      tp = recSimColl[track];
      if(tp.size() != 0) {
	is_gen_matched_[np_reco_] = 1;
	np_reco_toGen_++;
	TrackingParticleRef matchedTrackingParticle = tp.begin()->first;
	reco_nrSharedHits_[np_reco_] = tp.begin()->second;
	reco_nrSimHits_[np_reco_] = (matchedTrackingParticle->trackPSimHit(DetId::Tracker) ).size();
      }
    }
	//	  if( (track->numberOfValidHits() < tp.begin()->second) || ( reco_nrSimHits_[np_reco_] < tp.begin()->second) ){
	//    std::cout<<"RECOTOSIM NR_RECO < NR_SHARED"<<"nr reco track hits = "<<track->numberOfValidHits()<<" and nr shared hits = "<<tp.begin()->second<<", nr sim track hits = "<<reco_nrSimHits_[np_reco_]<<std::endl;
	//  }
	  
	
	//  std::cout<<"Found a reco track with matched sim track with:"<<"nr sim track hits = "<<reco_nrSimHits_[np_reco_]<<" and nr shared hits = "<<tp.begin()->second<<", nr reco track hits = "<<track->numberOfValidHits()<<", eta = "<<track->eta()<<", phi = "<<track->phi()<<std::endl;
	
    else{  // if fake
      fake_eta_[np_fake_] = track->eta();
      fake_pt_[np_fake_] = track->pt();
      fake_phi_[np_fake_] = track->phi();
	
      np_fake_++;
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

  trackValTree_->Branch("run_nr", &run_nr_, "run_nr/I");
  trackValTree_->Branch("evt_nr", &evt_nr_, "evt_nr/I");
  trackValTree_->Branch("lumi_nr", &lumi_nr_, "lumi_nr/I");
 
  trackValTree_->Branch("np_reco", &np_reco_, "np_reco/I");
  trackValTree_->Branch("np_reco_toGen", &np_reco_toGen_, "np_reco_toGen/I");
  trackValTree_->Branch("np_fake", &np_fake_, "np_fake/I");

  trackValTree_->Branch("reco_gen_matched", is_gen_matched_, "reco_gen_matched[np_reco]/I");
  trackValTree_->Branch("reco_pt", reco_pt_, "reco_pt[np_reco]/D");
  trackValTree_->Branch("reco_eta", reco_eta_, "reco_eta[np_reco]/D");
  trackValTree_->Branch("reco_phi", reco_phi_, "reco_phi[np_reco]/D");
  trackValTree_->Branch("reco_nrSimHits", reco_nrSimHits_, "reco_nrSimHits[np_reco]/I");
  trackValTree_->Branch("reco_nrRecoHits", reco_nrRecoHits_, "reco_nrRecoHits[np_reco]/I");
  trackValTree_->Branch("reco_nrSharedHits", reco_nrSharedHits_, "reco_nrSharedHits[np_reco]/D");

  trackValTree_->Branch("fake_pt", fake_pt_, "fake_pt[np_fake]/D");
  trackValTree_->Branch("fake_eta", fake_eta_, "fake_eta[np_fake]/D");
  trackValTree_->Branch("fake_phi", fake_phi_, "fake_phi[np_fake]/D");

  //-----------------------------------------------------------------------

  trackValTree_->Branch("np_gen",&np_gen_,"np_gen/I"); // needs to be filled in order to fill x[np]
  trackValTree_->Branch("np_gen_toReco", &np_gen_toReco_, "np_gen_ToReco/I");

  trackValTree_->Branch("gen_reco_matched", is_reco_matched_, "gen_reco_matched[np_gen]/I");
  trackValTree_->Branch("gen_pdgId", gen_pdgId_, "gen_pdgId[np_gen]/I");
  trackValTree_->Branch("gen_pt", gen_pt_, "gen_pt[np_gen]/D");
  trackValTree_->Branch("gen_eta", gen_eta_, "gen_eta[np_gen]/D");
  trackValTree_->Branch("gen_phi", gen_phi_, "gen_phi[np_gen]/D");
  trackValTree_->Branch("gen_dxy", gen_dxy_, "gen_dxy[np_gen]/D");
  trackValTree_->Branch("gen_dz", gen_dz_, "gen_dz[np_gen]/D");
  trackValTree_->Branch("gen_nrSimHits", gen_nrSimHits_, "gen_nrSimHits[np_gen]/I");
  trackValTree_->Branch("gen_nrRecoHits", gen_nrRecoHits_, "gen_nrRecoHits[np_gen]/I");
  trackValTree_->Branch("gen_nrSharedHits", gen_nrSharedHits_, "gen_nrSharedHits[np_gen]/D");  

  trackValTree_->Branch("gen_matched_pt", gen_matched_pt_, "gen_matched_pt[np_gen]/D");
  trackValTree_->Branch("gen_matched_eta", gen_matched_eta_, "gen_matched_eta[np_gen]/D");
  trackValTree_->Branch("gen_matched_phi", gen_matched_phi_, "gen_matched_phi[np_gen]/D");

  trackValTree_->Branch("gen_hit_pt", gen_hit_pt_, "gen_hit_pt[np_gen][50]/F");
  trackValTree_->Branch("gen_hit_eta", gen_hit_eta_, "gen_hit_eta[np_gen][50]/F");
  trackValTree_->Branch("gen_hit_phi", gen_hit_phi_, "gen_hit_phi[np_gen][50]/F");
  trackValTree_->Branch("gen_hit_layer", gen_hit_layer_, "gen_hit_layer[np_gen][50]/I");
  trackValTree_->Branch("gen_hit_subdetector", gen_hit_subdetector_, "gen_hit_subdetector[np_gen][50]/I");
}

//--------------------------functions----------------------------------
/*int MakeTrackValTree::hitPosition(const DetId& detId) const
{
  std::cout<<"function works"<<std::endl;
  return 1;
  }*/

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



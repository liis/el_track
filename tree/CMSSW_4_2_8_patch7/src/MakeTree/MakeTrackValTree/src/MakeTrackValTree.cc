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
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

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
      double getHitDistance(std::vector<PSimHit>::const_iterator, edm::ESHandle<TrackerGeometry>, bool);
      std::vector<PSimHit> getGoodHits(std::vector<PSimHit>, edm::ESHandle<TrackerGeometry>, bool);
      template<typename iter> std::vector<SimHitIdpr> getMatchedTPIds( iter, iter, const edm::Event&, const edm::ParameterSet&, bool);
      const TrackingRecHit* getHitPtr(edm::OwnVector<TrackingRecHit>::const_iterator);
      const TrackingRecHit* getHitPtr(trackingRecHit_iterator);
      int getNrSharedHits(std::vector<SimHitIdpr>&, TrackingParticle*);
    
  //---failed try to implement object sort-----
  //bool cmp(const edm::View<reco::ElectronSeed>::const_iterator, const edm::View<reco::ElectronSeed>::const_iterator);
  //     bool test_cmp(int, int);

  // ----------member data ---------------------------
  double ptMinTP, minRapidityTP, maxRapidityTP, tipTP, lipTP, signalOnlyTP, chargedOnlyTP, stableOnlyTP, minAbsEtaTP, maxAbsEtaTP;
  std::vector<int> pdgIdTP;
  int minHitTP;

  TTree *trackValTree_;

  int np_gen_, np_gen_toReco_, np_reco_, np_reco_toGen_, np_fake_, run_nr_, evt_nr_, lumi_nr_;
  int is_reco_matched_[MAXPART], is_gen_matched_[MAXPART], is_ecalDrivenSeed_[MAXPART], is_trackerDrivenSeed_[MAXPART];
  
  float gen_hit_pt_[MAXPART][MAXHIT], gen_hit_eta_[MAXPART][MAXHIT], gen_hit_phi_[MAXPART][MAXHIT];
  int gen_hit_subdetector_[MAXPART][MAXHIT], gen_hit_layer_[MAXPART][MAXHIT], gen_matchedRecHit_subdetector_[MAXPART][MAXHIT], gen_matchedRecHit_layer_[MAXPART][MAXHIT], gen_badRecHit_subdetector_[MAXPART][MAXHIT], gen_badRecHit_layer_[MAXPART][MAXHIT], gen_missedRecHit_subdetector_[MAXPART][MAXHIT], gen_missedRecHit_layer_[MAXPART][MAXHIT];
  int rec_hit_isMatched_[MAXPART][MAXHIT], rec_hit_subdetector_[MAXPART][MAXHIT], rec_hit_layer_[MAXPART][MAXHIT];


  double gen_pt_[MAXPART], gen_eta_[MAXPART], gen_phi_[MAXPART], gen_matched_pt_[MAXPART], gen_matched_eta_[MAXPART], gen_matched_phi_[MAXPART], gen_dxy_[MAXPART], gen_dz_[MAXPART], gen_ptAtLast_[MAXPART], gen_bremFraction_[MAXPART], gen_matched_seed_quality_[MAXPART];
  double reco_pt_[MAXPART], reco_eta_[MAXPART], reco_phi_[MAXPART], fake_pt_[MAXPART], fake_eta_[MAXPART], fake_phi_[MAXPART];
  int gen_pdgId_[MAXPART], gen_nrSimHits_[MAXPART], gen_nrUniqueSimHits_[MAXPART], gen_nrRecoHits_[MAXPART], gen_nrMatchedRecHits_[MAXPART], gen_nrSpuriousRecHits_[MAXPART], gen_nrLostSimHits_[MAXPART], gen_nr_matched_seed_hits_[MAXPART];
  int reco_nrRecoHits_[MAXPART], reco_nrSimHits_[MAXPART];
  double reco_nrSharedHits_[MAXPART], gen_nrSharedHits_[MAXPART];
  
  bool debug_, is_gsf_, hitdebug_;
  edm::InputTag track_label_gsf_, track_label_, el_seed_label_;

  TrackerHitAssociator* hitAssociator;

  const edm::ParameterSet& conf_;

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
MakeTrackValTree::MakeTrackValTree(const edm::ParameterSet& iConfig):
  conf_(iConfig) //to be able to pass the full set of configuration parameters
{
   //now do what ever initialization is needed
  //--------------get cut thresholds for sim tracks----------------
  //  conf_ = iConfig;

  debug_ = iConfig.getParameter<bool>("debug");
  hitdebug_ = iConfig.getParameter<bool>("hitdebug");
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

double MakeTrackValTree::getHitDistance(std::vector<PSimHit>::const_iterator it_hit, edm::ESHandle<TrackerGeometry> tracker, bool hitdebug = false){
  DetId dId = DetId(it_hit->detUnitId() );
  Local3DPoint local_p = it_hit->localPosition();
  const GeomDetUnit* detunit = tracker->idToDetUnit(dId.rawId());
  GlobalPoint global_p = detunit->toGlobal(local_p);
  if(hitdebug)
    std::cout<<"Hit global position x = "<<global_p.x()<<", y = "<<global_p.y()<<", z = "<<global_p.z()<<", r = "<<global_p.mag()<<std::endl;
  
  double distance = static_cast<double>(global_p.mag());
  return distance;
}

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
    layerNumber = static_cast<int>(pxbid.layer());
    if(hitdebug)
      std::cout<<"Hit at pixel barrel layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }

  else if(subdetId == PixelSubdetector::PixelEndcap){
    PXFDetId pxfid(dId.rawId());
    layerNumber = static_cast<int>(pxfid.disk());
    if(hitdebug)
   std::cout<<"Hit at pixel endcap layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
      
  else if( subdetId == StripSubdetector::TIB){
    TIBDetId tibid(dId.rawId());
    layerNumber = static_cast<int>(tibid.layer());
    if(hitdebug)
      std::cout<<"Hit at TIB layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
    
  else if( subdetId == StripSubdetector::TOB){
    TOBDetId tobid(dId.rawId());
    layerNumber = static_cast<int>(tobid.layer());
    if(hitdebug)
      std::cout<<"Hit at TOB layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
  else if( subdetId == StripSubdetector::TID){
    TIDDetId tidid(dId.rawId());
    layerNumber = static_cast<int>(tidid.wheel());
    if(hitdebug)
      std::cout<<"Hit at TID layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }
  else if( subdetId == StripSubdetector::TEC ){
    TECDetId tecid(dId.rawId());
    layerNumber = static_cast<int>(tecid.wheel());
    if(hitdebug)
      std::cout<<"Hit at TEC layer = "<<layerNumber<<", corresponding to subdetId = "<<subdetId<<std::endl;
  }

  std::vector<int> hit_position;
  hit_position.push_back(subdetId);
  hit_position.push_back(static_cast<int>(layerNumber) );

  return hit_position;
}

std::vector<PSimHit>  MakeTrackValTree::getGoodHits(std::vector<PSimHit> simhits, edm::ESHandle<TrackerGeometry> tracker, bool hitdebug = false) 
//loop over all hits in the event and omit subsequent hit in the same subdetector layer 
{  
  std::vector<PSimHit> good_hits;
  for(std::vector<PSimHit>::const_iterator it_hit = simhits.begin(); it_hit != simhits.end(); it_hit++ ){
    DetId dId = DetId(it_hit->detUnitId() ); 
    std::vector<int> hitposition = getHitPosition(dId, hitdebug);
    double distance = getHitDistance(it_hit, tracker, hitdebug_);

    int hit_subdetector = hitposition.at(0);
    int hit_layer = hitposition.at(1);
      
    if( !good_hits.size() )
      good_hits.push_back(*it_hit);
    else{ //compare to the previous hit position
      DetId dId_good = DetId(good_hits.back().detUnitId() );

      int good_hit_subdetector = (getHitPosition(dId_good)).at(0);
      int good_hit_layer = (getHitPosition(dId_good)).at(1);
      std::vector<PSimHit>::const_iterator last_good_hit = good_hits.end();
      --last_good_hit;
      
      //---check whether the distance from the primary vertex increases------
      if( distance < getHitDistance(last_good_hit, tracker, hitdebug_) )
      	break;
      // ---- check for the same subdetector hits in the same layer--------
      if( hit_subdetector == good_hit_subdetector && hit_layer == good_hit_layer){
	good_hits.pop_back(); //replace the previous hit with a later one
	good_hits.push_back(*it_hit);
      }
      else
	good_hits.push_back(*it_hit);
      
    }
  } // end loop over simhits
  
  if(hitdebug)
    std::cout<<"good hits size = "<<good_hits.size()<<std::endl;
  
  return good_hits;
}

template<typename iter>
std::vector<SimHitIdpr> MakeTrackValTree::getMatchedTPIds( iter begin, iter end, const edm::Event& iEvent, const edm::ParameterSet& conf, bool debug = false)
// loop over rec hits between 'iter begin' and 'iter end' and return a vector of associated tracking particle IDs
{
  hitAssociator = new TrackerHitAssociator(iEvent, conf); // ?? not needed 

  std::vector<SimHitIdpr> matchedSimIdsTot;
  matchedSimIdsTot.clear(); // clean up the matched Ids for each seed

  int ri = 1; // rec hit count
  for( iter it = begin; it != end; it++, ri++){
    const TrackingRecHit *hit = getHitPtr(it);

    std::vector<SimHitIdpr> matchedSimIds = hitAssociator->associateHitId(*hit); //find the IDs of tracking particles associated to the rec hit

    for( unsigned int i = 0; i < matchedSimIds.size(); i++ ) // save all matched IDs of a seed to a single vector
      matchedSimIdsTot.push_back(matchedSimIds[i]);

    //////////debug/////////////
    if(debug){
      std::cout<<" Seed rec hit # " << ri << " valid=" << it->isValid() << " det id = " << it->geographicalId().rawId();
      if(matchedSimIds.size()){
	for( unsigned int i = 0; i < matchedSimIds.size(); i++)      
	  std::cout<< ", associated Sim Track Id " << matchedSimIds[i].first<<std::endl; 
      }
      else
	std::cout<<", No matched tracking particle found"<<std::endl;
    }
    /////////////////////////
  }

  return matchedSimIdsTot;
}

const TrackingRecHit* MakeTrackValTree::getHitPtr(edm::OwnVector<TrackingRecHit>::const_iterator iter) {return &*iter;}
const TrackingRecHit* MakeTrackValTree::getHitPtr(trackingRecHit_iterator iter) {return &**iter;}

int MakeTrackValTree::getNrSharedHits(std::vector<SimHitIdpr>& recHitMatchedTPIds, TrackingParticle* tp) {
  int nshared = 0;  
  
  for( unsigned int i = 0; i<recHitMatchedTPIds.size(); i++){   
    for( TrackingParticle::g4t_iterator g4T = tp->g4Track_begin(); g4T != tp->g4Track_end(); g4T++){
      if( (*g4T).trackId() == recHitMatchedTPIds[i].first && tp->eventId() == recHitMatchedTPIds[i].second )
	nshared++;
	  
    }
  }
  return nshared;
}

/* //failed try to implement object sort
bool MakeTrackValTree::cmp(const edm::View<reco::ElectronSeed>::const_iterator it_s1,  const edm::View<reco::ElectronSeed>::const_iterator it_s2){
  return it_s1->nHits() < it_s2->nHits() ;
}

bool MakeTrackValTree::test_cmp(int i, int j) {
  return (i < j);
}
*/
// ------------ method called for each event  ------------
void
MakeTrackValTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  run_nr_ = iEvent.id().run();
  lumi_nr_ = iEvent.id().luminosityBlock();
  evt_nr_ = iEvent.id().event();

  //reinitialize at each event
  for (unsigned int i = 0; i < MAXPART; i++){
    is_reco_matched_[i] = -10;
    is_ecalDrivenSeed_[i] = -10;
    is_trackerDrivenSeed_[i] = -10;
    gen_eta_[i] = -999;
    gen_phi_[i] = -999;
    gen_pt_[i] = -999;
    gen_ptAtLast_[i] = -999;
    gen_bremFraction_[i] = -10;
    gen_pdgId_[i] = -999;
    gen_dxy_[i] = -10;
    gen_dz_[i] = -10;
    gen_nrSimHits_[i] = -10;
    gen_nrMatchedRecHits_[i] = -10;
    gen_nrSpuriousRecHits_[i] = -10;
    gen_nrLostSimHits_[i] = -10;

    gen_nrUniqueSimHits_[i] = -10;
    gen_nrRecoHits_[i] = -10;
    gen_nrSharedHits_[i] = -10;
    
    gen_matched_eta_[i] = -999;
    gen_matched_pt_[i] = -999;
    gen_matched_phi_[i] = -999;

    gen_nr_matched_seed_hits_[i] = -10;
    gen_matched_seed_quality_[i] = -10;

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
      gen_matchedRecHit_subdetector_[i][h] = -10;
      gen_matchedRecHit_layer_[i][h] = -10;
      gen_missedRecHit_subdetector_[i][h] = -10;
      gen_missedRecHit_layer_[i][h] = -10;
      gen_badRecHit_subdetector_[i][h] = -10;
      gen_badRecHit_layer_[i][h] = -10;

      rec_hit_isMatched_[i][h] = -1;
      rec_hit_subdetector_[i][h] = -10;
      rec_hit_layer_[i][h] = -10;
    }
  }

  edm::InputTag track_label; 
  if(is_gsf_)
    track_label = track_label_gsf_;
  else
    track_label = track_label_;

  edm::Handle<edm::View<reco::Track> > trackCollection; //reconstructed tracks
  iEvent.getByLabel(track_label, trackCollection);
  //iEvent.getByLabel("elGsfTracksWithQuality", trackCollection); //produced collection in aod file -- find out the definition

  edm::Handle<edm::View<reco::ElectronSeed> > elSeedCollection; 
  iEvent.getByLabel(el_seed_label_, elSeedCollection);

  //edm::Handle<edm::View<TrajectorySeed> > seedCollection; // generic track seeds
  //iEvent.getByLabel("electronMergedSeeds", seedCollection);

  if(debug_)
    std::cout<<"Reco track label = "<<track_label<<", electron seed label = "<<el_seed_label_<<std::endl;

  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ; //simulated tracks
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionHeff);

  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());
  TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP,minAbsEtaTP,maxAbsEtaTP);
       
  //------------ get track associators ----------------------
  edm::ESHandle<TrackAssociatorBase> myAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", myAssociator);

  //  std::cout<<"Getting associateSimToReco directly: "<<std::endl;
  reco::SimToRecoCollection simRecColl = myAssociator->associateSimToReco(trackCollection, TPCollectionHeff, &iEvent);

  // std::cout<<"Getting associateRecoToSim directly: "<<std::endl;
  reco::RecoToSimCollection recSimColl = myAssociator->associateRecoToSim(trackCollection, TPCollectionHeff, &iEvent);

  hitAssociator = new TrackerHitAssociator(iEvent, conf_);

  //------------------------------------------------------------

  edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP; //?
  iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP", parametersDefinerTP);   

  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);

  //---------------------Loop over electron track seeds---------------
  /*
  std::cout<<"Looping over electron seed collection of size "<<elSeedCollection->size()<<std::endl;
   
  for( edm::View<reco::ElectronSeed>::const_iterator it_seed = elSeedCollection->begin(); it_seed != elSeedCollection->end(); it_seed++){

    std::cout<<"Electron seed with numbr of hits = "<<it_seed->nHits()<<std::endl;
    TrajectorySeed::range rechits = it_seed->recHits(); //get recHits associated to the seed

    std::vector<SimHitIdpr> matchedSimIdsTot = getMatchedTPIds<edm::OwnVector<TrackingRecHit>::const_iterator>(it_seed->recHits().first, it_seed->recHits().second, iEvent, conf_, false); //get a vector of TP IDs associated to the recHits of the seed


    std::cout<<"simIds.size = "<<matchedSimIdsTot.size()<<std::endl;
    if( hitdebug_){
      std::cout<<"Matched TP IDs for seed rec hits: ";
      for(unsigned int i = 0; i < matchedSimIdsTot.size(); i++)
	std::cout<<matchedSimIdsTot[i].first<<", ";
      std::cout<<std::endl;
    }
      
  }
  */
  
  //----------------------Loop over tracking particles--------------------------------
  np_gen_ = 0;
  np_gen_toReco_ = 0;

  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

    if( !tpSelector(*tp) ) continue;
    if( tp->pt() < 1 ) continue; // such tracks are irrelevant for electrons

    ParticleBase::Point vertex = parametersDefinerTP->vertex(iEvent,iSetup,*tp);
    gen_dxy_[np_gen_] = -vertex.x()*sin(tp->phi())+vertex.y()*cos(tp->phi()); 
    gen_dz_[np_gen_] = vertex.z()-(vertex.x()*tp->px()+vertex.y()*tp->py())/tp->pt()*tp->pz()/tp->pt();

    gen_eta_[np_gen_] = tp->eta();
    gen_phi_[np_gen_] = tp->phi();
    gen_pt_[np_gen_] = tp->pt();
    gen_pdgId_[np_gen_] = tp->pdgId();

    std::vector<PSimHit> simhits=tp->trackPSimHit(DetId::Tracker);
    int nr_simhits = simhits.size();
    std::vector<PSimHit> goodhits=getGoodHits(simhits, tracker, false); //hitdebug_); //omit subsequent hits in the same subdetector layer
    int nr_goodhits = goodhits.size();

    gen_nrSimHits_[np_gen_] =  nr_simhits;
    gen_nrUniqueSimHits_[np_gen_] = nr_goodhits;

    std::vector<unsigned int> tpids; //For each tracking particle save a vector of corresponding track ID-s
    for (TrackingParticle::g4t_iterator g4T=tp->g4Track_begin(); g4T!=tp->g4Track_end(); g4T++)
      tpids.push_back(g4T->trackId());

    if(debug_){
      std::cout<<"----------------------------------------"<<std::endl;
      std::cout<<"Found Tracking particle with pt = "<<tp->pt()<<", sim hits size = "<<simhits.size()<<", Nr unique hits = "<<nr_goodhits<<", track ID-s size = "<<tpids.size()<<std::endl;
    }
    //---------------------- check mathed seeds ---------------------
    if(debug_)
      std::cout<<"Looping over electron seed collection of size "<<elSeedCollection->size()<<std::endl;

    int sn = 1;
    std::vector<double> seedMatchQuality;
    seedMatchQuality.clear();
    std::vector<int> seedNrMatchedHits;
    seedNrMatchedHits.clear();
    std::vector<edm::View<reco::ElectronSeed>::const_iterator> bestSeed;
    bestSeed.clear();

    for( edm::View<reco::ElectronSeed>::const_iterator it_seed = elSeedCollection->begin(); it_seed != elSeedCollection->end(); it_seed++, sn++){
      
      if(debug_)
	std::cout<<"Electron seed # "<<sn<<" with numbr of hits = "<<it_seed->nHits()<<std::endl;

      TrajectorySeed::range rechits = it_seed->recHits(); //get recHits associated to the seed
      std::vector<SimHitIdpr> matchedSimIdsTot = getMatchedTPIds<edm::OwnVector<TrackingRecHit>::const_iterator>(it_seed->recHits().first, it_seed->recHits().second, iEvent, conf_, false); //get a vector of TP IDs associated to the recHits of the seed
      int nshared = getNrSharedHits(matchedSimIdsTot, tp);

      if(debug_){
	std::cout<<"simIds.size = "<<matchedSimIdsTot.size();
	if(it_seed->isEcalDriven())
	  std::cout<<", IsEcalDriven";
	if(it_seed->isTrackerDriven())
	  std::cout<<", IsTrackerDriven";
	std::cout<<std::endl;
	
	std::cout<<"TP shared hits with a seed = "<<nshared<<std::endl;
      }

      if( hitdebug_){
	std::cout<<"Matched TP IDs for seed rec hits: ";
	for(unsigned int i = 0; i < matchedSimIdsTot.size(); i++)
	  std::cout<<matchedSimIdsTot[i].first<<", ";
	std::cout<<std::endl;
      }
	
      double match_quality = (double)nshared/(double)(it_seed->nHits());

      if(seedMatchQuality.size() == 0)
	bestSeed.push_back(it_seed); //fill the first seed to the vector
      else{
	if( match_quality > seedMatchQuality.at(0) ){ // compare to the best quality existing seed, if better replace
	  bestSeed.pop_back(); 
	  bestSeed.push_back(it_seed);
	}
      }

      seedNrMatchedHits.push_back(nshared);
      if(nshared > 1)
	seedMatchQuality.push_back( match_quality );
      else
	seedMatchQuality.push_back(0);

      seedMatchQuality.push_back( match_quality );
      seedNrMatchedHits.push_back(nshared);

      std::sort(seedMatchQuality.begin(), seedMatchQuality.end(), std::greater<double>() );
      std::sort(seedNrMatchedHits.begin(), seedNrMatchedHits.end(), std::greater<double>() );
    }
      

    if( bestSeed.size() ){ //save ecalDriven/tracker driven information for the best matched seed
      if( bestSeed.at(0)->isEcalDriven() )
	is_ecalDrivenSeed_[np_gen_] = 1;
      else
	is_ecalDrivenSeed_[np_gen_] = 0;
    
      if(bestSeed.at(0)->isTrackerDriven())
	is_trackerDrivenSeed_[np_gen_] = 1;
      else
	is_trackerDrivenSeed_[np_gen_] = 0;
    }
    
    if(seedMatchQuality.size() && seedNrMatchedHits.size()){
      gen_nr_matched_seed_hits_[np_gen_] = seedNrMatchedHits.at(0);
      gen_matched_seed_quality_[np_gen_] = seedMatchQuality.at(0);
      if(debug_){
	std::cout<<"Match quality with seed = "<< seedMatchQuality.at(0)<<std::endl;
	std::cout<<"Nr matched seed hits = "<< seedNrMatchedHits.at(0)<<std::endl;
      }
    }
    
    
    //--------------------------------------------------------------------

    int hit_count = 0;
    if(nr_goodhits > 13 && abs(tp->eta()) < 0.9)
      std::cout<<"Found track with nr hits = "<<nr_goodhits<<", pt = "<<tp->pt()<<", eta = "<<tp->eta()<<", ecalDriven = "<<is_ecalDrivenSeed_[np_gen_]<<", trackerDriven = "<<is_trackerDrivenSeed_[np_gen_]<<std::endl;

    for(std::vector<PSimHit>::const_iterator it_hit = goodhits.begin(); it_hit != goodhits.end(); it_hit++){

      GlobalVector global_p = getHitMomentum(it_hit,tracker, hitdebug_);      
      float hit_eta_global = global_p.eta();
      float hit_phi_global = global_p.phi();
      float pt_at_entry = global_p.perp();
      float hit_pt_eff = pt_at_entry/tp->pt();
      DetId dId = DetId(it_hit->detUnitId() ); 
      std::vector<int> hit_position= getHitPosition(dId, hitdebug_); // (subdetector, layer)
        
      int subdetector = hit_position.at(0);
      int layerNumber = hit_position.at(1);


      if(nr_goodhits > 13 && abs(tp->eta()) < 0.9){
	std::cout<<" subdetector = "<<hit_position.at(0)<<std::endl;
	std::cout<<"layernumber = "<<hit_position.at(1)<<std::endl;
	std::cout<<" distance = "<< getHitDistance(it_hit, tracker, hitdebug_)<<std::endl;
      }

      //-----------fill tree entries for each hit----------
      gen_hit_pt_[np_gen_][hit_count] = pt_at_entry;
      gen_hit_eta_[np_gen_][hit_count] = hit_eta_global;
      gen_hit_phi_[np_gen_][hit_count] = hit_phi_global;

      gen_hit_layer_[np_gen_][hit_count] = layerNumber;
      gen_hit_subdetector_[np_gen_][hit_count] = subdetector;        
      //---------------------------------------------------
      if(hit_count+1 == nr_goodhits){
	gen_ptAtLast_[np_gen_] = pt_at_entry;
	gen_bremFraction_[np_gen_] = 1 - pt_at_entry/tp->pt(); //fraction of pT lost via brehmstrahlung
	//	if(debug_)
	if(nr_goodhits > 15 && tp->eta() < 0.9)
	  std::cout<<"pt at last hit = "<<pt_at_entry<<", brem fraction = "<<gen_bremFraction_[np_gen_]<<std::endl;
      }

      if(hitdebug_)
	std::cout<<"hit efficiency at hit "<<hit_count<<" = "<<hit_pt_eff<<std::endl;
      
      hit_count++;
    }

    //-------------manually look for a best matched reco track-----------------------
    std::vector<SimHitIdpr> matchedSimIds;
    std::vector<edm::RefToBase<reco::Track> > bestMatchRecoTrack;        
 
    std::vector<const TrackingRecHit*> matchedRecHits_best;
    std::vector<const TrackingRecHit*> badRecHits_best;

    for( edm::View<reco::Track>::size_type i=0; i != trackCollection->size(); i++){ // loop over reco tracks and check for matched hits
      edm::RefToBase<reco::Track> track(trackCollection, i);      
      int nr_invalid_hits = 0;

      std::vector<const TrackingRecHit*> matchedRecHits;
      std::vector<const TrackingRecHit*> badRecHits;

      for(trackingRecHit_iterator rechit_it = track->recHitsBegin(); rechit_it != track->recHitsEnd(); rechit_it++){ //loop over rec hits	
	matchedSimIds = hitAssociator->associateHitId(**rechit_it); // find simulated tracks causing the reconstructed hit
	
	const TrackingRecHit* rechit = &**rechit_it;	

	bool goodhit = false;
	for( unsigned int iSimId = 0; iSimId < matchedSimIds.size(); iSimId++){
	  for( unsigned int iTpId = 0; iTpId< tpids.size(); iTpId++){
	    if(matchedSimIds[iSimId].first == tpids[iTpId]){ 
	      goodhit = true;
	      break;
	    }
	  }	    
	}

	if(matchedSimIds.size() == 0)
	  nr_invalid_hits++;
	else if(goodhit)
	  matchedRecHits.push_back(rechit);
	else
	  badRecHits.push_back(rechit);
      } // end for rechit


      if( bestMatchRecoTrack.size() == 0){ //choose the reco track with best match to tracking particle
	//	std::cout<<"Saving best track with pt = "<<track->pt()<<", nr_good_hits = "<<matchedRecHits.size()<<std::endl;
        bestMatchRecoTrack.push_back(track);
	
	matchedRecHits_best = matchedRecHits;
	badRecHits_best = badRecHits;
      }
      else if( matchedRecHits.size() > matchedRecHits_best.size() ){
	//	std::cout<<"Overwrite with next track with nr_good_hits = "<<matchedRecHits.size()<<std::endl;
        bestMatchRecoTrack.pop_back();
        bestMatchRecoTrack.push_back( track );

	matchedRecHits_best = matchedRecHits;
	badRecHits_best = badRecHits;
      }
    } // end of reco tracks
    gen_nrMatchedRecHits_[np_gen_] = matchedRecHits_best.size();
    gen_nrSpuriousRecHits_[np_gen_] = badRecHits_best.size();

    for(std::vector<const TrackingRecHit*>::size_type ihit=0; ihit !=matchedRecHits_best.size(); ihit++){//write mached rec hit parameters to the tree
      DetId dId_rec = matchedRecHits_best.at(ihit)->geographicalId();

      std::vector<int> matchedhit_pos = getHitPosition(dId_rec, false);
      if(matchedhit_pos.size()){
	gen_matchedRecHit_subdetector_[np_gen_][ihit] = matchedhit_pos.at(0);
	gen_matchedRecHit_layer_[np_gen_][ihit] = matchedhit_pos.at(1);
      }
      else
	std::cout<<"Failed to get hit position information."<<std::endl;
    }

    for(std::vector<const TrackingRecHit*>::size_type ihit=0; ihit !=badRecHits_best.size(); ihit++){//write spurious rec hit parameters to the tree
      DetId dId_rec = badRecHits_best.at(ihit)->geographicalId();
      std::vector<int> badhit_pos = getHitPosition(dId_rec, false);

      if(badhit_pos.size()){
        gen_badRecHit_subdetector_[np_gen_][ihit] = badhit_pos.at(0);
        gen_badRecHit_layer_[np_gen_][ihit] = badhit_pos.at(1);
      } 
      else
	std::cout<<"Failed to get hit position information."<<std::endl;
    }

    if(nr_goodhits > 15 && abs(tp->eta()) < 0.9)
      std::cout<<"nr matched rec hits = "<<matchedRecHits_best.size()<<std::endl;
  
    if(debug_){
      if(bestMatchRecoTrack.size() )
	std::cout<<"Best matched reco track with nr good hits = "<<matchedRecHits_best.size()<<", nr bad hits = "<<badRecHits_best.size()<<", pt = "<<bestMatchRecoTrack.at(0)->pt()<<std::endl;
      else
	std::cout<<"No matching reco-track found"<<std::endl;
    }

    //---------------------Loop over simHits and check for missed hits-------------------------------- 
    int matched_hit_count = 0;
    int missed_hit_count = 0;
    int i_simhit = 0;
    for(std::vector<PSimHit>::const_iterator it_simhit = simhits.begin(); it_simhit != simhits.end(); it_simhit++, i_simhit++){
      DetId dId_sim = DetId(it_simhit->detUnitId() );
      bool matched = false;
      for(std::vector<const TrackingRecHit*>::size_type ihit=0; ihit !=matchedRecHits_best.size(); ihit++){
	DetId dId_reco = DetId(matchedRecHits_best.at(ihit)->geographicalId());
	if(dId_sim == dId_reco)
	  matched = true;
      }      
      if(matched == true)
	matched_hit_count++;
      else{
	missed_hit_count++;
	gen_missedRecHit_subdetector_[np_gen_][i_simhit] = getHitPosition(dId_sim, false).at(0);
	gen_missedRecHit_layer_[np_gen_][i_simhit] = getHitPosition(dId_sim, false).at(1);
      }
    }
    if(debug_)
      std::cout<<"Nr matched sim hits = "<<matched_hit_count<<", Nr missed sim hits = "<<missed_hit_count<<std::endl;

        
    //--------check for matched reco track defined by AssociatorByHits (efficiency and fake-rate plots)-----------
    const reco::Track* matchedTrackPointer=0;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
    int nSharedHits(0), nRecoTrackHits(0);

    if(simRecColl.find(tpr) != simRecColl.end()){
      rt = simRecColl[tpr]; // find sim-to-reco association
      if ( rt.size()!=0 ) {	
	matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track                                                 
	nSharedHits = rt.begin()->second;
	nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
	
	if(debug_)
	  std::cout<<"AssociatorByHits:associated reco track with:"<<"nr reco hits = "<<nRecoTrackHits<<" quality = "<<nSharedHits<<", pt = "<<matchedTrackPointer->pt()<<", eta = "<<matchedTrackPointer->eta()<<", phi = "<< matchedTrackPointer->phi()<<std::endl;
	
	is_reco_matched_[np_gen_] = 1;
	gen_nrSharedHits_[np_gen_] = nSharedHits;
	gen_nrRecoHits_[np_gen_] = nRecoTrackHits; 
	
	gen_matched_pt_[np_gen_] = tp->pt();
	gen_matched_eta_[np_gen_] = tp->eta();
	gen_matched_phi_[np_gen_] = tp->phi();

        np_gen_toReco_++; //count the number of simTracks that have a recoTrack associated
      }
    }
    else
      is_reco_matched_[np_gen_] = 0;

    np_gen_++; //count selected sim tracks passing the selection (important to keep in the end of the loop for tree items)      
  }
 
  
  //--------------------Loop over reco tracks-----------------------------------
  np_reco_ = 0;
  np_reco_toGen_ = 0;
  np_fake_ = 0;

  for( edm::View<reco::Track>::size_type i=0; i<trackCollection->size(); i++){
    edm::RefToBase<reco::Track> track(trackCollection, i);
    
    if(track->pt() < 1 ) continue; // reduce noise, irrelevant for electrons

    reco_pt_[np_reco_] = track->pt();
    reco_phi_[np_reco_] = track->phi();
    reco_eta_[np_reco_] = track->eta();
    reco_nrRecoHits_[np_reco_] = track->numberOfValidHits();

    //--------------check rec-hits one-by one--------------------------    

    hitAssociator = new TrackerHitAssociator(iEvent, conf_);
                                                                                       
    for(trackingRecHit_iterator rechit_it = track->recHitsBegin(); rechit_it != track->recHitsEnd(); rechit_it++){
      const TrackingRecHit* rechit = &**rechit_it; // def from /interface/TrackAssociatorByHits.h L97

      if(rechit->isValid()){
	DetId dId_rec = rechit->geographicalId();
	std::vector<int> rechit_pos = getHitPosition(dId_rec, false);

	std::vector<SimHitIdpr> matchedSimIds = hitAssociator->associateHitId(**rechit_it);
	
	for(unsigned int j=0; j < matchedSimIds.size(); j++){                                                                                                                                    
          //for(unsigned int 

	  //	  if(matchedSimIds[j].first ==

	  //std::cout<<"matched id = "<<matchedSimIds[j].first<<std::endl;                                                                                                                               
	}  

	//	std::vector<PSimHit> matchedSimHits = hitAssociator->associateHit(**rechit_it);

	//	for( unsigned int isim = 0; isim<matchedSimHits.size(); isim++){
	//	  unsigned tkId = matchedSimHits[isim].trackId();

	// std::cout<<"track id = "<<tkId<<std::endl;

      

	//	for(unsigned int j=0; j < matchedSimIds.size(); j++){
	  //std::cout<<"matched id = "<<matchedSimIds[j].first<<std::endl;
	//}

      }
    }
    
    
    //---------------------check association maps to sim-tracks---------------------
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
    //--------------------------------------------------------------------------------
 
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
  trackValTree_->Branch("is_ecalDrivenSeed", is_ecalDrivenSeed_, "is_ecalDrivenSeed[np_gen]/I");
  trackValTree_->Branch("is_trackerDrivenSeed", is_trackerDrivenSeed_, "is_trackerDrivenSeed[np_gen]/I");
  trackValTree_->Branch("gen_pdgId", gen_pdgId_, "gen_pdgId[np_gen]/I");
  trackValTree_->Branch("gen_pt", gen_pt_, "gen_pt[np_gen]/D");
  trackValTree_->Branch("gen_ptAtLast", gen_ptAtLast_, "gen_ptAtLast[np_gen]/D");
  trackValTree_->Branch("gen_bremFraction", gen_bremFraction_, "gen_bremFraction[np_gen]/D");
  trackValTree_->Branch("gen_eta", gen_eta_, "gen_eta[np_gen]/D");
  trackValTree_->Branch("gen_phi", gen_phi_, "gen_phi[np_gen]/D");
  trackValTree_->Branch("gen_dxy", gen_dxy_, "gen_dxy[np_gen]/D");
  trackValTree_->Branch("gen_dz", gen_dz_, "gen_dz[np_gen]/D");
  trackValTree_->Branch("gen_nrSimHits", gen_nrSimHits_, "gen_nrSimHits[np_gen]/I");
  trackValTree_->Branch("gen_nrUniqueSimHits", gen_nrUniqueSimHits_, "gen_nrUniqueSimHits[np_gen]/I");
  trackValTree_->Branch("gen_nrMatchedRecHits", gen_nrMatchedRecHits_, "gen_nrMatchedRecHits[np_gen]/I");
  trackValTree_->Branch("gen_nrSpuriousRecHits", gen_nrSpuriousRecHits_, "gen_nrSpuriousRecHits[np_gen]/I");
  trackValTree_->Branch("gen_nrLostSimHits", gen_nrLostSimHits_, "gen_nrLostSimHits[np_gen]/I");

  trackValTree_->Branch("gen_nrRecoHits", gen_nrRecoHits_, "gen_nrRecoHits[np_gen]/I");
  trackValTree_->Branch("gen_nrSharedHits", gen_nrSharedHits_, "gen_nrSharedHits[np_gen]/D");  
  trackValTree_->Branch("gen_nrMatchedSeedHits", gen_nr_matched_seed_hits_, "gen_nrMatchedSeedHits[np_gen]/I");
  trackValTree_->Branch("gen_matchedSeedQuality", gen_matched_seed_quality_, "gen_matchedSeedQuality[np_gen]/D");

  trackValTree_->Branch("gen_matched_pt", gen_matched_pt_, "gen_matched_pt[np_gen]/D");
  trackValTree_->Branch("gen_matched_eta", gen_matched_eta_, "gen_matched_eta[np_gen]/D");
  trackValTree_->Branch("gen_matched_phi", gen_matched_phi_, "gen_matched_phi[np_gen]/D");

  trackValTree_->Branch("gen_hit_pt", gen_hit_pt_, "gen_hit_pt[np_gen][50]/F");
  trackValTree_->Branch("gen_hit_eta", gen_hit_eta_, "gen_hit_eta[np_gen][50]/F");
  trackValTree_->Branch("gen_hit_phi", gen_hit_phi_, "gen_hit_phi[np_gen][50]/F");
  trackValTree_->Branch("gen_hit_layer", gen_hit_layer_, "gen_hit_layer[np_gen][50]/I");
  trackValTree_->Branch("gen_hit_subdetector", gen_hit_subdetector_, "gen_hit_subdetector[np_gen][50]/I");
  trackValTree_->Branch("gen_matchedRecHit_layer", gen_matchedRecHit_layer_, "gen_matchedRecHit_layer[np_gen][50]/I");
  trackValTree_->Branch("gen_matchedRecHit_subdetector", gen_matchedRecHit_subdetector_, "gen_matchedRecHit_subdetector[np_gen][50]/I");
  trackValTree_->Branch("gen_missedRecHit_layer", gen_missedRecHit_layer_, "gen_missedRecHit_layer[np_gen][50]/I");
  trackValTree_->Branch("gen_missedRecHit_subdetector", gen_missedRecHit_subdetector_, "gen_missedRecHit_subdetector[np_gen][50]/I");
  trackValTree_->Branch("gen_badRecHit_layer", gen_badRecHit_layer_, "gen_badRecHit_layer[np_gen][50]/I");
  trackValTree_->Branch("gen_badRecHit_subdetector", gen_badRecHit_subdetector_, "gen_badRecHit_subdetector[np_gen][50]/I");
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



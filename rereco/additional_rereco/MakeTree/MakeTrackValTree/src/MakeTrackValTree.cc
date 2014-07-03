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
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

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
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimGeneral/TrackingAnalysis/interface/SimHitTPAssociationProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
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
      const TrackingRecHit* getHitPtr(edm::OwnVector<TrackingRecHit>::const_iterator);
      const TrackingRecHit* getHitPtr(trackingRecHit_iterator);

      template<typename iter> bool isMatchedRecHitById( iter, TrackingParticle* );
      std::pair<std::vector<const reco::ElectronSeed*>, double> findMatchedSeed(edm::Handle<edm::View<reco::ElectronSeed>>, TrackingParticle*); 
      std::map<int, std::vector<PSimHit> > getSimHits(std::vector<std::string>, const edm::Event& );
      std::vector<int> getHitPosition(DetId, bool);
  //  bool comp_simHits(const PSimHit, const PSimHit);
  bool comparator(const std::pair<int,int>&, const std::pair<int,int>& );

  // ------------------ thresholds for TPselector ----------
  double ptMinTP, minRapidityTP, maxRapidityTP, tipTP, lipTP, signalOnlyTP, chargedOnlyTP, stableOnlyTP; 
  std::vector<int> pdgIdTP;
  int minHitTP;
  // -------------------- for tree -------------
  TTree *trackValTree_;
  int run_nr_, evt_nr_, lumi_nr_;
  int np_gen_, np_gen_toReco_, np_reco_, np_reco_toGen_,np_fake_;
  int is_reco_matched_[MAXPART], is_gen_matched_[MAXPART];
  int gen_pdgId_[MAXPART], gen_nr_simhits_[MAXPART], gen_matched_seed_nshared_[MAXPART], gen_matched_seed_okCharge_[MAXPART], is_ecalDrivenSeed_[MAXPART], is_trackerDrivenSeed_[MAXPART];

  double reco_pt_[MAXPART], reco_eta_[MAXPART], reco_phi_[MAXPART], fake_pt_[MAXPART], fake_eta_[MAXPART], fake_phi_[MAXPART];

  double gen_pt_[MAXPART],  gen_eta_[MAXPART], gen_phi_[MAXPART], gen_matched_pt_[MAXPART], gen_matched_eta_[MAXPART], gen_matched_phi_[MAXPART], gen_matched_qoverp_[MAXPART], gen_matched_cotth_[MAXPART], gen_matched_theta_[MAXPART], gen_matched_seed_quality_[MAXPART];

  double gen_matched_rec_eta_[MAXPART], gen_matched_rec_theta_[MAXPART], gen_matched_rec_pt_[MAXPART], gen_matched_rec_qoverp_[MAXPART], gen_matched_rec_cotth_[MAXPART], gen_matched_rec_phi_[MAXPART];

  bool is_gsf_;
  edm::InputTag track_label_gsf_, track_label_, el_seed_label_;

  std::vector<std::string> trackerContainers;
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

  trackerContainers = iConfig.getParameter<std::vector<std::string> >("ROUList");
}


MakeTrackValTree::~MakeTrackValTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

const TrackingRecHit* MakeTrackValTree::getHitPtr(edm::OwnVector<TrackingRecHit>::const_iterator iter) {return &*iter;}
const TrackingRecHit* MakeTrackValTree::getHitPtr(trackingRecHit_iterator iter) {return &**iter;}


template<typename iter>
bool MakeTrackValTree::isMatchedRecHitById(iter it, TrackingParticle* tp){

  std::vector<unsigned int> tpids; //Save a vector of track ID-s, related to the tracking particle
  for (TrackingParticle::g4t_iterator g4T=tp->g4Track_begin(); g4T!=tp->g4Track_end(); g4T++)
    tpids.push_back(g4T->trackId());  
    
  const TrackingRecHit *hit = getHitPtr(it); //find the IDs of Sim. tracks, associated to the rec hit
  std::vector<SimHitIdpr> matchedSimIds = hitAssociator->associateHitId(*hit); 

  bool goodhit = false;
    
  for( unsigned int iRhId = 0; iRhId < matchedSimIds.size(); iRhId++) //loop over sim IDs associated to rec hit
    for(unsigned int iTpId = 0; iTpId < tpids.size(); iTpId++){ //loop over sim IDs associated to a TP
      if( matchedSimIds[iRhId].first == tpids[iTpId] && matchedSimIds[iRhId].second == tp->eventId() ) {
  	goodhit = true;
  	break;
      }
    }

  return goodhit;
}

std::map<int, std::vector<PSimHit> > MakeTrackValTree::getSimHits(std::vector<std::string>  trackerContainers, const edm::Event& iEvent){
  // get all sim hits in the event and save them in vectors by trackId of simTracks (later match to TPs)

  std::map< int, std::vector<PSimHit> > simHitMap;

  for ( auto const& trackerContainer : trackerContainers){
    edm::Handle<std::vector<PSimHit> > simHits;
    edm::InputTag tag_hits("g4SimHits", trackerContainer);
    iEvent.getByLabel(tag_hits, simHits);

    
    for (std::vector<PSimHit>::const_iterator ihit = simHits->begin(); ihit != simHits->end(); ihit++) {
      simHitMap[(*ihit).trackId()].push_back((*ihit));
    }
  }
  return simHitMap;
}

std::pair<std::vector<const reco::ElectronSeed*>, double> MakeTrackValTree::findMatchedSeed(edm::Handle<edm::View<reco::ElectronSeed>> elSeedCollection, TrackingParticle* tp){
  int nr_matched_seeds = 0;

  std::vector<const reco::ElectronSeed*> bestSeed; // a reco seed that best matches the tracking particle
  bestSeed.clear();
  double best_seed_quality = 0;
  for( edm::View<reco::ElectronSeed>::const_iterator it_seed = elSeedCollection->begin(); it_seed != elSeedCollection->end(); it_seed++){
    
    const reco::ElectronSeed* elSeed = &*it_seed; // remove iterator for saving output
    TrajectorySeed::range rechits = it_seed->recHits(); //get recHits associated to the seed

    unsigned int matched_hits = 0;
    for( edm::OwnVector<TrackingRecHit>::const_iterator rechit_it = rechits.first; rechit_it != rechits.second; rechit_it++){
      if(isMatchedRecHitById(rechit_it, tp) )
      	matched_hits++; // count matched hits
    }
    double match_quality = (double)matched_hits/(double)(it_seed->nHits());
    
    if( match_quality > best_seed_quality ){ // compare to the best quality seed, if better replace, initial best_seed_quality=0
      nr_matched_seeds++;
      if( bestSeed.size() )
	bestSeed.pop_back(); 
      bestSeed.push_back(elSeed);
      best_seed_quality = match_quality;
    }
    
    //    std::cout<<"Matched seed nr. "<<nr_matched_seeds+1<<"Matched hit quality = " << (double)(matched_hits/it_seed->nHits()) << std::endl;     
    //vstd::cout<<"Numbr of rechits = "<<it_seed->nHits()<<", subdet1 = "<<it_seed->subDet1()<<", subdet2 = "<<it_seed->subDet2()<<", tracker_dr//      iven = "<<it_seed->isTrackerDriven()<<", ecal driven = "<<it_seed->isEcalDriven()<<", charge = "<<(int)it_seed->getCharge()<<std::endl; 

   } // <-- end loop over electron seeds

  // std::cout<<"number of matched seeds related to the TP = "<< nr_matched_seeds << std::endl; // Undersatnd why so many (multiple objects for one seed)
  
  std::pair <std::vector<const reco::ElectronSeed*>, double> best_seed_info;
  best_seed_info = std::make_pair( bestSeed, best_seed_quality);


  //  if( (best_seed_info.first).size() > 0)
  //  std::cout<<"best seed charge "<<best_seed_info.first[0]->getCharge()<<std::endl;

  return best_seed_info; 
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

std::vector<PSimHit>


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
     gen_nr_simhits_[i] = -99;

     gen_matched_pt_[i] = -99;
     gen_matched_eta_[i] = -99;
     gen_matched_phi_[i] = -99;
     gen_matched_qoverp_[i] = -99;
     gen_matched_theta_[i] = -99;
     gen_matched_cotth_[i] = -99;

     //-- seeds matched to TPs ---
     is_ecalDrivenSeed_[i] = -10;
     is_trackerDrivenSeed_[i] = -10;
     gen_matched_seed_okCharge_[i] = -10;
     gen_matched_seed_nshared_[i] = -10;
     gen_matched_seed_quality_[i] = -10;

     //------ reco tracks --
     is_gen_matched_[i] = -10;
     reco_pt_[i] = -99;
     reco_eta_[i] = -99;
     reco_phi_[i] = -99;
     fake_pt_[i] = -99;
     fake_eta_[i] = -99;
     fake_phi_[i] = -99;

     gen_matched_rec_pt_[i] = -99;
     gen_matched_rec_eta_[i] = -99;
     gen_matched_rec_theta_[i] = -99;
     gen_matched_rec_qoverp_[i] = -99;
     gen_matched_rec_cotth_[i] = -99;
     gen_matched_rec_phi_[i] = -99;
     
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

   edm::Handle<TrackingVertexCollection> tvH; //For checking PU vertices
   iEvent.getByLabel("mix","MergedTrackTruth",tvH);

   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get(tracker);

   TrackingVertexCollection tv = *tvH;

   const TrackingParticleCollection tPCeff = *(TPCollection.product());
   TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP);

   edm::ESHandle<TrackAssociatorBase> myAssociator;
   iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", myAssociator);
   reco::SimToRecoCollection simRecColl = myAssociator->associateSimToReco(trackCollection, TPCollection, &iEvent, &iSetup);
   reco::RecoToSimCollection recSimColl = myAssociator->associateRecoToSim(trackCollection, TPCollection, &iEvent, &iSetup);
   hitAssociator = new TrackerHitAssociator(iEvent, conf_); //to access functions from the hitAsssociator code

   edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP;
   iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP", parametersDefinerTP);
   
   np_gen_ = 0;
   np_gen_toReco_ = 0;

   std::map<int, std::vector<PSimHit>> simHitMap = getSimHits(trackerContainers, iEvent); // get all simHits in the event

   for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
     TrackingParticleRef tpr(TPCollection, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
     
     TrackingParticle::Point vertexTP;
     TrackingParticle::Vector momentumTP;

     TrackingParticle::Point vertex = parametersDefinerTP->vertex(iEvent,iSetup,tpr);
     TrackingParticle::Vector momentum = parametersDefinerTP->momentum(iEvent,iSetup,tpr);


     if( !tpSelector(*tp) ) continue;
     if( tp->pt() < 1 ) continue; // such tracks are irrelevant for electrons

     
     momentumTP = tp->momentum();
     vertexTP = tp->vertex(); 

     double dxySim = -vertex.x()*sin(momentum.phi()) + vertex.y()*cos(momentum.phi());
     double dzSim = vertex.z() - (vertex.x()*momentum.x()+vertex.y()*momentum.y())/sqrt(momentum.perp2())* momentum.z()/sqrt(momentum.perp2());

     std::cout<<"dxySim = " << dxySim << ", dzSim = " << dzSim<< std::endl;
     
     gen_eta_[np_gen_] = tp->eta();
     gen_phi_[np_gen_] = tp->phi();
     gen_pt_[np_gen_] = tp->pt();
     gen_pdgId_[np_gen_] = tp->pdgId();
     
     //--------check for matched reco track defined by AssociatorByHits (efficiency and fake-rate plots)-----------
     const reco::Track* matchedTrackPointer=0;
     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
     //     int nSharedHits = 0; int nRecoTrackHits=0;

     if(simRecColl.find(tpr) != simRecColl.end()){
       rt = simRecColl[tpr]; // find sim-to-reco association
       if ( rt.size()!=0 ) {
	 matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track                                                 
	 //	 nSharedHits = rt.begin()->second;
	 //	 nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
	 is_reco_matched_[np_gen_] = 1;

	 //-------------- efficiencies-----------------
	 gen_matched_eta_[np_gen_] = tp->eta();
	 gen_matched_pt_[np_gen_] = tp->pt();

	 //------------- for resolutions, additional gen vars  -------------------
	 gen_matched_phi_[np_gen_] = tp->phi();
	 gen_matched_qoverp_[np_gen_] = tp->charge()/(tp->px()*tp->px() + tp->py()*tp->py() + tp->pz()*tp->pz());
	 gen_matched_theta_[np_gen_] = tp->theta();
	 gen_matched_cotth_[np_gen_] = 1./tan(tp->theta());
	 //------------ for resolutions, gen matched reco vars --------------------
	 gen_matched_rec_pt_[np_gen_] = matchedTrackPointer->pt();
	 gen_matched_rec_theta_[np_gen_] = matchedTrackPointer->theta();
	 gen_matched_rec_eta_[np_gen_] = matchedTrackPointer->eta();
	 gen_matched_rec_qoverp_[np_gen_] = matchedTrackPointer->qoverp();
	 gen_matched_rec_cotth_[np_gen_] = 1./(tan(matchedTrackPointer->theta()) );
	 gen_matched_rec_phi_[np_gen_] = matchedTrackPointer->phi();
	 //-------------------------------------------------------------------------

	 np_gen_toReco_++; //count reco matched TPs
	 //	 std::cout<<"TP matched with RECO track"<<std::endl;
       }
     }
     else
       is_reco_matched_[np_gen_] = 0;

     if( findMatchedSeed(elSeedCollection, tp).first.size() > 0.5 && findMatchedSeed(elSeedCollection, tp).second >0.5 ){
       const reco::ElectronSeed* bestSeed = findMatchedSeed(elSeedCollection, tp).first.at(0);

       gen_matched_seed_quality_[np_gen_] = findMatchedSeed(elSeedCollection,tp).second;

       gen_matched_seed_nshared_[np_gen_] = (int)bestSeed->nHits();
       gen_matched_seed_okCharge_[np_gen_] = (int)bestSeed->getCharge()*(int)tp->charge(); //Check whether seed charge matches TP charge

       is_ecalDrivenSeed_[np_gen_] = 0;
       is_trackerDrivenSeed_[np_gen_] = 0;
       if( bestSeed->isEcalDriven() )
	 is_ecalDrivenSeed_[np_gen_] = 1;    
       if(bestSeed->isTrackerDriven())
	 is_trackerDrivenSeed_[np_gen_] = 1;
 
     }
     gen_nr_simhits_[np_gen_] = tp->numberOfTrackerHits();
     
     //check whether TP is from PU vertices
     /*     double vtx_z_PU = vertexTP.z();
     for (size_t j = 0; j < tv.size(); j++) {
       if (tp->eventId().event() == tv[j].eventId().event()) {
	 vtx_z_PU = tv[j].position().z();
	 //	 std::cout<<"PU vertex -- with position: " << vtx_z_PU << std::endl; 
       }
     }
     */

     //     std::cout<<"Matched seed charge = "<<gen_matched_seed_okCharge_[np_gen_]<<", match quality = "<< gen_matched_seed_quality_[np_gen_]<<", nr matched seed hits = "<< gen_matched_seed_nshared_[np_gen_]<<std::endl;     

     std::cout<<"Found tracking particle with pt = "<<tp->pt()<<", eta = "<<tp->eta()<<std::endl; 
     //     std::cout<<"nr shared hits = "<<nSharedHits<<", nr reco hits "<<nRecoTrackHits<<" matched reco pt = "<<gen_matched_rec_pt_[np_gen_]<<", gen. qoverp = "<<gen_matched_qoverp_[np_gen_]<<", matched reco qoverp = "<<gen_matched_rec_qoverp_[np_gen_]<<std::endl;


  
 

     //-------------------------- Get sim hits of the TP ----------------------------------------
     std::vector<PSimHit> simhits;
     std::vector<PSimHit> simhitsTP;
     //std::vector< std::vector<PSimHit>::const_iterator > simhits_sorted;
     simhitsTP.clear(); //vector of simhits for sorting for each TP

     std::pair<double, int> hit_dist_idx; //for sorting list of hits by distance
     std::vector< std::pair<double,int> > hit_dist_idx_vec;

     int hit_idx=0; //keep track of indexes in simhits_sorted
     for( TrackingParticle::g4t_iterator g4T = tp->g4Track_begin(); g4T != tp->g4Track_end(); g4T++){ 
       std::map<int, std::vector<PSimHit> >::iterator simHitMap_simTrack = simHitMap.find(g4T->trackId()); // get simhits associated to simTrack

       simhits.clear();
       if( simHitMap_simTrack != simHitMap.end() ){
       	 simhits = simHitMap_simTrack->second;

	 for ( std::vector<PSimHit>::const_iterator it_hit = simhits.begin(); it_hit != simhits.end(); ++it_hit){
	   const PSimHit ihit = *it_hit;   

	   DetId dId = DetId(it_hit->detUnitId() );
	   const GeomDetUnit* detunit = tracker->idToDetUnit(dId.rawId());

	   //	   LocalVector local_p = it_hit->momentumAtEntry(); //get hit momentum
	   //GlobalVector global_p = detunit->toGlobal(local_p);
	   
	   Local3DPoint local_x = it_hit->localPosition(); //get hit distance
	   GlobalPoint global_x = detunit->toGlobal(local_x);

	   hit_dist_idx = std::make_pair( (double)global_x.mag(), hit_idx);
	   hit_dist_idx_vec.push_back(hit_dist_idx);

	   simhitsTP.push_back(ihit); 	  
	   hit_idx++; //index simhits in tracking particle
	 } //<--- end loop over simhits in simTrack
       }
     } // <-- end loop over g4T simTracks of TP

     
     std::vector<PSimHit> sorted_simhits;
     std::vector<PSimHit> cleaned_simhits;

     std::vector<std::vector<int> > simhit_position;
     simhit_position.clear();
     cleaned_simhits.clear();
     sorted_simhits.reserve(simhitsTP.size());
     

     if ( hit_dist_idx_vec.size() && simhitsTP.size() ){
       std::sort(hit_dist_idx_vec.begin(), hit_dist_idx_vec.end()); // sort global positions along with simhits_sorted indices

       for (unsigned int i = 0; i != hit_dist_idx_vec.size(); i++){ // reorder simhitsTP
	 sorted_simhits.push_back(simhitsTP.at(hit_dist_idx_vec.at(i).second ) );

	 DetId dId = DetId(sorted_simhits[i].detUnitId() );
	 simhit_position.push_back( getHitPosition(dId) );
	 
       }

       for (unsigned int i = 0; i != sorted_simhits.size(); i++){ // clean up sorted sim hits

	 DetId dId = DetId(sorted_simhits[i].detUnitId() );
	 const GeomDetUnit* detunit = tracker->idToDetUnit(dId.rawId());

	 Local3DPoint local_x = sorted_simhits[i].localPosition(); //get hit distance
	 GlobalPoint global_x = detunit->toGlobal(local_x);

	 LocalVector local_p = sorted_simhits[i].momentumAtEntry(); //get hit momentum
 	 GlobalVector global_p = detunit->toGlobal(local_p);
	 
	 std::cout<<"---------------"<<std::endl;
	 std::vector<int> hitposition = getHitPosition(dId, true); //get hit position wrt subdetectors                                                                                                    
	 std::cout<<"hit pt = "<<global_p.perp()<<", eta = "<<global_p.eta()<<", phi = "<<global_p.phi()<<std::endl;
	 std::cout<<"dist = "<<global_x.mag()<<std::endl;


	 if (global_p.perp() < 0.5 ) //we dont consider TP-s < 1 GeV at this point, this is noise (no hope to match to ECAL)
	   // || (global_p.perp() < 0.02 && 
	   //	     ((i>0 && simhit_position[i][0] == simhit_position[i-1][0] && simhit_position[i][1] == simhit_position[i-1][1] ) || 
	   //	      (i<sorted_simhits.size()-1 && simhit_position[i][0] == simhit_position[i+1][0] && simhit_position[i][1] == simhit_position[i+1][1] )
	   //)))
	   {
	     std::cout<<"REMOVE THIS HIT" << std::endl;
	     //	     sorted_simhits.erase(sorted_simhits.begin() + i); // remove repeated hit
	   }
	 else
	   cleaned_simhits.push_back(sorted_simhits[i]);

       }
       std::cout<<"Finish loop over sorted simhits"<<std::endl;
     }

    
     std::cout<<"Nr simhits from TP = "<< gen_nr_simhits_[np_gen_]<<std::endl;
     std::cout<<"counted hits in simTracks from TP = "<<hit_idx<<std::endl;
     std::cout<<"cleaned_simhits.size() = " <<cleaned_simhits.size()<<std::endl;

     np_gen_++; // count tracking particles that pass the standard quality cuts
   } // <-- end loop over tracking particles

   //   std::cout<<"nr. sel TP-s = " <<np_gen_<<std::endl;

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

   //   for( edm::View<reco::ElectronSeed>::const_iterator it_seed = elSeedCollection->begin(); it_seed != elSeedCollection->end(); it_seed++){
   //  TrajectorySeed::range rechits = it_seed->recHits(); //get recHits associated to the seed
   // } //FIXME add fake seeds


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
  trackValTree_->Branch("np_gen_toReco", &np_gen_toReco_, "np_gen_toReco/I");

  trackValTree_->Branch("gen_reco_matched", is_reco_matched_, "gen_reco_matched[np_gen]/I");

  trackValTree_->Branch("gen_pdgId", gen_pdgId_, "gen_pdgId[np_gen]/I");
  trackValTree_->Branch("gen_pt", gen_pt_, "gen_pt[np_gen]/D");
  trackValTree_->Branch("gen_eta", gen_eta_, "gen_eta[np_gen]/D");
  trackValTree_->Branch("gen_phi", gen_phi_, "gen_phi[np_gen]/D");
  trackValTree_->Branch("gen_nr_simhits", gen_nr_simhits_, "gen_nr_simhits_[np_gen]/I");

  //----------------- TP simToReco matching-----------------
  trackValTree_->Branch("gen_matched_pt", gen_matched_pt_, "gen_matched_pt[np_gen]/D");
  trackValTree_->Branch("gen_matched_qoverp", gen_matched_qoverp_, "gen_matched_qoverp[np_gen]/D");
  trackValTree_->Branch("gen_matched_eta", gen_matched_eta_, "gen_matched_eta[np_gen]/D");
  trackValTree_->Branch("gen_matched_theta", gen_matched_theta_, "gen_matched_theta[np_gen]/D");
  trackValTree_->Branch("gen_matched_cotth", gen_matched_cotth_, "gen_matched_cotth[np_gen]/D");
  trackValTree_->Branch("gen_matched_phi", gen_matched_phi_, "gen_matched_phi[np_gen]/D");

  trackValTree_->Branch("gen_matched_rec_pt", gen_matched_rec_pt_, "gen_matched_rec_pt[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_eta", gen_matched_rec_eta_, "gen_matched_rec_eta[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_qoverp", gen_matched_rec_qoverp_, "gen_matched_rec_qoverp[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_theta", gen_matched_rec_theta_, "gen_matched_rec_theta[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_cotth", gen_matched_rec_cotth_, "gen_matched_rec_cotth[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_phi", gen_matched_rec_phi_, "gen_matched_rec_phi[np_gen]/D");
  //----------------- TP seed matching-------------------------
  trackValTree_->Branch("gen_matchedSeedOkCharge", gen_matched_seed_okCharge_, "gen_matchedSeedOkCharge[np_gen]/I");
  trackValTree_->Branch("gen_matchedSeedQuality", gen_matched_seed_quality_, "gen_matchedSeedQuality[np_gen]/D");
  trackValTree_->Branch("gen_nrMatchedSeedHits", gen_matched_seed_nshared_, "gen_nrMatchedSeedHits[np_gen]/I");
  trackValTree_->Branch("is_ecalDrivenSeed", is_ecalDrivenSeed_, "is_ecalDrivenSeed[np_gen]/I");
  trackValTree_->Branch("is_trackerDrivenSeed", is_trackerDrivenSeed_, "is_trackerDrivenSeed[np_gen]/I");

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

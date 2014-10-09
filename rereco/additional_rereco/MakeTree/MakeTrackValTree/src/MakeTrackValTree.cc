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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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

#include "SimDataFormats/Track/interface/SimTrack.h"
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
      GlobalVector getHitMomentum(std::vector<PSimHit>::const_iterator, edm::ESHandle<TrackerGeometry>, bool);
      std::vector<PSimHit> getSimHitsTP(TrackingParticle*, std::map<int, std::vector<PSimHit> >, edm::ESHandle<TrackerGeometry>,  bool, bool, bool);
      bool isGoodTrack(edm::RefToBase<reco::Track>, math::XYZPoint, std::vector<math::XYZPoint>, bool); // use instead of AnalyticalTrackSelector
  math::XYZPoint getClosestPVpoint( edm::RefToBase<reco::Track>, std::vector<math::XYZPoint> ); 


  // ------------------ thresholds for TPselector ----------
  double ptMinTP, minRapidityTP, maxRapidityTP, tipTP, lipTP, signalOnlyTP, chargedOnlyTP, stableOnlyTP; 
  std::vector<int> pdgIdTP;
  int minHitTP;
  // -------------------- for tree -------------
  TTree *trackValTree_;
  int run_nr_, evt_nr_, lumi_nr_;
  int np_gen_, np_gen_toReco_, np_reco_, np_reco_toGen_,np_fake_;
  int is_reco_matched_[MAXPART], is_reco_matched_sim_[MAXPART], is_gen_matched_[MAXPART], is_charge_matched_[MAXPART];
  double simrec_simD_quality_[MAXPART], simrec_recD_quality_[MAXPART];
  int gen_pdgId_[MAXPART], gen_nr_simhits_[MAXPART], gen_matched_seed_nshared_[MAXPART], gen_matched_seed_okCharge_[MAXPART], is_ecalDrivenSeed_[MAXPART], is_trackerDrivenSeed_[MAXPART];

  double reco_pt_[MAXPART], reco_eta_[MAXPART], reco_phi_[MAXPART], fake_pt_[MAXPART], fake_eta_[MAXPART], fake_phi_[MAXPART], fake_dxy_[MAXPART], fake_dz_[MAXPART], fake_nr_rechits_[MAXPART];

  // parameters of TPs and reco matched TPs
  double gen_pt_[MAXPART], gen_eta_[MAXPART], gen_phi_[MAXPART], gen_matched_pt_[MAXPART], gen_matched_qoverp_[MAXPART], gen_matched_cotth_[MAXPART], gen_matched_eta_[MAXPART], gen_matched_theta_[MAXPART], gen_matched_phi_[MAXPART], gen_matched_dz_[MAXPART], gen_matched_dxy_[MAXPART], gen_dxy_[MAXPART], gen_dz_[MAXPART], gen_ptAtLast_[MAXPART], gen_bremFraction_[MAXPART], gen_matched_seed_quality_[MAXPART];  

  // parameters of reco tracks that have been matched to TPs
  double gen_matched_rec_eta_[MAXPART], gen_matched_rec_theta_[MAXPART], gen_matched_rec_pt_[MAXPART], gen_matched_rec_qoverp_[MAXPART], gen_matched_rec_cotth_[MAXPART], gen_matched_rec_phi_[MAXPART], gen_matched_rec_dxy_[MAXPART], gen_matched_rec_dz_[MAXPART];
  int gen_matched_rec_nhits_[MAXPART];
  double pt_pull_[MAXPART], theta_pull_[MAXPART], phi_pull_[MAXPART], dxy_pull_[MAXPART], dz_pull_[MAXPART], qoverp_pull_[MAXPART]; 

  bool is_gsf_, is_single_part_;
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
  is_single_part_ = iConfig.getParameter<bool>("isSinglePart");
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

  // --- properly initialize the map ---- //FIXME!!!
  // for 
  // simHitMap.insert(std::pair<int, std::vector<PSimHit> >(  simHits);
  
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

std::vector<PSimHit> MakeTrackValTree::getSimHitsTP(TrackingParticle* tp, std::map<int, std::vector<PSimHit>> simHitMap,  edm::ESHandle<TrackerGeometry> tracker, bool sort = true, bool cleanup = true, bool debug = false ){
  // -- get sim hits in a tracking particle, optionally order by distance from PV and clean from noise and low pt hits
  // input arguments: TP; map of all simHits in event (mapped by trackId); tracker_geometry

  std::vector<PSimHit> simhits; // simhits from a simTrack (by trackID)
  std::vector<PSimHit> simhitsTP; //all simhits in a TP
  std::vector<PSimHit> simhits_sorted; //simhits sorted by distance from pv
  std::vector<PSimHit> simhits_cleaned; //simhits cleaned up from noise hits
  simhitsTP.clear();
  simhits_sorted.clear();
  simhits_cleaned.clear();

  std::pair<double, int> hit_dist_idx; //for sorting list of hits by distance
  std::vector< std::pair<double,int> > hit_dist_idx_vec;

  std::vector<std::vector<int> > simhit_position;
  
  int hit_idx=0; //keep track of indexes in simhits_sorted
  int nr_simtracks =0;
  for( TrackingParticle::g4t_iterator g4T = tp->g4Track_begin(); g4T != tp->g4Track_end(); g4T++){  //loop over sim tracks in TP
    nr_simtracks++;
    std::map<int, std::vector<PSimHit> >::iterator simHitMap_simTrack = simHitMap.find(g4T->trackId()); // get simhits associated to simTrack
    simhits.clear();
    if( simHitMap_simTrack != simHitMap.end() ){
      simhits = simHitMap_simTrack->second;
      for ( std::vector<PSimHit>::const_iterator it_hit = simhits.begin(); it_hit != simhits.end(); ++it_hit){
	const PSimHit ihit = *it_hit;

	DetId dId = DetId(it_hit->detUnitId() );
	const GeomDetUnit* detunit = tracker->idToDetUnit(dId.rawId());
	
	Local3DPoint local_x = it_hit->localPosition(); //get hit distance from PV (to sort in time order)
	GlobalPoint global_x = detunit->toGlobal(local_x);

	hit_dist_idx = std::make_pair( (double)global_x.mag(), hit_idx); //save hit distance with index for sorting
	hit_dist_idx_vec.push_back(hit_dist_idx);

	simhitsTP.push_back(ihit);
	hit_idx++; //index simhits in tracking particle
      } //<--- end loop over simhits in simTrack
    }
	//    std::cout<<"sim track with nrhits = " << simhits.size() << std::endl;
  } // <-- end loop over g4T simTracks of TP

  if(debug){
    std::cout<< "--------------------" << std::endl;
    std::cout<<"nr sim tracks = " << nr_simtracks << std::endl;
  }
  if ( hit_dist_idx_vec.size() && simhitsTP.size() ){
    std::sort(hit_dist_idx_vec.begin(), hit_dist_idx_vec.end()); // sort global positions along with simhits_sorted indices
  
    for (unsigned int i = 0; i != hit_dist_idx_vec.size(); i++){ // reorder simhitsTP
      simhits_sorted.push_back(simhitsTP.at(hit_dist_idx_vec.at(i).second ) );

      DetId dId = DetId(simhits_sorted[i].detUnitId() );
      simhit_position.push_back( getHitPosition(dId) );
    }

    

    for (unsigned int i = 0; i != simhits_sorted.size(); i++){ // clean up sorted sim hits
      DetId dId = DetId(simhits_sorted[i].detUnitId() );
      const GeomDetUnit* detunit = tracker->idToDetUnit(dId.rawId());

      Local3DPoint local_x = simhits_sorted[i].localPosition(); //get hit distance
      GlobalPoint global_x = detunit->toGlobal(local_x);

      LocalVector local_p = simhits_sorted[i].momentumAtEntry(); //get hit momentum
      GlobalVector global_p = detunit->toGlobal(local_p);

      if (global_p.perp() > 0.9 )  //0.9 the rest is noise or impossible to reach ECal
		simhits_cleaned.push_back(simhits_sorted[i]);
      else if (debug)
	std::cout<<"REMOVE THIS HIT" << std::endl;

      if(debug){
	std::vector<int> hitposition = getHitPosition(dId, false); //get hit position wrt subdetectors
	std::cout<<"hit pt = "<<global_p.perp()<<", eta = "<<global_p.eta()<<", phi = "<<global_p.phi()<<std::endl;
	std::cout<<"dist = "<<global_x.mag()<<std::endl;
	std::cout<<"---------------"<<std::endl;
      }
    } // <-- end loop over sorted simhits collection for cleaning
  } // <-- end if simhitsTP.size()

  if(debug){
    std::cout<<"Nr simhits from TP = "<< gen_nr_simhits_[np_gen_]<<std::endl;
    std::cout<<"counted hits in simTracks from TP = "<<hit_idx<<std::endl;
    std::cout<<"simhits_cleaned.size() = " <<simhits_cleaned.size()<<std::endl;
  }

  return simhits_cleaned;
}


math::XYZPoint MakeTrackValTree::getClosestPVpoint(edm::RefToBase<reco::Track> trk, std::vector<math::XYZPoint> points){

  math::XYZPoint PV;
  PV.SetXYZ(0,0,0);
  float dz_max = 1000;

  for (std::vector<math::XYZPoint>::const_iterator point = points.begin(); point != points.end(); ++point) { 
    double dz_PV = trk->dz(*point);

    //    std::cout<<"Checking vertex with dxy = "<<trk->dxy(*point) << ", dz = " << dz_PV << std::endl;

    if( fabs(dz_PV) < fabs(dz_max) ){
      PV = *point;
      dz_max = dz_PV;
    }
      
  }
  return PV;
}

bool MakeTrackValTree::isGoodTrack(edm::RefToBase<reco::Track> trk, math::XYZPoint bsPosition, std::vector<math::XYZPoint> points, bool debug = false){
  // implementation of/RecoTracker/FinalTrackSelectors/python/selectHighPurity_cfi.py for reco GSF tracks
  
  double max_dxy = 3;
  double max_dz = 30;
  double res_par[2] = {0.003, 0.01};
  double d0_par1[2] = {0.3, 4.0};
  double dz_par1[2] = {0.35, 4.0};
  double d0_par2[2] = {0.4, 4.0};
  double dz_par2[2] = {0.4, 4.0};
  
  int min_layers = 3;
  int min_3Dlayers = 3;
  int max_lostLayers = 2; 
  double chi2n_par = 0.7;

  int nlayers = trk->hitPattern().trackerLayersWithMeasurement();
  int nlayers3D = trk->hitPattern().pixelLayersWithMeasurement() + trk->hitPattern().numberOfValidStripLayersWithMonoAndStereo();
  int nlayersLost = trk->hitPattern().trackerLayersWithoutMeasurement();

  double d0 = trk->dxy(bsPosition), d0E = trk->d0Error();
  double dz = trk->dz(bsPosition), dzE = trk->dzError();
  
  if( debug ){
	std::cout << "trk pt = " << trk->pt() << ", dxy = " << d0 << ", dz = " << dz << ", ndof = "<< trk->ndof() << std::endl;
	std::cout << "nlayers = " << nlayers << ", nlayers3D = " << nlayers3D << ", nlayersLost = " << nlayersLost <<std::endl;
  }
  //--------- distance from PV / beam spot -----------

  // parametrized d0 resolution for the track pt
  double nomd0E = sqrt(res_par[0]*res_par[0]+(res_par[1]/std::max(trk->pt(),1e-9))*(res_par[1]/std::max(trk->pt(),1e-9)));
  // parametrized z0 resolution for the track pt and eta
  double nomdzE = nomd0E*(std::cosh(trk->eta()) );
  double dzCut = std::min( pow(dz_par1[0]*nlayers,dz_par1[1])*nomdzE,
					  pow(dz_par2[0]*nlayers,dz_par2[1])*dzE );
  double d0Cut = std::min( pow(d0_par1[0]*nlayers,d0_par1[1])*nomd0E,
					  pow(d0_par2[0]*nlayers,d0_par2[1])*d0E );

  // ---- PrimaryVertex compatibility cut
  bool primaryVertexZCompatibility(false);
  bool primaryVertexD0Compatibility(false);

  
  for (std::vector<math::XYZPoint>::const_iterator point = points.begin(); point != points.end(); ++point) { //at least one primary vertex required
	double dzPV = trk->dz(*point); //re-evaluate the dz with respect to the vertex position
	double d0PV = trk->dxy(*point); //re-evaluate the dxy with respect to the vertex position
	if (abs(dzPV) < dzCut)  primaryVertexZCompatibility = true;
	if (abs(d0PV) < d0Cut) primaryVertexD0Compatibility = true;

	if(primaryVertexZCompatibility && primaryVertexD0Compatibility) {
	  if ( debug )
		std::cout << "distances " << dzPV << " " << d0PV << " vs " << dzCut << " " << d0Cut << std::endl;
	  break;
	}
  }
  if (points.size() < 0.5){ // for all single particle events as vertex ndof=0
	if (abs(dz) < max_dz )
	  primaryVertexZCompatibility = true;
	if (abs(d0) < max_dxy )
	  primaryVertexD0Compatibility = true;
  }
  

  if( !(primaryVertexZCompatibility && primaryVertexD0Compatibility) ){
	if( debug )
	  std::cout<<"didnt pass the ip cuts"<<std::endl;
	return false;
  }
  //--------- hit pattern quality ------------
  
  if(nlayers < min_layers) return false;
  if(nlayers3D < min_3Dlayers) return false;
  if(nlayersLost > max_lostLayers) return false;
  
  //--------------- chi2 ---------------

  int count1dhits = 0;
  double chi2n = trk->normalizedChi2();
  for (trackingRecHit_iterator ith = trk->recHitsBegin(); ith != trk->recHitsEnd(); ++ith) {
	const TrackingRecHit * hit = ith->get();
	if (hit->isValid()) {
	  if (typeid(*hit) == typeid(SiStripRecHit1D))
		++count1dhits;
	}
  }
  if (count1dhits > 0) {
	double chi2 = trk->chi2();
	chi2n = (chi2+count1dhits)/double(trk->ndof()+count1dhits);
  }

  if( debug )
	std::cout<<"chi2 = " << chi2n << std::endl;

  if (chi2n > chi2n_par*nlayers) return false;
  
  return true;
}


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
     is_reco_matched_sim_[i] = -10;
     is_charge_matched_[i] = -10;
     
     simrec_simD_quality_[i] = -10;
     simrec_recD_quality_[i] = -10;
     
     gen_eta_[i] = -99;
     gen_phi_[i] = -99;
     gen_pt_[i] = -99;
     gen_pdgId_[i] = -99;
     gen_nr_simhits_[i] = -99;

     gen_eta_[i] = -999;
     gen_phi_[i] = -999;
     gen_pt_[i] = -999;
     gen_ptAtLast_[i] = -999;
     gen_bremFraction_[i] = -10;
     gen_pdgId_[i] = -999;
     gen_dxy_[i] = -10;
     gen_dz_[i] = -10;
     gen_nr_simhits_[i] = -10;
     
     //------for resolutions-------
     gen_matched_pt_[i] = -99;
     gen_matched_eta_[i] = -99;
     gen_matched_phi_[i] = -99;
     gen_matched_qoverp_[i] = -99;
     gen_matched_theta_[i] = -99;
     gen_matched_cotth_[i] = -99;
     gen_matched_dxy_[i] = -999;
     gen_matched_dz_[i] = -999;
     
     //----- pulls -----------
     pt_pull_[i] = -99;
     qoverp_pull_[i] = -99;
     theta_pull_[i] = -99;
     phi_pull_[i] = -99;
     dxy_pull_[i] = -99;
     dz_pull_[i] = -99;
     
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
     fake_dxy_[i] = -99;
     fake_dz_[i] = -99;
     fake_nr_rechits_[i] = -99;
	 
     gen_matched_rec_pt_[i] = -99;
     gen_matched_rec_eta_[i] = -99;
     gen_matched_rec_theta_[i] = -99;
     gen_matched_rec_qoverp_[i] = -99;
     gen_matched_rec_cotth_[i] = -99;
     gen_matched_rec_phi_[i] = -99;
     gen_matched_rec_dxy_[i] = -99;
     gen_matched_rec_dz_[i] = -99;
     gen_matched_rec_nhits_[i] = -99;
   }

   edm::InputTag track_label; 
   if(is_gsf_)
     track_label = track_label_gsf_;
   else
     track_label = track_label_;

   //   std::cout<<"RECO track label = "<<track_label<<std::endl;
   
   edm::Handle<edm::View<reco::Track> > trackCollection; //reconstructed tracks
   iEvent.getByLabel(track_label, trackCollection);
   
   edm::Handle<edm::View<reco::ElectronSeed> > elSeedCollection; 
   iEvent.getByLabel(el_seed_label_, elSeedCollection);

   // Tracking particles
   edm::Handle<TrackingParticleCollection>  TPCollection ; //simulated tracks
   iEvent.getByLabel("mix","MergedTrackTruth",TPCollection);
   //   std::cout<<"Number of TPs in the event = " << TPCollection->size() << std::endl;
   
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get(tracker);

   edm::Handle<reco::BeamSpot> recoBeamSpotHandle; //get beam spot position
   iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle);
   reco::BeamSpot bs = *recoBeamSpotHandle;
   math::XYZPoint bsPosition = bs.position();
   
   const TrackingParticleCollection tPCeff = *(TPCollection.product());
   TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP);

   edm::ESHandle<TrackAssociatorBase> myAssociator;
   edm::ESHandle<TrackAssociatorBase> myAssociatorRecoDenom;
   edm::ESHandle<TrackAssociatorBase> myAssociatorSimDenom;
   iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", myAssociator);
   iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHitsRecoDenom", myAssociatorRecoDenom);
   iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHitsSimDenom", myAssociatorSimDenom);

   reco::SimToRecoCollection simRecCollRecoDenom = myAssociatorRecoDenom->associateSimToReco(trackCollection,TPCollection, &iEvent, &iSetup);
   reco::SimToRecoCollection simRecCollSimDenom = myAssociatorSimDenom->associateSimToReco(trackCollection,TPCollection, &iEvent, &iSetup);

   //   reco::SimToRecoCollection simRecColl = myAssociator->associateSimToReco(trackCollection, TPCollection, &iEvent, &iSetup);
   reco::RecoToSimCollection recSimColl = myAssociator->associateRecoToSim(trackCollection, TPCollection, &iEvent, &iSetup);
   hitAssociator = new TrackerHitAssociator(iEvent, conf_); //to access functions from the hitAsssociator code

   edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP;
   iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP", parametersDefinerTP);
   

   edm::Handle<TrackingVertexCollection> tvH; //For checking PU vertices
   iEvent.getByLabel("mix","MergedTrackTruth",tvH);
   TrackingVertexCollection tv = *tvH;
   //   std::cout<<"nr tracking vertices = "<<tv.size()<<std::endl;

   float primary_tv_dz = 1000; 
   float max_nr_daughters = 0;

   for (size_t j = 0; j < tv.size(); j++) { //loop over trackingVertices to find a leading-one
     if ( tv[j].sourceTracks().size() == 0 ){ // if no mother track (why > 1 present in case of Zee?)
       //       std::cout<<" Tracking vertex wih nr. source tracks = "<<tv[j].sourceTracks().size()<<", dxy = "<< sqrt( pow(tv[j].position().X(),2) + pow(tv[j].position().Y(), 2) ) <<", dz = "<< tv[j].position().Z() << ", involume = "<<tv[j].inVolume()<<", nr daughters = "<< tv[j].nDaughterTracks()<<std::endl; 

       
       if( tv[j].nDaughterTracks() > max_nr_daughters){ //pick the leading one (contains max number of daughter-tracks)
	 max_nr_daughters = tv[j].nDaughterTracks();
	 primary_tv_dz = tv[j].position().Z();
       }
     }

     if ( !(tv[j].eventId().bunchCrossing()== 0 && tv[j].eventId().event() == 0) ) //check whether there are PU vertices
       std::cout<< "Found PU vertex with bc = "<< tv[j].eventId().bunchCrossing() << std::endl;
   }
   
   // Select good primary vertices for use in subsequent track selection
   edm::Handle<reco::VertexCollection> hVtx;
   iEvent.getByLabel("pixelVertices", hVtx);
   //   std::cout<<"Nr primary vertices = " << hVtx->size()<<std::endl;

   std::vector<math::XYZPoint> PV_points;
   math::XYZPoint leading_PV_point;
   int leading_pv_index = 0;

   PV_points.clear();
   int vertex_count = 0;
   double dz_min = 1000;

   for (reco::VertexCollection::const_iterator it = hVtx->begin(); it != hVtx->end(); it++) {
     reco::Vertex vtx = *it;
     if( !(it->isFake()) && it->ndof() >=2 ){
       PV_points.push_back(it->position());
       //       std::cout<<"reco primary vertex with dz(leading tp, vertex) = "<< fabs(it->position().Z() - primary_tv_dz) <<", vertex nr = "<<vertex_count<<std::endl;
       
       if( fabs(it->position().Z() - primary_tv_dz) < dz_min){
	 dz_min = fabs(it->position().Z() - primary_tv_dz); //choose the reco PV that is closest to leading TrackingVertex
	 leading_pv_index = vertex_count;
       }
       
       vertex_count++;
     }
   }

   //   std::cout<<"dz min = "<<dz_min<<std::endl;
  
   if( PV_points.size() ) //write vertex that is closest in dz to the leading TrackingVertex
     leading_PV_point = PV_points[leading_pv_index];
   else if(is_single_part_) // accept all vertices for singleParticle samples later, here use placeholder
     leading_PV_point.SetXYZ(0,0,0);
   else{
     leading_PV_point.SetXYZ(1000,1000,1000);
     std::cout<< "WARNING! No primary vertex found in the event! Setting leading_PV_Point to a placeholder value. PV mathing will be skipped. " <<std::endl;
   }
     

   edm::Handle<SimTrackContainer>  simTrackCollection ; //simulated tracks
   iEvent.getByLabel("g4SimHits", simTrackCollection);

   SimTrackContainer simTC = *(simTrackCollection.product() );

   for (size_t j = 0; j <simTC.size(); j++){
     EncodedEventId simTrackEvt = simTC[j].eventId();
     if ( !(simTrackEvt.event() == 0 && simTrackEvt.bunchCrossing() == 0) )
       std::cout<<"Sim track evt Id = " <<simTrackEvt.event()<< ", bunchCrossing = " << simTrackEvt.bunchCrossing() <<std::endl;
   }

   std::map<int, std::vector<PSimHit>> simHitMap = getSimHits(trackerContainers, iEvent); // get all simHits in the event
   //   std::cout<<"nr tps = "<<tPCeff.size() << std::endl;

   np_gen_ = 0;
   np_gen_toReco_ = 0;
   for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
     TrackingParticleRef tpr(TPCollection, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
     
     TrackingParticle::Point vertexTP;
     TrackingParticle::Vector momentumTP;

     TrackingParticle::Point vertex = parametersDefinerTP->vertex(iEvent,iSetup,tpr);
     TrackingParticle::Vector momentum = parametersDefinerTP->momentum(iEvent,iSetup,tpr);
     
     if( !tpSelector(*tp) ) continue; // be careful with not adding preselections -- memory issues
     if( tp->pt() < 2 ) continue; // such tracks are irrelevant for electrons
     
     std::vector<PSimHit> simhits_TP = getSimHitsTP(tp, simHitMap, tracker, true, true, false); //last bool is for debug 
     if ( !simhits_TP.size() ) continue; // if no simhit with pt > 0.9 (FIXME: maybe check the last hit momentum too)
	 
     momentumTP = tp->momentum();
     vertexTP = tp->vertex(); 
	 
     gen_eta_[np_gen_] = tp->eta();
     gen_phi_[np_gen_] = tp->phi();
     gen_pt_[np_gen_] = tp->pt();
     gen_pdgId_[np_gen_] = tp->pdgId();
     gen_nr_simhits_[np_gen_] = tp->numberOfTrackerHits();
     gen_dxy_[np_gen_] = -vertex.x()*sin(momentum.phi()) + vertex.y()*cos(momentum.phi());
     gen_dz_[np_gen_] = vertex.z() - (vertex.x()*momentum.x()+vertex.y()*momentum.y())/sqrt(momentum.perp2())* momentum.z()/sqrt(momentum.perp2());	 


     if ( !(tp->eventId().bunchCrossing()== 0 && tp->eventId().event() == 0) )
       std::cout<<"found PU particle!"<<std::endl;
	 
     
     //	 std::cout<<"bunch crossing = "<<tp->eventId().bunchCrossing()<<", event = "<<tp->eventId().event() << std::endl; //0, 0 indicates primary vertex, the rest is PU
     //     std::cout<<"found tracking particle nr "<<np_gen_+1<<"with pt = "<<tp->pt() <<", eta = "<<tp->eta() << ", phi = "<<tp->phi() << "nr hits = " <<tp->numberOfTrackerHits() << ", charge = " << tp->charge()<<std::endl;
	 //	 if( abs(gen_dxy_[np_gen_]) < 0.5 )
	 //	 std::cout<<"pt TP = " << tp->pt() << ", dxy TP = "<<gen_dxy_[np_gen_] << "TP parent vertex z = "<<vertexTP.z()<<". vertex z = " << vertex.z()<<std::endl;

     
     for(std::vector<PSimHit>::const_iterator it_hit = simhits_TP.begin(); it_hit != simhits_TP.end(); it_hit++){ //get the fraction of brehmstrahlung at last hit
       if( it_hit+1 == simhits_TP.end() ){
	 GlobalVector global_p = getHitMomentum(it_hit,tracker, false); 
	 gen_ptAtLast_[np_gen_] = global_p.perp();
	 gen_bremFraction_[np_gen_] = 1 - global_p.perp()/tp->pt(); //fraction of pT lost via brehmstrahlung
	 //		 std::cout<<"pt at last = "<<global_p.perp()<<"TP pt = "<<tp->pt()<<std::endl;
       }
     }
	 
     //--------check for matched reco track defined by AssociatorByHits (efficiency and fake-rate plots)-----------
     const reco::Track* matchedTrackPointer=0;
     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;

     if(simRecCollSimDenom.find(tpr) != simRecCollSimDenom.end()){
       rt = simRecCollSimDenom[tpr];
       if (rt.size() != 0){
	 is_reco_matched_sim_[np_gen_] = 1;
	 simrec_simD_quality_[np_gen_] = rt.begin()->second;
       }
     }
     else{
       is_reco_matched_sim_[np_gen_] = 0;
       simrec_simD_quality_[np_gen_] = 0;
     }
     
     if(simRecCollRecoDenom.find(tpr) != simRecCollRecoDenom.end()){
       rt = simRecCollRecoDenom[tpr]; // find sim-to-reco association
       if ( rt.size()!=0 ) {
	 is_reco_matched_[np_gen_] = 1;
	 simrec_recD_quality_[np_gen_] = rt.begin()->second;

	 matchedTrackPointer = rt.begin()->first.get(); //pointer to corresponding reco track                                                 
	 	 
	 //	 std::cout<<"TP mathced with reco denom. nr shared hits = "<<nSharedHits<<", nr reco hits = "<<nRecoTrackHits <<std::endl;
	 //		 std::cout<<"matched with reco track with pt = " << matchedTrackPointer->pt() <<", eta = "<<matchedTrackPointer->eta() << ", phi = " << matchedTrackPointer->phi() << std::endl;

	 is_charge_matched_[np_gen_] = tp->charge()*matchedTrackPointer->charge();

	 //-------------- efficiencies-----------------
	 gen_matched_eta_[np_gen_] = tp->eta();
	 gen_matched_pt_[np_gen_] = tp->pt();
	 
	 //------------- for resolutions, additional gen vars  -------------------
	 gen_matched_phi_[np_gen_] = tp->phi();
	 gen_matched_qoverp_[np_gen_] = tp->charge()/(tp->px()*tp->px() + tp->py()*tp->py() + tp->pz()*tp->pz());
	 gen_matched_theta_[np_gen_] = tp->theta();
	 gen_matched_cotth_[np_gen_] = 1./tan(tp->theta());
	 gen_matched_dz_[np_gen_] = gen_dz_[np_gen_];
	 gen_matched_dxy_[np_gen_] =  gen_dxy_[np_gen_];
	 
	 //------------ for resolutions, gen matched reco vars --------------------
	 gen_matched_rec_nhits_[np_gen_] = matchedTrackPointer->numberOfValidHits();
	 gen_matched_rec_pt_[np_gen_] = matchedTrackPointer->pt();
	 gen_matched_rec_theta_[np_gen_] = matchedTrackPointer->theta();
	 gen_matched_rec_eta_[np_gen_] = matchedTrackPointer->eta();
	 gen_matched_rec_qoverp_[np_gen_] = matchedTrackPointer->qoverp();
	 gen_matched_rec_cotth_[np_gen_] = 1./(tan(matchedTrackPointer->theta()) );
	 gen_matched_rec_phi_[np_gen_] = matchedTrackPointer->phi();
	 
	 gen_matched_rec_dz_[np_gen_] = matchedTrackPointer->dz(bsPosition);
	 gen_matched_rec_dxy_[np_gen_] = matchedTrackPointer->dxy(bsPosition);
	 
	 //------------------- pulls ----------------------------------------
	 pt_pull_[np_gen_] = (matchedTrackPointer->pt() - tp->pt())/matchedTrackPointer->ptError();
	 qoverp_pull_[np_gen_] = (matchedTrackPointer->qoverp() - gen_matched_qoverp_[np_gen_])/matchedTrackPointer->qoverpError();
	 theta_pull_[np_gen_] = (matchedTrackPointer->theta() - tp->theta())/matchedTrackPointer->thetaError();
	 phi_pull_[np_gen_] = (matchedTrackPointer->phi() - tp->phi())/matchedTrackPointer->phiError();
	 dz_pull_[np_gen_] = (matchedTrackPointer->dz(bsPosition) - gen_dz_[np_gen_])/matchedTrackPointer->dzError();
	 dxy_pull_[np_gen_] = (matchedTrackPointer->dxy(bsPosition) - gen_dxy_[np_gen_])/matchedTrackPointer->dxyError();
	 
	 np_gen_toReco_++; //count reco matched TPs
	 //	 std::cout<<"TP matched with RECO track"<<std::endl;
       }
     }
     else{
       is_reco_matched_[np_gen_] = 0;
       simrec_recD_quality_[np_gen_] = 0;
     }

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
     
     //     std::cout<<"Matched seed charge = "<<gen_matched_seed_okCharge_[np_gen_]<<", match quality = "<< gen_matched_seed_quality_[np_gen_]<<", nr matched seed hits = "<< gen_matched_seed_nshared_[np_gen_]<<std::endl;     
     
     //	 std::cout<<"Found tracking particle with pt = "<<tp->pt()<<", eta = "<<tp->eta()<<", phi = " << tp->phi()<<std::endl; 
     // std::cout<<"nr shared hits = "<<nSharedHits<<", nr reco hits "<<nRecoTrackHits<<" matched reco pt = "<<gen_matched_rec_pt_[np_gen_]<< std::endl;
     //", gen. qoverp = "<<gen_matched_qoverp_[np_gen_]<<", matched reco qoverp = "<<gen_matched_rec_qoverp_[np_gen_]<<std::endl; 
     
     np_gen_++; // count tracking particles that pass the standard quality cuts
   } // <-- end loop over tracking particles

   //   std::cout<<"NR. SEL TP-s = " <<np_gen_<<std::endl;

   //------------------------ loop over reco tracks --------------------------
   np_reco_ = 0;
   np_reco_toGen_ = 0;
   np_fake_ = 0;

   for( edm::View<reco::Track>::size_type i=0; i<trackCollection->size(); i++){
     edm::RefToBase<reco::Track> track(trackCollection, i);
	 
     if( track->pt() < 2 ) continue; // reduce noise, irrelevant for electrons
     if (!isGoodTrack(track, bsPosition, PV_points, false) ) continue;
     
     math::XYZPoint PV_point = getClosestPVpoint(track, PV_points); //pick the vertex with closest dz from the track
     // std::cout<<"LEADING PV point dz distance to track = " <<track->dz(leading_PV_point) <<std::endl;
     //std::cout<<"Closest PV point dz distance to track = " <<track->dz(PV_point) <<std::endl;

     if( is_single_part_ || fabs( track->dz(leading_PV_point) - track->dz(PV_point) ) < 0.001 ){ //choose tracks from leading vertex (workaround for missing PU tracking particles)

       //std::cout<<"Found reco track with pt = " << track->pt() << " eta = " << track->eta() << ", dz leadPV = "<<track->dz(leading_PV_point) <<std::endl;
     
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

	   // if ( fabs( track->dz(leading_PV_point) - track->dz(PV_point) ) >= 0.001  ) //FIXME understand why these happen 
	   //  std::cout << "NO PV MATCH at dz_lead = " << track->dz(leading_PV_point) << ", dz closest = " << track->dz(PV_point) <<std::endl;
	   //  std::cout<<"MATCHED!"<<std::endl;
	 }
       }

       else{  // if fake
	 is_gen_matched_[np_reco_] = 0;
	 
	 fake_pt_[np_fake_] = track->pt();
	 fake_phi_[np_fake_] = track->phi();
	 fake_nr_rechits_[np_fake_] = track->numberOfValidHits();
	 fake_dxy_[np_fake_] = track->dxy(bsPosition);
	 fake_dz_[np_fake_] = track->dz(bsPosition);
       	 fake_eta_[np_fake_] = track->eta();
	 
	 np_fake_++;
       }
       np_reco_++;     
     } // <-- end if reco track is matched to leading vertex

     //     np_reco_++;     //moved within if for consistent indexing
   } // <-- end loop over reco tracks

   //   for( edm::View<reco::ElectronSeed>::const_iterator it_seed = elSeedCollection->begin(); it_seed != elSeedCollection->end(); it_seed++){
   //  TrajectorySeed::range rechits = it_seed->recHits(); //get recHits associated to the seed
   // } //FIXME add fake seeds


   trackValTree_->Fill();
   simHitMap.clear();	 
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

  trackValTree_->Branch("fake_dxy", fake_dxy_, "fake_dxy[np_fake]/D");
  trackValTree_->Branch("fake_dz", fake_dz_, "fake_dz[np_fake]/D");
  trackValTree_->Branch("fake_nr_rechits", fake_nr_rechits_, "fake_nr_rechits[np_fake]/D");

  // ----------------- TPs ----------------------------
  trackValTree_->Branch("np_gen", &np_gen_, "np_gen/I");
  trackValTree_->Branch("np_gen_toReco", &np_gen_toReco_, "np_gen_toReco/I");

  trackValTree_->Branch("gen_reco_matched", is_reco_matched_, "gen_reco_matched[np_gen]/I");
  trackValTree_->Branch("gen_reco_matched_sim", is_reco_matched_sim_, "gen_reco_matched_sim[np_gen]/I");
  trackValTree_->Branch("charge_matched", is_charge_matched_, "charge_matched[np_gen]/I");

  trackValTree_->Branch("gen_pdgId", gen_pdgId_, "gen_pdgId[np_gen]/I");
  trackValTree_->Branch("gen_pt", gen_pt_, "gen_pt[np_gen]/D");
  trackValTree_->Branch("gen_eta", gen_eta_, "gen_eta[np_gen]/D");
  trackValTree_->Branch("gen_phi", gen_phi_, "gen_phi[np_gen]/D");
  trackValTree_->Branch("gen_nr_simhits", gen_nr_simhits_, "gen_nr_simhits_[np_gen]/I");
  trackValTree_->Branch("gen_dxy", gen_dxy_, "gen_dxy[np_gen]/D");
  trackValTree_->Branch("gen_dz", gen_dz_, "gen_dz[np_gen]/D");
  trackValTree_->Branch("gen_ptAtLast", gen_ptAtLast_, "gen_ptAtLast[np_gen]/D");
  trackValTree_->Branch("gen_bremFraction", gen_bremFraction_, "gen_bremFraction[np_gen]/D");
  
  //----------------- TP simToReco matching-----------------
  trackValTree_->Branch("gen_matched_pt", gen_matched_pt_, "gen_matched_pt[np_gen]/D");
  trackValTree_->Branch("gen_matched_qoverp", gen_matched_qoverp_, "gen_matched_qoverp[np_gen]/D");
  trackValTree_->Branch("gen_matched_eta", gen_matched_eta_, "gen_matched_eta[np_gen]/D");
  trackValTree_->Branch("gen_matched_theta", gen_matched_theta_, "gen_matched_theta[np_gen]/D");
  trackValTree_->Branch("gen_matched_cotth", gen_matched_cotth_, "gen_matched_cotth[np_gen]/D");
  trackValTree_->Branch("gen_matched_phi", gen_matched_phi_, "gen_matched_phi[np_gen]/D");
  trackValTree_->Branch("gen_matched_dxy", gen_matched_dxy_, "gen_matched_dxy[np_gen]/D");
  trackValTree_->Branch("gen_matched_dz", gen_matched_dz_, "gen_matched_dz[np_gen]/D");
  
  trackValTree_->Branch("gen_matched_rec_pt", gen_matched_rec_pt_, "gen_matched_rec_pt[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_eta", gen_matched_rec_eta_, "gen_matched_rec_eta[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_qoverp", gen_matched_rec_qoverp_, "gen_matched_rec_qoverp[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_theta", gen_matched_rec_theta_, "gen_matched_rec_theta[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_cotth", gen_matched_rec_cotth_, "gen_matched_rec_cotth[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_phi", gen_matched_rec_phi_, "gen_matched_rec_phi[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_dxy", gen_matched_rec_dxy_, "gen_matched_rec_dxy[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_dz", gen_matched_rec_dz_, "gen_matched_rec_dz[np_gen]/D");
  trackValTree_->Branch("gen_matched_rec_nhits", gen_matched_rec_nhits_, "gen_matched_rec_nhits[np_gen]/I");

  // ------------------- pulls ---------------------------
  trackValTree_->Branch("pt_pull", pt_pull_, "pt_pull[np_gen]/D");
  trackValTree_->Branch("qoverp_pull", qoverp_pull_, "qoverp_pull[np_gen]/D");
  trackValTree_->Branch("theta_pull", theta_pull_, "theta_pull[np_gen]/D");
  trackValTree_->Branch("phi_pull", phi_pull_, "phi_pull[np_gen]/D");
  trackValTree_->Branch("dxy_pull", dxy_pull_, "dxy_pull[np_gen]/D");
  trackValTree_->Branch("dz_pull", dz_pull_, "dz_pull[np_gen]/D");
  
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

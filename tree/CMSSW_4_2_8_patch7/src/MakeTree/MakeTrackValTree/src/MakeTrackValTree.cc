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

      // ----------member data ---------------------------
  double ptMinTP, minRapidityTP, maxRapidityTP, tipTP, lipTP, signalOnlyTP, chargedOnlyTP, stableOnlyTP, minAbsEtaTP, maxAbsEtaTP;
  std::vector<int> pdgIdTP;
  int minHitTP;

  TTree *trackValTree_;

  int np_gen_, np_gen_toReco_, np_reco_, np_reco_toGen_, np_fake_;
  int is_reco_matched_[MAXPART], is_gen_matched_[MAXPART];
  int gen_passhit075_[MAXPART][MAXHIT], gen_passhit080_[MAXPART][MAXHIT];
  float gen_hit_pt_[MAXPART][MAXHIT];

  double gen_pt_[MAXPART], gen_eta_[MAXPART], gen_phi_[MAXPART], gen_matched_pt_[MAXPART], gen_matched_eta_[MAXPART], gen_matched_phi_[MAXPART], gen_dxy_[MAXPART], gen_dz_[MAXPART], gen_passhit3_pt_[MAXPART], gen_passhit3_eta_[MAXPART], gen_passlast_pt_[MAXPART], gen_passlast_eta_[MAXPART];
  double reco_pt_[MAXPART], reco_eta_[MAXPART], reco_phi_[MAXPART], fake_pt_[MAXPART], fake_eta_[MAXPART], fake_phi_[MAXPART];
  int gen_pdgId_[MAXPART], gen_nrSimHits_[MAXPART], gen_nrRecoHits_[MAXPART], reco_nrRecoHits_[MAXPART], reco_nrSimHits_[MAXPART];
  double reco_nrSharedHits_[MAXPART], gen_nrSharedHits_[MAXPART];
  
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
    
    gen_passhit3_pt_[i] = -999;
    gen_passhit3_eta_[i] = -999;
    gen_passlast_pt_[i] = -999;
    gen_passlast_eta_[i] = -999;

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
      gen_passhit075_[i][h] = -1;
      gen_passhit080_[i][h] = -1;
    }
  }

  edm::Handle<edm::View<reco::Track> >   trackCollection; //reconstructed tracks
  iEvent.getByLabel("cutsRecoTracksHp", trackCollection);
  //  iEvent.getByLabel("generalTracks", trackCollection);

  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ; //simulated tracks
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());
  if(debug_)
    std::cout <<"Number of simulated tracks = "<<tPCeff.size() << std::endl;

  edm::ESHandle<TrackAssociatorBase> theAssociator; //create track associators for MC truth matching
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits" ,theAssociator);

  reco::SimToRecoCollection simRecColl = theAssociator->associateSimToReco(trackCollection, TPCollectionHeff, &iEvent);
  reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(trackCollection, TPCollectionHeff, &iEvent);

  edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
  iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP", parametersDefinerTP);   

  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);


  /*
  edm::Handle<edm::PSimHitContainer> hitsPxlBrlLow;
  iEvent.getByLabel("g4SimHits","TrackerHitsPixelBarrelLowTof", hitsPxlBrlLow);
  
    if (!hitsPxlBrlLow.isValid()) {
    edm::LogError("MakeTrackValTree::analyze")
      << "Unable to find TrackerHitsPixelBarrelLowTof in event!";
    return;                                             
  }
  */

  TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,               
								 signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP,minAbsEtaTP,maxAbsEtaTP);  
  
  //---------------------Loop over simulated hits------------------------------
  edm::PSimHitContainer::const_iterator itHit;

  /*
  //pixel barrel low container
  for (itHit = hitsPxlBrlLow->begin(); itHit != hitsPxlBrlLow->end(); ++itHit){
    DetId detid=DetId(itHit->detUnitId());
    std::cout<<"DetId = "<<detid<<std::endl;
  }
  */


  //----------------------Loop over sim tracks--------------------------------
  np_gen_ = 0;
  np_gen_toReco_ = 0;


  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

    if( !tpSelector(*tp) ) continue;
    if(tp->pt() < 1) continue; // reduce noise. such tracks are irrelevant for electrons - never makes it to ECAL

    ParticleBase::Point vertex = parametersDefinerTP->vertex(iEvent,iSetup,*tp);
    gen_dxy_[np_gen_] = -vertex.x()*sin(tp->phi())+vertex.y()*cos(tp->phi()); 
    gen_dz_[np_gen_] = vertex.z()-(vertex.x()*tp->px()+vertex.y()*tp->py())/tp->pt()*tp->pz()/tp->pt();

    gen_eta_[np_gen_] = tp->eta();
    gen_phi_[np_gen_] = tp->phi();
    gen_pt_[np_gen_] = tp->pt();
    gen_pdgId_[np_gen_] = tp->pdgId();

    std::vector<PSimHit> simhits=tp->trackPSimHit(DetId::Tracker);        
    double nr_simhits = simhits.size();
    gen_nrSimHits_[np_gen_] =  nr_simhits;

  

    
    bool hitdebug = false;
    if(hitdebug)
      std::cout<<"TRACKING PARTICLE with pt = "<<tp->pt()<<std::endl;

    unsigned int hit_count = 0;
    for(std::vector<PSimHit>::const_iterator it_hit = simhits.begin(); it_hit != simhits.end(); it_hit++){

      DetId dId = DetId(it_hit->detUnitId() );
      unsigned int subdetId = static_cast<unsigned int>(dId.subdetId());
      int layerNumber = 0;
      if(subdetId == PixelSubdetector::PixelBarrel){
	PXBDetId pxbid(dId.rawId());
	layerNumber = pxbid.layer();
	if(hitdebug)
	  std::cout<<"Hit at pixel barrel layer = "<<layerNumber<<std::endl;
      }
      else if(subdetId == PixelSubdetector::PixelEndcap){
	PXFDetId pxfid(dId.rawId());
	layerNumber = pxfid.disk();
	if(hitdebug)
	  std::cout<<"Hit at pixel endcap layer = "<<layerNumber<<std::endl;
      }
      
      else if( subdetId == StripSubdetector::TIB){
	TIBDetId tibid(dId.rawId());
	layerNumber = tibid.layer();
	if(hitdebug)
	  std::cout<<"Hit at strip TIB layer = "<<layerNumber<<std::endl;
      }
    
      else if( subdetId == StripSubdetector::TOB){
	TOBDetId tobid(dId.rawId());
	layerNumber = tobid.layer();
	if(hitdebug)
	  std::cout<<"Hit at strip TOB layer = "<<layerNumber<<std::endl;
      }
      else if( subdetId == StripSubdetector::TID){
	TIDDetId tidid(dId.rawId());
	layerNumber = tidid.wheel();
	if(hitdebug)
	  std::cout<<"Hit at strip TID layer = "<<layerNumber<<std::endl;
      }
      else if( subdetId == StripSubdetector::TEC ){
	TECDetId tecid(dId.rawId());
	layerNumber = tecid.wheel();
	if(hitdebug)
	  std::cout<<"Hit at strip TEC layer = "<<layerNumber<<std::endl;
      }

      LocalVector local_p = it_hit->momentumAtEntry();
      const GeomDetUnit* detunit = tracker->idToDetUnit(dId.rawId());
      GlobalVector global_p = detunit->toGlobal(local_p);
      float pt_at_entry = global_p.perp();
      if(hitdebug)
	std::cout<<"momentum at entry = "<<pt_at_entry<<std::endl;

      float hit_pt_eff = pt_at_entry/tp->pt();

      for(unsigned int hit_nr = 0; hit_nr !=nr_simhits; hit_nr++){
	if( hit_nr == hit_count) { //check the track pt efficiency at each track hit
	  pt_at_entry = gen_hit_pt_[np_gen_][hit_nr];

	  //-----------efficiencies, possibly unnecessary----------
	  if( hit_pt_eff > 0.75 )
	    gen_passhit075_[np_gen_][hit_nr] = 1; //pass	  
	  else
	    gen_passhit075_[np_gen_][hit_nr] = 0; //hit exists, but fails efficiency
	  
	  if( hit_pt_eff > 0.8 )
	    gen_passhit080_[np_gen_][hit_nr] = 1;
	  else
	    gen_passhit080_[np_gen_][hit_nr] = 0;
	  //-----------------------------------------------
	  break; //exit the loop after finding the hit of interest
	}
      }

      std::cout<<"hit efficiency at hit "<<hit_count<<" = "<<hit_pt_eff<<std::endl;
      hit_count++;
    }
    std::cout<<"Found track with pt"<<tp->pt()<<std::endl;
    std::cout<<"hit 0 = "<<gen_hit_pt_[np_gen_][0]<<std::endl;
    std::cout<<"hit 1 = "<<gen_hit_pt_[np_gen_][1]<<std::endl;
    

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
    
      //--------------------check the hit pattern-----------------
      /*if(true){ //add condition for fake tracks later
	const reco::HitPattern& p = track->hitPattern();
	if(debug_)
	  std::cout<<"-----------------trak hit pattern information for track nr "<<i<<"--------------"<<std::endl;      
	for(int i=0; i<p.numberOfHits(); i++){
	  int hit = p.getHitPattern(i);
	  if(p.validHitFilter(hit) && p.pixelHitFilter(hit) ){
	    if(debug_)
	      std::cout<<"found hit in pixel detector layer "<<p.getLayer(hit)<<" from "<<hit<<std::endl;
	  }
	}
	if(debug_)
	  std::cout<<"-----------------------------------------------------"<<std::endl;
      }*/

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

  //  trackValTree_->Branch("gen_passhit3_pt", gen_passhit3_pt_, "gen_passhit3_pt[np_gen]/D");
  //trackValTree_->Branch("gen_passhit3_eta", gen_passhit3_eta_, "gen_passhit3_eta[np_gen]/D");
  //trackValTree_->Branch("gen_passlast_pt", gen_passlast_pt_, "gen_passlast_pt[np_gen]/D");
  //trackValTree_->Branch("gen_passlast_eta", gen_passlast_eta_, "gen_passlast_eta[np_gen]/D");
  
  trackValTree_->Branch("gen_matched_pt", gen_matched_pt_, "gen_matched_pt[np_gen]/D");
  trackValTree_->Branch("gen_matched_eta", gen_matched_eta_, "gen_matched_eta[np_gen]/D");
  trackValTree_->Branch("gen_matched_phi", gen_matched_phi_, "gen_matched_phi[np_gen]/D");

  trackValTree_->Branch("gen_hit_pt", gen_hit_pt_, "gen_hit_pt[np_gen][50]/F");
  trackValTree_->Branch("gen_passhit075", gen_passhit075_, "gen_passhit075_hit1[np_gen][50]/I"); //only one variable dimension allowed, sync with MAXHIT
  trackValTree_->Branch("gen_passhit080", gen_passhit080_, "gen_passhit080_hit1[np_gen][50]/I");
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



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

#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"


//
// class declaration
//

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

   edm::InputTag track_label; 
   if(is_gsf_)
     track_label = track_label_gsf_;
   else
     track_label = track_label_;

   edm::Handle<edm::View<reco::Track> > trackCollection; //reconstructed tracks
   iEvent.getByLabel(track_label, trackCollection);

   edm::Handle<edm::View<reco::Track> > selTrackCollection; //reconstructed tracks with some preselection requirements
   iEvent.getByLabel("elGsfTracksWithQuality", selTrackCollection);
   
   edm::Handle<edm::View<reco::ElectronSeed> > elSeedCollection; 
   iEvent.getByLabel(el_seed_label_, elSeedCollection);

   edm::Handle<TrackingParticleCollection>  TPCollection ; //simulated tracks
   iEvent.getByLabel("mix","MergedTrackTruth",TPCollection);


   std::cout<<" nr. electron seeds: " << elSeedCollection->size()<<std::endl;
   std::cout<<" nr. tracks: " << trackCollection->size()<<std::endl;
   std::cout<<" nr. tracking particles: " << TPCollection->size() <<std::endl;

   //   std::cout<<" nr. GSF tracks: " << selTrackCollection->size()<<std::endl;

   /*   int sn = 0;
   for( edm::View<reco::ElectronSeed>::const_iterator it_seed = elSeedCollection->begin(); it_seed != elSeedCollection->end(); it_seed++, sn++){
      
     std::cout<<"Electron seed # "<<sn<<" with numbr of hits = "<<it_seed->nHits()<<std::endl;

     }*/

   //   const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());
   // TrackingParticleSelector tpSelector; // = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP,minAbsEtaTP,maxAbsEtaTP);
   
   /*   for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
     TrackingParticleRef tpr(TPCollectionHeff, i);

     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

     //     if( !tpSelector(*tp) ) continue;
     if( tp->pt() < 1 ) continue;

     std::cout<<"Found generated track with pt = "<<tp->pt()<<", eta = " << tp->eta()<<std::endl;
     }*/

}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTrackValTree::beginJob()
{
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

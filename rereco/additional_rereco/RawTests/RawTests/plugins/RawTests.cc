// -*- C++ -*-
//
// Package:    RawTests/RawTests
// Class:      RawTests
// 
/**\class RawTests RawTests.cc RawTests/RawTests/plugins/RawTests.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Liis Rebane (ETHZ) [liis]
//         Created:  Mon, 29 Sep 2014 12:39:58 GMT
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
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"


//
// class declaration
//

class RawTests : public edm::EDAnalyzer {
   public:
      explicit RawTests(const edm::ParameterSet&);
      ~RawTests();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
RawTests::RawTests(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


RawTests::~RawTests()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RawTests::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   // ------ check whether SimTracks from Pileup events are present in the events ----------
   edm::Handle<SimTrackContainer>  simTrackCollection ; //simulated tracks                                                      
   iEvent.getByLabel("g4SimHits", simTrackCollection);
   SimTrackContainer simTC = *(simTrackCollection.product() );
   
   // for (size_t j = 0; j <simTC.size(); j++){
     //   EncodedEventId simTrackEvt = simTC[j].eventId();
     //     if ( !(simTrackEvt.event() == 0 && simTrackEvt.bunchCrossing() == 0) )
     //     std::cout<<"Sim track evt Id = " <<simTrackEvt.event()<< ", bunchCrossing = " << simTrackEvt.bunchCrossing() <<std::endl;
   //}

   

   edm::Handle<TrackingParticleCollection>  TPCollection ; //simulated tracks
   iEvent.getByLabel("mix","MergedTrackTruth",TPCollection);
   const TrackingParticleCollection tPCeff = *(TPCollection.product());

   for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from simulated  tracks   
     TrackingParticleRef tpr(TPCollection, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

     //     if ( (tp->eventId().bunchCrossing()== 0 && tp->eventId().event() == 0) )
       std::cout<<"found PU particle with evtID = "<<tp->eventId().event()<<", bc = "<<tp->eventId().bunchCrossing()<<std::endl;

   }


   edm::Handle<TrackingVertexCollection> tvH; //For checking PU vertices
   iEvent.getByLabel("mix","MergedTrackTruth",tvH);
   TrackingVertexCollection tv = *tvH;

   /*
   for (size_t j = 0; j < tv.size(); j++) { //loop over trackingVertices to find a leading-one
     std::cout<< "Found PU vertex with evt = " << tv[j].eventId().event()<<", bc = "<< tv[j].eventId().bunchCrossing() << std::endl;
   }
   */

   


}


// ------------ method called once each job just before starting event loop  ------------
void 
RawTests::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RawTests::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
RawTests::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
RawTests::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
RawTests::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
RawTests::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RawTests::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RawTests);

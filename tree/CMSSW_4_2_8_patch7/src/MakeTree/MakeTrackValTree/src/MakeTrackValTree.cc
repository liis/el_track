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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
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

  TTree *trackValTree_;

  int nr_gen_trk_, nr_reco_trk_;
  double gen_pt_, gen_eta_;

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

  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ; //simulated tracks
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());

  nr_gen_trk_ = tPCeff.size();
  
  if(debug_)
    std::cout <<"Number of simulated tracks = "<<tPCeff.size() << std::endl;

  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from get tracks                                      
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

    gen_pt_ = tp->pt();
    gen_eta_ = tp->eta();

    //    std::cout<<"gen track pT = "<<(tp->pt())<<", eta = "<< tp->eta()<<std::endl;                                                           
    //if( tpSelector(*tp) ) nr_sel_simtracks++;
    
    trackValTree_->Fill();

  }



}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTrackValTree::beginJob()
{
  edm::Service<TFileService> fs;
  trackValTree_ = fs->make<TTree>("trackValTree","trackValTree");
  
  trackValTree_->Branch("nr_reco_trk", &nr_reco_trk_, "nr_reco_trk/I");

  trackValTree_->Branch("nr_gen_trk", &nr_gen_trk_, "nr_gen_trk/I");
  trackValTree_->Branch("gen_pt", &gen_pt_, "gen_pt_/D");
  trackValTree_->Branch("gen_eta", &gen_eta_, "gen_eta_/D");
  
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

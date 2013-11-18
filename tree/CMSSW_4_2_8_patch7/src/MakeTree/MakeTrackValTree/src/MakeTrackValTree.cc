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

#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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

  int np_gen_;
  int nr_gen_trk_, nr_reco_trk_;
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
  }


  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ; //simulated tracks
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());


  
  if(debug_)
    std::cout <<"Number of simulated tracks = "<<tPCeff.size() << std::endl;

  TrackingParticleSelector tpSelector = TrackingParticleSelector(ptMinTP, minRapidityTP, maxRapidityTP,tipTP, lipTP, minHitTP,               
								 signalOnlyTP, chargedOnlyTP, stableOnlyTP,pdgIdTP,minAbsEtaTP,maxAbsEtaTP);  
  np_gen_ = 0;
  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ // get information from get tracks   
    TrackingParticleRef tpr(TPCollectionHeff, i);
    TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

    if( tpSelector(*tp) ){
      gen_eta_[np_gen_] = tp->eta();
      gen_pt_[np_gen_] = tp->pt();

      np_gen_++;
    }

  }
  
  if(debug_)
    std::cout<<"nr of selected gen tracks = "<<np_gen_<<std::endl;
  
  trackValTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTrackValTree::beginJob()
{
  edm::Service<TFileService> fs;
  trackValTree_ = fs->make<TTree>("trackValTree","trackValTree");
 
  trackValTree_->Branch("np_reco_", &nr_reco_trk_, "np_reco_/I");

  trackValTree_->Branch("np_gen",&np_gen_,"np_gen/I"); // needs to be filled in order to fill x[np]
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



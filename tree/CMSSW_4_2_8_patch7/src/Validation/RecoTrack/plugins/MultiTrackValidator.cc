#include "Validation/RecoTrack/interface/MultiTrackValidator.h"
#include "DQMServices/ClientConfig/interface/FitSlicesYTool.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/plugins/CosmicParametersDefinerForTPESProducer.h"
#include "Validation/RecoTrack/interface/MTVHistoProducerAlgoFactory.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Ref.h"

#include "TMath.h"
#include <TF1.h>

using namespace std;
using namespace edm;

MultiTrackValidator::MultiTrackValidator(const edm::ParameterSet& pset):MultiTrackValidatorBase(pset){
  //theExtractor = IsoDepositExtractorFactory::get()->create( extractorName, extractorPSet);

  ParameterSet psetForHistoProducerAlgo = pset.getParameter<ParameterSet>("histoProducerAlgoBlock");
  string histoProducerAlgoName = psetForHistoProducerAlgo.getParameter<string>("ComponentName");
  histoProducerAlgo_ = MTVHistoProducerAlgoFactory::get()->create(histoProducerAlgoName ,psetForHistoProducerAlgo);
  histoProducerAlgo_->setDQMStore(dbe_);

  dirName_ = pset.getParameter<std::string>("dirName");
  associatormap = pset.getParameter< edm::InputTag >("associatormap");
  UseAssociators = pset.getParameter< bool >("UseAssociators");

  m_dEdx1Tag = pset.getParameter< edm::InputTag >("dEdx1Tag");
  m_dEdx2Tag = pset.getParameter< edm::InputTag >("dEdx2Tag");


  tpSelector = TrackingParticleSelector(pset.getParameter<double>("ptMinTP"),
					pset.getParameter<double>("minRapidityTP"),
					pset.getParameter<double>("maxRapidityTP"),
					pset.getParameter<double>("tipTP"),
					pset.getParameter<double>("lipTP"),
					pset.getParameter<int>("minHitTP"),
					pset.getParameter<bool>("signalOnlyTP"),
					pset.getParameter<bool>("chargedOnlyTP"),
					pset.getParameter<bool>("stableOnlyTP"),
					pset.getParameter<std::vector<int> >("pdgIdTP"),
					pset.getParameter<double>("minAbsEtaTP"),
					pset.getParameter<double>("maxAbsEtaTP"));
					



  cosmictpSelector = CosmicTrackingParticleSelector(pset.getParameter<double>("ptMinTP"),
						    pset.getParameter<double>("minRapidityTP"),
						    pset.getParameter<double>("maxRapidityTP"),
						    pset.getParameter<double>("tipTP"),
						    pset.getParameter<double>("lipTP"),
						    pset.getParameter<int>("minHitTP"),
						    pset.getParameter<bool>("chargedOnlyTP"),
						    pset.getParameter<std::vector<int> >("pdgIdTP"));

  runStandalone = pset.getParameter<bool>("runStandalone");

    
  if (!UseAssociators) {
    associators.clear();
    associators.push_back(associatormap.label());
  }
}


MultiTrackValidator::~MultiTrackValidator(){delete histoProducerAlgo_;}

void MultiTrackValidator::beginRun(Run const&, EventSetup const& setup) {
  //  dbe_->showDirStructure();

  //int j=0;  //is This Necessary ???
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){
      dbe_->cd();
      InputTag algo = label[www];
      string dirName=dirName_;
      if (algo.process()!="")
	dirName+=algo.process()+"_";
      if(algo.label()!="")
	dirName+=algo.label()+"_";
      if(algo.instance()!="")
	dirName+=algo.instance()+"_";      
      if (dirName.find("Tracks")<dirName.length()){
	dirName.replace(dirName.find("Tracks"),6,"");
      }
      string assoc= associators[ww];
      if (assoc.find("Track")<assoc.length()){
	assoc.replace(assoc.find("Track"),5,"");
      }
      dirName+=assoc;
      std::replace(dirName.begin(), dirName.end(), ':', '_');

      dbe_->setCurrentFolder(dirName.c_str());
      
      // vector of vector initialization
      histoProducerAlgo_->initialize(); //TO BE FIXED. I'D LIKE TO AVOID THIS CALL

      dbe_->goUp(); //Is this really necessary ???
      string subDirName = dirName + "/simulation";
      dbe_->setCurrentFolder(subDirName.c_str());

      //Booking histograms concerning with simulated tracks
      histoProducerAlgo_->bookSimHistos();

      dbe_->cd();
      dbe_->setCurrentFolder(dirName.c_str());

      //Booking histograms concerning with reconstructed tracks
      histoProducerAlgo_->bookRecoHistos();
      if (runStandalone) histoProducerAlgo_->bookRecoHistosForStandaloneRunning();

      if (UseAssociators) {
	edm::ESHandle<TrackAssociatorBase> theAssociator;
	for (unsigned int w=0;w<associators.size();w++) {
	  setup.get<TrackAssociatorRecord>().get(associators[w],theAssociator);
	  associator.push_back( theAssociator.product() );
	}//end loop w
      }
    }//end loop www
  }// end loop ww
}

void MultiTrackValidator::analyze(const edm::Event& event, const edm::EventSetup& setup){
  using namespace reco;
  
    

  edm::LogInfo("TrackValidator") << "\n====================================================" << "\n"
				 << "Analyzing new event" << "\n"
				 << "====================================================\n" << "\n";
  edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
  setup.get<TrackAssociatorRecord>().get(parametersDefiner,parametersDefinerTP);    
  
  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ;
  event.getByLabel(label_tp_effic,TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());
  
  edm::Handle<TrackingParticleCollection>  TPCollectionHfake ;
  event.getByLabel(label_tp_fake,TPCollectionHfake);
  const TrackingParticleCollection tPCfake = *(TPCollectionHfake.product());
  
  //if (tPCeff.size()==0) {edm::LogInfo("TrackValidator") 
  //<< "TP Collection for efficiency studies has size = 0! Skipping Event." ; return;}
  //if (tPCfake.size()==0) {edm::LogInfo("TrackValidator") 
  //<< "TP Collection for fake rate studies has size = 0! Skipping Event." ; return;}
  
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  event.getByLabel(bsSrc,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;      
  
  int w=0; //counter counting the number of sets of histograms
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){

      //
      //get collections from the event
      //
      edm::Handle<View<Track> >  trackCollection;
      if(!event.getByLabel(label[www], trackCollection)&&ignoremissingtkcollection_)continue;
      //if (trackCollection->size()==0) {
      //edm::LogInfo("TrackValidator") << "TrackCollection size = 0!" ; 
      //continue;
      //}
      reco::RecoToSimCollection recSimColl;
      reco::SimToRecoCollection simRecColl;
      //associate tracks
      if(UseAssociators){
	edm::LogVerbatim("TrackValidator") << "Analyzing " 
					   << label[www].process()<<":"
					   << label[www].label()<<":"
					   << label[www].instance()<<" with "
					   << associators[ww].c_str() <<"\n";
	
	LogTrace("TrackValidator") << "Calling associateRecoToSim method" << "\n";
	recSimColl=associator[ww]->associateRecoToSim(trackCollection,
						      TPCollectionHfake,
						      &event);
	LogTrace("TrackValidator") << "Calling associateSimToReco method" << "\n";
	simRecColl=associator[ww]->associateSimToReco(trackCollection,
						      TPCollectionHeff, 
						      &event);
      }
      else{
	edm::LogVerbatim("TrackValidator") << "Analyzing " 
					   << label[www].process()<<":"
					   << label[www].label()<<":"
					   << label[www].instance()<<" with "
					   << associatormap.process()<<":"
					   << associatormap.label()<<":"
					   << associatormap.instance()<<"\n";
	
	Handle<reco::SimToRecoCollection > simtorecoCollectionH;
	event.getByLabel(associatormap,simtorecoCollectionH);
	simRecColl= *(simtorecoCollectionH.product()); 
	
	Handle<reco::RecoToSimCollection > recotosimCollectionH;
	event.getByLabel(associatormap,recotosimCollectionH);
	recSimColl= *(recotosimCollectionH.product()); 
      }


      // ########################################################
      // fill simulation histograms (LOOP OVER TRACKINGPARTICLES)
      // ########################################################

      //compute number of tracks per eta interval
      //
      edm::LogVerbatim("TrackValidator") << "\n# of TrackingParticles: " << tPCeff.size() << "\n";
      int ats(0);  //This counter counts the number of simTracks that are "associated" to recoTracks
      int st(0);   //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )
      for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){ //loop over TPs collection for tracking efficiency
	TrackingParticleRef tpr(TPCollectionHeff, i);
	TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
	ParticleBase::Vector momentumTP; 
	ParticleBase::Point vertexTP;
	double dxySim(0);
	double dzSim(0); 

	
	//---------- THIS PART HAS TO BE CLEANED UP. THE PARAMETER DEFINER WAS NOT MEANT TO BE USED IN THIS WAY ----------
	//If the TrackingParticle is collison like, get the momentum and vertex at production state
	if(parametersDefiner=="LhcParametersDefinerForTP")
	  {
	    if(! tpSelector(*tp)) continue;
	    momentumTP = tp->momentum();
	    vertexTP = tp->vertex();
	    //Calcualte the impact parameters w.r.t. PCA
	    ParticleBase::Vector momentum = parametersDefinerTP->momentum(event,setup,*tp);
	    ParticleBase::Point vertex = parametersDefinerTP->vertex(event,setup,*tp);
	    dxySim = (-vertex.x()*sin(momentum.phi())+vertex.y()*cos(momentum.phi()));
	    dzSim = vertex.z() - (vertex.x()*momentum.x()+vertex.y()*momentum.y())/sqrt(momentum.perp2()) 
	      * momentum.z()/sqrt(momentum.perp2());
	  }
	//If the TrackingParticle is comics, get the momentum and vertex at PCA
	if(parametersDefiner=="CosmicParametersDefinerForTP")
	  {
	    if(! cosmictpSelector(*tp,&bs,event,setup)) continue;	
	    momentumTP = parametersDefinerTP->momentum(event,setup,*tp);
	    vertexTP = parametersDefinerTP->vertex(event,setup,*tp);
	    dxySim = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
	    dzSim = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y())/sqrt(momentumTP.perp2()) 
	      * momentumTP.z()/sqrt(momentumTP.perp2());
	  }
	//---------- THE PART ABOVE HAS TO BE CLEANED UP. THE PARAMETER DEFINER WAS NOT MEANT TO BE USED IN THIS WAY ----------

	st++;   //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )

	// in the coming lines, histos are filled using as input
	// - momentumTP 
	// - vertexTP 
	// - dxySim
	// - dzSim

	histoProducerAlgo_->fill_generic_simTrack_histos(w,momentumTP,vertexTP);


	// ##############################################
	// fill RecoAssociated SimTracks' histograms
	// ##############################################
	bool isRecoMatched(false);
	const reco::Track* matchedTrackPointer=0;
	std::vector<std::pair<RefToBase<Track>, double> > rt;
	int nSharedHits(0),nRecoTrackHits(0);
	if(simRecColl.find(tpr) != simRecColl.end()){
	  rt = (std::vector<std::pair<RefToBase<Track>, double> >) simRecColl[tpr];
	  if (rt.size()!=0) {
	    ats++; //This counter counts the number of simTracks that have a recoTrack associated
	    isRecoMatched = true; matchedTrackPointer = rt.begin()->first.get();
	    nSharedHits = (int) rt.begin()->second;
	    nRecoTrackHits = matchedTrackPointer->numberOfValidHits();
	    edm::LogVerbatim("TrackValidator") << "TrackingParticle #" << st 
					       << " with pt=" << sqrt(momentumTP.perp2()) 
					       << " associated with quality:" << rt.begin()->second <<"\n";
	  }
	}else{
	  edm::LogVerbatim("TrackValidator") 
	    << "TrackingParticle #" << st
	    << " with pt,eta,phi: " 
	    << sqrt(momentumTP.perp2()) << " , "
	    << momentumTP.eta() << " , "
	    << momentumTP.phi() << " , "
	    << " NOT associated to any reco::Track" << "\n";
	}
	

	

	std::vector<PSimHit> simhits=tp->trackPSimHit(DetId::Tracker);	
	//for(std::vector<PSimHit>::const_iterator it=simhits.begin(); it!=simhits.end(); it++){
	//  std::cout << "simHit p: " << it->pabs() << std::endl;
	//}

        int nSimHits = simhits.size();
	
	//WARNING: it actually has to be called also for simTrack that are not matched to recoTrack
	histoProducerAlgo_->fill_recoAssociated_simTrack_histos(w,*tp,momentumTP,vertexTP,dxySim,dzSim,nSimHits,matchedTrackPointer);


	if ( isRecoMatched ) {
	  ParticleBase::Vector momTP = parametersDefinerTP->momentum(event,setup,*(tpr.get()));
	  ParticleBase::Point vtxTP = parametersDefinerTP->vertex(event,setup,*(tpr.get()));		 	 
	  int chargeTP = tpr->charge();

	  histoProducerAlgo_->fill_ResoAndPull_recoTrack_histos(w,*tp,momTP,vtxTP,chargeTP,
								*matchedTrackPointer,bs.position());
	  histoProducerAlgo_->fill_hitFinding_histos(w,*tp,momTP,vtxTP,chargeTP,
						     nSimHits,nSharedHits,nRecoTrackHits,
						     *matchedTrackPointer);
	}
	

      } // End  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){

      //if (st!=0) h_tracksSIM[w]->Fill(st);  // TO BE FIXED
      
      
      // ##############################################
      // fill recoTracks histograms (LOOP OVER TRACKS)
      // ##############################################
      edm::LogVerbatim("TrackValidator") << "\n# of reco::Tracks with "
					 << label[www].process()<<":"
					 << label[www].label()<<":"
					 << label[www].instance()
					 << ": " << trackCollection->size() << "\n";

      int at(0); //This counter counts the number of recoTracks that are associated to SimTracks
      int rT(0); //This counter counts the number of recoTracks in general


      // dE/dx
      // at some point this could be generalized, with a vector of tags and a corresponding vector of Handles
      // I'm writing the interface such to take vectors of ValueMaps
      edm::Handle<edm::ValueMap<reco::DeDxData> > dEdx1Handle;
      edm::Handle<edm::ValueMap<reco::DeDxData> > dEdx2Handle;
      std::vector<edm::ValueMap<reco::DeDxData> > v_dEdx;
      v_dEdx.clear();
      //std::cout << "PIPPO: label is " << label[www] << std::endl;
      if (label[www].label()=="generalTracks") {
	try {
	  event.getByLabel(m_dEdx1Tag, dEdx1Handle);
	  const edm::ValueMap<reco::DeDxData> dEdx1 = *dEdx1Handle.product();
	  event.getByLabel(m_dEdx2Tag, dEdx2Handle);
	  const edm::ValueMap<reco::DeDxData> dEdx2 = *dEdx2Handle.product();
	  v_dEdx.push_back(dEdx1);
	  v_dEdx.push_back(dEdx2);
	} catch (cms::Exception e){
	  LogTrace("TrackValidator") << "exception found: " << e.what() << "\n";
	}
      }
      //end dE/dx

      for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
	RefToBase<Track> track(trackCollection, i);
	rT++;
	
	bool isSimMatched(false);
	std::vector<std::pair<TrackingParticleRef, double> > tp;
	if(recSimColl.find(track) != recSimColl.end()){
	  tp = recSimColl[track];
	  if (tp.size()!=0) {
	    isSimMatched = true;
	    at++;
	    edm::LogVerbatim("TrackValidator") << "reco::Track #" << rT << " with pt=" << track->pt() 
					       << " associated with quality:" << tp.begin()->second <<"\n";
	  }
	} else {
	  edm::LogVerbatim("TrackValidator") << "reco::Track #" << rT << " with pt=" << track->pt()
					     << " NOT associated to any TrackingParticle" << "\n";		  
	}
	

	histoProducerAlgo_->fill_generic_recoTrack_histos(w,*track,bs.position(),isSimMatched);

	// dE/dx
	//	reco::TrackRef track2  = reco::TrackRef( trackCollection, i );
	if (v_dEdx.size() > 0) histoProducerAlgo_->fill_dedx_recoTrack_histos(w,track, v_dEdx);
	//if (v_dEdx.size() > 0) histoProducerAlgo_->fill_dedx_recoTrack_histos(track2, v_dEdx);


	//Fill other histos
 	//try{ //Is this really necessary ????

	if (tp.size()==0) continue;	

	histoProducerAlgo_->fill_simAssociated_recoTrack_histos(w,*track);

	/*
	  } // End of try{
	  catch (cms::Exception e){
	  LogTrace("TrackValidator") << "exception found: " << e.what() << "\n";
	  }
	*/
	
      } // End of for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){

      //TO BE FIXED
      //if (at!=0) h_tracks[w]->Fill(at);
      //h_fakes[w]->Fill(rT-at);
      //nrec_vs_nsim[w]->Fill(rT,st);

      edm::LogVerbatim("TrackValidator") << "Total Simulated: " << st << "\n"
					 << "Total Associated (simToReco): " << ats << "\n"
					 << "Total Reconstructed: " << rT << "\n"
					 << "Total Associated (recoToSim): " << at << "\n"
					 << "Total Fakes: " << rT-at << "\n";

      w++;
    } // End of  for (unsigned int www=0;www<label.size();www++){
  } //END of for (unsigned int ww=0;ww<associators.size();ww++){
}

void MultiTrackValidator::endRun(Run const&, EventSetup const&) {
  int w=0;
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){
      if(!skipHistoFit && runStandalone)	histoProducerAlgo_->finalHistoFits(w);
      if (runStandalone) histoProducerAlgo_->fillProfileHistosFromVectors(w);
      histoProducerAlgo_->fillHistosFromVectors(w);
      w++;
    }    
  }
  if ( out.size() != 0 && dbe_ ) dbe_->save(out);
}




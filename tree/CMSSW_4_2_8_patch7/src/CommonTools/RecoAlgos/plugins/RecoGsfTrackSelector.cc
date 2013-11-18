/** \class reco::RecoGsfTrackSelector
 *
 * Filter to select tracks according to pt, rapidity, tip, lip, number of hits, chi2
 *
 * \author Giuseppe Cerati, INFN
 *
 *  $Date: 2009/03/04 13:11:30 $
 *  $Revision: 1.1 $
 *
 */
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/RecoAlgos/interface/TrackSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoGsfTrackSelector.h"

namespace reco {
  typedef ObjectSelector<RecoGsfTrackSelector> RecoGsfTrackSelector;
  DEFINE_FWK_MODULE(RecoGsfTrackSelector);
}

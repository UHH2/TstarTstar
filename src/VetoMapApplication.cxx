#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/TstarTstar/include/PuppiCHSMatching.h"


#include "UHH2/TstarTstar/include/VetoMapApplication.h"

using namespace std;
using namespace uhh2;

VetoMapApplicator::VetoMapApplicator(Context& ctx) {
    
    AK4cleaner.reset( new JetCleaner(ctx, HotZoneVetoId() ) );
    HOTVRcleaner.reset( new TopJetCleaner(ctx, HotZoneVetoId() ) );

    PuppiCHSMatcher.reset(new PuppiCHSMatching(ctx, "CHS_matched"));

    h_is_btagevent = ctx.get_handle<bool>("is_btagevent");
    h_ST_AK4 = ctx.get_handle<double>("ST_AK4");
    h_ST_HOTVR = ctx.get_handle<double>("ST_HOTVR");
    h_CHS_matched = ctx.get_handle<vector<Jet>>("CHS_matched");

}

bool VetoMapApplicator::process(Event& event){

    // clean jet collections based on veto ID
    if(!(AK4cleaner->process(event))) return false;
    if(!(HOTVRcleaner->process(event))) return false;

    // need to repeat AK4 CHS matching
    PuppiCHSMatcher->process(event);

    // repeat N(AK4)
    bool pass_njet = (event.jets->size()>3);
    if(!pass_njet) return false;
    
    // repeat N(HOTVR)
    bool pass_fat_njet = (event.topjets->size()>0);
    if(!pass_fat_njet) return false;

    // repeat the b-tagging cut
    BTag bJetID = BTag(BTag::algo::DEEPJET, BTag::wp::WP_MEDIUM);
    bool pass_btagcut = false;
    for (const auto & jet : event.get(h_CHS_matched)) {if(bJetID(jet, event)) pass_btagcut = true;}
    event.set(h_is_btagevent, pass_btagcut);

    // recalculate ST
    // st calculation, for usage in plotting later on
    double st = 0.;
    for(const auto & lepton : *event.electrons) st += lepton.pt();
    for(const auto & lepton : *event.muons) st += lepton.pt();
    st += event.met->pt();
    double stHOTVR = st;

    for(const auto & jet : *event.jets) st += jet.pt();
    for(const auto & jet : *event.topjets) stHOTVR += jet.pt();

    event.set(h_ST_AK4, st);
    event.set(h_ST_HOTVR, stHOTVR);

    // repeat the ST > 500 cut
    if(stHOTVR < 500) return false;

    return true;
}

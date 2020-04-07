// -*- C++ -*-
//
// Package:    RPC-DTTrigger/DTRPCTimingUpdate
// Class:      DTRPCTimingUpdate
// 
/**\class DTRPCTimingUpdate DTRPCTimingUpdate.cc RPC-DTTrigger/DTRPCTimingUpdate/plugins/DTRPCTimingUpdate.cc
*/
// Original Author:  Jiwon Park
//         Created:  Thu, 09 Jan 2020 04:46:39 GMT

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

using namespace edm;

class DTRPCTimingUpdate : public edm::stream::EDProducer<> {
  public:
    explicit DTRPCTimingUpdate(const edm::ParameterSet&);
    ~DTRPCTimingUpdate();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    edm::ESHandle<DTGeometry> dtGeo;
    edm::ESHandle<RPCGeometry> rpcGeo;

    edm::EDGetTokenT<RPCRecHitCollection> src_;
    edm::EDGetTokenT<DTRecSegment4DCollection> dt4DSegments_;

    //GEANT4 simhits
    edm::EDGetTokenT<edm::PSimHitContainer> DTsimHitToken_;
    edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> RPCdigisimlinkToken_;

    int label_;
    int numDigi_switch;
};

/* Member Functions */

DTRPCTimingUpdate::DTRPCTimingUpdate(const edm::ParameterSet& iConfig):
  src_(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("src") )),
  dt4DSegments_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dt4DSegments") )),
  DTsimHitToken_(consumes<PSimHitContainer>(iConfig.getUntrackedParameter<edm::InputTag>("DTsimHitLabel", edm::InputTag("g4SimHits:MuonDTHits")) )),
  RPCdigisimlinkToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(iConfig.getParameter<edm::InputTag>("rpcSimLinkLabel") )),
  label_(iConfig.getUntrackedParameter<int>("label"))
{

  produces<RPCRecHitCollection>().setBranchAlias("out");
  produces<RPCRecHitCollection>("nocorr").setBranchAlias("out_nocorr");
  produces<RPCRecHitCollection>("all").setBranchAlias("out_all");

}


DTRPCTimingUpdate::~DTRPCTimingUpdate()
{
}

void
DTRPCTimingUpdate::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  iSetup.get<MuonGeometryRecord>().get( dtGeo );    
  iSetup.get<MuonGeometryRecord>().get( rpcGeo );

  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(src_, rpcRecHits);
  edm::Handle<DTRecSegment4DCollection> all4DSegments;
  iEvent.getByToken(dt4DSegments_, all4DSegments);
  edm::Handle<PSimHitContainer> DTsimHit;
  iEvent.getByToken(DTsimHitToken_, DTsimHit);
  edm::Handle<edm::DetSetVector<RPCDigiSimLink>> thelinkDigis;
  iEvent.getByToken(RPCdigisimlinkToken_, thelinkDigis);

  if (!all4DSegments.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find DTRecSegment4D with label "<< all4DSegments << std::endl;
    return;
  }
  if (!rpcRecHits.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find RPCRecHitDigiCollection with label "<< rpcRecHits << std::endl;
    return;
  }

  std::unique_ptr<RPCRecHitCollection> out(new RPCRecHitCollection());
  std::vector<RPCRecHit> vec_hits;
  std::unique_ptr<RPCRecHitCollection> out_nocorr(new RPCRecHitCollection());
  std::vector<RPCRecHit> vec_hits_nocorr;
  std::unique_ptr<RPCRecHitCollection> out_all(new RPCRecHitCollection());
  std::vector<RPCRecHit> vec_hits_all;


  // RPC RecHit Loop
  RPCDetId tmp_id;
  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

    RPCRecHit aRecHit = *rpcIt;
    LocalError tmp_err1(-99.,-99.,-99.);//Set dummy number for later
    aRecHit.setError(tmp_err1);

    RPCRecHit aNocorrRecHit = *rpcIt;
    LocalError tmp_err2(-99.,-99.,-99.);//Set dummy number for later
    aNocorrRecHit.setError(tmp_err2);

    RPCRecHit allRecHit = *rpcIt;
    LocalError tmp_err3(-99.,-99.,-99.);//Set dummy number for later
    aRecHit.setError(tmp_err3);

    RPCDetId rpc_id = (RPCDetId)(*rpcIt).rpcId();
    if (rpc_id.region() != 0) continue; //skip the barrels


    // DT Loop
    //Ref of loop: https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/RecoLocalMuon/DTSegment/src/DTSegment4DT0Corrector.cc
    DTRecSegment4DCollection::id_iterator dtChamberId;
    for (dtChamberId = all4DSegments->id_begin(); dtChamberId != all4DSegments->id_end(); ++dtChamberId) {
  
      numDigi_switch = label_;
      DTRecSegment4DCollection::range dtRange = all4DSegments->get(*dtChamberId);
      if (numDigi_switch == 1){
        dtRange.first++;
        if (dtRange.first != dtRange.second) continue; //Check only 1 digi(segment) or not
        dtRange.first--;
      }
  
      for (DTRecSegment4DCollection::const_iterator dtSegment = dtRange.first; dtSegment != dtRange.second; ++dtSegment) {
  
        DTRecSegment4D tmpseg = *dtSegment;
        if (!tmpseg.hasZed()) continue; //Looking at the segment that has Z projection FIXME Check needed
        if (dtSegment->dimension() != 4) continue; //Check needed //MB4 doesn't have theta superlayer

        const DTChamberId dt_id = dtSegment->chamberId();
        if (dt_id.station() != rpc_id.station()) continue; //Check only nearest station
 
        const LocalPoint& dt_lp = tmpseg.localPosition();
        const LocalVector& dt_dir = tmpseg.localDirection();
        const LocalError& dt_posErr = tmpseg.localPositionError();
  
        const GeomDet* gdet=dtGeo->idToDet(tmpseg.geographicalId());
        const BoundPlane &DTSurface = gdet->surface();
  
        float Pos_xx = dt_posErr.xx();
        float Pos_xy = dt_posErr.xy();
        float Pos_yy = dt_posErr.yy();
  
        float Xo = dt_lp.x();
        float Yo = dt_lp.y();
        float Zo = dt_lp.z();
        float dx = dt_dir.x();
        float dy = dt_dir.y();
        float dz = dt_dir.z();
  
        int cptype = 0;
        bool dt_simmatched = false;
  
        for (PSimHitContainer::const_iterator DTsimIt = DTsimHit->begin(); DTsimIt != DTsimHit->end(); DTsimIt++) {
          cptype = DTsimIt->particleType();
          DTChamberId dtsim_id = DTsimIt->detUnitId();
  
          if (abs(cptype) != 13) continue;
          if (dt_id.wheel() == dtsim_id.wheel() && dt_id.station() == dtsim_id.station() && dt_id.sector() == dtsim_id.sector()){
            DTSuperLayerId dt_slID(DTsimIt->detUnitId());
            int superlayer = dt_slID.superlayer();
  
            dt_simmatched = true; 
          }
        }
  
        if (!dt_simmatched) continue;

        if (rpcIt->BunchX() == 0){
          LocalError error_all(allRecHit.localPositionError().xx(), 1.0, aRecHit.localPositionError().yy());
          allRecHit.setError(error_all);
        }

        const BoundPlane& RPCSurface = rpcGeo->idToDet(rpc_id)->surface();
        GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));
        LocalPoint CenterRollinDTFrame = DTSurface.toLocal(CenterPointRollGlobal);
        float D = CenterRollinDTFrame.z();
        if (abs(D) > 200) continue;

        LocalPoint lp_dtframe(Xo + dx*D/dz, Yo + dy*D/dz, D);
        GlobalPoint dt_gp_ext = dtGeo->idToDet(dt_id)->surface().toGlobal(lp_dtframe);
        LocalPoint lp_extrapol = rpcGeo->idToDet(rpc_id)->surface().toLocal(dt_gp_ext);

        //local distance
        LocalPoint lp_rpc(0.0,0.0,0.0);
        lp_rpc = aRecHit.localPosition();
        const auto& rpc_bound = RPCSurface.bounds(); //Check if extrapolated segment is within RPC chamber
        //std::cout << rpc_bound.length() << " " << rpc_bound.width() << std::endl;
        //https://github.com/cms-sw/cmssw/blob/04a85943355965447837cff7da77bde47a88e176/RecoLocalMuon/RPCRecHit/src/TracktoRPC.cc#L33
        if (not( fabs(lp_extrapol.x()) < rpc_bound.width()/2 && fabs(lp_extrapol.y()) < rpc_bound.length()/2 )) continue;

        float Dx = abs(lp_rpc.x() - lp_extrapol.x());
        float Dy = abs(lp_rpc.y() - lp_extrapol.y());

        if (Dx < 2){//arbitrary cutoff here FIXME
          //std::cout << lp_extrapol << " / " << lp_rpc << std::endl;
          //https://github.com/cms-sw/cmssw/blob/9a33fb13bee1a546877a4b581fa63876043f38f0/SimMuon/MCTruth/src/RPCHitAssociator.cc#L74-L92
          int fstrip = rpcIt->firstClusterStrip();
          int cls = rpcIt->clusterSize();
          int bx = rpcIt->BunchX();

          std::vector<SimHitIdpr> matched;
          std::set<RPCDigiSimLink> links;
          for (int i = fstrip; i < fstrip + cls; ++i) {

            for (edm::DetSetVector<RPCDigiSimLink>::const_iterator itlink = thelinkDigis->begin(); itlink != thelinkDigis->end(); itlink++) {
              for (edm::DetSet<RPCDigiSimLink>::const_iterator digi_iter = itlink->data.begin(); digi_iter != itlink->data.end(); ++digi_iter) {
                uint32_t detid = digi_iter->getDetUnitId();
                int str = digi_iter->getStrip();
                int bunchx = digi_iter->getBx();

                if (detid == rpc_id && str == i && bunchx == bx) {
                  links.insert(*digi_iter);
                }
              }
            }

            if (links.empty()) std::cout << "Unmatched simHit!" << std::endl;

            for (std::set<RPCDigiSimLink>::iterator itlink = links.begin(); itlink != links.end(); ++itlink) {
              SimHitIdpr currentId(itlink->getTrackId(), itlink->getEventId());
              if (find(matched.begin(), matched.end(), currentId) == matched.end()) matched.push_back(currentId);
            }
            //std::cout << matched.size() << " ";
          }
          if (!links.empty()) {
            //https://github.com/cms-sw/cmssw/blob/master/SimMuon/RPCDigitizer/src/RPCSynchronizer.cc#L107
            //https://github.com/cms-sw/cmssw/blob/master/SimMuon/RPCDigitizer/python/muonRPCDigis_cfi.py#L39
            double c = 299792458;  // [m/s]
            //light speed in [cm/ns]
            double cspeed = c * 1e+2 * 1e-9;
            //signal propagation speed [cm/ns]
            double sspeed = 0.66 * cspeed;

            //https://github.com/cms-sw/cmssw/blob/04a85943355965447837cff7da77bde47a88e176/RecoLocalMuon/RPCRecHit/src/TracktoRPC.cc
            const GeomDet *geomDet2 = rpcGeo->idToDet(rpc_id);
            const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(geomDet2);
            const RectangularStripTopology *top_ = dynamic_cast<const RectangularStripTopology *>(&(aroll->topology()));
            //LocalPoint xmin = top_->localPosition(0.);
            //LocalPoint xmax = top_->localPosition((float)aroll->nstrips());
            //float rsize = fabs(xmax.x() - xmin.x());
            float stripLenHalf = top_->stripLength() / 2;
            //region == 0: distanceFromEdge = half_stripL + simHitPos.y() (RPCSynchronizer.cc)
            //This suggests readout at -y end!
            float prop_length = stripLenHalf - lp_extrapol.y(); //This is ref
            //float prop_length = stripLenHalf + lp_extrapol.y(); //sgnPlus


            //std::cout << links.begin()->getTimeOfFlight() << " ";
            if (abs(links.begin()->getParticleType()) == 13 and links.begin()->getBx() == 0 and links.begin()->getMomentumAtEntry().perp() > 15) {
              if (aRecHit.corrTime() != 0.) continue;
              //aRecHit.setCorrTime(links.begin()->getTimeOfFlight() - (rpcIt->time() - stripLenHalf/sspeed + prop_length/sspeed));
              //aRecHit.setCorrTime(rpcIt->time() - stripLenHalf/sspeed + prop_length/sspeed);
              aRecHit.setTimeAndError(rpcIt->time() - stripLenHalf/sspeed + prop_length/sspeed, rpcIt->timeError());
              aRecHit.setCorrTime(rpcIt->time());//Trick! swap time and corrTime
              LocalError error(lp_extrapol.y(), aRecHit.localPositionError().xy(), aRecHit.localPositionError().yy());
              aRecHit.setError(error);

              aNocorrRecHit.setTimeAndError(rpcIt->time(), rpcIt->timeError());
              aNocorrRecHit.setCorrTime(rpcIt->time());
              LocalError error2(lp_extrapol.y(), aNocorrRecHit.localPositionError().xy(), aNocorrRecHit.localPositionError().yy());
              aNocorrRecHit.setError(error2);

              allRecHit.setTimeAndError(rpcIt->time() - stripLenHalf/sspeed + prop_length/sspeed, rpcIt->timeError());
              allRecHit.setCorrTime(rpcIt->time());
              LocalError error3(lp_extrapol.y(), allRecHit.localPositionError().xy(), allRecHit.localPositionError().yy());
              allRecHit.setError(error3);
            }
            //else aRecHit.setCorrTime(-99);
            else continue;

            if(tmp_id != rpc_id){
              if(!vec_hits.empty()) out->put(rpc_id, vec_hits.begin(), vec_hits.end());
              tmp_id = rpc_id;
              vec_hits.clear();
              vec_hits.push_back(aRecHit);

              if(!vec_hits_nocorr.empty()) out_nocorr->put(rpc_id, vec_hits_nocorr.begin(), vec_hits_nocorr.end());
              tmp_id = rpc_id;
              vec_hits_nocorr.clear();
              vec_hits_nocorr.push_back(aNocorrRecHit);
            }
            else{
              vec_hits.push_back(aRecHit);
              vec_hits_nocorr.push_back(aNocorrRecHit);
            }

//            if(tmp_id != rpc_id){
//              if(!vec_hits_nocorr.empty()) out_nocorr->put(rpc_id, vec_hits_nocorr.begin(), vec_hits_nocorr.end());
//              tmp_id = rpc_id;
//              vec_hits_nocorr.clear();
//              vec_hits_nocorr.push_back(aNocorrRecHit);
//            }
//            else vec_hits_nocorr.push_back(aNocorrRecHit);
          }

        }//Window cut
      }
    }

    //test if recHits are correctly filled in  loop
    if(tmp_id != rpc_id){
      if(!vec_hits_all.empty()) out_all->put(rpc_id, vec_hits_all.begin(), vec_hits_all.end());
      tmp_id = rpc_id;
      vec_hits_all.clear();
      if(allRecHit.localPositionError().xy() > 0){ //positive means passed pre-selection upto bx=0)
        vec_hits_all.push_back(allRecHit);
      }
    }
    else vec_hits_all.push_back(allRecHit);

  }//Rechit loop
  iEvent.put(std::move(out));
  iEvent.put(std::move(out_nocorr), "nocorr");
  iEvent.put(std::move(out_all), "all");

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DTRPCTimingUpdate::beginStream(edm::StreamID)
{
}
// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DTRPCTimingUpdate::endStream() {
}
// ------------ method called when starting to processes a run  ------------
/*
void
DTRPCTimingUpdate::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DTRPCTimingUpdate::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DTRPCTimingUpdate::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DTRPCTimingUpdate::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DTRPCTimingUpdate::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DTRPCTimingUpdate);

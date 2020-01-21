// -*- C++ -*-
//
// Package:    RPC-DTTrigger/DTRPCTimingUpdate
// Class:      DTRPCTimingUpdate
// 
/**\class DTRPCTimingUpdate DTRPCTimingUpdate.cc RPC-DTTrigger/DTRPCTimingUpdate/plugins/DTRPCTimingUpdate.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jiwon Park
//         Created:  Thu, 09 Jan 2020 04:46:39 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

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

    //edm::EDGetTokenT<MuonDigiCollection<RPCDetId,RPCDigi>> src_;
    edm::EDGetTokenT<RPCRecHitCollection> src_;

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
DTRPCTimingUpdate::DTRPCTimingUpdate(const edm::ParameterSet& iConfig):
  src_(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("src") ))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
  produces<RPCRecHitCollection>(); 
}


DTRPCTimingUpdate::~DTRPCTimingUpdate()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DTRPCTimingUpdate::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   iEvent.put(std::make_unique<ExampleData2>(*pIn));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(src_, rpcRecHits);
  std::unique_ptr<RPCRecHitCollection> out(new RPCRecHitCollection());
 
  RPCDetId tmp_id;
  std::vector<RPCRecHit> vec_hits;
  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {
    RPCRecHit aRecHit = *rpcIt;
    RPCDetId rpc_id = (RPCDetId)(*rpcIt).rpcId();
    aRecHit.setCorrTime(999);

    if(tmp_id != rpc_id){
      if(!vec_hits.empty()) out->put(rpc_id, vec_hits.begin(), vec_hits.end());
      tmp_id = rpc_id;
      vec_hits.clear();
      vec_hits.push_back(aRecHit);
    }
    else vec_hits.push_back(aRecHit);
  }
  iEvent.put(std::move(out));

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

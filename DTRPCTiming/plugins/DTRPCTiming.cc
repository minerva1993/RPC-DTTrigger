// -*- C++ -*-
//
// Package:    RPC+CSCTrigger/DTRPCTiming
// Class:      DTRPCTiming

//
// Original Author:  Ms. Choi Jieun
// Modified to DT study :  Mr. Jiwon Park

// system include files
#include <memory>
#include <iostream>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

using namespace edm;
using namespace std;

class DTRPCTiming : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit DTRPCTiming(const edm::ParameterSet&);
    ~DTRPCTiming();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    //https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1TMuon/interface/GeometryTranslator.h
    const DTGeometry& getDTGeometry() const { return *dtGeo; }

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    TTree *tree;
    TH1D *EventInfo;
/*
    TH1D *h_ME31NDigis;
    TH1D *h_ME41NDigis;
    TH1D *h_ME31NDigis0;
    TH1D *h_ME41NDigis0;

    TH1D *NRecHits;
    TH1D *h_RE31NRecHits;
    TH1D *h_RE41NRecHits;
    TH1D *h_RE31NRecHits0;
    TH1D *h_RE41NRecHits0;

    TH1D *h_S3NDigis;
    TH1D *h_S4NDigis;
    TH1D *h_S3NRecHits;
    TH1D *h_S4NRecHits;

    TH1D *h_cscME31NSimHits;
    TH1D *h_cscME41NSimHits;
    TH1D *h_rpcME31NSimHits;
    TH1D *h_rpcME41NSimHits;
*/
    unsigned int b_EVENT, b_RUN, b_LUMI;
/*
    double b_ME31NDigis;
    double b_ME41NDigis;

    unsigned int bx_ME31NDigis;
    unsigned int bx_ME41NDigis;
    double b_ME31NDigis_Total;
    double b_ME41NDigis_Total;

    double pure_ME31NDigis_Total;
    double pure_ME41NDigis_Total;

    unsigned int b_S3NRecHits;
    unsigned int b_S4NRecHits;
    unsigned int b_RE31NRecHits;
    unsigned int b_RE41NRecHits;
    unsigned int bx_RE31NRecHits;
    unsigned int bx_RE41NRecHits;

    unsigned int b_S3NDigis;
    unsigned int b_S4NDigis;

    unsigned int b_ME31NSimHits;
    unsigned int b_ME41NSimHits;
    unsigned int b_RE31NSimHits;
    unsigned int b_RE41NSimHits;

    int b_ptype;
    TH1D *h_ptype;

    int nRPC;
    int nCSC;
    int b_rpcBX;
    int b_cscBX;

    TH1D *h_xNMatchedME31;
    TH1D *h_xNMatchedME41;

    TH1D *h_yNMatchedME31;
    TH1D *h_yNMatchedME41;

    TH1D *h_simValidME31x;
    TH1D *h_simValidME31y;
    TH1D *h_simValidME41x;
    TH1D *h_simValidME41y;
  
    TH2D *h_MatchedME31;
    TH2D *h_MatchedME41;

    TH2D *h_RatioME31;
    TH2D *h_RatioME41;

    double sME31x[100];
    double sME31y[100];
    double sME41x[100];
    double sME41y[100];

    bool isValidME31x[100];
    bool isValidME31y[100];
    bool isValidME41x[100];
    bool isValidME41y[100];

    double ME31[25][25];
    double ME41[25][25];

    bool isMatchME31[25][25];
    bool isMatchME41[25][25];
*/
    int EventNum;
//    int label_;
//    int numDigi_switch;

    edm::ESHandle<DTGeometry> dtGeo;
    edm::ESHandle<RPCGeometry> rpcGeo;

    edm::Handle<DTRecSegment4DCollection> all4DSegments;
    edm::EDGetTokenT<DTRecSegment4DCollection> dt4DSegments;
    //std::vector<DTRecSegment4D> dtseg;

    edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
    edm::Handle<RPCRecHitCollection> rpcRecHits;

    //GEANT4 simhits
    edm::EDGetTokenT<edm::PSimHitContainer> DTsimHitToken;

    //GlobalPoint getDTGlobalPosition(DTChamberId rawId, const DTRecSegment4DCollection& dt4Dseg) const;
    GlobalPoint getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const;
    //std::pair<RPCRecHitCollection::const_iterator, float*> matchingRPC(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct, float dx_cutoff, float dy_cutoff) const;
};

/*
GlobalPoint
DTRPCTiming::getDTGlobalPosition(DTChamberId rawId, const DTRecSegment4DCollection& dt4Dseg) const{

  DTChamberId cscId = CSCDetId(rawId);
  DTChamberId key_id(cscId.endcap(), cscId.station(), cscId.ring(), cscId.chamber(), CSCConstants::KEY_CLCT_LAYER);

  const auto& dtChamber = getDTGeometry().chamber(cscId);
  float fractional_strip = lct.getFractionalStrip();
  //CSCs have 6 layers. The key (refernce) layer is the third layer (CSCConstant)
  const auto& layer_geo = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();

  // LCT::getKeyWG() also starts from 0
  float wire = layer_geo->middleWireOfGroup(lct.getKeyWG() + 1);
  const LocalPoint& csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  const GlobalPoint& csc_gp = cscGeo->idToDet(key_id)->surface().toGlobal(csc_intersect);

  return csc_gp;

}
*/
GlobalPoint
DTRPCTiming::getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const{

  RPCDetId rpcid = RPCDetId(rpcId);
  const LocalPoint& rpc_lp = rpcIt.localPosition();
  const GlobalPoint& rpc_gp = rpcGeo->idToDet(rpcid)->surface().toGlobal(rpc_lp);
 
  return rpc_gp;

}
/*
//RPCRecHitCollection::const_iterator
std::pair<RPCRecHitCollection::const_iterator, float*>
DTRPCTiming::matchingRPC(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct, float dx_cutoff, float dy_cutoff) const{
//const float & makes this faster?

  GlobalPoint gp_cscint(0.0,0.0,0.0);
  gp_cscint = getCSCGlobalPosition(rawId, lct);

  float min_distance = std::numeric_limits<float>::max();
  RPCRecHitCollection::const_iterator rpc_hit_matched;
  float min_DxDy[2] = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  
  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {
    RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();

    GlobalPoint gp_rpc(0.0,0.0,0.0);
    gp_rpc = getRPCGlobalPosition(rpcid, *rpcIt);

    if (rpcid.region() == 0) continue; //skip the barrels

    float Dx = abs(gp_rpc.x()-gp_cscint.x());
    float Dy = abs(gp_rpc.y()-gp_cscint.y());
    float distance = sqrt(Dx*Dx + Dy*Dy);

    if (rawId.station() == 3 && rawId.ring() == 1 && rpcid.station() == 3 && rpcid.ring() == 1){
      if (Dx < dx_cutoff && Dy < dy_cutoff && distance < min_distance){
        min_DxDy[0]=Dx;
        min_DxDy[1]=Dy;
        min_distance = distance;
        rpc_hit_matched = rpcIt;
      }
    }
  }
//  return rpc_hit_matched;
  return std::make_pair(rpc_hit_matched,min_DxDy);
}
*/
DTRPCTiming::DTRPCTiming(const edm::ParameterSet& iConfig)
{
  dt4DSegments = consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dt4DSegments"));
  auto RPCDigiLabel = iConfig.getParameter<edm::InputTag>("simMuonRPCDigis");
  rpcRecHitsToken_ = consumes<RPCRecHitCollection>(edm::InputTag(RPCDigiLabel.label(), "" ));

  //GEANT4 simhit
  DTsimHitToken = consumes<PSimHitContainer>(iConfig.getUntrackedParameter<edm::InputTag>("DTsimHitLabel", edm::InputTag("g4SimHits:MuonDTHits")));

  //numDigi label
  //label_ = iConfig.getUntrackedParameter<int>("label");

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "Tree for RPC+CSCTrigger");

  EventInfo = fs->make<TH1D>("EventInfo","Event Information",2,0,2);
  EventInfo->GetXaxis()->SetBinLabel(1,"Total Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Selected Number of Events");
/*
  h_ME31NDigis = fs->make<TH1D>("h_ME31NDigis", "number of digi per chamber (ME3/1)", 10, 0, 10);
  h_ME31NDigis->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME31NDigis->GetYaxis()->SetTitle("Number of chamber");

  h_ME41NDigis = fs->make<TH1D>("h_ME41NDigis", "number of digi per chamber (ME4/1)", 10, 0, 10);
  h_ME41NDigis->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME41NDigis->GetYaxis()->SetTitle("Number of chamber");

  h_ME31NDigis0 = fs->make<TH1D>("h_ME31NDigis_withBX0", "number of digi per chamber (ME3/1)", 10, 0, 10);
  h_ME31NDigis0->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME31NDigis0->GetYaxis()->SetTitle("Number of chamber");

  h_ME41NDigis0 = fs->make<TH1D>("h_ME41NDigis_withBX0", "number of digi per chamber (ME4/1)", 10, 0, 10);
  h_ME41NDigis0->GetXaxis()->SetTitle("Number of digi per chamber");
  h_ME41NDigis0->GetYaxis()->SetTitle("Number of chamber");

  NRecHits = fs->make<TH1D>("NRecHits", "", 10, 0, 10);
  NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  NRecHits->GetYaxis()->SetTitle("Number of chamber");

  h_RE31NRecHits = fs->make<TH1D>("h_RE31NRecHits", "number of rechit per chamber (RE3/1)", 20, 0, 20);
  h_RE31NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE31NRecHits->GetYaxis()->SetTitle("Number of chamber");
   
  h_RE41NRecHits = fs->make<TH1D>("h_RE41NRecHits", "number of rechit per chamber (RE4/1)", 20, 0, 20);
  h_RE41NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE41NRecHits->GetYaxis()->SetTitle("Number of chamber");

  h_RE31NRecHits0 = fs->make<TH1D>("h_RE31NRecHits_withBX0", "number of rechit per chamber (RE3/1)", 20, 0, 20);
  h_RE31NRecHits0->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE31NRecHits0->GetYaxis()->SetTitle("Number of chamber");
   
  h_RE41NRecHits0 = fs->make<TH1D>("h_RE41NRecHits_withBX0", "number of rechit per chamber (RE4/1)", 20, 0, 20);
  h_RE41NRecHits0->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RE41NRecHits0->GetYaxis()->SetTitle("Number of chamber");

  h_S3NDigis = fs->make<TH1D>("h_S3NDigis", "number of digi in station 3", 10, 0, 10);
  h_S3NDigis->GetXaxis()->SetTitle("Number of digi");
  h_S3NDigis->GetYaxis()->SetTitle("Number of chamber");
  
  h_S4NDigis = fs->make<TH1D>("h_S4NDigis", "number of digi in station 4", 10, 0, 10);
  h_S4NDigis->GetXaxis()->SetTitle("Number of digi");
  h_S4NDigis->GetYaxis()->SetTitle("Number of chamber");

  h_S3NRecHits = fs->make<TH1D>("h_S3NRecHits", "number of rechit in station 3", 20, 0, 20);
  h_S3NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_S3NRecHits->GetYaxis()->SetTitle("Number of chamber");
  
  h_S4NRecHits = fs->make<TH1D>("h_S4NRecHits", "number of rechit in station 4", 20, 0, 20);
  h_S4NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_S4NRecHits->GetYaxis()->SetTitle("Number of chamber");

  h_cscME31NSimHits = fs->make<TH1D>("h_cscME31NSimHits", "number of cscsimhit per chamber (ME3/1)", 10, 0, 10);
  h_cscME31NSimHits->GetXaxis()->SetTitle("Number of simhit per chamber");
  h_cscME31NSimHits->GetYaxis()->SetTitle("Number of chamber");

  h_cscME41NSimHits = fs->make<TH1D>("h_cscME41NSimHits", "number of cscsimhit per chamber (ME4/1)", 10, 0, 10);
  h_cscME41NSimHits->GetXaxis()->SetTitle("Number of simhit per chamber");
  h_cscME41NSimHits->GetYaxis()->SetTitle("Number of chamber");

  h_rpcME31NSimHits = fs->make<TH1D>("h_rpcME31NSimHits", "number of rpcsimhit per chamber (ME3/1)", 10, 0, 10);
  h_rpcME31NSimHits->GetXaxis()->SetTitle("Number of simhit per chamber");
  h_rpcME31NSimHits->GetYaxis()->SetTitle("Number of chamber");

  h_rpcME41NSimHits = fs->make<TH1D>("h_rpcME41NSimHits", "number of rpcsimhit per chamber (ME4/1)", 10, 0, 10);
  h_rpcME41NSimHits->GetXaxis()->SetTitle("Number of simhit per chamber");
  h_rpcME41NSimHits->GetYaxis()->SetTitle("Number of chamber");

  h_xNMatchedME31 = fs->make<TH1D>("h_xNMatchedME31", "Matching Efficiency in ME3/1", 25, 0, 25);
  h_xNMatchedME31->GetXaxis()->SetTitle("X cutoff (cm)");
  h_xNMatchedME31->GetYaxis()->SetTitle("Matched (%)");

  h_xNMatchedME41 = fs->make<TH1D>("h_xNMatchedME41", "Matching Efficiency in ME4/1", 25, 0, 25);
  h_xNMatchedME41->GetXaxis()->SetTitle("X cutoff (cm)");
  h_xNMatchedME41->GetYaxis()->SetTitle("Matched (%)");

  h_yNMatchedME31 = fs->make<TH1D>("h_yNMatchedME31", "Matching Efficiency in ME3/1", 25, 0, 25);
  h_yNMatchedME31->GetXaxis()->SetTitle("Y cutoff (cm)");
  h_yNMatchedME31->GetYaxis()->SetTitle("Matched (%)");

  h_yNMatchedME41 = fs->make<TH1D>("h_yNMatchedME41", "Matching Efficiency in ME4/1", 25, 0, 25);
  h_yNMatchedME41->GetXaxis()->SetTitle("Y cutoff (cm)");
  h_yNMatchedME41->GetYaxis()->SetTitle("Matched (%)");

  h_MatchedME31 = fs->make<TH2D>("h_MatchedME31", "Matching efficiency in ME31", 25, 0, 25, 25, 0, 25);
  h_MatchedME31->GetXaxis()->SetTitle("X cutoff (cm)");
  h_MatchedME31->GetYaxis()->SetTitle("Y cutoff (cm)");

  h_MatchedME41 = fs->make<TH2D>("h_MatchedME41", "Matching efficiency in ME41", 25, 0, 25, 25, 0, 25);
  h_MatchedME41->GetXaxis()->SetTitle("X cutoff (cm)");
  h_MatchedME41->GetYaxis()->SetTitle("Y cutoff (cm)");

  h_RatioME31 = fs->make<TH2D>("h_RatioME31", "Matching efficiency in ME31", 25, 0, 25, 25, 0, 25);
  h_RatioME31->GetXaxis()->SetTitle("X cutoff (cm)");
  h_RatioME31->GetYaxis()->SetTitle("Y cutoff (cm)");

  h_RatioME41 = fs->make<TH2D>("h_RatioME41", "Matching efficiency in ME41", 25, 0, 25, 25, 0, 25);
  h_RatioME41->GetXaxis()->SetTitle("X cutoff (cm)");
  h_RatioME41->GetYaxis()->SetTitle("Y cutoff (cm)");

  h_ptype = fs->make<TH1D>("h_ptype", "", 4,0,4);
  h_ptype->GetXaxis()->SetBinLabel(1,"Electron");
  h_ptype->GetXaxis()->SetBinLabel(2,"Muon");
  h_ptype->GetXaxis()->SetBinLabel(3,"Pion");
  h_ptype->GetXaxis()->SetBinLabel(4,"Proton");
  h_ptype->GetXaxis()->SetTitle("Particle type");
  h_ptype->GetYaxis()->SetTitle("Entries");

  h_simValidME31x = fs->make<TH1D>("h_simValidME31x", "Validation Percentage in ME3/1", 100, 0, 100);
  h_simValidME31x->GetXaxis()->SetTitle("X cutoff (mm)");
  h_simValidME31x->GetYaxis()->SetTitle("Matched (%)");

  h_simValidME31y = fs->make<TH1D>("h_simValidME31y", "Validation Percentage in ME3/1", 100, 0, 100);
  h_simValidME31y->GetXaxis()->SetTitle("Y cutoff (mm)");
  h_simValidME31y->GetYaxis()->SetTitle("Matched (%)");

  h_simValidME41x = fs->make<TH1D>("h_simValidME41x", "Validation Percentage in ME4/1", 100, 0, 100);
  h_simValidME41x->GetXaxis()->SetTitle("X cutoff (mm)");
  h_simValidME41x->GetYaxis()->SetTitle("Matched (%)");

  h_simValidME41y = fs->make<TH1D>("h_simValidME41y", "Validation Percentage in ME4/1", 100, 0, 100);
  h_simValidME41y->GetXaxis()->SetTitle("Y cutoff (mm)");
  h_simValidME41y->GetYaxis()->SetTitle("Matched (%)");
*/
}

DTRPCTiming::~DTRPCTiming()
{

}

void
DTRPCTiming::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  iSetup.get<MuonGeometryRecord>().get( dtGeo );     
  iSetup.get<MuonGeometryRecord>().get( rpcGeo );     

  EventInfo->Fill(0.5);

  iEvent.getByToken(dt4DSegments, all4DSegments);
  iEvent.getByToken(rpcRecHitsToken_, rpcRecHits);
/*
  //GEANT4 simhits
  edm::Handle<PSimHitContainer> CSCsimHit;
  iEvent.getByToken(CSCsimHitToken, CSCsimHit);

  PSimHitContainer::const_iterator CSCsimIt;

  if (!corrlcts.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find CSCCorrleatedLCTDigiCollection with label "<< corrlcts << std::endl;
    return;
  }
  if (!rpcRecHits.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find RPCRecHitDigiCollection with label "<< rpcRecHits << std::endl;
    return;
  }
*/
  cout << "\nNew event" << endl;
  b_EVENT = b_RUN = b_LUMI = 0;

  b_EVENT  = iEvent.id().event();
  b_RUN    = iEvent.id().run();
  b_LUMI   = iEvent.id().luminosityBlock();

  EventNum++;
  cout << "Event " << EventNum << endl;
/*
  nRPC = nCSC = 0;
  b_S3NRecHits = b_S4NRecHits = 0;
  bx_RE31NRecHits = bx_RE41NRecHits =0;
  b_RE31NRecHits = b_RE41NRecHits = 0;
  b_ME31NSimHits = b_ME41NSimHits = b_RE31NSimHits = b_RE41NSimHits = 0;

  //to check rechit info
  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

    nRPC++;
    RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();

    if(rpcid.station() == 3) b_S3NRecHits++;
    if(rpcid.station() == 4) b_S4NRecHits++;

    if(rpcid.station() == 3 && rpcid.ring() == 1) b_RE31NRecHits++;
    if(rpcid.station() == 4 && rpcid.ring() == 1) b_RE41NRecHits++;
    if((*rpcIt).BunchX() == 0 && rpcid.station() == 3 && rpcid.ring() == 1) bx_RE31NRecHits++;
    if((*rpcIt).BunchX() == 0 && rpcid.station() == 4 && rpcid.ring() == 1) bx_RE41NRecHits++;

  }

  //to check cscsimhit Info
  for (CSCsimIt = CSCsimHit->begin(); CSCsimIt != CSCsimHit->end(); CSCsimIt++) {
    CSCDetId cscsimhit_id(CSCsimIt->detUnitId());
    b_ptype = CSCsimIt->particleType();
//    cout << "ptype :: " << b_ptype << endl;
    if ( abs(b_ptype) != 11 && abs(b_ptype) != 13 && abs(b_ptype) != 2212 && abs(b_ptype) != 211 ) cout << "NONEPARTICLE TYPE :: " << abs(b_ptype) << endl;
    
    if (abs(b_ptype) == 11) h_ptype->Fill(0.5);
    if (abs(b_ptype) == 13) h_ptype->Fill(1.5);
    if (abs(b_ptype) == 211) h_ptype->Fill(2.5);
    if (abs(b_ptype) == 2212) h_ptype->Fill(3.5);
    
    if (cscsimhit_id.station() == 3 && cscsimhit_id.ring() == 1) b_ME31NSimHits++;
    if (cscsimhit_id.station() == 4 && cscsimhit_id.ring() == 1) b_ME41NSimHits++;
  }

  //to check rpcsimhit Info

  for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++){  

    nCSC++;
    CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
    b_ME31NDigis = b_ME41NDigis = 0;

    bx_ME31NDigis = bx_ME41NDigis =0;
    b_cscBX = b_rpcBX = 0;

    b_S3NDigis = b_S4NDigis = 0;

    GlobalPoint gp_cscint;

    numDigi_switch = label_;
    if (numDigi_switch == 1){
      range1.first++;
      if (range1.first != range1.second) continue; // check that there are two digis in the chamber, there is probably a better way but it works...
      range1.first--;
    }
    else if (numDigi_switch == 4){
      range1.first++;
      range1.first++;
      range1.first++;
      range1.first++;
      if (range1.first != range1.second) continue; // check that there are two digis in the chamber, there is probably a better way but it works...
      range1.first--;
      range1.first--;
      range1.first--;
      range1.first--;
    }
//    else cout << "option is not decleared!" << endl;
   
    for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++){
      const CSCDetId csc_id((*csc).first.rawId());

      b_cscBX = lct->getBX();

      gp_cscint = GlobalPoint(0.0,0.0,0.0);
      gp_cscint = getCSCGlobalPosition(csc_id, *lct);
      if (abs(gp_cscint.eta()) < 1.9) continue;

      if (csc_id.station() == 3 && csc_id.ring() == 1) pure_ME31NDigis_Total++;
      if (csc_id.station() == 4 && csc_id.ring() == 1) pure_ME41NDigis_Total++;

      const auto cscCham = getCSCGeometry().chamber(csc_id);
      float fractional_strip = lct->getFractionalStrip();
      //CSCs have 6 layers. The key (refernce) layer is the third layer (CSCConstant)
      const auto layer_geo = cscCham->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();

      // LCT::getKeyWG() also starts from 0
      float wire = layer_geo->middleWireOfGroup(lct->getKeyWG() + 1);
      const LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);

      int cptype = 0;
      bool csc_simmatched = false;

      for (int i = 0; i < 100; i++){
        isValidME31x[i] = false;
        isValidME31y[i] = false;
        isValidME41x[i] = false;
        isValidME41y[i] = false;
      }

      for (CSCsimIt = CSCsimHit->begin(); CSCsimIt != CSCsimHit->end(); CSCsimIt++) {

        cptype = CSCsimIt->particleType();
        //const GlobalPoint sim_gp = cscGeo->idToDet(csc_id)->surface().toGlobal(CSCsimIt->localPosition());
        const LocalPoint lp_cscsim = CSCsimIt->localPosition();
        CSCDetId cscsim_id(CSCsimIt->detUnitId());

        //sim validation with window
        if (abs(cptype) != 13) continue;
        if (csc_id.endcap() == cscsim_id.endcap() && csc_id.station() == cscsim_id.station() && csc_id.ring() == cscsim_id.ring() && csc_id.chamber() == cscsim_id.chamber()){

          cout << "I'm the same" << endl;

          if (sqrt(csc_intersect.x()-lp_cscsim.x())*(csc_intersect.x()-lp_cscsim.x())+(csc_intersect.y()-lp_cscsim.y())*(csc_intersect.y()-lp_cscsim.y()) < 0.5) csc_simmatched = true;
  
          float sDx = abs(csc_intersect.x() - lp_cscsim.x());
          float sDy = abs(csc_intersect.y() - lp_cscsim.y());
  
          //ME31
          if (csc_id.station() == 3 && csc_id.ring() == 1) {
            for (int i = 0; i < 100; i++) {
              if (sDx < i/100.) isValidME31x[i] = true;
              if (sDy < i/100.) isValidME31y[i] = true;
            }
          }
  
          //ME41
          if (csc_id.station() == 4 && csc_id.ring() == 1) {
            for (int i = 0; i < 100; i++) {
              if (sDx < i/100.) isValidME41x[i] = true;
              if (sDy < i/100.) isValidME41y[i] = true;
            }
          }

        }
      }

      for (int i = 0; i < 100; i++){
        if (isValidME31x[i]) sME31x[i]++;
        if (isValidME31y[i]) sME31y[i]++;
        if (isValidME41x[i]) sME41x[i]++;
        if (isValidME41y[i]) sME41y[i]++;
      }

      if (!csc_simmatched && abs(cptype) != 13) continue;

      double xslope = gp_cscint.x()/gp_cscint.z();
      double yslope = gp_cscint.y()/gp_cscint.z();

      if (csc_id.station() == 3) b_S3NDigis++;
      if (csc_id.station() == 4) b_S4NDigis++;

      if (csc_id.station() == 3 && csc_id.ring() == 1){
        b_ME31NDigis = b_ME31NDigis + 1;
        b_ME31NDigis_Total = b_ME31NDigis_Total + 1;
        if (b_cscBX == 6) bx_ME31NDigis++;
      }
      else if (csc_id.station() == 4 && csc_id.ring() == 1){
        b_ME41NDigis = b_ME41NDigis + 1;
        b_ME41NDigis_Total = b_ME41NDigis_Total + 1;
        if (b_cscBX == 6) bx_ME41NDigis++;
      }

      for (int i=0; i<25; i++){
        for (int j=0; j<25; j++){
          isMatchME31[i][j] = false;
          isMatchME41[i][j] = false;
        }
      }

      for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

        RPCDetId rpcid = (RPCDetId)(*rpcIt).rpcId();
        if (rpcid.region() == 0) continue; //skip the barrels

        GlobalPoint gp_rpc(0.0,0.0,0.0);
        gp_rpc = getRPCGlobalPosition(rpcid, *rpcIt);

        if (gp_rpc.z() * gp_cscint.z() < 0 ) continue;

        double dz = gp_rpc.z() - gp_cscint.z();
        double dx = dz*xslope;
        double dy = dz*yslope;

        GlobalPoint gp_transcsc(gp_cscint.x()+dx, gp_cscint.y()+dy, gp_rpc.z());
        LocalPoint lp_extrapol = rpcGeo->idToDet(rpcid)->surface().toLocal(gp_transcsc);

        //local distance
        LocalPoint lp_rpc(0.0,0.0,0.0);
        lp_rpc = (*rpcIt).localPosition();
        float Dx = abs(lp_rpc.x() - lp_extrapol.x());
        float Dy = abs(lp_rpc.y() - lp_extrapol.y());

        //global distance
//        float Dx = abs(gp_rpc.x()-gp_cscint.x());
//        float Dy = abs(gp_rpc.y()-gp_cscint.y());

        if (csc_id.station() == 3 && csc_id.ring() == 1 && rpcid.station() == 3 && rpcid.ring() == 1){
          for (int i = 0; i < 25; i++){
            for (int j = 0; j < 25; j++){
              if (Dx < i+1 && Dy < j+1) isMatchME31[i][j] = true;
            }
          }
        }
        if (csc_id.station() == 4 && csc_id.ring() == 1 && rpcid.station() == 4 && rpcid.ring() == 1){
          for (int i = 0; i < 25; i++){
            for (int j = 0; j < 25; j++){
              if (Dx < i+1 && Dy < j+1) isMatchME41[i][j] = true;
            }
          }
        }

      }//RPCRecHit loop

      for (int i = 0; i < 25; i++){
        for (int j = 0; j < 25; j++){
         
          if (isMatchME31[i][j]) ME31[i][j]++;
          if (isMatchME41[i][j]) ME41[i][j]++;
          
        }
      }

    }//CSCLCT loop

    if (b_ME31NDigis != 0) h_ME31NDigis->Fill(b_ME31NDigis);
    if (b_ME41NDigis != 0) h_ME41NDigis->Fill(b_ME41NDigis);

    if (b_ME31NSimHits != 0) h_cscME31NSimHits->Fill(b_ME31NSimHits);
    if (b_ME41NSimHits != 0) h_cscME41NSimHits->Fill(b_ME41NSimHits);
    if (b_RE31NSimHits != 0) h_rpcME31NSimHits->Fill(b_RE31NSimHits);
    if (b_RE41NSimHits != 0) h_rpcME41NSimHits->Fill(b_RE41NSimHits);

    if (bx_ME31NDigis != 0) h_ME31NDigis0->Fill(bx_ME31NDigis);
    if (bx_ME41NDigis != 0) h_ME41NDigis0->Fill(bx_ME41NDigis);

    if (b_S3NDigis != 0) h_S3NDigis->Fill(b_S3NDigis);
    if (b_S4NDigis != 0) h_S4NDigis->Fill(b_S4NDigis);

  }//CSCChamber loop

  NRecHits->Fill(nRPC);
  h_S3NRecHits->Fill(b_S3NRecHits);
  h_S4NRecHits->Fill(b_S4NRecHits);
  h_RE31NRecHits->Fill(b_RE31NRecHits);
  h_RE41NRecHits->Fill(b_RE41NRecHits);
  h_RE31NRecHits0->Fill(bx_RE31NRecHits);
  h_RE41NRecHits0->Fill(bx_RE41NRecHits);
*/
  tree->Fill();
  EventInfo->Fill(1.5);

}

void 
DTRPCTiming::beginJob()
{

  EventNum = 0;
/*
  tree->Branch("EVENT", &b_EVENT, "EVENT/i");
  tree->Branch("RUN"  , &b_RUN  , "RUN/i");
  tree->Branch("LUMI" , &b_LUMI , "LUMI/i");

  tree->Branch("ME31NDigis" , &b_ME31NDigis , "ME31NDigis/i");
  tree->Branch("ME41NDigis" , &b_ME41NDigis , "ME41NDigis/i");

  tree->Branch("RE31NRecHits" , &b_RE31NRecHits , "RE31NRecHits/i");
  tree->Branch("RE41NRecHits" , &b_RE41NRecHits , "RE41NRecHits/i");

  tree->Branch("ME31NSimHits" , &b_ME31NSimHits , "ME31NSimHits/i");
  tree->Branch("ME41NSimHits" , &b_ME41NSimHits , "ME41NSimHits/i");

  tree->Branch("RE31NSimHits" , &b_RE31NSimHits , "RE31NSimHits/i");
  tree->Branch("RE41NSimHits" , &b_RE41NSimHits , "RE41NSimHits/i");

  tree->Branch("cscBX" , &b_cscBX , "cscBX/i");
  tree->Branch("rpcBX" , &b_rpcBX , "rpcBX/i");

  tree->Branch("ptype" , &b_ptype , "ptype/i");

  b_ME31NDigis_Total = b_ME41NDigis_Total = 0;
  pure_ME31NDigis_Total = pure_ME41NDigis_Total = 0;

  for (int i=0; i<25; i++){
    for (int j=0; j<25; j++){
      ME31[i][j] = 0;
      ME41[i][j] = 0;
    }
  }

  for (int i=0; i<100; i++){
      sME31x[i] = 0;
      sME31y[i] = 0;
      sME41x[i] = 0;
      sME41y[i] = 0;
  }
*/
}

void 
DTRPCTiming::endJob() 
{
/*
  if (b_ME31NDigis_Total != 0){
    for (int i=0; i<25; i++){
      h_xNMatchedME31->SetBinContent(i+1, ME31[i][24]/b_ME31NDigis_Total*100);
      h_xNMatchedME41->SetBinContent(i+1, ME41[i][24]/b_ME41NDigis_Total*100);
    }
  }
  if (b_ME41NDigis_Total != 0){
    for (int j=0; j< 25; j++){
      h_yNMatchedME31->SetBinContent(j+1, ME31[24][j]/b_ME31NDigis_Total*100);
      h_yNMatchedME41->SetBinContent(j+1, ME41[24][j]/b_ME41NDigis_Total*100);
    }
  }
    
  if (b_ME31NDigis_Total != 0 && b_ME41NDigis_Total != 0){
    for (int i=0; i<25; i++){
      for (int j=0; j<25; j++){
        h_MatchedME31->SetBinContent(i+1,j+1,ME31[i][j]/b_ME31NDigis_Total*100);
        h_MatchedME41->SetBinContent(i+1,j+1,ME41[i][j]/b_ME41NDigis_Total*100);  
      }
    } 
  }

  tree->Fill();  

  cout << "Including simhit" << endl;
  cout << "matched ratio at 13cm * 9cm ME31 " << ME31[12][8]/b_ME31NDigis_Total*100 << endl;
  cout << "matched ratio at 13cm * 9cm ME41 " << ME41[12][8]/b_ME41NDigis_Total*100 << endl;

  cout << "\nPure Total # of Digis: ME31 " << pure_ME31NDigis_Total << ":: ME41 " << pure_ME41NDigis_Total << endl;
  cout << "simmatching Total # of Digis: ME31 " << b_ME31NDigis_Total << ":: ME41 " << b_ME41NDigis_Total << endl;

  if (pure_ME31NDigis_Total != 0 && pure_ME41NDigis_Total != 0){
    for (int i=0; i<25; i++){
      for (int j=0; j<25; j++){
        h_RatioME31->SetBinContent(i+1,j+1,ME31[i][j]/pure_ME31NDigis_Total*100);
        h_RatioME41->SetBinContent(i+1,j+1,ME41[i][j]/pure_ME41NDigis_Total*100);  
      }
    } 
  }

  if (pure_ME31NDigis_Total != 0 && pure_ME41NDigis_Total != 0){
    for (int i=0; i<100; i++){
      //sim validation
      h_simValidME31x->SetBinContent(i+1,sME31x[i]/pure_ME31NDigis_Total*100);
      h_simValidME31y->SetBinContent(i+1,sME31y[i]/pure_ME31NDigis_Total*100);
      h_simValidME41x->SetBinContent(i+1,sME41x[i]/pure_ME31NDigis_Total*100);
      h_simValidME41y->SetBinContent(i+1,sME41y[i]/pure_ME31NDigis_Total*100);
    }
  }
*/
}

void
DTRPCTiming::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(DTRPCTiming);

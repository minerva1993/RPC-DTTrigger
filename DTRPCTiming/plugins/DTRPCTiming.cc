// -*- C++ -*-
//
// Package:    RPC+CSCTrigger/DTRPCTiming
// Class:      DTRPCTiming
//
// IMPORTANT COMMENT: 'Digi' is for DT, even though currently we are using segment
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

    //DT Digi (using segments, instead)
    TH1D *h_MBNDigis[4]; //Station
    TH1D *h_WNDigis[5]; //Wheel
    TH2D *h_SWNDigis; //Station - wheel 

    //RPC
    TH1D *h_NRecHits;
    TH1D *h_RBNRecHits[6]; //Station, separating in/out in RB1,2
    TH1D *h_WNRecHits[5]; //Wheel
    TH2D *h_SWNRecHits;

    //No Bx=0 histos / simHit multiplicity now

    unsigned int b_EVENT, b_RUN, b_LUMI;

    //DT
    double b_DTNDigis[4][5]; //0~3: station 1~4, 0~4: wheel -2~2
    //Not looking at BX now
    //unsigned int bx_DTNDigis
    double b_DTNDigis_Total[4][5];
    double pure_DTNDigis_Total[4][5];

    //RPC
    unsigned int b_RPCNRecHits[6][5]; //0~5: station 1~2(in/out)~3~4, 0~4: wheel -2~2
    //Not looking at BX now
    //unsigned int bx_RPCNRecHits

/*  No simhit and matching for now
    unsigned int b_DTNSimHits;
    unsigned int b_RPCNSimHits;

    int b_ptype;
    TH1D *h_ptype;
*/
    int nRPC;
    int nDT;
/*
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

    edm::EDGetTokenT<DTRecSegment4DCollection> dt4DSegments;
    edm::Handle<DTRecSegment4DCollection> all4DSegments;

    edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
    edm::Handle<RPCRecHitCollection> rpcRecHits;

    //GEANT4 simhits
    edm::EDGetTokenT<edm::PSimHitContainer> DTsimHitToken;

    GlobalPoint getDTGlobalPosition(DTChamberId rawId, const DTRecSegment4D& dt4DIt) const;
    GlobalPoint getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const;
    //std::pair<RPCRecHitCollection::const_iterator, float*> matchingRPC(CSCDetId rawId, const CSCCorrelatedLCTDigi& lct, float dx_cutoff, float dy_cutoff) const;
};

GlobalPoint
DTRPCTiming::getDTGlobalPosition(DTChamberId rawId, const DTRecSegment4D& dt4DIt) const{

  DTChamberId dtId = DTChamberId(rawId);
  const LocalPoint& dt_lp = dt4DIt.localPosition();
  const GlobalPoint& dt_gp = dtGeo->idToDet(dtId)->surface().toGlobal(dt_lp);

  return dt_gp;

}

GlobalPoint
DTRPCTiming::getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const{

  RPCDetId rpc_id = RPCDetId(rpcId);
  const LocalPoint& rpc_lp = rpcIt.localPosition();
  const GlobalPoint& rpc_gp = rpcGeo->idToDet(rpc_id)->surface().toGlobal(rpc_lp);
 
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
    RPCDetId rpc_id = (RPCDetId)(*rpcIt).rpcId();

    GlobalPoint gp_rpc(0.0,0.0,0.0);
    gp_rpc = getRPCGlobalPosition(rpc_id, *rpcIt);

    if (rpc_id.region() == 0) continue; //skip the barrels

    float Dx = abs(gp_rpc.x()-gp_cscint.x());
    float Dy = abs(gp_rpc.y()-gp_cscint.y());
    float distance = sqrt(Dx*Dx + Dy*Dy);

    if (rawId.station() == 3 && rawId.ring() == 1 && rpc_id.station() == 3 && rpc_id.ring() == 1){
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

  //DT
  for(int i=0; i<4; i++){
    h_MBNDigis[i] = fs->make<TH1D>(Form("h_MB%iNDigis",i+1), Form("Number of digi per chamber (MB%i)",i+1), 10, 0, 10);
    h_MBNDigis[i]->GetXaxis()->SetTitle("Number of digi per chamber");
    h_MBNDigis[i]->GetYaxis()->SetTitle("Number of chamber");
  }
  for(int i=0; i<5; i++){
    h_WNDigis[i] = fs->make<TH1D>(Form("h_W%iNDigis",i-2), Form("Number of digi per chamber (W%i)",i-2), 10, 0, 10);
    h_WNDigis[i]->GetXaxis()->SetTitle("Number of digi per chamber");
    h_WNDigis[i]->GetYaxis()->SetTitle("Number of chamber");
  }

  h_SWNDigis = fs->make<TH2D>("h_SWNDigis", "Number of digi per chamber", 4, 0.5, 4.5, 5, -2.5, 2.5);
  h_SWNDigis->GetXaxis()->SetTitle("Station");
  h_SWNDigis->GetYaxis()->SetTitle("Wheel");
  h_SWNDigis->GetXaxis()->SetBinLabel(1,"MB1");
  h_SWNDigis->GetXaxis()->SetBinLabel(2,"MB2");
  h_SWNDigis->GetXaxis()->SetBinLabel(3,"MB3");
  h_SWNDigis->GetXaxis()->SetBinLabel(4,"MB4");
  h_SWNDigis->GetYaxis()->SetBinLabel(1,"W-2");
  h_SWNDigis->GetYaxis()->SetBinLabel(2,"W-1");
  h_SWNDigis->GetYaxis()->SetBinLabel(3,"W0");
  h_SWNDigis->GetYaxis()->SetBinLabel(4,"W+1");
  h_SWNDigis->GetYaxis()->SetBinLabel(5,"W+2");

  //RPC
  h_NRecHits = fs->make<TH1D>("h_NRecHits", "", 10, 0, 10);
  h_NRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_NRecHits->GetYaxis()->SetTitle("Number of chamber");

  for(int i=0; i<6; i++){
    const char* inout = "";
    int st = 0;
    if(i < 4){
      if(i==0 or i==2) inout = "in";
      else             inout = "out";
      if(i<2) st = 1;
      else    st = 2;
    }
    else st = i-1;

    h_RBNRecHits[i] = fs->make<TH1D>(Form("h_RB%i%sNRecHits",st,inout), Form("Number of rechit per chamber (RB%i%s)",st,inout), 10, 0, 10);
    h_RBNRecHits[i]->GetXaxis()->SetTitle("Number of rechit per chamber");
    h_RBNRecHits[i]->GetYaxis()->SetTitle("Number of chamber");
  }
  for(int i=0; i<5; i++){
    h_WNRecHits[i] = fs->make<TH1D>(Form("h_W%iNRecHits",i-2), Form("Number of rechit per chamber (W%i)",i-2), 10, 0, 10);
    h_WNRecHits[i]->GetXaxis()->SetTitle("Number of rechit per chamber");
    h_WNRecHits[i]->GetYaxis()->SetTitle("Number of chamber");
  }
  h_SWNRecHits = fs->make<TH2D>("h_SWNRecHits", "Number of rechit per chamber", 6, 0.5, 6.5, 5, -2.5, 2.5);
  h_SWNRecHits->GetXaxis()->SetTitle("Station");
  h_SWNRecHits->GetYaxis()->SetTitle("Wheel");
  h_SWNRecHits->GetXaxis()->SetBinLabel(1,"RB1in");
  h_SWNRecHits->GetXaxis()->SetBinLabel(2,"RB1out");
  h_SWNRecHits->GetXaxis()->SetBinLabel(3,"RB2in");
  h_SWNRecHits->GetXaxis()->SetBinLabel(4,"RB2out");
  h_SWNRecHits->GetXaxis()->SetBinLabel(5,"RB3");
  h_SWNRecHits->GetXaxis()->SetBinLabel(6,"RB4");
  h_SWNRecHits->GetYaxis()->SetBinLabel(1,"W-2");
  h_SWNRecHits->GetYaxis()->SetBinLabel(2,"W-1");
  h_SWNRecHits->GetYaxis()->SetBinLabel(3,"W0");
  h_SWNRecHits->GetYaxis()->SetBinLabel(4,"W+1");
  h_SWNRecHits->GetYaxis()->SetBinLabel(5,"W+2");

/*
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
*/
  if (!all4DSegments.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find DTRecSegment4D with label "<< all4DSegments << std::endl;
    return;
  }
  if (!rpcRecHits.isValid()) {
    edm::LogInfo("DataNotFound") << "can't find RPCRecHitDigiCollection with label "<< rpcRecHits << std::endl;
    return;
  }

  //cout << "\nNew event" << endl;
  b_EVENT = b_RUN = b_LUMI = 0;

  b_EVENT  = iEvent.id().event();
  b_RUN    = iEvent.id().run();
  b_LUMI   = iEvent.id().luminosityBlock();

  EventNum++;
  //cout << "Event " << EventNum << endl;

  nRPC = nDT = 0;
  for(int i=0; i<6; i++){
    for(int j=0; j<5; j++) b_RPCNRecHits[i][j]=0;
  } 

  //b_ME31NSimHits = b_ME41NSimHits = b_RE31NSimHits = b_RE41NSimHits = 0;

  //to check rechit info
  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

    /// Ring id: Wheel number in Barrel (from -2 to +2) Ring Number in Endcap (from 1 to 3)
    /// Station id : For Barrel: the four groups of chambers at same r (distance from beam axis) and increasing phi
    ///              For Endcap: the three groups of chambers at same z (distance from interaction point), i.e. the disk
    /// Layer id: each station can have two layers of chambers: layer 1 is the inner chamber and layer 2 is the outer chamber (when present)
    /// Only in Barrel: RB1 and RB2.
    RPCDetId rpc_id = (RPCDetId)(*rpcIt).rpcId();
    if (rpc_id.region() != 0) continue; //skip the endcap
    nRPC++;
    int idxRPCRing = rpc_id.ring() + 2; //0~4
    int idxRPCStation = -1;
    if(rpc_id.station() <= 2) idxRPCStation = (rpc_id.station()-1)*(rpc_id.station()) + (rpc_id.layer()-1); //0~3
    else idxRPCStation = rpc_id.station() + 1; 

    h_SWNRecHits->Fill(idxRPCStation+1, rpc_id.ring(), 1); //Want to draw ring -2~2
    b_RPCNRecHits[idxRPCStation][idxRPCRing]++;
    //if((*rpcIt).BunchX() == 0 && rpc_id.station() == 3 && rpc_id.ring() == 1) bx_RE31NRecHits++;
    //if((*rpcIt).BunchX() == 0 && rpc_id.station() == 4 && rpc_id.ring() == 1) bx_RE41NRecHits++;
  }
  /*
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
  */
  //to check rpcsimhit Info

  //std::cout << "\t Number of DT Segments in this event = " << all4DSegments->size() << std::endl;

  //Ref of loop: https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/RecoLocalMuon/DTSegment/src/DTSegment4DT0Corrector.cc
  DTRecSegment4DCollection::id_iterator dtChamberId;
  for(dtChamberId = all4DSegments->id_begin(); dtChamberId != all4DSegments->id_end(); ++dtChamberId) {

    nDT++;
    for(int i=0; i<4; i++){
      for(int j=0; j<5; j++) b_DTNDigis[i][j]=0;
    }
    //bx_ME31NDigis = bx_ME41NDigis =0;
    //b_cscBX = b_rpcBX = 0;

    GlobalPoint gp_dt;
    /*
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
    */

    DTRecSegment4DCollection::range dtRange = all4DSegments->get(*dtChamberId);
    for(DTRecSegment4DCollection::const_iterator dtSegment = dtRange.first; dtSegment != dtRange.second; ++dtSegment){

      DTRecSegment4D tmpseg = *dtSegment;
      if(!tmpseg.hasZed()) continue; //Looking at the segment that has Z projection FIXME Check needed
      if(dtSegment->dimension() != 4) continue; //Check needed //MB4 doesn't have theta superlayer
      const DTChamberId dt_id = dtSegment->chamberId();
      int idxDTStation = dt_id.station()-1;
      int idxDTWheel = dt_id.wheel() + 2;
      //b_cscBX = lct->getBX();

      gp_dt = GlobalPoint(0.0,0.0,0.0);
      gp_dt = getDTGlobalPosition(dt_id, tmpseg);
      if (abs(gp_dt.eta()) > 1.2) continue; //FIXME No RPC, DT

      float Zo = gp_dt.z();
      //cout << Zo << " / ";

      pure_DTNDigis_Total[idxDTStation][idxDTWheel]++;
      /*
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
      */

      //No sim match, should be same as pure_*, no bx for now
      h_SWNDigis->Fill(dt_id.station(), dt_id.wheel(), 1);
      b_DTNDigis_Total[idxDTStation][idxDTWheel]++;
      b_DTNDigis[idxDTStation][idxDTWheel]++;
      /*
      for (int i=0; i<25; i++){
        for (int j=0; j<25; j++){
          isMatchME31[i][j] = false;
          isMatchME41[i][j] = false;
        }
      }

      for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {

        RPCDetId rpc_id = (RPCDetId)(*rpcIt).rpcId();
        if (rpc_id.region() == 0) continue; //skip the barrels

        GlobalPoint gp_rpc(0.0,0.0,0.0);
        gp_rpc = getRPCGlobalPosition(rpc_id, *rpcIt);

        if (gp_rpc.z() * gp_cscint.z() < 0 ) continue;

        double dz = gp_rpc.z() - gp_cscint.z();
        double dx = dz*xslope;
        double dy = dz*yslope;

        GlobalPoint gp_transcsc(gp_cscint.x()+dx, gp_cscint.y()+dy, gp_rpc.z());
        LocalPoint lp_extrapol = rpcGeo->idToDet(rpc_id)->surface().toLocal(gp_transcsc);

        //local distance
        LocalPoint lp_rpc(0.0,0.0,0.0);
        lp_rpc = (*rpcIt).localPosition();
        float Dx = abs(lp_rpc.x() - lp_extrapol.x());
        float Dy = abs(lp_rpc.y() - lp_extrapol.y());

        //global distance
//        float Dx = abs(gp_rpc.x()-gp_cscint.x());
//        float Dy = abs(gp_rpc.y()-gp_cscint.y());

        if (csc_id.station() == 3 && csc_id.ring() == 1 && rpc_id.station() == 3 && rpc_id.ring() == 1){
          for (int i = 0; i < 25; i++){
            for (int j = 0; j < 25; j++){
              if (Dx < i+1 && Dy < j+1) isMatchME31[i][j] = true;
            }
          }
        }
        if (csc_id.station() == 4 && csc_id.ring() == 1 && rpc_id.station() == 4 && rpc_id.ring() == 1){
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
      */
    }//CSCLCT loop

    for(int i=0; i<4; i++)
      h_MBNDigis[i]->Fill(b_DTNDigis[i][0]+b_DTNDigis[i][1]+b_DTNDigis[i][2]+b_DTNDigis[i][3]+b_DTNDigis[i][4]);
    for(int j=0; j<5; j++)
      h_WNDigis[j]->Fill(b_DTNDigis[0][j]+b_DTNDigis[1][j]+b_DTNDigis[2][j]+b_DTNDigis[3][j]);
    /*
    if (b_ME31NSimHits != 0) h_cscME31NSimHits->Fill(b_ME31NSimHits);
    if (b_ME41NSimHits != 0) h_cscME41NSimHits->Fill(b_ME41NSimHits);
    if (b_RE31NSimHits != 0) h_rpcME31NSimHits->Fill(b_RE31NSimHits);
    if (b_RE41NSimHits != 0) h_rpcME41NSimHits->Fill(b_RE41NSimHits);

    if (bx_ME31NDigis != 0) h_ME31NDigis0->Fill(bx_ME31NDigis);
    if (bx_ME41NDigis != 0) h_ME41NDigis0->Fill(bx_ME41NDigis);
    */

  }//CSCChamber loop

  h_NRecHits->Fill(nRPC);
  for(int i=0; i<6; i++)
    h_RBNRecHits[i]->Fill(b_RPCNRecHits[i][0]+b_RPCNRecHits[i][1]+b_RPCNRecHits[i][2]+b_RPCNRecHits[i][3]+b_RPCNRecHits[i][4]);
  for(int j=0; j<5; j++)
    h_WNRecHits[j]->Fill(b_RPCNRecHits[0][j]+b_RPCNRecHits[1][j]+b_RPCNRecHits[2][j]+b_RPCNRecHits[3][j]+b_RPCNRecHits[4][j]+b_RPCNRecHits[5][j]);

  tree->Fill();
  EventInfo->Fill(1.5);

}

void 
DTRPCTiming::beginJob()
{

  EventNum = 0;

  tree->Branch("EVENT", &b_EVENT, "EVENT/i");
  tree->Branch("RUN"  , &b_RUN  , "RUN/i");
  tree->Branch("LUMI" , &b_LUMI , "LUMI/i");
/*
  tree->Branch("MB1NDigis", &b_MB1NDigis, "MB1NDigis/i");
  tree->Branch("MB2NDigis", &b_MB2NDigis, "MB2NDigis/i");
  tree->Branch("MB3NDigis", &b_MB3NDigis, "MB3NDigis/i");
  tree->Branch("MB4NDigis", &b_MB4NDigis, "MB4NDigis/i");
  tree->Branch("WNeg2NDigis", &b_WNeg2NDigis, "WNeg2NDigis/i");
  tree->Branch("WNeg1NDigis", &b_WNeg1NDigis, "WNeg1NDigis/i");
  tree->Branch("W0NDigis",    &b_W0NDigis,    "W0NDigis/i");
  tree->Branch("WPos1NDigis", &b_WPos1NDigis, "WPos1NDigis/i");
  tree->Branch("WPos2NDigis", &b_WPos2NDigis, "WPos2NDigis/i");

  tree->Branch("MB1NRecHits", &b_RB1NRecHits, "MB1NRecHits/i");
  tree->Branch("MB2NRecHits", &b_RB2NRecHits, "MB2NRecHits/i");
  tree->Branch("MB3NRecHits", &b_RB3NRecHits, "MB3NRecHits/i");
  tree->Branch("MB4NRecHits", &b_RB4NRecHits, "MB4NRecHits/i");
  tree->Branch("WNeg2NRecHits", &b_WNeg2NRecHits, "WNeg2NRecHits/i");
  tree->Branch("WNeg1NRecHits", &b_WNeg1NRecHits, "WNeg1NRecHits/i");
  tree->Branch("W0NRecHits",    &b_W0NRecHits,    "W0NRecHits/i");
  tree->Branch("WPos1NRecHits", &b_WPos1NRecHits, "WPos1NRecHits/i");
  tree->Branch("WPos2NRecHits", &b_WPos2NRecHits, "WPos2NRecHits/i");

  tree->Branch("ME31NSimHits" , &b_ME31NSimHits , "ME31NSimHits/i");
  tree->Branch("ME41NSimHits" , &b_ME41NSimHits , "ME41NSimHits/i");

  tree->Branch("RE31NSimHits" , &b_RE31NSimHits , "RE31NSimHits/i");
  tree->Branch("RE41NSimHits" , &b_RE41NSimHits , "RE41NSimHits/i");

  tree->Branch("cscBX" , &b_cscBX , "cscBX/i");
  tree->Branch("rpcBX" , &b_rpcBX , "rpcBX/i");

  tree->Branch("ptype" , &b_ptype , "ptype/i");
*/
    for(int i=0; i<4; i++){
      for(int j=0; j<5; j++) b_DTNDigis_Total[i][j]=0;
    }
  //pure_ME31NDigis_Total = pure_ME41NDigis_Total = 0;
/*
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

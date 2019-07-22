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
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

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
    TH1D *h_DTBNDigis[4]; //Station
    TH1D *h_DTWNDigis[5]; //Wheel
    TH2D *h_DTSWNDigis; //Station - wheel
    TH1D *h_DTBNSimHits[4];
    TH1D *h_DTWNSimHits[5];
    TH2D *h_DTSWNSimHits;

    //RPC
    TH1D *h_RPCNRecHits;
    TH1D *h_RPCBNRecHits[6]; //Station, separating in/out in RB1,2
    TH1D *h_RPCWNRecHits[5]; //Wheel
    TH2D *h_RPCSWNRecHits;
    TH1D *h_RPCTimeRes;
    TH1D *h_RPCNonMuTimeRes;

    //N simHit per digi; in 1 digi case, count numbers in same chamber
    TH1D *h_NSPD;

    //No Bx=0 histos / simHit multiplicity now

    TH1D *h_DTBsimValidx[4];
    TH1D *h_DTBsimValidy[4];
    TH1D *h_DTWsimValidx[5];
    TH1D *h_DTWsimValidy[5];
/*
    TH1D *h_xNMatchedME31;
    TH1D *h_xNMatchedME41;

    TH1D *h_yNMatchedME31;
    TH1D *h_yNMatchedME41;
    TH2D *h_MatchedME31;
    TH2D *h_MatchedME41;

    TH2D *h_RatioME31;
    TH2D *h_RatioME41;
*/

    //Declare variables
    unsigned int b_EVENT, b_RUN, b_LUMI;

    //DT
    double b_DTNDigis[4][5]; //0~3: station 1~4, 0~4: wheel -2~2
    //Not looking at BX now
    //unsigned int bx_DTNDigis
    double b_DTNDigis_Total[4][5];
    double pure_DTNDigis_Total[4][5];
    unsigned int b_DTNSimHits[4][5];

    //RPC
    unsigned int b_RPCNRecHits[6][5]; //0~5: station 1~2(in/out)~3~4, 0~4: wheel -2~2
    //Not looking at BX now
    //unsigned int bx_RPCNRecHits
    unsigned int b_RPCNSimHits;

    int b_ptype;
    TH1D *h_ptype;

    int nRPC;
    int nDT;
/*
    int b_rpcBX;
    int b_cscBX;
*/
    double sDTx[4][5][100];//station-wheel-distance
    double sDTy[4][5][100];
    bool DTisValidx[4][5][100];
    bool DTisValidy[4][5][100];
/*
    double ME31[25][25];
    double ME41[25][25];

    bool isMatchME31[25][25];
    bool isMatchME41[25][25];
*/
    int EventNum;
    int label_;
    int numDigi_switch;

    edm::ESHandle<DTGeometry> dtGeo;
    edm::ESHandle<RPCGeometry> rpcGeo;

    edm::EDGetTokenT<DTRecSegment4DCollection> dt4DSegments_;
    edm::Handle<DTRecSegment4DCollection> all4DSegments;

    edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
    edm::Handle<RPCRecHitCollection> rpcRecHits;

    //GEANT4 simhits
    edm::EDGetTokenT<edm::PSimHitContainer> DTsimHitToken_;
    edm::Handle<PSimHitContainer> DTsimHit;
    edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> RPCdigisimlinkToken_;
    edm::Handle<edm::DetSetVector<RPCDigiSimLink>> thelinkDigis;

    GlobalPoint getDTGlobalPosition(DTChamberId rawId, const DTRecSegment4D& dt4DIt) const;
    GlobalPoint getRPCGlobalPosition(RPCDetId rpcId, const RPCRecHit& rpcIt) const;
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

DTRPCTiming::DTRPCTiming(const edm::ParameterSet& iConfig)
{
  dt4DSegments_ = consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dt4DSegments"));
  auto RPCDigiLabel = iConfig.getParameter<edm::InputTag>("simMuonRPCDigis");
  rpcRecHitsToken_ = consumes<RPCRecHitCollection>(edm::InputTag(RPCDigiLabel.label(), "" ));
  DTsimHitToken_ = consumes<PSimHitContainer>(iConfig.getUntrackedParameter<edm::InputTag>("DTsimHitLabel", edm::InputTag("g4SimHits:MuonDTHits")));
  RPCdigisimlinkToken_ = consumes<edm::DetSetVector<RPCDigiSimLink>>(iConfig.getParameter<edm::InputTag>("rpcSimLinkLabel"));

  //numDigi label
  label_ = iConfig.getUntrackedParameter<int>("label");

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "Tree for RPC+CSCTrigger");

  EventInfo = fs->make<TH1D>("EventInfo","Event Information",2,0,2);
  EventInfo->GetXaxis()->SetBinLabel(1,"Total Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Selected Number of Events");

  //DT
  for(int i=0; i<4; i++){
    h_DTBNDigis[i] = fs->make<TH1D>(Form("h_DTMB%iNDigis",i+1), Form("Number of digi per chamber (MB%i)",i+1), 10, 0, 10);
    h_DTBNDigis[i]->GetXaxis()->SetTitle("Number of digi per chamber");
    h_DTBNDigis[i]->GetYaxis()->SetTitle("Number of chamber");

    h_DTBNSimHits[i] = fs->make<TH1D>(Form("h_DTMB%iNSimHits",i+1), Form("Number of DT simhit per chamber (MB%i)",i+1), 10, 0, 10);
    h_DTBNSimHits[i]->GetXaxis()->SetTitle("Number of simhit per chamber");
    h_DTBNSimHits[i]->GetYaxis()->SetTitle("Number of chamber");

    h_DTBsimValidx[i] = fs->make<TH1D>(Form("h_DTMB%isimValidx",i+1), Form("Validation Percentage in MB%i",i+1), 100, 0, 100);
    h_DTBsimValidx[i]->GetXaxis()->SetTitle("X cutoff (mm)");
    h_DTBsimValidx[i]->GetYaxis()->SetTitle("Matched (%)");

    h_DTBsimValidy[i] = fs->make<TH1D>(Form("h_DTMB%isimValidy",i+1), Form("Validation Percentage in MB%i",i+1), 100, 0, 100);
    h_DTBsimValidy[i]->GetXaxis()->SetTitle("Y cutoff (mm)");
    h_DTBsimValidy[i]->GetYaxis()->SetTitle("Matched (%)");
  }
  for(int i=0; i<5; i++){
    h_DTWNDigis[i] = fs->make<TH1D>(Form("h_DTW%iNDigis",i-2), Form("Number of digi per chamber (W%i)",i-2), 10, 0, 10);
    h_DTWNDigis[i]->GetXaxis()->SetTitle("Number of digi per chamber");
    h_DTWNDigis[i]->GetYaxis()->SetTitle("Number of chamber");

    h_DTWNSimHits[i] = fs->make<TH1D>(Form("h_DTW%iNSimHits",i-2), Form("Number of DT simhit per chamber (W%i)",i-2), 10, 0, 10);
    h_DTWNSimHits[i]->GetXaxis()->SetTitle("Number of simhit per chamber");
    h_DTWNSimHits[i]->GetYaxis()->SetTitle("Number of chamber");

    h_DTWsimValidx[i] = fs->make<TH1D>(Form("h_DTW%isimValidx",i-2), Form("Validation Percentage in W%i",i-2), 100, 0, 100);
    h_DTWsimValidx[i]->GetXaxis()->SetTitle("X cutoff (mm)");
    h_DTWsimValidx[i]->GetYaxis()->SetTitle("Matched (%)");

    h_DTWsimValidy[i] = fs->make<TH1D>(Form("h_DTW%isimValidWy",i-2), Form("Validation Percentage in W%i",i-2), 100, 0, 100);
    h_DTWsimValidy[i]->GetXaxis()->SetTitle("Y cutoff (mm)");
    h_DTWsimValidy[i]->GetYaxis()->SetTitle("Matched (%)");
  }

  h_DTSWNDigis = fs->make<TH2D>("h_DTSWNDigis", "Number of digi per chamber", 4, 0.5, 4.5, 5, -2.5, 2.5);
  h_DTSWNDigis->GetXaxis()->SetTitle("Station");
  h_DTSWNDigis->GetYaxis()->SetTitle("Wheel");
  h_DTSWNDigis->GetXaxis()->SetBinLabel(1,"MB1");
  h_DTSWNDigis->GetXaxis()->SetBinLabel(2,"MB2");
  h_DTSWNDigis->GetXaxis()->SetBinLabel(3,"MB3");
  h_DTSWNDigis->GetXaxis()->SetBinLabel(4,"MB4");
  h_DTSWNDigis->GetYaxis()->SetBinLabel(1,"W-2");
  h_DTSWNDigis->GetYaxis()->SetBinLabel(2,"W-1");
  h_DTSWNDigis->GetYaxis()->SetBinLabel(3,"W0");
  h_DTSWNDigis->GetYaxis()->SetBinLabel(4,"W+1");
  h_DTSWNDigis->GetYaxis()->SetBinLabel(5,"W+2");

  //RPC
  h_RPCTimeRes = fs->make<TH1D>("h_RPCTimeRes", "RPC time residual (TOF-time)", 50, 0, 50);
  h_RPCTimeRes->GetXaxis()->SetTitle("Time residual");
  h_RPCTimeRes->GetYaxis()->SetTitle("Number of rechits");

  h_RPCNonMuTimeRes = fs->make<TH1D>("h_RPCNonMuTimeRes", "RPC time residual for non-muon simhits(TOF-time)", 50, 0, 50);
  h_RPCNonMuTimeRes->GetXaxis()->SetTitle("Time residual");
  h_RPCNonMuTimeRes->GetYaxis()->SetTitle("Number of rechits");

  h_RPCNRecHits = fs->make<TH1D>("h_RPCNRecHits", "", 10, 0, 10);
  h_RPCNRecHits->GetXaxis()->SetTitle("Number of rechit per chamber");
  h_RPCNRecHits->GetYaxis()->SetTitle("Number of chamber");

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

    h_RPCBNRecHits[i] = fs->make<TH1D>(Form("h_RPCB%i%sNRecHits",st,inout), Form("Number of rechit per chamber (RB%i%s)",st,inout), 10, 0, 10);
    h_RPCBNRecHits[i]->GetXaxis()->SetTitle("Number of rechit per chamber");
    h_RPCBNRecHits[i]->GetYaxis()->SetTitle("Number of chamber");
  }
  for(int i=0; i<5; i++){
    h_RPCWNRecHits[i] = fs->make<TH1D>(Form("h_RPCW%iNRecHits",i-2), Form("Number of rechit per chamber (W%i)",i-2), 10, 0, 10);
    h_RPCWNRecHits[i]->GetXaxis()->SetTitle("Number of rechit per chamber");
    h_RPCWNRecHits[i]->GetYaxis()->SetTitle("Number of chamber");
  }
  h_RPCSWNRecHits = fs->make<TH2D>("h_RPCSWNRecHits", "Number of rechit per chamber", 6, 0.5, 6.5, 5, -2.5, 2.5);
  h_RPCSWNRecHits->GetXaxis()->SetTitle("Station");
  h_RPCSWNRecHits->GetYaxis()->SetTitle("Wheel");
  h_RPCSWNRecHits->GetXaxis()->SetBinLabel(1,"RB1in");
  h_RPCSWNRecHits->GetXaxis()->SetBinLabel(2,"RB1out");
  h_RPCSWNRecHits->GetXaxis()->SetBinLabel(3,"RB2in");
  h_RPCSWNRecHits->GetXaxis()->SetBinLabel(4,"RB2out");
  h_RPCSWNRecHits->GetXaxis()->SetBinLabel(5,"RB3");
  h_RPCSWNRecHits->GetXaxis()->SetBinLabel(6,"RB4");
  h_RPCSWNRecHits->GetYaxis()->SetBinLabel(1,"W-2");
  h_RPCSWNRecHits->GetYaxis()->SetBinLabel(2,"W-1");
  h_RPCSWNRecHits->GetYaxis()->SetBinLabel(3,"W0");
  h_RPCSWNRecHits->GetYaxis()->SetBinLabel(4,"W+1");
  h_RPCSWNRecHits->GetYaxis()->SetBinLabel(5,"W+2");

  h_NSPD = fs->make<TH1D>("h_NSPD", "number of simHit per digi", 20, 0, 20);
  h_NSPD->GetXaxis()->SetTitle("Number of simHit");
  h_NSPD->GetYaxis()->SetTitle("Number of digi (segment)");

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
*/
  h_ptype = fs->make<TH1D>("h_ptype", "", 4,0,4);
  h_ptype->GetXaxis()->SetBinLabel(1,"Electron");
  h_ptype->GetXaxis()->SetBinLabel(2,"Muon");
  h_ptype->GetXaxis()->SetBinLabel(3,"Pion");
  h_ptype->GetXaxis()->SetBinLabel(4,"Proton");
  h_ptype->GetXaxis()->SetTitle("Particle type");
  h_ptype->GetYaxis()->SetTitle("Entries");
}

DTRPCTiming::~DTRPCTiming()
{

}

void
DTRPCTiming::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  EventInfo->Fill(0.5);

  iSetup.get<MuonGeometryRecord>().get( dtGeo );     
  iSetup.get<MuonGeometryRecord>().get( rpcGeo );     

  iEvent.getByToken(dt4DSegments_, all4DSegments);
  iEvent.getByToken(DTsimHitToken_, DTsimHit);
  iEvent.getByToken(rpcRecHitsToken_, rpcRecHits);
  iEvent.getByToken(RPCdigisimlinkToken_, thelinkDigis);

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
  for(int i=0; i<4; i++){
    for(int j=0; j<5; j++) b_DTNSimHits[i][j]=0;
  }

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

    h_RPCSWNRecHits->Fill(idxRPCStation+1, rpc_id.ring(), 1); //Want to draw ring -2~2
    b_RPCNRecHits[idxRPCStation][idxRPCRing]++;
    //if((*rpcIt).BunchX() == 0 && rpc_id.station() == 3 && rpc_id.ring() == 1) bx_RE31NRecHits++;
    //if((*rpcIt).BunchX() == 0 && rpc_id.station() == 4 && rpc_id.ring() == 1) bx_RE41NRecHits++;

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

      if(links.empty()) cout << "Unmatched simHit!" << endl;

      for(std::set<RPCDigiSimLink>::iterator itlink = links.begin(); itlink != links.end(); ++itlink){
        SimHitIdpr currentId(itlink->getTrackId(), itlink->getEventId());
        if(find(matched.begin(), matched.end(), currentId) == matched.end()) matched.push_back(currentId);
      }
      //cout << matched.size() << " ";
    }
    if(!links.empty()){
      //cout << links.begin()->getTimeOfFlight() << " ";
      if(abs(links.begin()->getParticleType()) != 13){
        h_RPCNonMuTimeRes->Fill(links.begin()->getTimeOfFlight() - rpcIt->time());
      }
      else{
        //https://github.com/cms-sw/cmssw/blob/dc968f763c733222c535f6fe1de69c0e2082ad7c/TrackPropagation/RungeKutta/src/PathToPlane2Order.cc#L51
        if(links.begin()->getBx() != 0 or links.begin()->getMomentumAtEntry().perp() >15) continue;
        h_RPCTimeRes->Fill(links.begin()->getTimeOfFlight() - rpcIt->time());
      }
    }
  }

  //GEANT4 simhits
  PSimHitContainer::const_iterator DTsimIt;
  for (DTsimIt = DTsimHit->begin(); DTsimIt != DTsimHit->end(); DTsimIt++) {
    DTChamberId dt_simId = DTsimIt->detUnitId();
    if(dt_simId.station() == 4) continue; //No theta layer in MB4
    //cout << dt_simId << "/";
    b_ptype = DTsimIt->particleType();
    //cout << "ptype :: " << b_ptype << endl;
    //if ( abs(b_ptype) != 11 && abs(b_ptype) != 13 && abs(b_ptype) != 2212 && abs(b_ptype) != 211 ) cout << "NONEPARTICLE TYPE :: " << abs(b_ptype) << endl;
    
    if (abs(b_ptype) == 11) h_ptype->Fill(0.5);
    if (abs(b_ptype) == 13) h_ptype->Fill(1.5);
    if (abs(b_ptype) == 211) h_ptype->Fill(2.5);
    if (abs(b_ptype) == 2212) h_ptype->Fill(3.5);

    b_DTNSimHits[dt_simId.station()-1][dt_simId.wheel()+2]++;
  }
 
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

    numDigi_switch = label_;
    DTRecSegment4DCollection::range dtRange = all4DSegments->get(*dtChamberId);
    if (numDigi_switch == 1){
      dtRange.first++;
      if (dtRange.first != dtRange.second) continue; //Check only 1 digi(segment) or not
      dtRange.first--;
    }

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
      if(abs(gp_dt.eta()) > 1.2) continue; //FIXME No RPC, DT

      float Zo = gp_dt.z();

      pure_DTNDigis_Total[idxDTStation][idxDTWheel]++;

      const LocalPoint& dt_lp = tmpseg.localPosition();
      int cptype = 0;
      bool dt_simmatched = false;

      for(int i=0; i<4; i++){
        for(int j=0; j < 5; j++){
          for(int k=0; k < 100; k++){
            DTisValidx[i][j][k] = false;
            DTisValidy[i][j][k] = false;
          }
        }
      }

      int NSPD = 0;
      for (DTsimIt = DTsimHit->begin(); DTsimIt != DTsimHit->end(); DTsimIt++) {

        cptype = DTsimIt->particleType();
        //const GlobalPoint sim_gp = cscGeo->idToDet(csc_id)->surface().toGlobal(CSCsimIt->localPosition());
        const LocalPoint lp_dtsim = DTsimIt->localPosition();
        DTChamberId dtsim_id = DTsimIt->detUnitId();

        //sim validation with window
        if (abs(cptype) != 13) continue;
        if (dt_id.wheel() == dtsim_id.wheel() && dt_id.station() == dtsim_id.station() && dt_id.sector() == dtsim_id.sector()){

          //cout << "I'm the same" << endl;
          NSPD++;
          dt_simmatched = true; //Check only if the segment has a simHit in the same chamber
          //if (sqrt(dt_lp.x()-lp_dtsim.x())*(dt_lp.x()-lp_dtsim.x())+(dt_lp.y()-lp_dtsim.y())*(dt_lp.y()-lp_dtsim.y()) < 0.5) dt_simmatched = true;
  
          float sDx = abs(dt_lp.x() - lp_dtsim.x());
          float sDy = abs(dt_lp.y() - lp_dtsim.y());
          for (int k=0; k<100; k++) {
            if (sDx < k/100.) DTisValidx[idxDTStation][idxDTWheel][k] = true;
            if (sDy < k/100.) DTisValidy[idxDTStation][idxDTWheel][k] = true;
          }
        }
      }

      h_NSPD->Fill(NSPD);

      for(int i=0; i<4; i++){
        for(int j=0; j < 5; j++){
          for(int k=0; k < 100; k++){
            if(DTisValidx[i][j][k]) sDTx[i][j][k]++;
            if(DTisValidy[i][j][k]) sDTy[i][j][k]++;
          }
        }
      }

      //if (!dt_simmatched && abs(cptype) != 13) continue;
      if (!dt_simmatched) continue;
      /*
      double xslope = gp_cscint.x()/gp_cscint.z();
      double yslope = gp_cscint.y()/gp_cscint.z();
      */

      //No sim match, should be same as pure_*, no bx for now
      h_DTSWNDigis->Fill(dt_id.station(), dt_id.wheel(), 1);
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

    for(int i=0; i<4; i++){
      double tmp1 = 0; unsigned int tmp2 = 0;
      for(int j=0; j<5; j++){
        tmp1 += b_DTNDigis[i][j];
        tmp2 += b_DTNSimHits[i][j];
      }
      h_DTBNDigis[i]->Fill(tmp1);
      h_DTBNSimHits[i]->Fill(tmp2);
    }
    for(int j=0; j<5; j++){
      double tmp1 = 0; unsigned int tmp2 = 0;
      for(int i=0; i<4; i++){
        tmp1 += b_DTNDigis[i][j];
        tmp2 += b_DTNSimHits[i][j];
      }
      h_DTWNDigis[j]->Fill(tmp1);
      h_DTWNSimHits[j]->Fill(tmp2);
    }
    /*
    if (bx_ME31NDigis != 0) h_ME31NDigis0->Fill(bx_ME31NDigis);
    if (bx_ME41NDigis != 0) h_ME41NDigis0->Fill(bx_ME41NDigis);
    */

  }//CSCChamber loop

  h_RPCNRecHits->Fill(nRPC);
  for(int i=0; i<6; i++){
    unsigned int tmp = 0;
    for(int j=0; j<5; j++) tmp += b_RPCNRecHits[i][j];
    h_RPCBNRecHits[i]->Fill(tmp);
  }
  for(int j=0; j<5; j++){
    unsigned int tmp = 0;
    for(int i=0; i<6; i++) tmp += b_RPCNRecHits[i][j];
    h_RPCWNRecHits[j]->Fill(tmp);
  }

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
    for(int j=0; j<5; j++){
      b_DTNDigis_Total[i][j]=0;
      pure_DTNDigis_Total[i][j]=0;

      for(int k=0; k < 100; k++){
        sDTx[i][j][k] = 0;
        sDTy[i][j][k] = 0;
      }
    }
  }
/*
  for (int i=0; i<25; i++){
    for (int j=0; j<25; j++){
      ME31[i][j] = 0;
      ME41[i][j] = 0;
    }
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
*/
  for (int k=0; k<100; k++){
    //sim validation
    for(int i=0; i<4; i++){
      double tmp1 = 0; double tmp2 = 0; double tmp3 = 0;
      for(int j=0; j<5; j++){
        tmp1 += sDTx[i][j][k];
        tmp2 += sDTy[i][j][k];
        tmp3 += pure_DTNDigis_Total[i][j];
      }
      h_DTBsimValidx[i]->SetBinContent(k+1,tmp1/tmp3*100);
      h_DTBsimValidy[i]->SetBinContent(k+1,tmp2/tmp3*100);
    }
    for(int j=0; j<5; j++){
      double tmp1 = 0; double tmp2 = 0; double tmp3 = 0;
      for(int i=0; i<4; i++){
        tmp1 += sDTx[i][j][k];
        tmp2 += sDTy[i][j][k];
        tmp3 += pure_DTNDigis_Total[i][j];
      }
      h_DTWsimValidx[j]->SetBinContent(k+1,tmp1/tmp3*100);
      h_DTWsimValidy[j]->SetBinContent(k+1,tmp2/tmp3*100);
    }
  }
}

void
DTRPCTiming::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(DTRPCTiming);

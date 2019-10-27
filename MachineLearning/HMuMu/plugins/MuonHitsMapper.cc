// -*- C++ -*-
//
// Package:    MachineLearning/MuonHitsMapper
// Class:      MuonHitsMapper
// 
/**\class MuonHitsMapper MuonHitsMapper.cc MachineLearning/MuonHitsMapper/plugins/MuonHitsMapper.cc

 Description: Information about muons and hits around them

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Adish Pradeep Vartak
//         Created:  Tue, 22 Oct 2019 12:34:12 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>

// ROOT include files
#include <TTree.h>
#include <TVector3.h>

// CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetEnumerators.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/Records/interface/TkStripCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h" 

class MuonHitsMapper : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit MuonHitsMapper(const edm::ParameterSet&);
        ~MuonHitsMapper();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        bool addEventInfo;

        const edm::EDGetTokenT<std::vector<reco::Muon> >              muonsToken;
        const edm::EDGetTokenT<std::vector<reco::GenParticle> >       gensToken;
        const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelClustersToken;
        const edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster> > stripClustersToken;

        TTree* tree;

        unsigned long run;
        unsigned long lumi;
        unsigned long long event;

        std::vector<TVector3>               genmup;
        std::vector<char>                   genmuq;
        std::vector<TVector3>               muonp;
        std::vector<char>                   muonq;
        std::vector<TVector3>               mutkp;
        std::vector<char>                   mutkq;

        std::vector<TVector3>               hits;
        std::vector<std::vector<int> >      muidx;
        std::vector<char>                   hittype;
        std::vector<double>                 hiterrxx, hiterrxy, hiterryy;

        std::vector<std::vector<TVector3> > trackpos;
};

MuonHitsMapper::MuonHitsMapper(const edm::ParameterSet& iConfig):
    addEventInfo(iConfig.existsAs<bool>("addEventInfo") ? iConfig.getParameter<bool>("addEventInfo") : false),
    muonsToken        (consumes<std::vector<reco::Muon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    gensToken         (consumes<std::vector<reco::GenParticle> >      (iConfig.getParameter<edm::InputTag>("gens"))),
    pixelClustersToken(consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
    stripClustersToken(consumes<edmNew::DetSetVector<SiStripCluster> >(iConfig.getParameter<edm::InputTag>("stripClusters")))
{
   usesResource("TFileService");
}


MuonHitsMapper::~MuonHitsMapper() {
}

void MuonHitsMapper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

    Handle<std::vector<reco::Muon> > muonsH;
    iEvent.getByToken(muonsToken, muonsH);
    
    Handle<std::vector<reco::GenParticle> > gensH;
    iEvent.getByToken(gensToken, gensH);
 
    Handle<edmNew::DetSetVector<SiPixelCluster> > pixelClustersH;
    iEvent.getByToken(pixelClustersToken, pixelClustersH);
  
    Handle<edmNew::DetSetVector<SiStripCluster> > stripClustersH;
    iEvent.getByToken(stripClustersToken, stripClustersH);
  
    ESHandle<TrackerGeometry> trackerGeometry;
    iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometry);

    ESHandle<PixelClusterParameterEstimator> pixelCPEH;
    iSetup.get<TkPixelCPERecord>().get("PixelCPEGeneric", pixelCPEH);
    const PixelCPEBase* pix_cpe = dynamic_cast<const PixelCPEBase*>(&(*pixelCPEH));
 
    ESHandle<StripClusterParameterEstimator> stripCPEH;
    iSetup.get<TkStripCPERecord>().get("StripCPEfromTrackAngle", stripCPEH);
    const StripClusterParameterEstimator* str_cpe = &(*stripCPEH);

    ESHandle<TransientTrackBuilder> ttBuilderH;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilderH);
 
    muonp.clear();
    muonq.clear();
    mutkp.clear();
    mutkq.clear();

    hits.clear();
    muidx.clear();
    trackpos.clear();
    hittype.clear();
    hiterrxx.clear();
    hiterrxy.clear();
    hiterryy.clear();

    typedef edm::Ref<edmNew::DetSetVector<SiStripCluster>, SiStripCluster> ClusterStripRef;
    typedef edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> ClusterPixelRef;

    std::vector<std::vector<std::pair<ClusterPixelRef, DetId> > > pix_clusters;
    std::vector<std::vector<std::pair<ClusterStripRef, DetId> > > str_clusters;
    std::vector<reco::MuonRef> murefs;

    for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
        if (muons_iter->innerTrack().isNull()) continue;

        TVector3 muonv3;
        muonv3.SetXYZ(muons_iter->px(), muons_iter->py(), muons_iter->pz());
        muonp.push_back(muonv3);
        muonq.push_back(char(muons_iter->charge()));

        TVector3 mutkv3;
        mutkv3.SetXYZ(muons_iter->innerTrack()->px(), muons_iter->innerTrack()->py(), muons_iter->innerTrack()->pz());
        mutkq.push_back(char(muons_iter->innerTrack()->charge()));

        std::vector<std::pair<ClusterPixelRef, DetId> > pixc;
        std::vector<std::pair<ClusterStripRef, DetId> > strc;
        reco::TrackRef track = muons_iter->innerTrack();
        for (trackingRecHit_iterator tkhits_iter = track->recHitsBegin(); tkhits_iter != track->recHitsEnd(); ++tkhits_iter) {
            if (!(*tkhits_iter)->isValid()) continue;
            
            const BaseTrackerRecHit* tracker_hit = dynamic_cast<const BaseTrackerRecHit*>(*tkhits_iter);
            ClusterPixelRef pix_cluster = tracker_hit->firstClusterRef().cluster_pixel();
            ClusterStripRef str_cluster = tracker_hit->firstClusterRef().cluster_strip();
            if (pix_cluster.isNonnull()) pixc.push_back(std::pair<ClusterPixelRef, DetId>(pix_cluster, tracker_hit->geographicalId()));
            if (str_cluster.isNonnull()) strc.push_back(std::pair<ClusterStripRef, DetId>(str_cluster, tracker_hit->geographicalId()));
        }
        pix_clusters.push_back(pixc);
        str_clusters.push_back(strc);
        reco::MuonRef mref(muonsH, muons_iter - muonsH->begin());        
        murefs.push_back(mref);
    }

    for (auto clusterSet = pixelClustersH->begin(); clusterSet != pixelClustersH->end(); ++clusterSet) {
        DetId detId(clusterSet->detId());
        auto detUnit = dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry->idToDet(detId));
        auto pixelTopology = dynamic_cast<const PixelTopology*>(&(detUnit->specificTopology()));
   
        for(auto cluster = clusterSet->begin(); cluster != clusterSet->end(); cluster++) {
            LocalPoint localPoint = pixelTopology->localPosition(MeasurementPoint(cluster->x(), cluster->y()));
            GlobalPoint globalPoint = detUnit->surface().toGlobal(localPoint);
            TVector3 hit;
            hit.SetXYZ(globalPoint.x(), globalPoint.y(), globalPoint.z());
    
            std::tuple<LocalPoint, LocalError, SiPixelRecHitQuality::QualWordType> cpe_tuple = pix_cpe->getParameters(*cluster, *detUnit);
            LocalError pos_err = std::get<1>(cpe_tuple);
            double hitexx = pos_err.xx();
            double hitexy = pos_err.xy();
            double hiteyy = pos_err.yy();
    
            char hitt = 0;
            GeomDetEnumerators::SubDetector subDet = detUnit->specificType().subDetector();
            if (subDet == GeomDetEnumerators::P1PXB)  hitt = 1;
            if (subDet == GeomDetEnumerators::P1PXEC) hitt = 2;

            std::vector<int> vmuidx; 
            std::vector<TVector3> vtpos; 
            ClusterPixelRef cref = edmNew::makeRefTo(pixelClustersH, cluster);
            for (std::size_t i = 0; i < muonp.size(); i++) {
                bool matchedToMu = false;
                bool add_cluster = false;
                for (std::size_t j = 0; j < pix_clusters[i].size(); j++) {
                    if (cref == pix_clusters[i][j].first) matchedToMu = true;
                }
                for (std::size_t j = 0; j < pix_clusters[i].size(); j++) {
                    auto du = dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry->idToDet(pix_clusters[i][j].second));
                    auto pixTop = dynamic_cast<const PixelTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = pixTop->localPosition(MeasurementPoint(pix_clusters[i][j].first->x(), pix_clusters[i][j].first->y()));
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    if (hv.DeltaR(hit) < 0.1) add_cluster = true;
                }
                for (std::size_t j = 0; j < str_clusters[i].size(); j++) {
                    auto du = dynamic_cast<const StripGeomDetUnit*>(trackerGeometry->idToDet(str_clusters[i][j].second));
                    auto strTop = dynamic_cast<const StripTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = strTop->localPosition(str_clusters[i][j].first->barycenter());
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    if (hv.DeltaR(hit) < 0.1) add_cluster = true;
                }
                if (add_cluster) {
                    int midx = int(i);
                    if (matchedToMu) midx *= -1;
                    vmuidx.push_back(midx);
                    TVector3 tpos;
                    GlobalPoint trackpoint = ttBuilderH->build(&(*(murefs[i]->innerTrack()))).trajectoryStateClosestToPoint(globalPoint).position();
                    tpos.SetXYZ(trackpoint.x(), trackpoint.y(), trackpoint.z());
                }
            }
            if (vmuidx.size() == 0) continue;
    
            hits.push_back(hit);
            muidx.push_back(vmuidx);
            trackpos.push_back(vtpos);
            hittype.push_back(hitt);
            hiterrxx.push_back(hitexx);
            hiterrxy.push_back(hitexy);
            hiterryy.push_back(hiteyy);
        }
    }

    for (auto clusterSet = stripClustersH->begin(); clusterSet != stripClustersH->end(); ++clusterSet) {
        DetId detId(clusterSet->detId());
        auto detUnit = dynamic_cast<const StripGeomDetUnit*>(trackerGeometry->idToDet(detId));
        auto stripTopology = dynamic_cast<const StripTopology*>(&(detUnit->specificTopology()));
    
        for(auto cluster = clusterSet->begin(); cluster != clusterSet->end(); cluster++) {
            LocalPoint localPoint = stripTopology->localPosition(cluster->barycenter());
            GlobalPoint globalPoint = detUnit->surface().toGlobal(localPoint);
            TVector3 hit;
            hit.SetXYZ(globalPoint.x(), globalPoint.y(), globalPoint.z());
    
            std::pair<LocalPoint, LocalError> cpe_pair = str_cpe->localParameters(*cluster, *detUnit);
            LocalError pos_err = cpe_pair.second;
            double hitexx = pos_err.xx();
            double hitexy = pos_err.xy();
            double hiteyy = pos_err.yy();
    
            char hitt = 0;
            GeomDetEnumerators::SubDetector subDet = detUnit->specificType().subDetector();
            if (subDet == GeomDetEnumerators::TIB) hitt = 3;
            if (subDet == GeomDetEnumerators::TOB) hitt = 4;
            if (subDet == GeomDetEnumerators::TEC) hitt = 5;

            std::vector<int> vmuidx;
            std::vector<TVector3> vtpos; 
            ClusterStripRef cref = edmNew::makeRefTo(stripClustersH, cluster);
            for (std::size_t i = 0; i < muonp.size(); i++) {
                bool matchedToMu = false;
                bool add_cluster = false;
                for (std::size_t j = 0; j < str_clusters[i].size(); j++) {
                    if (cref == str_clusters[i][j].first) matchedToMu = true;
                }
                for (std::size_t j = 0; j < pix_clusters[i].size(); j++) {
                    auto du = dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry->idToDet(pix_clusters[i][j].second));
                    auto pixTop = dynamic_cast<const PixelTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = pixTop->localPosition(MeasurementPoint(pix_clusters[i][j].first->x(), pix_clusters[i][j].first->y()));
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    if (hv.DeltaR(hit) < 0.1) add_cluster = true;
                }
                for (std::size_t j = 0; j < str_clusters[i].size(); j++) {
                    auto du = dynamic_cast<const StripGeomDetUnit*>(trackerGeometry->idToDet(str_clusters[i][j].second));
                    auto strTop = dynamic_cast<const StripTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = strTop->localPosition(str_clusters[i][j].first->barycenter());
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    if (hv.DeltaR(hit) < 0.1) add_cluster = true;
                }
                if (add_cluster) {
                    int midx = int(i);
                    if (matchedToMu) midx *= -1;
                    vmuidx.push_back(midx);
                    TVector3 tpos;
                    GlobalPoint trackpoint = ttBuilderH->build(&(*(murefs[i]->innerTrack()))).trajectoryStateClosestToPoint(globalPoint).position();
                    tpos.SetXYZ(trackpoint.x(), trackpoint.y(), trackpoint.z());
                }
            }
            if (vmuidx.size() == 0) continue;
 
            hits.push_back(hit);
            muidx.push_back(vmuidx);
            trackpos.push_back(vtpos);
            hittype.push_back(hitt);
            hiterrxx.push_back(hitexx);
            hiterrxy.push_back(hitexy);
            hiterryy.push_back(hiteyy);
        }
    }

    if (gensH.isValid()) {
        for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if (abs(gens_iter->pdgId()) != 13 || gens_iter->status() != 1) continue;
            TVector3 genv3;
            genv3.SetXYZ(gens_iter->px(), gens_iter->py(), gens_iter->pz());
            genmup.push_back(genv3);
            genmuq.push_back(char(gens_iter->charge()));
        }
    }

    tree->Fill();
}


void MuonHitsMapper::beginJob() {
    // Access the TFileService
    edm::Service<TFileService> fs;

    // Create the TTree
    tree = fs->make<TTree>("tree"       , "tree");

    // Event coordinates
    if (addEventInfo) {
    tree->Branch("event"                , &event                               , "event/i");
    tree->Branch("run"                  , &run                                 , "run/i");
    tree->Branch("lumi"                 , &lumi                                , "lumi/i");
    }

    // Generator-level info
    tree->Branch("genmup"               , "std::vector<TVector3>"              , &genmup    , 32000, 0);
    tree->Branch("genmuq"               , "std::vector<char>"                  , &genmuq    , 32000, 0);

    // Muon info
    tree->Branch("muonp"                , "std::vector<TVector3>"              , &muonp     , 32000, 0);
    tree->Branch("muonq"                , "std::vector<char>"                  , &muonq     , 32000, 0);
    tree->Branch("mutkp"                , "std::vector<TVector3>"              , &mutkp     , 32000, 0);
    tree->Branch("mutkq"                , "std::vector<char>"                  , &mutkq     , 32000, 0);

    // Hits info
    tree->Branch("hits"                 , "std::vector<TVector3>"              , &hits      , 32000, 0);
    tree->Branch("muidx"                , "std::vector<std::vector<int> >"     , &muidx     , 32000, 0);
    tree->Branch("hittype"              , "std::vector<char>"                  , &hittype   , 32000, 0);
    tree->Branch("hiterrxx"             , "std::vector<double>"                , &hiterrxx  , 32000, 0);
    tree->Branch("hiterrxy"             , "std::vector<double>"                , &hiterrxy  , 32000, 0);
    tree->Branch("hiterryy"             , "std::vector<double>"                , &hiterryy  , 32000, 0);

    // Track info
    tree->Branch("trackpos"             , "std::vector<std::vector<TVector3> >", &trackpos  , 32000, 0);
}

void MuonHitsMapper::endJob() {
}

void MuonHitsMapper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MuonHitsMapper);

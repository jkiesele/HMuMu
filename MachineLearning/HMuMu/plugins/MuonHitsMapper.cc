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
#include <string>

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

#include "DataFormats/Math/interface/deltaPhi.h"

template<typename _RandomAccessIterator>
inline std::vector<size_t>
retsort(_RandomAccessIterator __first, _RandomAccessIterator __last){
    typedef typename std::iterator_traits<_RandomAccessIterator>::value_type
            _ValueType;
    std::vector<_ValueType> copy(__first,__last); //copy
    std::vector<size_t> sortedilo;
    std::sort(copy.begin(),copy.end());

    for(_RandomAccessIterator it=copy.begin();it!=copy.end();++it){
        //get the position in input
        size_t pos=std::find(__first,__last,*it)-__first;
        while(std::find(sortedilo.begin(),sortedilo.end(),pos)!=sortedilo.end())
            pos=std::find(__first+pos+1,__last,*it)-__first;
        sortedilo.push_back(pos);
    }
    return sortedilo;
}

template<class T>
std::vector<T> sortByIndex(const std::vector<T> & in, const std::vector<size_t>& sorting){
    std::vector<T> out;
    for(size_t i=0;i<sorting.size();i++)
        out.push_back( in.at(sorting.at(i)) );
    return out;
}



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
        double radiusForNearbyHits;
        std::string pixelCPETag;
        std::string stripCPETag;
        std::string ttkBuildTag;

        const edm::EDGetTokenT<std::vector<reco::Muon> >              muonsToken;
        const edm::EDGetTokenT<std::vector<reco::GenParticle> >       gensToken;
        const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelClustersToken;
        const edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster> > stripClustersToken;

        TTree* tree;

        unsigned long run;
        unsigned long lumi;
        unsigned long long event;

        //TVector3 genmup;
        std::vector<float> genmup;
        char     genmuq;
        //TVector3 muonp;
        std::vector<float> muonp;
        char     muonq;
        //TVector3 mutkp;
        std::vector<float> mutkp;
        char     mutkq;

        //other muon properties
        std::vector<float> muonprop;

        //std::vector<TVector3> hits;
        std::vector<std::vector<float> > hits;
        std::vector<char>     hittype;
        std::vector<float> hitmatched;
        std::vector<float>   hiterrxx, hiterrxy, hiterryy;
        //std::vector<TVector3> trackpos;
        std::vector<std::vector<float> > trackpos;
        std::vector<float> hitDRdisttotrack,hiteucldisttotrack;

        // muon chambers segments and muon track branches
        std::vector<float> trackmuposx, trackmuposy,  trackmuposxerr, trackmuposyerr;
        std::vector<unsigned int> trackmupostation;
        std::vector<float>  segx, segy, segxerr, segyerr, segmudr, segmudrerr ;
        //utils

        int muoncounter;


};

MuonHitsMapper::MuonHitsMapper(const edm::ParameterSet& iConfig):
    addEventInfo       (iConfig.existsAs<bool>       ("addEventInfo")    ? iConfig.getParameter<bool>       ("addEventInfo")    : false),
    radiusForNearbyHits(iConfig.existsAs<double>     ("hitSearchRadius") ? iConfig.getParameter<double>     ("hitSearchRadius") : 1.0),
    pixelCPETag        (iConfig.existsAs<std::string>("pixelCPE")        ? iConfig.getParameter<std::string>("pixelCPE")        : "PixelCPEGeneric"),
    stripCPETag        (iConfig.existsAs<std::string>("stripCPE")        ? iConfig.getParameter<std::string>("stripCPE")        : "StripCPEfromTrackAngle"),
    ttkBuildTag        (iConfig.existsAs<std::string>("ttkBuilder")      ? iConfig.getParameter<std::string>("ttkBuilder")      : "TransientTrackBuilder"),
    muonsToken        (consumes<std::vector<reco::Muon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    gensToken         (consumes<std::vector<reco::GenParticle> >      (iConfig.getParameter<edm::InputTag>("gens"))),
    pixelClustersToken(consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
    stripClustersToken(consumes<edmNew::DetSetVector<SiStripCluster> >(iConfig.getParameter<edm::InputTag>("stripClusters")))
{
   usesResource("TFileService");
   muoncounter=0;
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
    iSetup.get<TkPixelCPERecord>().get(pixelCPETag, pixelCPEH);
    const PixelCPEBase* pix_cpe = dynamic_cast<const PixelCPEBase*>(&(*pixelCPEH));
 
    ESHandle<StripClusterParameterEstimator> stripCPEH;
    iSetup.get<TkStripCPERecord>().get(stripCPETag, stripCPEH);
    const StripClusterParameterEstimator* str_cpe = &(*stripCPEH);

    ESHandle<TransientTrackBuilder> ttBuilderH;
    iSetup.get<TransientTrackRecord>().get(ttkBuildTag, ttBuilderH);
 
    typedef edm::Ref<edmNew::DetSetVector<SiStripCluster>, SiStripCluster> ClusterStripRef;
    typedef edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> ClusterPixelRef;



    for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {



        genmup.clear();
        muonp.clear();
        mutkp.clear();
        muonprop.clear();

        TVector3 gmu;
        TVector3 rmu;

        if (muons_iter->innerTrack().isNull()) continue;
        //muonp.SetXYZ(muons_iter->px(), muons_iter->py(), muons_iter->pz());
        rmu.SetXYZ(muons_iter->px(), muons_iter->py(), muons_iter->pz());

        double gen_dR = 1e10;
        if (gensH.isValid()) {
            auto use_gens_iter = gensH->end();
            for (auto gens_iter = gensH->begin(); gens_iter <= gensH->end(); ++gens_iter) {
                if(gens_iter == gensH->end()) break;
                if (abs(gens_iter->pdgId()) != 13 || gens_iter->status() != 1) continue;
                //genmup.SetXYZ(gens_iter->px(), gens_iter->py(), gens_iter->pz());
                gmu.SetXYZ(gens_iter->px(), gens_iter->py(), gens_iter->pz());
                if (gmu.DeltaR(rmu) > gen_dR) continue;
                gen_dR = gmu.DeltaR(rmu);
                use_gens_iter=gens_iter;
            }
            if(use_gens_iter != gensH->end()){
                genmup.push_back(use_gens_iter->px());
                genmup.push_back(use_gens_iter->py());
                genmup.push_back(use_gens_iter->pz());
                genmup.push_back(use_gens_iter->eta());
                genmup.push_back(use_gens_iter->phi());
                genmuq = char(use_gens_iter->charge());
            }
        }
        if (gen_dR > 0.1) continue;


        muoncounter++;

        //if((muoncounter<17925-2000))continue;

        muonp.push_back(muons_iter->px());
        muonp.push_back(muons_iter->py());
        muonp.push_back(muons_iter->pz());
        muonp.push_back(muons_iter->eta());
        muonp.push_back(muons_iter->phi());
        muonq = char(muons_iter->charge());

        muonprop.push_back(muons_iter->dxyError());
        muonprop.push_back(muons_iter->innerTrack()->chi2());
        muonprop.push_back(muons_iter->innerTrack()->etaError());//indicative of rel pzerr
        muonprop.push_back(muons_iter->innerTrack()->normalizedChi2());
        muonprop.push_back(muons_iter->innerTrack()->outerPx());
        muonprop.push_back(muons_iter->innerTrack()->outerPy());
        muonprop.push_back(muons_iter->innerTrack()->outerPz());

        segmudr.clear();
        segmudrerr.clear();
        trackmuposx.clear();
        trackmuposy.clear();
        trackmuposxerr.clear();
        trackmuposyerr.clear();
        trackmupostation.clear();
        segx.clear();
        segy.clear();
        segxerr.clear();
        segyerr.clear();

        for (const auto &ch : muons_iter->matches()) {  // loop over matched chambers
            for(std::vector<reco::MuonSegmentMatch>::const_iterator matseg = ch.segmentMatches.begin(); matseg != ch.segmentMatches.end(); matseg++) {
                segx.push_back(matseg->x);
                segy.push_back(matseg->y);
                segxerr.push_back(matseg->xErr);
                segyerr.push_back(matseg->yErr);
                float dr_ = sqrt( pow(matseg->x - ch.x,2) + pow(matseg->y - ch.y,2 ) );
                segmudr.push_back(sqrt(matseg->x*matseg->x + matseg->y*matseg->y));
                segmudrerr.push_back(dr_ * sqrt( (matseg->xErr/matseg->x)*(matseg->xErr/matseg->x) + (matseg->yErr/matseg->y)*(matseg->yErr/matseg->y) ) );
                trackmuposx.push_back(ch.x);
                trackmuposy.push_back(ch.y);
                trackmuposxerr.push_back(ch.xErr);
                trackmuposyerr.push_back(ch.yErr);
                trackmupostation.push_back(ch.station());
            }
            if(ch.segmentMatches.size()<1){
                segx.push_back(0);
                segy.push_back(0);
                segxerr.push_back(0);
                segyerr.push_back(0);
                segmudr.push_back(0);
                segmudrerr.push_back(0);
                trackmuposx.push_back(ch.x);
                trackmuposy.push_back(ch.y);
                trackmuposxerr.push_back(ch.xErr);
                trackmuposyerr.push_back(ch.yErr);
                trackmupostation.push_back(ch.station());
            }

        }




        //mutkp.SetXYZ(muons_iter->innerTrack()->px(), muons_iter->innerTrack()->py(), muons_iter->innerTrack()->pz());
        mutkp.push_back(muons_iter->innerTrack()->p());
        mutkp.push_back(muons_iter->innerTrack()->eta());
        mutkp.push_back(muons_iter->innerTrack()->phi());
        mutkq = char(muons_iter->innerTrack()->charge());

        hits.clear();
        hittype.clear();
        hitmatched.clear();
        hiterrxx.clear();
        hiterrxy.clear();
        hiterryy.clear();
        trackpos.clear();
        hitDRdisttotrack.clear();
        hiteucldisttotrack.clear();

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
        
                ClusterPixelRef cref = edmNew::makeRefTo(pixelClustersH, cluster);
                bool matchedToMu = false;
                bool add_cluster = false;
                for (std::size_t j = 0; j < pixc.size(); j++) {
                    if (cref == pixc[j].first) matchedToMu = true;
                }
                for (std::size_t j = 0; j < pixc.size(); j++) {
                    auto du = dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry->idToDet(pixc[j].second));
                    auto pixTop = dynamic_cast<const PixelTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = pixTop->localPosition(MeasurementPoint(pixc[j].first->x(), pixc[j].first->y()));
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    double rad = (hv-hit).Mag();
                    if (rad < radiusForNearbyHits) add_cluster = true;
                }
                for (std::size_t j = 0; j < strc.size(); j++) {
                    auto du = dynamic_cast<const StripGeomDetUnit*>(trackerGeometry->idToDet(strc[j].second));
                    auto strTop = dynamic_cast<const StripTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = strTop->localPosition(strc[j].first->barycenter());
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    double rad = (hv-hit).Mag();
                    if (rad < radiusForNearbyHits) add_cluster = true;
                }
                if (add_cluster) {
                    if (matchedToMu) hitt *= -1;
                    TVector3 tpos;
                    GlobalPoint trackpoint = ttBuilderH->build(&(*(muons_iter->innerTrack()))).trajectoryStateClosestToPoint(globalPoint).position();
                    tpos.SetXYZ(trackpoint.x(), trackpoint.y(), trackpoint.z());

                    std::vector<float> vhit;
                    vhit.push_back(hit.Mag());
                    vhit.push_back(hit.Eta());
                    vhit.push_back(reco::deltaPhi(hit.Phi(), muons_iter->phi()));
                    vhit.push_back(sqrt(hit.x()*hit.x() + hit.y()*hit.y()));
                    vhit.push_back(hit.x());
                    vhit.push_back(hit.y());
                    vhit.push_back(hit.z());

                    std::vector<float> vtpos;
                    vtpos.push_back(tpos.Mag());
                    vtpos.push_back(tpos.Eta());
                    vtpos.push_back(reco::deltaPhi(tpos.Phi(), muons_iter->phi()));
                    vtpos.push_back(sqrt(tpos.x()*tpos.x() + tpos.y()*tpos.y()));
                    vtpos.push_back(tpos.x());
                    vtpos.push_back(tpos.y());
                    vtpos.push_back(tpos.z());

                    auto distance = tpos.DeltaR(hit);
                    hitDRdisttotrack.push_back( distance);
                    hiteucldisttotrack.push_back((tpos - hit).Mag());

                    //hits.push_back(hit);
                    hits.push_back(vhit);
                    hittype.push_back(hitt);
                    if (matchedToMu)
                        hitmatched.push_back(1);
                    else
                        hitmatched.push_back(0);
                    hiterrxx.push_back(hitexx);
                    hiterrxy.push_back(hitexy);
                    hiterryy.push_back(hiteyy);
                    //trackpos.push_back(tpos);
                    trackpos.push_back(vtpos);
                }
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
        
                ClusterStripRef cref = edmNew::makeRefTo(stripClustersH, cluster);
                bool matchedToMu = false;
                bool add_cluster = false;
                for (std::size_t j = 0; j < strc.size(); j++) {
                    if (cref == strc[j].first) matchedToMu = true;
                }
                for (std::size_t j = 0; j < pixc.size(); j++) {
                    auto du = dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry->idToDet(pixc[j].second));
                    auto pixTop = dynamic_cast<const PixelTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = pixTop->localPosition(MeasurementPoint(pixc[j].first->x(), pixc[j].first->y()));
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    double rad = (hv-hit).Mag();
                    if (rad < radiusForNearbyHits) add_cluster = true;
                }
                for (std::size_t j = 0; j < strc.size(); j++) {
                    auto du = dynamic_cast<const StripGeomDetUnit*>(trackerGeometry->idToDet(strc[j].second));
                    auto strTop = dynamic_cast<const StripTopology*>(&(du->specificTopology()));
                    LocalPoint  lp = strTop->localPosition(strc[j].first->barycenter());
                    GlobalPoint gp = du->surface().toGlobal(lp);
                    TVector3 hv;
                    hv.SetXYZ(gp.x(), gp.y(), gp.z());
                    double rad = (hv-hit).Mag();
                    if (rad < radiusForNearbyHits) add_cluster = true;
                }
                if (add_cluster) {
                    if (matchedToMu) hitt *= -1;
                    TVector3 tpos;
                    GlobalPoint trackpoint = ttBuilderH->build(&(*(muons_iter->innerTrack()))).trajectoryStateClosestToPoint(globalPoint).position();
                    tpos.SetXYZ(trackpoint.x(), trackpoint.y(), trackpoint.z());

                    std::vector<float> vhit;
                    vhit.push_back(hit.Mag());
                    vhit.push_back(hit.Eta());
                    vhit.push_back(reco::deltaPhi(hit.Phi(), muons_iter->phi()));
                    vhit.push_back(sqrt(hit.x()*hit.x() + hit.y()*hit.y()));
                    vhit.push_back(hit.x());
                    vhit.push_back(hit.y());
                    vhit.push_back(hit.z());

                    std::vector<float> vtpos;
                    vtpos.push_back(tpos.Mag());
                    vtpos.push_back(tpos.Eta());
                    vtpos.push_back(reco::deltaPhi(tpos.Phi(), muons_iter->phi()));
                    vtpos.push_back(sqrt(tpos.x()*tpos.x() + tpos.y()*tpos.y()));
                    vtpos.push_back(tpos.x());
                    vtpos.push_back(tpos.y());
                    vtpos.push_back(tpos.z());


                    auto distance = tpos.DeltaR(hit);
                    hitDRdisttotrack.push_back( distance);
                    hiteucldisttotrack.push_back((tpos - hit).Mag());

                    //hits.push_back(hit);
                    hits.push_back(vhit);
                    hittype.push_back(hitt);
                    if (matchedToMu)
                        hitmatched.push_back(1);
                    else
                        hitmatched.push_back(0);
                    hiterrxx.push_back(hitexx);
                    hiterrxy.push_back(hitexy);
                    hiterryy.push_back(hiteyy);
                    //trackpos.push_back(tpos);
                    trackpos.push_back(vtpos);
                }
            }
        }
        //order all tracker hits by distance to track for truncation later
        //way faster here than in python later

        auto sortindices  =  retsort(hitDRdisttotrack.begin(),hitDRdisttotrack.end());
        hits = sortByIndex(hits, sortindices);
        hittype = sortByIndex(hittype, sortindices);
        hitmatched = sortByIndex(hitmatched, sortindices);
        hiterrxx = sortByIndex(hiterrxx, sortindices);
        hiterrxy = sortByIndex(hiterrxy, sortindices);
        hiterryy = sortByIndex(hiterryy, sortindices);
        trackpos = sortByIndex(trackpos, sortindices);

        hitDRdisttotrack = sortByIndex(hitDRdisttotrack, sortindices);
        hiteucldisttotrack = sortByIndex(hiteucldisttotrack, sortindices);

        tree->Fill();
    }
}


void MuonHitsMapper::beginJob() {
    // Access the TFileService
    edm::Service<TFileService> fs;

    // Create the TTree
    tree = fs->make<TTree>("tree"       , "tree");

    // Event coordinates
    if (addEventInfo) {
    tree->Branch("event"                , &event                              , "event/i" );
    tree->Branch("run"                  , &run                                , "run/i"   );
    tree->Branch("lumi"                 , &lumi                               , "lumi/i"  );
    }

    // Muon info
    tree->Branch("genmuq"               , &genmuq                             , "genmuq/b");
    tree->Branch("muonq"                , &muonq                              , "muonq/b" );
    tree->Branch("mutkq"                , &mutkq                              , "mutkq/b" );

    //tree->Branch("genmup"               , "TVector3"                          , &genmup    , 32000, 0);
    //tree->Branch("muonp"                , "TVector3"                          , &muonp     , 32000, 0);
    //tree->Branch("mutkp"                , "TVector3"                          , &mutkp     , 32000, 0);

    tree->Branch("genmup"               , "std::vector<float>"               , &genmup    , 32000, 0);
    tree->Branch("muonp"                , "std::vector<float>"               , &muonp     , 32000, 0);
    tree->Branch("mutkp"                , "std::vector<float>"               , &mutkp     , 32000, 0);
    tree->Branch("muonprop"                , "std::vector<float>"               , &muonprop     , 32000, 0);



    // Hits info
    //tree->Branch("hits"                 , "std::vector<TVector3>"             , &hits      , 32000, 0);
    tree->Branch("hits"                 , "std::vector<std::vector<float> >" , &hits      , 32000, 0);
    tree->Branch("hittype"              , "std::vector<char>"                 , &hittype   , 32000, 0);
    tree->Branch("hitmatched"              , "std::vector<float>"                 , &hitmatched   , 32000, 0);
    tree->Branch("hiterrxx"             , "std::vector<float>"               , &hiterrxx  , 32000, 0);
    tree->Branch("hiterrxy"             , "std::vector<float>"               , &hiterrxy  , 32000, 0);
    tree->Branch("hiterryy"             , "std::vector<float>"               , &hiterryy  , 32000, 0);
    //tree->Branch("trackpos"             , "std::vector<TVector3>" , &trackpos  , 32000, 0);
    tree->Branch("trackpos"             , "std::vector<std::vector<float> >" , &trackpos  , 32000, 0);
    tree->Branch("hitDRdisttotrack"              , "std::vector<float>"                 , &hitDRdisttotrack   , 32000, 0);
    tree->Branch("hiteucldisttotrack"              , "std::vector<float>"                 , &hiteucldisttotrack   , 32000, 0);



    // muon xy distance from segment 
    tree->Branch("segmudr"             , "std::vector<float> "  , &segmudr     , 32000, 0);
    tree->Branch("segmudrerr"          , "std::vector<float> "  , &segmudrerr  , 32000, 0);
 
    // muon segments info
    tree->Branch("segx"                , "std::vector<float> "  , &segx        , 32000, 0);
    tree->Branch("segy"                , "std::vector<float> "  , &segy        , 32000, 0);
    tree->Branch("segxerr"             , "std::vector<float> "  , &segxerr     , 32000, 0);
    tree->Branch("segyerr"             , "std::vector<float> "  , &segyerr     , 32000, 0);
 
    // muon track position in muon chambers
    tree->Branch("trackmuposx"         , "std::vector<float>"                , &trackmuposx     , 32000, 0);
    tree->Branch("trackmuposy"         , "std::vector<float>"                , &trackmuposy     , 32000, 0);
    tree->Branch("trackmuposxerr"      , "std::vector<float>"                , &trackmuposxerr     , 32000, 0);
    tree->Branch("trackmuposyerr"      , "std::vector<float>"                , &trackmuposyerr     , 32000, 0);
    tree->Branch("trackmupostation"    , "std::vector<unsigned int>"         , &trackmupostation   , 32000, 0);
}

void MuonHitsMapper::endJob() {
}

void MuonHitsMapper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MuonHitsMapper);

// -*- C++ -*-
//
// Package:    Validation/MuonIdentification
// Class:      MuonIdVal
//
/*

 Description:  Makes and fills lots of histograms using the various reco::Muon
               methods.


*/
//
// Original Author:  Jacob Ribnik
//         Created:  Wed Apr 18 13:48:08 CDT 2007
//
//

#ifndef Validation_MuonIdentification_MuonIdVal_h
#define Validation_MuonIdentification_MuonIdVal_h

// system include files
#include <string>

// user include files
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

// class MuonIdVal : public edm::EDAnalyzer {
class MuonIdVal : public DQMEDAnalyzer {
public:
  explicit MuonIdVal(const edm::ParameterSet &);
  ~MuonIdVal() override;

private:
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void Fill(MonitorElement *, float);

  template <class TSegmentCollectionConstIterator>
  bool matchSegment(TSegmentCollectionConstIterator,
                    const reco::MuonSegmentMatch &);


  template <class TSegmentCollectionHandle>
  void fillSegmentPosition2D(const TSegmentCollectionHandle &,
                       const edm::Handle<reco::MuonCollection> &,
                       const edm::ESHandle<GlobalTrackingGeometry> &);

  edm::ParameterSet iConfig;
  edm::ParameterSet parameters_;
  std::string eventInfoFolder_;
  std::string subsystemname_;

  // ----------member data ---------------------------
  edm::InputTag inputMuonCollection_;
  edm::InputTag inputDTRecSegment4DCollection_;
  edm::InputTag inputCSCSegmentCollection_;
  edm::InputTag inputGEMSegmentCollection_;
  edm::InputTag inputMuonTimeExtraValueMap_;
  edm::InputTag inputMuonCosmicCompatibilityValueMap_;
  edm::InputTag inputMuonShowerInformationValueMap_;
  edm::EDGetTokenT<reco::MuonCollection> inputMuonCollectionToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> inputDTRecSegment4DCollectionToken_;
  edm::EDGetTokenT<CSCSegmentCollection> inputCSCSegmentCollectionToken_;
  edm::EDGetTokenT<GEMSegmentCollection> inputGEMSegmentCollectionToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> inputMuonTimeExtraValueMapCombToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> inputMuonTimeExtraValueMapDTToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> inputMuonTimeExtraValueMapCSCToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonCosmicCompatibility>> inputMuonCosmicCompatibilityValueMapToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonShower>> inputMuonShowerInformationValueMapToken_;
  bool useTrackerMuons_;
  bool useGlobalMuons_;
  bool useTrackerMuonsNotGlobalMuons_;
  bool useGlobalMuonsNotTrackerMuons_;
  bool useGEM_;
  bool makeEnergyPlots_;
  bool makeTimePlots_;
  bool make2DPlots_;
  bool makeAllChamberPlots_;
  bool makeCosmicCompatibilityPlots_;
  bool makeShowerInformationPlots_;
  std::string baseFolder_;

  edm::Handle<reco::MuonCollection> muonCollectionH_;
  edm::Handle<DTRecSegment4DCollection> dtSegmentCollectionH_;
  edm::Handle<CSCSegmentCollection> cscSegmentCollectionH_;
  edm::Handle<GEMSegmentCollection> gemSegmentCollectionH_;
  edm::Handle<reco::MuonTimeExtraMap> combinedMuonTimeExtraValueMapH_;
  edm::Handle<reco::MuonTimeExtraMap> cscMuonTimeExtraValueMapH_;
  edm::Handle<reco::MuonTimeExtraMap> dtMuonTimeExtraValueMapH_;
  edm::Handle<edm::ValueMap<reco::MuonCosmicCompatibility>> muonCosmicCompatibilityValueMapH_;
  edm::Handle<edm::ValueMap<reco::MuonShower>> muonShowerInformationValueMapH_;
  edm::ESHandle<GlobalTrackingGeometry> geometry_;

  // trackerMuon == 0; globalMuon == 1
  // energy deposits
  MonitorElement *hEnergyEMBarrel[4];
  MonitorElement *hEnergyHABarrel[4];
  MonitorElement *hEnergyHO[4];
  MonitorElement *hEnergyEMEndcap[4];
  MonitorElement *hEnergyHAEndcap[4];

  // time information
  MonitorElement *hMuonTimeNDOF[4];
  MonitorElement *hMuonTimeTimeAtIpInOut[4];
  MonitorElement *hMuonTimeTimeAtIpInOutErr[4];
  MonitorElement *hMuonTimeTimeAtIpOutIn[4];
  MonitorElement *hMuonTimeTimeAtIpOutInErr[4];
  MonitorElement *hMuonTimeExtraCombinedNDOF[4];
  MonitorElement *hMuonTimeExtraCombinedTimeAtIpInOut[4];
  MonitorElement *hMuonTimeExtraCombinedTimeAtIpInOutErr[4];
  MonitorElement *hMuonTimeExtraCombinedTimeAtIpOutIn[4];
  MonitorElement *hMuonTimeExtraCombinedTimeAtIpOutInErr[4];
  MonitorElement *hMuonTimeExtraCSCNDOF[4];
  MonitorElement *hMuonTimeExtraCSCTimeAtIpInOut[4];
  MonitorElement *hMuonTimeExtraCSCTimeAtIpInOutErr[4];
  MonitorElement *hMuonTimeExtraCSCTimeAtIpOutIn[4];
  MonitorElement *hMuonTimeExtraCSCTimeAtIpOutInErr[4];
  MonitorElement *hMuonTimeExtraDTNDOF[4];
  MonitorElement *hMuonTimeExtraDTTimeAtIpInOut[4];
  MonitorElement *hMuonTimeExtraDTTimeAtIpInOutErr[4];
  MonitorElement *hMuonTimeExtraDTTimeAtIpOutIn[4];
  MonitorElement *hMuonTimeExtraDTTimeAtIpOutInErr[4];

  // muonid
  MonitorElement *hCaloCompat[4];
  MonitorElement *hSegmentCompat[4];
  MonitorElement *hCaloSegmentCompat[4];
  MonitorElement *hMuonQualityTrkRelChi2[4];
  MonitorElement *hMuonQualityStaRelChi2[4];
  MonitorElement *hMuonQualityTrkKink[4];
  MonitorElement *hGlobalMuonPromptTightBool[4];
  MonitorElement *hTMLastStationLooseBool[4];
  MonitorElement *hTMLastStationTightBool[4];
  MonitorElement *hTM2DCompatibilityLooseBool[4];
  MonitorElement *hTM2DCompatibilityTightBool[4];
  MonitorElement *hTMOneStationLooseBool[4];
  MonitorElement *hTMOneStationTightBool[4];
  MonitorElement *hTMLastStationOptimizedLowPtLooseBool[4];
  MonitorElement *hTMLastStationOptimizedLowPtTightBool[4];
  MonitorElement *hGMTkChiCompatibilityBool[4];
  MonitorElement *hGMStaChiCompatibilityBool[4];
  MonitorElement *hGMTkKinkTightBool[4];
  MonitorElement *hTMLastStationAngLooseBool[4];
  MonitorElement *hTMLastStationAngTightBool[4];
  MonitorElement *hTMOneStationAngLooseBool[4];
  MonitorElement *hTMOneStationAngTightBool[4];
  MonitorElement *hTMLastStationOptimizedBarrelLowPtLooseBool[4];
  MonitorElement *hTMLastStationOptimizedBarrelLowPtTightBool[4];

  // cosmic compatibilities
  MonitorElement *hCombinedCosmicCompat[4];
  MonitorElement *hTimeCosmicCompat[4];
  MonitorElement *hB2BCosmicCompat[4];
  MonitorElement *hOverlapCosmicCompat[4];

  // by station

  // shower information
  MonitorElement *hMuonShowerSizeT[4][4];
  MonitorElement *hMuonShowerDeltaR[4][4];
  MonitorElement *hMuonAllHits[4][4];
  MonitorElement *hMuonHitsFromSegments[4][4];
  MonitorElement *hMuonUncorrelatedHits[4][4];

  MonitorElement *hDTPullxPropErr[4][4];
  MonitorElement *hDTPulldXdZPropErr[4][4];
  MonitorElement *hDTPullyPropErr[4][3];
  MonitorElement *hDTPulldYdZPropErr[4][3];
  MonitorElement *hDTDistWithSegment[4][4];
  MonitorElement *hDTDistWithNoSegment[4][4];
  MonitorElement *hDTPullDistWithSegment[4][4];
  MonitorElement *hDTPullDistWithNoSegment[4][4];

  MonitorElement *hCSCPullxPropErr[4][4];
  MonitorElement *hCSCPulldXdZPropErr[4][4];
  MonitorElement *hCSCPullyPropErr[4][4];
  MonitorElement *hCSCPulldYdZPropErr[4][4];
  MonitorElement *hCSCDistWithSegment[4][4];
  MonitorElement *hCSCDistWithNoSegment[4][4];
  MonitorElement *hCSCPullDistWithSegment[4][4];
  MonitorElement *hCSCPullDistWithNoSegment[4][4];

  // key: (muon_type_idx, station)
  std::map<std::tuple<int, int>, MonitorElement *> hGEMPullxPropErr;
  std::map<std::tuple<int, int>, MonitorElement *> hGEMPulldXdZPropErr;
  std::map<std::tuple<int, int>, MonitorElement *> hGEMPullyPropErr;
  std::map<std::tuple<int, int>, MonitorElement *> hGEMPulldYdZPropErr;
  std::map<std::tuple<int, int>, MonitorElement *> hGEMDistWithSegment;
  std::map<std::tuple<int, int>, MonitorElement *> hGEMDistWithNoSegment;
  std::map<std::tuple<int, int>, MonitorElement *> hGEMPullDistWithSegment;
  std::map<std::tuple<int, int>, MonitorElement *> hGEMPullDistWithNoSegment;

  // by chamber, trackerMuons only
  // DT:  [station][wheel][sector]
  // CSC: [endcap][station][ring][chamber]
  // GEM: [region][station][chamber]
  MonitorElement *hDTChamberDx[4][5][14];
  MonitorElement *hDTChamberDy[3][5][14];
  MonitorElement *hDTChamberEdgeXWithSegment[4][5][14];
  MonitorElement *hDTChamberEdgeXWithNoSegment[4][5][14];
  MonitorElement *hDTChamberEdgeYWithSegment[4][5][14];
  MonitorElement *hDTChamberEdgeYWithNoSegment[4][5][14];

  MonitorElement *hCSCChamberDx[2][4][4][36];
  MonitorElement *hCSCChamberDy[2][4][4][36];
  MonitorElement *hCSCChamberEdgeXWithSegment[2][4][4][36];
  MonitorElement *hCSCChamberEdgeXWithNoSegment[2][4][4][36];
  MonitorElement *hCSCChamberEdgeYWithSegment[2][4][4][36];
  MonitorElement *hCSCChamberEdgeYWithNoSegment[2][4][4][36];

  // key is (region, station, chamber)
  std::map<std::tuple<int, int, int>, MonitorElement *> hGEMChamberDx;
  std::map<std::tuple<int, int, int>, MonitorElement *> hGEMChamberDy;
  std::map<std::tuple<int, int, int>, MonitorElement *> hGEMChamberEdgeXWithSegment;
  std::map<std::tuple<int, int, int>, MonitorElement *> hGEMChamberEdgeXWithNoSegment;
  std::map<std::tuple<int, int, int>, MonitorElement *> hGEMChamberEdgeYWithSegment;
  std::map<std::tuple<int, int, int>, MonitorElement *> hGEMChamberEdgeYWithNoSegment;

  // segment matching "efficiency"
  MonitorElement *hSegmentIsAssociatedRZ;
  MonitorElement *hSegmentIsAssociatedXY;
  MonitorElement *hSegmentIsNotAssociatedRZ;
  MonitorElement *hSegmentIsNotAssociatedXY;
  MonitorElement *hSegmentIsBestDrAssociatedRZ;
  MonitorElement *hSegmentIsBestDrAssociatedXY;
  MonitorElement *hSegmentIsBestDrNotAssociatedRZ;
  MonitorElement *hSegmentIsBestDrNotAssociatedXY;
};


template <class TSegmentCollectionConstIterato>
bool MuonIdVal::matchSegment(TSegmentCollectionConstIterato segment,
                             const reco::MuonSegmentMatch & segmentMatch) {
  LocalPoint segmentLocalPosition = segment->localPosition();
  LocalVector segmentLocalDirection = segment->localDirection();
  LocalError segmentLocalPositionError = segment->localPositionError();
  LocalError segmentLocalDirectionError = segment->localDirectionError();

  if (fabs(segmentMatch.x - segmentLocalPosition.x()) > 1E-6) return false;
  if (fabs(segmentMatch.y - segmentLocalPosition.y()) > 1E-6) return false;
  if (fabs(segmentMatch.dXdZ - segmentLocalDirection.x() / segmentLocalDirection.z()) > 1E-6) return false;
  if (fabs(segmentMatch.dYdZ - segmentLocalDirection.y() / segmentLocalDirection.z()) > 1E-6) return false;
  if (fabs(segmentMatch.xErr - sqrt(segmentLocalPositionError.xx())) > 1E-6) return false;
  if (fabs(segmentMatch.yErr - sqrt(segmentLocalPositionError.yy())) > 1E-6) return false;
  if (fabs(segmentMatch.dXdZErr - sqrt(segmentLocalDirectionError.xx())) > 1E-6) return false;
  if (fabs(segmentMatch.dYdZErr - sqrt(segmentLocalDirectionError.yy())) > 1E-6) return false;

  return true;
}


template <class TSegmentCollectionHandle>
void MuonIdVal::fillSegmentPosition2D(const TSegmentCollectionHandle & segment_collection,
                                      const edm::Handle<reco::MuonCollection> & muon_collection,
                                      const edm::ESHandle<GlobalTrackingGeometry> & geometry) {
                                
  for (auto segment = segment_collection->begin(); segment != segment_collection->end(); ++segment) {
    const GeomDet *geom_det = geometry->idToDet(segment->geographicalId());
    GlobalPoint segmentGlobalPosition = geom_det->toGlobal(segment->localPosition());

    bool segmentFound = false;
    bool segmentBestDrFound = false;
    for (auto muon = muon_collection->begin(); muon != muon_collection->end(); ++muon) {
      if (not muon->isMatchesValid())
        continue;

      for (const auto & chamberMatch : muon->matches()) {
        for (const auto & segmentMatch : chamberMatch.segmentMatches) {
          if (matchSegment(segment, segmentMatch)) {
            segmentFound = true;
            if (segmentMatch.isMask(reco::MuonSegmentMatch::BestInStationByDR))
              segmentBestDrFound = true;
            break;
          }
        }  // segmentMatch

        if (segmentFound)
          break;
      }  // chamberMatch

      if (segmentFound)
        break;
    }  // muon

    if (segmentFound) {
      hSegmentIsAssociatedRZ->Fill(segmentGlobalPosition.z(), segmentGlobalPosition.perp());
      hSegmentIsAssociatedXY->Fill(segmentGlobalPosition.x(), segmentGlobalPosition.y());

      if (segmentBestDrFound) {
        hSegmentIsBestDrAssociatedRZ->Fill(segmentGlobalPosition.z(), segmentGlobalPosition.perp());
        hSegmentIsBestDrAssociatedXY->Fill(segmentGlobalPosition.x(), segmentGlobalPosition.y());
      }
    } else {
      hSegmentIsNotAssociatedRZ->Fill(segmentGlobalPosition.z(), segmentGlobalPosition.perp());
      hSegmentIsNotAssociatedXY->Fill(segmentGlobalPosition.x(), segmentGlobalPosition.y());
      hSegmentIsBestDrNotAssociatedRZ->Fill(segmentGlobalPosition.z(), segmentGlobalPosition.perp());
      hSegmentIsBestDrNotAssociatedXY->Fill(segmentGlobalPosition.x(), segmentGlobalPosition.y());
    }
  }  // gem segment

}


#endif

#ifndef DQMOffline_Muon_GEMEfficiencyAnalyzer_h
#define DQMOffline_Muon_GEMEfficiencyAnalyzer_h

#include "DQMOffline/Muon/interface/GEMOfflineDQMBase.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


class GEMEfficiencyAnalyzer : public GEMOfflineDQMBase {
public:
  explicit GEMEfficiencyAnalyzer(const edm::ParameterSet &);
  ~GEMEfficiencyAnalyzer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

protected:
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event &event, const edm::EventSetup &eventSetup) override;

private:
  struct GEMLayerData {
    GEMLayerData(Disk::DiskPointer surface,
                 std::vector<const GEMChamber*> chambers,
                 int region, int station, int layer)
        : surface(surface), chambers(chambers), region(region),
          station(station), layer(layer) {}

    Disk::DiskPointer surface;
    std::vector<const GEMChamber*> chambers;
    int region, station, layer;
  };

  void bookDetectorOccupancy(
      DQMStore::IBooker &, const GEMStation *, const GEMDetId&, const TString &, const TString &);
  void bookOccupancy(DQMStore::IBooker &, const GEMDetId&, const TString &, const TString &);
  void bookChamberOccupancy(DQMStore::IBooker &, const int, const GEMDetId&, const TString &, const TString &);
  void bookResolution(DQMStore::IBooker &, const GEMDetId&, const TString &, const TString &);

  inline bool isInsideOut(const reco::Track&);

  std::vector<GEMLayerData> buildGEMLayers(const edm::ESHandle<GEMGeometry>&);
  const reco::Track* getTrack(const reco::Muon&);
  std::pair<TrajectoryStateOnSurface, DetId> getStartingState(const reco::TransientTrack&, const GEMLayerData&, const edm::ESHandle<GlobalTrackingGeometry>&);
  std::pair<TrajectoryStateOnSurface, DetId> findStartingState(const reco::TransientTrack&, const GEMLayerData&, const edm::ESHandle<GlobalTrackingGeometry>&);
  bool skipLayer(const reco::Track*, const GEMLayerData&);
  bool checkBounds(const GlobalPoint&, const Plane&);
  const GEMEtaPartition* findEtaPartition(const GlobalPoint&, const std::vector<const GEMChamber*>&);
  const GEMRecHit *findMatchedHit(const GlobalPoint&, const GEMRecHitCollection::range &, const GEMEtaPartition*);

  edm::EDGetTokenT<GEMRecHitCollection> rechit_token_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muon_token_;

  MuonServiceProxy *muon_service_;

  bool is_cosmics_;
  bool use_global_muon_;
  float rdphi_cut_;
  std::vector<double> pt_binning_;
  int eta_nbins_;
  double eta_low_;
  double eta_up_;

  std::string folder_;

  TString title_;
  TString matched_title_;

  MEMap me_detector_;
  MEMap me_detector_matched_;

  MEMap me_muon_pt_; // region-station
  MEMap me_muon_eta_; // region-station
  MEMap me_muon_phi_; // region-station
  MEMap me_chamber_; // region-station-layer
  MEMap me_muon_pt_matched_;
  MEMap me_muon_eta_matched_;
  MEMap me_muon_phi_matched_;
  MEMap me_chamber_matched_;

  MEMap me_residual_rdphi_;    // global
  MEMap me_residual_y_;    // local
  MEMap me_pull_rdphi_;
  MEMap me_pull_y_;
};

inline bool GEMEfficiencyAnalyzer::isInsideOut(const reco::Track& track) {
  return track.innerPosition().mag2() > track.outerPosition().mag2();
}


#endif  // DQMOffline_Muon_GEMEfficiencyAnalyzer_h

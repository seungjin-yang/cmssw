#ifndef DQMOffline_Muon_GEMEfficiencyAnalyzer_h
#define DQMOffline_Muon_GEMEfficiencyAnalyzer_h

#include "DQMOffline/Muon/interface/GEMOfflineDQMBase.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

class GEMEfficiencyAnalyzer : public GEMOfflineDQMBase {
public:
  explicit GEMEfficiencyAnalyzer(const edm::ParameterSet &);
  ~GEMEfficiencyAnalyzer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

protected:
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event &event, const edm::EventSetup &eventSetup) override;

private:
  void bookDetectorOccupancy(
      DQMStore::IBooker &, const GEMStation *, const MEMapKey1 &, const TString &, const TString &);
  void bookOccupancy(DQMStore::IBooker &, const MEMapKey1 &, const TString &, const TString &);
  void bookResolution(DQMStore::IBooker &, const MEMapKey2 &, const TString &, const TString &);

  std::vector<std::pair<ReferenceCountingPointer<Disk>, std::vector<const GEMChamber*> > > buildSurfaces(const edm::ESHandle<GEMGeometry>&);
  const GEMRecHit *findMatchedHit(const GlobalPoint&, const GEMRecHitCollection::range &, const GEMEtaPartition*);
  const GEMEtaPartition* findEtaPartition(const GlobalPoint&, const std::vector<const GEMChamber*>&);
  bool checkBounds(const GlobalPoint&, const GeomDet*);

  edm::EDGetTokenT<GEMRecHitCollection> rechit_token_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muon_token_;

  MuonServiceProxy *muon_service_;

  bool use_global_muon_;
  float rdphi_cut_;

  std::vector<double> pt_binning_;
  int eta_nbins_;
  double eta_low_;
  double eta_up_;

  std::string folder_;

  TString title_;
  TString matched_title_;

  // MEMapRS, MEMapRSL, MEMapRSL, MEMapRSE, MEMapR

  MEMap1 me_detector_;
  MEMap1 me_detector_matched_;

  MEMap1 me_muon_pt_;
  MEMap1 me_muon_eta_;
  MEMap1 me_muon_phi_;
  MEMap1 me_muon_pt_matched_;
  MEMap1 me_muon_eta_matched_;
  MEMap1 me_muon_phi_matched_;

  MEMap2 me_residual_rdphi_;    // global
  MEMap2 me_residual_y_;    // local
  MEMap2 me_pull_rdphi_;
  MEMap2 me_pull_y_;
};

#endif  // DQMOffline_Muon_GEMEfficiencyAnalyzer_h

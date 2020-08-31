#ifndef DQMOffline_Muon_GEMEfficiencyAnalyzer_h
#define DQMOffline_Muon_GEMEfficiencyAnalyzer_h

#include "DQMOffline/Muon/interface/GEMOfflineDQMBase.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

class GEMEfficiencyAnalyzer : public GEMOfflineDQMBase {
public:
  explicit GEMEfficiencyAnalyzer(const edm::ParameterSet &);
  ~GEMEfficiencyAnalyzer() override;

protected:
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event &event, const edm::EventSetup &eventSetup) override;

private:
  // https://github.com/cms-sw/cmssw/blob/CMSSW_11_2_0_pre2/MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h#L56-L57
  enum class MuonType {kNone, kGlobal, kStandalone, kCosmic};

  void bookDetectorOccupancy(
      DQMStore::IBooker &, const GEMStation *, const MEMapKey1 &, const TString &, const TString &);
  void bookOccupancy(DQMStore::IBooker &, const MEMapKey2 &, const TString &, const TString &);
  void bookResolution(DQMStore::IBooker &, const MEMapKey3 &, const TString &, const TString &);

  const GEMEtaPartition* findEtaPartition(const GEMChamber*, const GlobalPoint&);
  const GEMRecHit *findMatchedHit(const float, const GEMRecHitCollection::range &);
  const int getCSCDetailBin(bool, bool, bool, bool);

  const std::string getDetName(const DetId&);
  const std::string getDetName(const unsigned int);
  const unsigned int getDetIdx(unsigned int, bool verbose=false);

  edm::EDGetTokenT<GEMRecHitCollection> rechit_token_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muon_token_;

  MuonServiceProxy *muon_service_;

  std::string muon_type_str_;
  MuonType muon_type_;
  float residual_x_cut_;
  bool use_fiducial_cut_;

  std::vector<double> pt_binning_;
  int eta_nbins_;
  double eta_low_;
  double eta_up_;

  std::string folder_;

  TString title_;
  TString matched_title_;

  // NOTE
  MEMap1 me_detector_;
  MEMap1 me_detector_matched_;

  MEMap1 me_debug_detector_incoming_;
  MEMap1 me_debug_detector_incoming_matched_;

  MEMap1 me_debug_detector_outgoing_;
  MEMap1 me_debug_detector_outgoing_matched_;

  MEMap2 me_muon_pt_;
  MEMap2 me_muon_eta_;
  MEMap2 me_muon_phi_;
  MEMap2 me_muon_pt_matched_;
  MEMap2 me_muon_eta_matched_;
  MEMap2 me_muon_phi_matched_;

  MEMap3 me_residual_x_;    // local
  MEMap3 me_residual_y_;    // local
  MEMap3 me_residual_phi_;  // global
  MEMap3 me_pull_x_;
  MEMap3 me_pull_y_;

  MonitorElement* me_csc_;
  MonitorElement* me_csc_detail_;

  MonitorElement* me_csc_matched_;
  MonitorElement* me_csc_detail_matched_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE
  //////////////////////////////////////////////////////////////////////////////
  MonitorElement* me_debug_normalized_chi2_;
  MonitorElement* me_debug_normalized_chi2_matched_;

  MonitorElement* me_debug_start_state_x_err_det_;
  MonitorElement* me_debug_start_state_x_err_det_matched_;

  MonitorElement* me_debug_error_propagation_;
  MonitorElement* me_debug_num_valid_chambers_per_layer_;

  MonitorElement* me_debug_in_out_det_;
  MonitorElement* me_debug_in_out_det_matched_;

  MonitorElement* me_debug_in_out_det_incoming_;
  MonitorElement* me_debug_in_out_det_incoming_matched_;

  MonitorElement* me_debug_in_out_det_outgoing_;
  MonitorElement* me_debug_in_out_det_outgoing_matched_;

  MonitorElement* me_debug_in_dest_det_;
  MonitorElement* me_debug_in_dest_det_matched_;

  MonitorElement* me_debug_in_dest_det_incoming_;
  MonitorElement* me_debug_in_dest_det_incoming_matched_;

  MonitorElement* me_debug_in_dest_det_outgoing_;
  MonitorElement* me_debug_in_dest_det_outgoing_matched_;

  MonitorElement* me_debug_unmatched_;
  MonitorElement* me_debug_unmatched_no_hit_;

  MonitorElement* me_debug_min_residual_x_;
  MonitorElement* me_debug_min_residual_x_inner_det_;

  MonitorElement* me_debug_ip_incoming_;
  MonitorElement* me_debug_ip_incoming_matched_;
  MonitorElement* me_debug_ip_outgoing_;
  MonitorElement* me_debug_ip_outgoing_matched_;

  unsigned int debug_start_det_idx_;

  const std::vector<std::string> det_labels = {
    "Tracker", "Ecal", "Hcal",
    "GE-1/1",
    "ME-1/1A", "ME-1/1B", "ME-1/2", "ME-1/3",
    "RE-1/2", "RE-1/3",
    "RE-2/2", "RE-2/3",
    "ME-2/1", "ME-2/2",
    "ME-3/1", "ME-3/2",
    "RE-3/2", "RE-3/3",
    "ME-4/1", "ME-4/2",
    "RE-4/2", "RE-4/3",
    "MB",
    "RB",
    "other",
  };

};

#endif  // DQMOffline_Muon_GEMEfficiencyAnalyzer_h

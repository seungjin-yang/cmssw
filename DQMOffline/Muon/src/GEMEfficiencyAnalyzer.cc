#include "DQMOffline/Muon/interface/GEMEfficiencyAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"


GEMEfficiencyAnalyzer::GEMEfficiencyAnalyzer(const edm::ParameterSet& pset) : GEMOfflineDQMBase(pset) {
  rechit_token_ = consumes<GEMRecHitCollection>(pset.getParameter<edm::InputTag>("recHitTag"));
  muon_token_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"));

  auto muon_service_parameter = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  muon_service_ = new MuonServiceProxy(muon_service_parameter, consumesCollector());

  muon_type_str_ = pset.getUntrackedParameter<std::string>("muonType");
  if (muon_type_str_.compare("Global") == 0)
    muon_type_ = MuonType::kGlobal;
  else if (muon_type_str_.compare("Standalone") == 0)
    muon_type_ = MuonType::kStandalone;
  else if (muon_type_str_.compare("Cosmic") == 0)
    muon_type_ = MuonType::kCosmic;
  else
    muon_type_ = MuonType::kNone;

  residual_x_cut_ = static_cast<float>(pset.getParameter<double>("residualXCut"));

  pt_binning_ = pset.getUntrackedParameter<std::vector<double> >("ptBinning");
  eta_nbins_ = pset.getUntrackedParameter<int>("etaNbins");
  eta_low_ = pset.getUntrackedParameter<double>("etaLow");
  eta_up_ = pset.getUntrackedParameter<double>("etaUp");

  folder_ = pset.getUntrackedParameter<std::string>("folder");

  title_ = TString::Format("%s Muon", muon_type_str_.c_str());
  matched_title_ = title_ + TString::Format(" (|x_{Muon} - x_{Hit}| < %.1f)", residual_x_cut_);
}

GEMEfficiencyAnalyzer::~GEMEfficiencyAnalyzer() {}

void GEMEfficiencyAnalyzer::bookHistograms(DQMStore::IBooker& ibooker,
                                           edm::Run const& run,
                                           edm::EventSetup const& isetup) {
  edm::ESHandle<GEMGeometry> gem;
  isetup.get<MuonGeometryRecord>().get(gem);
  if (not gem.isValid()) {
    edm::LogError(log_category_) << "GEMGeometry is invalid" << std::endl;
    return;
  }

  for (const GEMRegion* region : gem->regions()) {
    const int region_number = region->region();
    const char* region_sign = region_number > 0 ? "+" : "-";

    for (const GEMStation* station : region->stations()) {
      const int station_number = station->station();

      const MEMapKey1 key1{region_number, station_number};
      const auto&& station_name_suffix = TString::Format("_ge%s%d1", region_sign, station_number);
      const auto&& station_title_suffix = TString::Format(" : GE %s%d/1", region_sign, station_number);
      bookDetectorOccupancy(ibooker, station, key1, station_name_suffix, station_title_suffix);

      const int num_etas = getNumEtaPartitions(station);

      if (station_number == 1) {
        for (const bool is_odd : {true, false}) {
          std::tuple<int, int, bool> key2{region_number, station_number, is_odd};
          const TString&& parity_name_suffix = station_name_suffix + (is_odd ? "_odd" : "_even");
          const TString&& parity_title_suffix =
              station_title_suffix + (is_odd ? ", Odd Superchamber" : ", Even Superchamber");
          bookOccupancy(ibooker, key2, parity_name_suffix, parity_title_suffix);

          for (int ieta = 1; ieta <= num_etas; ieta++) {
            const TString&& ieta_name_suffix = parity_name_suffix + Form("_ieta%d", ieta);
            const TString&& ieta_title_suffix = parity_title_suffix + Form(", i#eta = %d", ieta);
            const MEMapKey3 key3{region_number, station_number, is_odd, ieta};
            bookResolution(ibooker, key3, ieta_name_suffix, ieta_title_suffix);
          }  // ieta
        }    // is_odd

      } else {
        std::tuple<int, int, bool> key2{region_number, station_number, false};
        bookOccupancy(ibooker, key2, station_name_suffix, station_title_suffix);

        for (int ieta = 1; ieta <= num_etas; ieta++) {
          const MEMapKey3 key3{region_number, station_number, false, ieta};
          const TString&& ieta_name_suffix = station_name_suffix + Form("_ieta%d", ieta);
          const TString&& ieta_title_suffix = station_title_suffix + Form(", i#eta = %d", ieta);
          bookResolution(ibooker, key3, ieta_name_suffix, ieta_title_suffix);
        }  // ieta
      }
    }  // station
  }    // region

  //////////////////////////////////////////////////////////////////////////////
  // NOTE debug - efficiency
  //////////////////////////////////////////////////////////////////////////////

  ibooker.setCurrentFolder(folder_ + "/Efficiency");

  TString csc_title = "the location of segments for CSC hitting at GEM (before matching with GEM rechits)";

  me_csc_ = ibooker.book1D("csc_seg", csc_title, 5, -0.5, 4.5);
  me_csc_matched_ = ibooker.book1D("csc_seg_matched", csc_title, 5, -0.5, 4.5);

  std::vector<std::string> csc_labels = {"None", "1", "2", "3", "4"};

  for (unsigned int idx = 0; idx < csc_labels.size(); idx++) {
    me_csc_->setBinLabel(idx  + 1, csc_labels[idx]);
    me_csc_matched_->setBinLabel(idx + 1, csc_labels[idx]);
  }  

  std::vector<std::string> csc_detail_labels = {
    "None",
    "1", "2", "3", "4",
    "1+2", "1+3", "1+4", "2+3", "2+4", "3+4",
    "1+2+3", "1+2+4", "1+3+4", "2+3+4",
    "All",
  };

  me_csc_detail_ = ibooker.book1D("csc_seg_detail", csc_title, csc_detail_labels.size(), 0.5, csc_detail_labels.size() + 0.5);
  me_csc_detail_matched_ = ibooker.book1D("csc_seg_detail_matched", csc_title, csc_detail_labels.size(), 0.5, csc_detail_labels.size() + 0.5);

  for (unsigned int idx = 0; idx < csc_detail_labels.size(); idx++) {
    me_csc_detail_->setBinLabel(idx + 1,  csc_detail_labels[idx]);
    me_csc_detail_matched_->setBinLabel(idx + 1,  csc_detail_labels[idx]);
  }

  me_debug_start_state_x_err_ = ibooker.book1D("start_state_x_err", "x error", 100, 0, 0.2);
  me_debug_start_state_x_err_matched_ = ibooker.book1D("start_state_x_err_matched", "x error_matched", 100, 0, 0.2);

  me_debug_dest_state_x_err_ = ibooker.book1D("dest_state_x_err", "x error", 100, 0, 10);
  me_debug_dest_state_x_err_matched_ = ibooker.book1D("dest_state_x_err_matched", "x error_matched", 100, 0, 10);

  me_debug_start_state_x_err_det_ = ibooker.book1D(
      "start_state_x_err_det", "x error",
      det_labels.size(), -0.5, det_labels.size() - 0.5);

  me_debug_start_state_x_err_det_matched_ = ibooker.book1D(
      "start_state_x_err_det_matched", "x error",
      det_labels.size(), -0.5, det_labels.size() - 0.5);

  me_debug_start_state_x_err_det_detail_ = ibooker.book2D(
      "start_state_x_err_det_detail", "x error",
      det_labels.size(), -0.5, det_labels.size() - 0.5,
      40, 0, 0.2);

  me_debug_start_state_x_err_det_detail_matched_ = ibooker.book2D(
      "start_state_x_err_det_detail_matched", "x error",
      det_labels.size(), -0.5, det_labels.size() - 0.5,
      40, 0, 0.2);

  me_debug_normalized_chi2_ = ibooker.book1D(
      "normalized_chi2", "#chi^2 / ndof",
      100, 0, 10);

  me_debug_normalized_chi2_matched_ = ibooker.book1D(
      "normalized_chi2_matched", "#chi^2 / ndof (matched)",
      100, 0, 10);


  me_debug_in_out_det_ = ibooker.book2D(
      "in_out_det", ";innerDet;outerDet;",
      det_labels.size(), -0.5, det_labels.size() - 0.5,
      det_labels.size(), -0.5, det_labels.size() - 0.5);

  me_debug_in_out_det_matched_ = ibooker.book2D(
      "in_out_det_matched", ";innerDet;outerDet;",
      det_labels.size(), -0.5, det_labels.size() - 0.5,
      det_labels.size(), -0.5, det_labels.size() - 0.5);

  for (unsigned int idx = 0; idx < det_labels.size(); idx++) {
    me_debug_in_out_det_->setBinLabel(idx + 1, det_labels[idx]);
    me_debug_in_out_det_matched_->setBinLabel(idx + 1, det_labels[idx]);

    me_debug_in_out_det_->setBinLabel(idx + 1, det_labels[idx], 2);
    me_debug_in_out_det_matched_->setBinLabel(idx + 1, det_labels[idx], 2);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  //////////////////////////////////////////////////////////////////////////////
  ibooker.setCurrentFolder(folder_ + "/Debug");
  me_debug_error_propagation_ = ibooker.book2D(
      "error_propagation",
      "Error;#sigma_{x}, Start;#sigma_{x}, Destination",
      40, 0, 2,
      100, 0, 20);

  me_debug_num_valid_chambers_per_layer_ = ibooker.book1D(
      "num_valid_chambers_per_layer", ";Number of chambers per layer;Entries",
      37, -0.5, 36.5);
}

void GEMEfficiencyAnalyzer::bookDetectorOccupancy(DQMStore::IBooker& ibooker,
                                                  const GEMStation* station,
                                                  const MEMapKey1& key,
                                                  const TString& name_suffix,
                                                  const TString& title_suffix) {
  ibooker.setCurrentFolder(folder_ + "/Efficiency");
  BookingHelper helper(ibooker, name_suffix, title_suffix);

  const auto&& superchambers = station->superChambers();
  if (not checkRefs(superchambers)) {
    edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
    return;
  }

  // the number of GEMChambers per GEMStation
  const int num_ch = superchambers.size() * superchambers.front()->nChambers();
  // the number of eta partitions per GEMChamber
  const int num_etas = getNumEtaPartitions(station);

  me_detector_[key] = helper.book2D("detector", title_, num_ch, 0.5, num_ch + 0.5, num_etas, 0.5, num_etas + 0.5);

  me_detector_matched_[key] =
      helper.book2D("detector_matched", matched_title_, num_ch, 0.5, num_ch + 0.5, num_etas, 0.5, num_etas + 0.5);

  setDetLabelsEta(me_detector_[key], station);
  setDetLabelsEta(me_detector_matched_[key], station);
}

void GEMEfficiencyAnalyzer::bookOccupancy(DQMStore::IBooker& ibooker,
                                          const MEMapKey2& key,
                                          const TString& name_suffix,
                                          const TString& title_suffix) {
  ibooker.setCurrentFolder(folder_ + "/Efficiency");
  BookingHelper helper(ibooker, name_suffix, title_suffix);

  me_muon_pt_[key] = helper.book1D("muon_pt", title_, pt_binning_, "Muon p_{T} [GeV]");
  me_muon_eta_[key] = helper.book1D("muon_eta", title_, eta_nbins_, eta_low_, eta_up_, "Muon |#eta|");
  me_muon_phi_[key] = helper.book1D("muon_phi", title_, 36, -M_PI, M_PI, "Muon #phi");

  me_muon_pt_matched_[key] = helper.book1D("muon_pt_matched", matched_title_, pt_binning_, "Muon p_{T} [GeV]");
  me_muon_eta_matched_[key] =
      helper.book1D("muon_eta_matched", matched_title_, eta_nbins_, eta_low_, eta_up_, "Muon |#eta|");
  me_muon_phi_matched_[key] =
      helper.book1D("muon_phi_matched", matched_title_, 36, -M_PI, M_PI, "Muon #phi");

}

void GEMEfficiencyAnalyzer::bookResolution(DQMStore::IBooker& ibooker,
                                           const MEMapKey3& key,
                                           const TString& name_suffix,
                                           const TString& title_suffix) {
  ibooker.setCurrentFolder(folder_ + "/Resolution");
  BookingHelper helper(ibooker, name_suffix, title_suffix);

  // NOTE Residual & Pull
  me_residual_x_[key] = helper.book1D("residual_x", title_, 50, -residual_x_cut_, residual_x_cut_, "Residual in Local X [cm]");
  me_residual_y_[key] = helper.book1D("residual_y", title_, 60, -12.0, 12.0, "Residual in Local Y [cm]");
  me_residual_phi_[key] = helper.book1D("residual_phi", title_, 80, -0.008, 0.008, "Residual in Global #phi [rad]");

  me_pull_x_[key] = helper.book1D("pull_x", title_, 60, -3.0, 3.0, "Pull in Local X");
  me_pull_y_[key] = helper.book1D("pull_y", title_, 60, -3.0, 3.0, "Pull in Local Y");
}

void GEMEfficiencyAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  // FIXME
  if (muon_type_ == MuonType::kNone) {
    edm::LogError(log_category_) << "wrong muon name" << std::endl;
    return;
  }


  edm::Handle<GEMRecHitCollection> rechit_collection;
  event.getByToken(rechit_token_, rechit_collection);
  if (not rechit_collection.isValid()) {
    edm::LogError(log_category_) << "GEMRecHitCollection is invalid" << std::endl;
    return;
  }

  edm::Handle<edm::View<reco::Muon> > muon_view;
  event.getByToken(muon_token_, muon_view);
  if (not muon_view.isValid()) {
    edm::LogError(log_category_) << "View<Muon> is invalid" << std::endl;
    return;
  }

  edm::ESHandle<GEMGeometry> gem;
  setup.get<MuonGeometryRecord>().get(gem);
  if (not gem.isValid()) {
    edm::LogError(log_category_) << "GEMGeometry is invalid" << std::endl;
    return;
  }

  edm::ESHandle<CSCGeometry> csc;
  setup.get<MuonGeometryRecord>().get(csc);
  if (not csc.isValid()) {
    edm::LogError(log_category_) << "CSCGeometry is invalid" << std::endl;
    return;
  }

  edm::ESHandle<TransientTrackBuilder> transient_track_builder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", transient_track_builder);
  if (not transient_track_builder.isValid()) {
    edm::LogError(log_category_) << "TransientTrackRecord is invalid" << std::endl;
    return;
  }

  muon_service_->update(setup);

  edm::ESHandle<Propagator>&& propagator_any = muon_service_->propagator("SteppingHelixPropagatorAny");
  if (not propagator_any.isValid()) {
    edm::LogError(log_category_) << "SteppingHelixPropagatorAny is invalid" << std::endl;
    return;
  }

  edm::ESHandle<Propagator>&& propagator_along = muon_service_->propagator("SteppingHelixPropagatorAlong");
  if (not propagator_along.isValid()) {
    edm::LogError(log_category_) << "SteppingHelixPropagatorAlong is invalid" << std::endl;
    return;
  }

  edm::ESHandle<Propagator>&& propagator_opposite = muon_service_->propagator("SteppingHelixPropagatorOpposite");
  if (not propagator_opposite.isValid()) {
    edm::LogError(log_category_) << "SteppingHelixPropagatorOpposite is invalid" << std::endl;
    return;
  }

  edm::ESHandle<MagneticField>&& magnetic_field = muon_service_->magneticField();

  std::map<PropagationDirection, edm::ESHandle<Propagator> > propagator_map;
  propagator_map.emplace(PropagationDirection::anyDirection, std::move(propagator_any));
  propagator_map.emplace(PropagationDirection::alongMomentum, std::move(propagator_along));
  propagator_map.emplace(PropagationDirection::oppositeToMomentum, std::move(propagator_opposite));

  for (const reco::Muon& muon : *muon_view) {
    const reco::Track* track = nullptr;
    if ((muon_type_ == MuonType::kGlobal) and muon.innerTrack().isNonnull()) {
      track = muon.innerTrack().get();

    } else if ((muon_type_ != MuonType::kGlobal) and muon.outerTrack().isNonnull()) {
      track = muon.outerTrack().get();
    }

    if (track == nullptr) {
      edm::LogError(log_category_) << "failed to get muon track" << std::endl;
      continue;
    }

    // NOTE TransientTrack
    const reco::TransientTrack&& transient_track = transient_track_builder->build(track);
    if (not transient_track.isValid()) {
      edm::LogError(log_category_) << "failed to build TransientTrack" << std::endl;
      continue;
    }

    const int num_csc_st1 = muon.numberOfSegments(1, MuonSubdetId::CSC);
    const int num_csc_st2 = muon.numberOfSegments(2, MuonSubdetId::CSC);
    const int num_csc_st3 = muon.numberOfSegments(3, MuonSubdetId::CSC);
    const int num_csc_st4 = muon.numberOfSegments(4, MuonSubdetId::CSC);

    const bool has_csc_st1 = num_csc_st1 > 0;
    const bool has_csc_st2 = num_csc_st2 > 0;
    const bool has_csc_st3 = num_csc_st3 > 0;
    const bool has_csc_st4 = num_csc_st4 > 0;

    const bool has_no_csc_st = not (has_csc_st1 or has_csc_st2 or has_csc_st3 or has_csc_st4);
    if (has_no_csc_st) {
      std::cout << "no csc Station" << std::endl;
      continue;
    }

    TrajectoryStateOnSurface start_state;
    PropagationDirection propagation_direction;
    bool okay = false;

    unsigned int start_raw_det_id = 0;
    unsigned int inner_det_id = 0;
    unsigned int outer_det_id = 0;
    bool is_from_outside_to_inside = false;

    switch (muon_type_) {
      case MuonType::kGlobal: {
        // start_state = transient_track.outermostMeasurementState().freeState();
        start_state = std::move(transient_track.outermostMeasurementState());
        // propagation_direction = PropagationDirection::alongMomentum;
        propagation_direction = PropagationDirection::anyDirection;
        okay = true;
        break;
      }

      case MuonType::kStandalone: {
        // start_state = transient_track.outermostMeasurementState().freeState();
        start_state = std::move(transient_track.outermostMeasurementState());
        // propagation_direction = PropagationDirection::oppositeToMomentum;
        propagation_direction = PropagationDirection::anyDirection;
        okay = true;
        break;
      }

      case MuonType::kCosmic: {
        float x2_in = std::hypot(track->innerPosition().x(), track->innerPosition().y());
        float x2_out = std::hypot(track->outerPosition().x(), track->outerPosition().y());

        float p2_in = track->innerMomentum().mag2();
        float p2_out = track->outerMomentum().mag2();

        if (x2_in < x2_out) {
          start_state = std::move(transient_track.outermostMeasurementState());
          start_raw_det_id = track->outerDetId();

          inner_det_id = track->innerDetId();
          outer_det_id = track->outerDetId();

          // start_state = std::move(transient_track.innermostMeasurementState());
          // start_raw_det_id = track->innerDetId();

          const std::string in_det_name = getDetName(track->innerDetId());
          const std::string out_det_name = getDetName(track->outerDetId());

        } else {
          start_state = std::move(transient_track.innermostMeasurementState());
          start_raw_det_id = track->innerDetId();

          inner_det_id = track->outerDetId();
          outer_det_id = track->innerDetId();

          // start_state = std::move(transient_track.outermostMeasurementState());
          // start_raw_det_id = track->outerDetId();

          const std::string in_det_name = getDetName(track->outerDetId());
          const std::string out_det_name = getDetName(track->innerDetId());

          std::swap(p2_in, p2_out);
        }

        propagation_direction = (p2_in > p2_out) ? PropagationDirection::oppositeToMomentum : PropagationDirection::alongMomentum;
        // propagation_direction = PropagationDirection::anyDirection;
        is_from_outside_to_inside = p2_out > p2_in;
        okay = true;
        break;
      }

      default: {
        propagation_direction = PropagationDirection::anyDirection;
        break;
      }
    }

    if (not okay) {
      edm::LogError(log_category_) << "sth wrong" << std::endl;
    }

    if (not start_state.isValid()) {
      edm::LogError(log_category_) << "A starting TrajectoryStateOnSurface is not valid!" << std::endl;
    }

    const unsigned int start_det_idx = getDetIdx(start_raw_det_id);
    const unsigned int inner_det_idx = getDetIdx(inner_det_id);
    const unsigned int outer_det_idx = getDetIdx(outer_det_id);
    const float start_state_x_err = std::sqrt(start_state.localError().positionError().xx());


    const int csc_detail_bin = getCSCDetailBin(has_csc_st1, has_csc_st2, has_csc_st3, has_csc_st4);

    // NOTE Currently, only GE11 is being considered.
    for (const GEMRegion* region : gem->regions()) {

      // NOTE
      bool is_opposite_region = (muon.eta() * region->region() < 0);

      if (muon_type_ == MuonType::kCosmic) {
        if (is_from_outside_to_inside xor is_opposite_region)
          continue;

      } else {
        // pp collision
        if (is_opposite_region)
          continue;
      }

      for (const GEMStation* station : region->stations()) {

        std::map<int, std::vector<const GEMChamber*> > chambers_per_layer;
        for (int layer_id : {1, 2}) {
          chambers_per_layer.emplace(layer_id, std::vector<const GEMChamber*>());
          chambers_per_layer[layer_id].reserve(station->superChambers().size());
        }

        for (const GEMSuperChamber* super_chamber : station->superChambers()) {
          for (const GEMChamber* chamber : super_chamber->chambers()) {
            chambers_per_layer[chamber->id().layer()].push_back(chamber);
          }
        }

        for (const auto [layer_id, layer] : chambers_per_layer) {
          int debug_num_valid_chambers_per_layer = 0;

          for (const GEMChamber* chamber : layer) {            
            const BoundPlane& bound_plane = chamber->surface();

            const auto&& tsos = propagator_map[propagation_direction]->propagate(start_state, bound_plane);
            if (not tsos.isValid()) {
              edm::LogInfo(log_category_) << "failed to propagate ," << chamber->id() << std::endl;
              continue;
            }

            const GlobalPoint&& tsos_global_pos = tsos.globalPosition();
            const LocalError&& tsos_err = tsos.localError().positionError();
            const GEMEtaPartition* eta_partition = findEtaPartition(chamber, tsos_global_pos);
            if (eta_partition == nullptr) {
              edm::LogInfo(log_category_) << "failed to find GEMEtaPartition" << std::endl; 
              continue;
            }

            auto debug_range = rechit_collection->get(chamber->id());
            auto num_rechits_on_chamber = std::distance(debug_range.first, debug_range.second);

            if (muon_type_ == MuonType::kCosmic) {
              std::cout << "GEM_DEBUG: " << num_rechits_on_chamber << std::endl;
            }
            me_debug_error_propagation_->Fill(start_state_x_err, std::sqrt(tsos_err.xx()));
            me_debug_in_out_det_->Fill(inner_det_idx, outer_det_idx);
            debug_num_valid_chambers_per_layer++;

            if (has_no_csc_st) me_csc_->Fill(0);
            if (has_csc_st1)   me_csc_->Fill(1);
            if (has_csc_st2)   me_csc_->Fill(2);
            if (has_csc_st3)   me_csc_->Fill(3);
            if (has_csc_st4)   me_csc_->Fill(4);
            me_csc_detail_->Fill(csc_detail_bin);

            const GEMDetId&& gem_id = eta_partition->id();

            bool is_odd = gem_id.station() == 1 ? (gem_id.chamber() % 2 == 1) : false;
            const std::tuple<int, int> key1{gem_id.region(), gem_id.station()};
            const std::tuple<int, int, bool> key2{gem_id.region(), gem_id.station(), is_odd};
            const std::tuple<int, int, bool, int> key3{gem_id.region(), gem_id.station(), is_odd, gem_id.roll()};

            const int chamber_bin = getDetOccXBin(gem_id, gem);

            fillME(me_detector_, key1, chamber_bin, gem_id.roll());
            fillME(me_muon_pt_, key2, muon.pt());
            fillME(me_muon_eta_, key2, std::fabs(muon.eta()));
            fillME(me_muon_phi_, key2, muon.phi());

            me_debug_start_state_x_err_->Fill(start_state_x_err);
            me_debug_start_state_x_err_det_->Fill(start_det_idx);
            me_debug_start_state_x_err_det_detail_->Fill(start_det_idx, start_state_x_err);

            me_debug_dest_state_x_err_->Fill(std::sqrt(tsos_err.xx()));

            me_debug_normalized_chi2_->Fill(track->normalizedChi2());

            const LocalPoint&& tsos_local_pos = tsos.localPosition();
            const GEMRecHit* matched_hit = findMatchedHit(tsos_local_pos.x(), rechit_collection->get(gem_id));
            if (matched_hit == nullptr) {
              continue;
            }

            if (has_no_csc_st) me_csc_matched_->Fill(0);
            if (has_csc_st1)   me_csc_matched_->Fill(1);
            if (has_csc_st2)   me_csc_matched_->Fill(2);
            if (has_csc_st3)   me_csc_matched_->Fill(3);
            if (has_csc_st4)   me_csc_matched_->Fill(4);
            me_csc_detail_matched_->Fill(csc_detail_bin);
            me_debug_in_out_det_matched_->Fill(inner_det_idx, outer_det_idx);


            fillME(me_detector_matched_, key1, chamber_bin, gem_id.roll());
            fillME(me_muon_pt_matched_, key2, muon.pt());
            fillME(me_muon_eta_matched_, key2, std::fabs(muon.eta()));
            fillME(me_muon_phi_matched_, key2, muon.phi());

            me_debug_start_state_x_err_matched_->Fill(start_state_x_err);
            me_debug_start_state_x_err_det_matched_->Fill(start_det_idx);
            me_debug_start_state_x_err_det_detail_matched_->Fill(start_det_idx, start_state_x_err);

            me_debug_dest_state_x_err_matched_->Fill(std::sqrt(tsos_err.xx()));

            me_debug_normalized_chi2_matched_->Fill(track->normalizedChi2());

            const LocalPoint&& hit_local_pos = matched_hit->localPosition();
            const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit_local_pos);

            const float residual_x = tsos_local_pos.x() - hit_local_pos.x();
            const float residual_y = tsos_local_pos.y() - hit_local_pos.y();
            const float residual_phi = reco::deltaPhi(tsos_global_pos.barePhi(), hit_global_pos.barePhi());

            const LocalError&& hit_err = matched_hit->localPositionError();

            const float pull_x = residual_x / std::sqrt(tsos_err.xx() + hit_err.xx());
            const float pull_y = residual_y / std::sqrt(tsos_err.yy() + hit_err.yy());

            fillME(me_residual_x_, key3, residual_x);
            fillME(me_residual_y_, key3, residual_y);
            fillME(me_residual_phi_, key3, residual_phi);

            fillME(me_pull_x_, key3, pull_x);
            fillME(me_pull_y_, key3, pull_y);
          } // chamber

          me_debug_num_valid_chambers_per_layer_->Fill(debug_num_valid_chambers_per_layer);

        } // layer
      } // GEMStation
    } // GEMRegion

  }    // Muon
}

const GEMEtaPartition* GEMEfficiencyAnalyzer::findEtaPartition(const GEMChamber* chamber, const GlobalPoint& global_point){
  for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
    const LocalPoint&& local_point = eta_partition->toLocal(global_point);
    const LocalPoint local_point_2d(local_point.x(), local_point.y(), 0.0f);
    if (eta_partition->surface().bounds().inside(local_point_2d)) 
      return eta_partition;
  }

  return nullptr;
}


const GEMRecHit* GEMEfficiencyAnalyzer::findMatchedHit(const float track_local_x,
                                                       const GEMRecHitCollection::range& range) {
  float min_residual_x{residual_x_cut_};
  const GEMRecHit* closest_hit = nullptr;

  for (auto hit = range.first; hit != range.second; ++hit) {
    float residual_x = std::fabs(track_local_x - hit->localPosition().x());
    if (residual_x <= min_residual_x) {
      min_residual_x = residual_x;
      closest_hit = &(*hit);
    }
  }

  return closest_hit;
}


const int GEMEfficiencyAnalyzer::getCSCDetailBin(bool has_st1, bool has_st2, bool has_st3, bool has_st4) {
  int xbin = 0;
  int num_stations = 0;

  if (has_st1) {
    xbin += 1;
    num_stations++;
  }

  if (has_st2) {
    xbin += 2;
    num_stations++;
  }

  if (has_st3) {
    xbin += 3;
    num_stations++;
  }

  if (has_st4) {
    xbin += 4;
    num_stations++;
  }

  if (num_stations == 0)
    xbin = 1;
  else if (num_stations == 1)
    xbin += 1;
  else if (num_stations == 2) 
    xbin += has_st1 ? 3 : 4;
  else
    xbin += 6;

  return xbin;
}


const std::string GEMEfficiencyAnalyzer::getDetName(const DetId& det_id) {
  std::string det_name;
  if (det_id.det() == DetId::Detector::Muon) {
    if (det_id.subdetId() == MuonSubdetId::DT) {
      const DTWireId dt_id{det_id};
      std::stringstream buffer;
      buffer << dt_id;
      det_name = "DT: " + buffer.str();

    } else if (det_id.subdetId() == MuonSubdetId::CSC) {
      const CSCDetId csc_id{det_id};
      std::stringstream buffer;
      buffer << csc_id;
      det_name = "CSC: " + buffer.str();

    } else if (det_id.subdetId() == MuonSubdetId::RPC) {
      const RPCDetId rpc_id{det_id};
      std::stringstream buffer;
      buffer << rpc_id;
      det_name = "RPC: " + buffer.str();

    } else if (det_id.subdetId() == MuonSubdetId::GEM) {
      const GEMDetId gem_id{det_id};
      std::stringstream buffer;
      buffer << gem_id;
      det_name = "GEM: " + buffer.str();

    } else {
      det_name = "wrong muon det";
    }

  } else { 
    det_name = "wront det";
  }
  return det_name;
}


const std::string GEMEfficiencyAnalyzer::getDetName(const unsigned int raw_det_id) {
  const DetId det_id{raw_det_id};
  return getDetName(det_id);
}


const unsigned int GEMEfficiencyAnalyzer::getDetIdx(unsigned int raw_det_id) {
  DetId det_id{raw_det_id};
  std::string label = "unknown";

  if (det_id.det() == DetId::Detector::Muon) {
    if (det_id.subdetId() == MuonSubdetId::DT) {
      const DTWireId dt_id{det_id};
      label = "MB";

    } else if (det_id.subdetId() == MuonSubdetId::CSC) {
      const CSCDetId csc_id{det_id};
      label = Form("ME%d/%d", csc_id.zendcap() * csc_id.station(), csc_id.ring());

    } else if (det_id.subdetId() == MuonSubdetId::RPC) {
      const RPCDetId rpc_id{det_id};

      if (rpc_id.region() == 0) 
        label = "RB";
      else
        label = Form("RE%d/%d", rpc_id.region() * rpc_id.station(), rpc_id.ring());

    } else if (det_id.subdetId() == MuonSubdetId::GEM) {
      const GEMDetId gem_id{det_id};
      label = Form("GE%d/1", gem_id.region() * gem_id.station());
    }
  }

  unsigned int idx = std::find(det_labels.begin(), det_labels.end(), label) - det_labels.begin();
  if (idx == det_labels.size())
    idx--;

  return idx;
}

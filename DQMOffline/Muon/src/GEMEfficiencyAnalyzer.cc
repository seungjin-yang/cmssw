#include "DQMOffline/Muon/interface/GEMEfficiencyAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Math/interface/deltaPhi.h"

GEMEfficiencyAnalyzer::GEMEfficiencyAnalyzer(const edm::ParameterSet& pset) : GEMOfflineDQMBase(pset) {
  rechit_token_ = consumes<GEMRecHitCollection>(pset.getParameter<edm::InputTag>("recHitTag"));
  muon_token_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"));

  auto muon_service_parameter = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  muon_service_ = new MuonServiceProxy(muon_service_parameter, consumesCollector());
  propagator_name_ = pset.getUntrackedParameter<std::string>("propagatorName");

  use_global_muon_ = pset.getUntrackedParameter<bool>("useGlobalMuon");

  residual_x_cut_ = static_cast<float>(pset.getParameter<double>("residualXCut"));

  pt_binning_ = pset.getUntrackedParameter<std::vector<double> >("ptBinning");
  eta_nbins_ = pset.getUntrackedParameter<int>("etaNbins");
  eta_low_ = pset.getUntrackedParameter<double>("etaLow");
  eta_up_ = pset.getUntrackedParameter<double>("etaUp");

  folder_ = pset.getUntrackedParameter<std::string>("folder");

  title_ = (use_global_muon_ ? "Global Muon" : "Standalone Muon");
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

  TString csc_title = "the location of segments for CSC hitting at GEM (before matching with GEM rechits)";

  me_csc_ = ibooker.book1D("zzz_csc_seg_st", csc_title, 4, 0.5, 4.5);
  for (int xbin = 1; xbin <= 4; xbin++) {
    me_csc_->setBinLabel(xbin, std::to_string(xbin));
  }  

  std::vector<std::string> csc_detail_labels = {
    "1", "2", "3", "4",
    "1+2", "1+3", "1+4", "2+3", "2+4", "3+4",
    "1+2+3", "1+2+4", "1+3+4", "2+3+4",
    "All", "None",
  };

  me_csc_detail_ = ibooker.book1D("zzz_csc_seg_st_detail", csc_title, csc_detail_labels.size(), 0.5, csc_detail_labels.size() + 0.5);
  for (unsigned int idx = 0; idx < csc_detail_labels.size(); idx++) {
    me_csc_detail_->setBinLabel(idx + 1,  csc_detail_labels[idx]);
  }

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
  }

  edm::ESHandle<GEMGeometry> gem;
  setup.get<MuonGeometryRecord>().get(gem);
  if (not gem.isValid()) {
    edm::LogError(log_category_) << "GEMGeometry is invalid" << std::endl;
    return;
  }

  edm::ESHandle<TransientTrackBuilder> transient_track_builder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", transient_track_builder);
  if (not transient_track_builder.isValid()) {
    edm::LogError(log_category_) << "TransientTrackRecord is invalid" << std::endl;
    return;
  }

  muon_service_->update(setup);
  edm::ESHandle<Propagator>&& propagator = muon_service_->propagator(propagator_name_);
  if (not propagator.isValid()) {
    edm::LogError(log_category_) << "Propagator is invalid" << std::endl;
    return;
  }


  for (const reco::Muon& muon : *muon_view) {


    
    const reco::Track* track = nullptr;
    if (use_global_muon_ and muon.innerTrack().isNonnull()) {
      track = muon.innerTrack().get();

    } else if ((not use_global_muon_) and muon.outerTrack().isNonnull()) {
      track = muon.outerTrack().get();
    }

    if (track == nullptr) {
      edm::LogError(log_category_) << "failed to get muon track" << std::endl;
      continue;
    }

    const reco::TransientTrack&& transient_track = transient_track_builder->build(track);
    if (not transient_track.isValid()) {
      edm::LogError(log_category_) << "failed to build TransientTrack" << std::endl;
      continue;
    }

    const auto&& start_tsos = use_global_muon_ ? transient_track.outermostMeasurementState() : transient_track.innermostMeasurementState();



    if (not start_tsos.isValid()) {
      edm::LogError(log_category_) << "failed to get MeasurementState" << std::endl;
    }

    // unsigned int station_mask = muon.stationMask();
    // const bool has_csc_st1 = (station_mask & 1 << 4) != 0;
    // const bool has_csc_st2 = (station_mask & 1 << 5) != 0;
    // const bool has_csc_st3 = (station_mask & 1 << 6) != 0;
    // const bool has_csc_st4 = (station_mask & 1 << 7) != 0;

    const int num_csc_st1 = muon.numberOfSegments(1, MuonSubdetId::CSC);
    const int num_csc_st2 = muon.numberOfSegments(2, MuonSubdetId::CSC);
    const int num_csc_st3 = muon.numberOfSegments(3, MuonSubdetId::CSC);
    const int num_csc_st4 = muon.numberOfSegments(4, MuonSubdetId::CSC);

    const bool has_csc_st1 = num_csc_st1 > 0;
    const bool has_csc_st2 = num_csc_st2 > 0;
    const bool has_csc_st3 = num_csc_st3 > 0;
    const bool has_csc_st4 = num_csc_st4 > 0;

    const int csc_detail_bin = getCSCDetailBin(has_csc_st1, has_csc_st2, has_csc_st3, has_csc_st4);

    if (log_category_.find("Cosmic") != std::string::npos) {
      std::cout << log_category_ << ": "
                << Form("%d, %d, %d, %d --> %d", num_csc_st1, num_csc_st2, num_csc_st3, num_csc_st4, csc_detail_bin)
                << std::endl;
    }

    for (const GEMChamber* chamber : gem->chambers()) {
      const BoundPlane& bound_plane = chamber->surface();

      const TrajectoryStateOnSurface&& tsos = propagator->propagate(start_tsos, bound_plane);
      if (not tsos.isValid()) {
        edm::LogInfo(log_category_) << "failed to propagate ," << chamber->id() << std::endl;
        continue;
      }

      const GlobalPoint&& tsos_global_pos = tsos.globalPosition();
      const GEMEtaPartition* eta_partition = findEtaPartition(chamber, tsos_global_pos);
      if (eta_partition == nullptr) {
        edm::LogInfo(log_category_) << "failed to find GEMEtaPartition" << std::endl; 
        continue;
      }

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

      const LocalPoint&& tsos_local_pos = tsos.localPosition();
      const GEMRecHit* matched_hit = findMatchedHit(tsos_local_pos.x(), rechit_collection->get(gem_id));
      if (matched_hit == nullptr) {
        continue;
      }

      fillME(me_detector_matched_, key1, chamber_bin, gem_id.roll());
      fillME(me_muon_pt_matched_, key2, muon.pt());
      fillME(me_muon_eta_matched_, key2, std::fabs(muon.eta()));
      fillME(me_muon_phi_matched_, key2, muon.phi());

      const LocalPoint&& hit_local_pos = matched_hit->localPosition();
      const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit_local_pos);

      const float residual_x = tsos_local_pos.x() - hit_local_pos.x();
      const float residual_y = tsos_local_pos.y() - hit_local_pos.y();
      const float residual_phi = reco::deltaPhi(tsos_global_pos.barePhi(), hit_global_pos.barePhi());

      const LocalError&& tsos_err = tsos.localError().positionError();
      const LocalError&& hit_err = matched_hit->localPositionError();

      const float pull_x = residual_x / std::sqrt(tsos_err.xx() + hit_err.xx());
      const float pull_y = residual_y / std::sqrt(tsos_err.yy() + hit_err.yy());

      fillME(me_residual_x_, key3, residual_x);
      fillME(me_residual_y_, key3, residual_y);
      fillME(me_residual_phi_, key3, residual_phi);

      fillME(me_pull_x_, key3, pull_x);
      fillME(me_pull_y_, key3, pull_y);
    }  // GEMChamber
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

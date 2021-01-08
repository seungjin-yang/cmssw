#include "DQMOffline/Muon/interface/GEMEfficiencyAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Validation/MuonGEMHits/interface/GEMValidationUtils.h"

#include "TVector2.h"


GEMEfficiencyAnalyzer::GEMEfficiencyAnalyzer(const edm::ParameterSet& pset) : GEMOfflineDQMBase(pset) {
  rechit_token_ = consumes<GEMRecHitCollection>(pset.getParameter<edm::InputTag>("recHitTag"));
  muon_token_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"));

  auto muon_service_parameter = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  muon_service_ = new MuonServiceProxy(muon_service_parameter, consumesCollector());

  is_cosmics_ = pset.getUntrackedParameter<bool>("isCosmics");
  use_global_muon_ = pset.getUntrackedParameter<bool>("useGlobalMuon");

  rdphi_cut_ = static_cast<float>(pset.getParameter<double>("RDeltaPhiCut"));

  pt_binning_ = pset.getUntrackedParameter<std::vector<double> >("ptBinning");
  eta_nbins_ = pset.getUntrackedParameter<int>("etaNbins");
  eta_low_ = pset.getUntrackedParameter<double>("etaLow");
  eta_up_ = pset.getUntrackedParameter<double>("etaUp");

  folder_ = pset.getUntrackedParameter<std::string>("folder");

  // FIXME looks ugly
  title_ = (use_global_muon_ ? "Global Muon" : "Standalone Muon");
  matched_title_ = title_ + TString::Format(" (R_{Muon} #times |#phi_{Muon} - #phi_{Hit}| < %.1f cm)", rdphi_cut_);
}

GEMEfficiencyAnalyzer::~GEMEfficiencyAnalyzer() {}

void GEMEfficiencyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("recHitTag", edm::InputTag("gemRecHits"));
  desc.add<edm::InputTag>("muonTag", edm::InputTag("muons")); // FIXME
  {
    edm::ParameterSetDescription psd0;
    psd0.setAllowAnything();
    desc.add<edm::ParameterSetDescription>("ServiceParameters", psd0);
  }
  desc.addUntracked<bool>("isCosmics", false);
  desc.addUntracked<bool>("useGlobalMuon", true); // FIXME looks ugly..
  desc.add<double>("RDeltaPhiCut", 2.0); // TODO need to be tuned
  desc.addUntracked<std::vector<double> >("ptBinning", {20. ,30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 200.});
  desc.addUntracked<int>("etaNbins", 7);
  desc.addUntracked<double>("etaLow", 1.5);
  desc.addUntracked<double>("etaUp", 2.2);
  desc.addUntracked<std::string>("folder", "GEM/GEMEfficiency/GEMEfficiencyAnalyzer"); // FIXME
  desc.addUntracked<std::string>("logCategory", "GEMEfficiencyAnalyzer"); // FIXME
  descriptions.add("gemEfficiencyAnalyzerDefault", desc);
}

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
    const int region_id = region->region();

    for (const GEMStation* station : region->stations()) {
      const int station_id = station->station();

      const GEMDetId&& rs_key = getReStKey(region_id, station_id);
      // FIXME waiting for yeckang:gemValidationTools PR
      // const auto&& station_name_suffix = GEMUtils::getSuffixName(region_id, station_id);
      // const auto&& station_title_suffix = GEMUtils::getSuffixTitle(region_id, station_id);
      const auto&& rs_name_suffix = TString::Format("_GE%+.2d", region_id * (station_id * 10 + 1));
      const auto&& rs_title_suffix = TString::Format(" : GE%+.2d", region_id * (station_id * 10 + 1));
      bookDetectorOccupancy(ibooker, station, rs_key, rs_name_suffix, rs_title_suffix);
      bookOccupancy(ibooker, rs_key, rs_name_suffix, rs_title_suffix);

      const auto&& superchambers = station->superChambers();
      if (not checkRefs(superchambers)) {
        edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
        return;
      }
      // the number of GEMChamber per each layer
      const int num_chambers_per_layer = superchambers.size();
      for (const GEMChamber* chamber : superchambers[0]->chambers()) {
        const int layer_id = chamber->id().layer();

        // FIXME waiting for yeckang:gemValidationTools PR
        // const auto&& rsl_name_suffix = GEMUtils::getSuffixName(region_id, station_id, layer_id);
        // const auto&& rsl_title_suffix = GEMUtils::getSuffixTitle(region_id, station_id, layer_id);
        const auto&& rsl_name_suffix = TString::Format("_GE%+.2d_L%d", region_id * (station_id * 10 + 1), layer_id);
        const auto&& rsl_title_suffix = TString::Format(" : GE%+.2d Layer %d", region_id * (station_id * 10 + 1), layer_id);
        const GEMDetId&& rsl_key = getReStLaKey(chamber->id());
        bookChamberOccupancy(ibooker, num_chambers_per_layer, rsl_key, rsl_name_suffix, rsl_title_suffix);
      }

      const auto& chambers = superchambers[0]->chambers();
      if (not checkRefs(chambers)) {
        edm::LogError(log_category_) << "failed to get a valid vector of GEMChamber ptrs" << std::endl;
        return;
      }

      for (const GEMEtaPartition* eta_partition : chambers[0]->etaPartitions()) {
        const int roll_id = eta_partition->id().roll();
        const GEMDetId&& rse_key = getReStEtKey(eta_partition->id());
        // FIXME waiting for yeckang:gemValidationTools PR
        // const TString&& rse_name_suffix = GEMUtils::getSuffixName(region_id, station_id) + TString::Format("_R%d", roll_id);
        // const TString&& rse_title_suffix = GEMUtils::getSuffixTitle(region_id, station_id) + TString::Format(" Roll %d", roll_id);
        const TString&& rse_name_suffix = TString::Format("_GE%+.2d_R%d", region_id * (station_id * 10 + 1), roll_id);
        const TString&& rse_title_suffix = TString::Format(" : GE%+.2d Roll %d", region_id * (station_id * 10 + 1), roll_id);
        bookResolution(ibooker, rse_key, rse_name_suffix, rse_title_suffix);
      }  // ieta
    }  // station
  }    // region
}

void GEMEfficiencyAnalyzer::bookDetectorOccupancy(DQMStore::IBooker& ibooker,
                                                  const GEMStation* station,
                                                  const GEMDetId& key,
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
                                          const GEMDetId& key,
                                          const TString& name_suffix,
                                          const TString& title_suffix) {
  ibooker.setCurrentFolder(folder_ + "/Efficiency");
  BookingHelper helper(ibooker, name_suffix, title_suffix);

  me_muon_pt_[key] = helper.book1D("muon_pt", title_, pt_binning_, "Muon p_{T} [GeV]");
  me_muon_eta_[key] = helper.book1D("muon_eta", title_, eta_nbins_, eta_low_, eta_up_, "Muon |#eta|");
  me_muon_phi_[key] = helper.book1D("muon_phi", title_, 108, -5., 355., "Muon #phi^{#circ}");

  me_muon_pt_matched_[key] = helper.book1D("muon_pt_matched", matched_title_, pt_binning_, "Muon p_{T} [GeV]");
  me_muon_eta_matched_[key] =
      helper.book1D("muon_eta_matched", matched_title_, eta_nbins_, eta_low_, eta_up_, "Muon |#eta|");
  me_muon_phi_matched_[key] =
      helper.book1D("muon_phi_matched", matched_title_, 108, -5., 355., "Muon #phi^{#circ}");

}

void GEMEfficiencyAnalyzer::bookChamberOccupancy(DQMStore::IBooker& ibooker,
                                                 const int num_chambers,
                                                 const GEMDetId& key,
                                                 const TString& name_suffix,
                                                 const TString& title_suffix) {
  ibooker.setCurrentFolder(folder_ + "/Efficiency");
  BookingHelper helper(ibooker, name_suffix, title_suffix);

  me_chamber_[key] = helper.book1D("muon_chamber", title_, num_chambers, 0.5, num_chambers + 0.5);
  me_chamber_matched_[key] = helper.book1D("muon_chamber_matched", title_, num_chambers, 0.5, num_chambers + 0.5);

}

void GEMEfficiencyAnalyzer::bookResolution(DQMStore::IBooker& ibooker,
                                           const GEMDetId& key,
                                           const TString& name_suffix,
                                           const TString& title_suffix) {
  ibooker.setCurrentFolder(folder_ + "/Resolution");
  BookingHelper helper(ibooker, name_suffix, title_suffix);

  me_residual_rdphi_[key] = helper.book1D("residual_rdphi", title_, 50, -5.0, 5.0, "R_{muon}#times(#phi_{muon}-#phi_{hit}) [cm #times rad.]");
  me_residual_y_[key] = helper.book1D("residual_y", title_, 60, -12.0, 12.0, "Residual in Local Y [cm]");

  // FIXME
  me_pull_rdphi_[key] = helper.book1D("pull_rdphi", title_, 60, -3.0, 3.0, "Pull R#Delta#phi");
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

  edm::ESHandle<GlobalTrackingGeometry> global_tracking_geometry;
  setup.get<GlobalTrackingGeometryRecord>().get(global_tracking_geometry);
  if (not global_tracking_geometry.isValid()) {
    edm::LogError(log_category_) << "GlobalTrackingGeometry is invalid" << std::endl;
    return;
  }

  edm::ESHandle<TransientTrackBuilder> transient_track_builder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", transient_track_builder);
  if (not transient_track_builder.isValid()) {
    edm::LogError(log_category_) << "TransientTrackRecord is invalid" << std::endl;
    return;
  }

  muon_service_->update(setup);
  edm::ESHandle<Propagator>&& propagator = muon_service_->propagator("SteppingHelixPropagatorAny");
  if (not propagator.isValid()) {
    edm::LogError(log_category_) << "Propagator is invalid" << std::endl;
    return;
  }

  const auto&& layer_vector = buildGEMLayers(gem);

  for (const reco::Muon& muon : *muon_view) {
    const reco::Track* track = getTrack(muon);
    if (track == nullptr) {
      edm::LogError(log_category_) << "failed to get muon track" << std::endl;
      continue;
    }

    const reco::TransientTrack&& transient_track = transient_track_builder->build(track);
    if (not transient_track.isValid()) {
      edm::LogError(log_category_) << "failed to build TransientTrack" << std::endl;
      continue;
    }

    for (const auto layer : layer_vector) {
      if (skipLayer(track, layer)) {
        // TODO LogInfo
        continue;
      }

      const auto&& [start_state, start_id] = getStartingState(transient_track, layer, global_tracking_geometry);
      if (not start_state.isValid()) {
        // TODO detail msg
        edm::LogInfo(log_category_) << "failed to get a starting state" << std::endl;
        continue;
      }

      // trajectory state on the destination surface
      const auto&& dest_state = propagator->propagate(start_state, *(layer.surface));
      if (not dest_state.isValid()) {
        auto msg = Form("failed to propagate a muon from DetId=(%d, %d) to GE%+.2d Layer %d",
                        static_cast<int>(start_id.det()),
                        start_id.subdetId(),
                        layer.region * (layer.station * 10 + 1),
                        layer.layer);
        edm::LogError(log_category_) << msg << std::endl;
        continue;
      }

      const GlobalPoint&& dest_global_pos = dest_state.globalPosition();
      if (not checkBounds(dest_global_pos, (*layer.surface))) {
        edm::LogInfo(log_category_) << "failed to pass checkBounds" << std::endl;
        continue;
      }

      const GEMEtaPartition* eta_partition = findEtaPartition(dest_global_pos, layer.chambers);
      if (eta_partition == nullptr) {
        edm::LogInfo(log_category_) << "failed to find an eta partition" << std::endl;
        continue;
      }

      const GEMDetId&& gem_id = eta_partition->id();

      const GEMDetId&& rs_key = getReStKey(gem_id);
      const GEMDetId&& rsl_key = getReStLaKey(gem_id);
      const GEMDetId&& rse_key = getReStEtKey(gem_id);

      const int chamber_bin = getDetOccXBin(gem_id, gem);

      const float muon_abs_eta = std::fabs(muon.eta());
      const float muon_phi_degree = toDegree(muon.phi());

      fillME(me_detector_, rs_key, chamber_bin, gem_id.roll());
      fillME(me_muon_pt_, rs_key, muon.pt());
      fillME(me_muon_eta_, rs_key, muon_abs_eta);
      fillME(me_muon_phi_, rs_key, muon_phi_degree);
      fillME(me_chamber_, rsl_key, gem_id.chamber());

      const GEMRecHit* matched_hit = findMatchedHit(
          dest_global_pos, rechit_collection->get(gem_id), eta_partition);

      if (matched_hit == nullptr) {
        continue;
      }

      fillME(me_detector_matched_, rs_key, chamber_bin, gem_id.roll());
      fillME(me_muon_pt_matched_, rs_key, muon.pt());
      fillME(me_muon_eta_matched_, rs_key, muon_abs_eta);
      fillME(me_muon_phi_matched_, rs_key, muon_phi_degree);
      fillME(me_chamber_matched_, rsl_key, gem_id.chamber());

      const LocalPoint&& dest_local_pos = eta_partition->toLocal(dest_global_pos);
      const LocalPoint&& hit_local_pos = matched_hit->localPosition();
      const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit_local_pos);

      // FIXME got the following error
      // const float dphi = reco::deltaPhi(dest_global_pos.phi(), hit_global_pos.phi());
      // /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_11_2_0_pre10/src/DataFormats/Math/interface/deltaPhi.h:19:17:
      // the type 'const Geom::Phi<float, Geom::MinusPiToPi>' of 'constexpr' variable 'o2pi' is not literal
      // TODO https://github.com/cms-sw/cmssw/issues/14862
      const float dphi = TVector2::Phi_mpi_pi(dest_global_pos.barePhi() - hit_global_pos.barePhi());
      const float residual_rdphi = dest_global_pos.perp() * dphi;
      const float residual_y = dest_local_pos.y() - hit_local_pos.y();

      const LocalError&& dest_local_err = dest_state.localError().positionError();
      const LocalError&& hit_local_err = matched_hit->localPositionError();

      const BoundPlane& bound_plane = eta_partition->surface();
      const GlobalError& dest_global_err = ErrorFrameTransformer().transform(dest_local_err, bound_plane);
      const GlobalError& hit_global_err = ErrorFrameTransformer().transform(hit_local_err, bound_plane);

      const float global_phi_err = std::hypot(dest_global_err.phierr(dest_global_pos), hit_global_err.phierr(hit_global_pos));

      // FIXME
      const float pull_rdphi = residual_rdphi / global_phi_err;
      const float pull_y = residual_y / std::sqrt(dest_local_err.yy() + hit_local_err.yy());

      fillME(me_residual_rdphi_, rse_key, residual_rdphi);
      fillME(me_residual_y_, rse_key, residual_y);

      fillME(me_pull_rdphi_, rse_key, pull_rdphi);
      fillME(me_pull_y_, rse_key, pull_y);
    }  // layer
  }    // Muon
}


std::vector<GEMEfficiencyAnalyzer::GEMLayerData>
GEMEfficiencyAnalyzer::buildGEMLayers(const edm::ESHandle<GEMGeometry>& gem) {
  std::vector<GEMLayerData> layer_vector;

  for (const GEMRegion* region : gem->regions()) {
    const int region_id = region->region();

    for (const GEMStation* station : region->stations()) {
      const int station_id = station->station();

      // layer_id - chambers
      std::map<int, std::vector<const GEMChamber*> > tmp_data;
      for (const GEMSuperChamber* super_chamber : station->superChambers()) {
        for (const GEMChamber* chamber : super_chamber->chambers()) {
          const int layer_id = chamber->id().layer();

          if (tmp_data.find(layer_id) == tmp_data.end())
            tmp_data.insert({layer_id, std::vector<const GEMChamber*>()});

          tmp_data[layer_id].push_back(chamber);
        } // GEMChamber
      } // GEMSuperChamber

      for (auto [layer_id, chamber_vector] : tmp_data) {
        // TODO checkRefs
        auto [rmin, rmax] = chamber_vector[0]->surface().rSpan();
        auto [zmin, zmax] = chamber_vector[0]->surface().zSpan();

        for (const GEMChamber* chamber : chamber_vector) {
          // the span of a bound surface in the global coordinates
          const auto [chamber_rmin, chamber_rmax] = chamber->surface().rSpan();
          const auto [chamber_zmin, chamber_zmax] = chamber->surface().zSpan();

          rmin = std::min(rmin, chamber_rmin);
          rmax = std::max(rmax, chamber_rmax);

          zmin = std::min(zmin, chamber_zmin);
          zmax = std::max(zmax, chamber_zmax);
        }

        // layer position and rotation
        const float layer_z = chamber_vector[0]->position().z();
        Surface::PositionType position(0.f, 0.f, layer_z);
        Surface::RotationType rotation;

        // the bounds from min and max R and Z in the local coordinates.
        // but disk seems to ignore z
        auto bounds = new SimpleDiskBounds(rmin, rmax, zmin - layer_z, zmax - layer_z);
        const Disk::DiskPointer&& layer = Disk::build(position, rotation, bounds);

        layer_vector.emplace_back(layer, chamber_vector, region_id, station_id, layer_id);

      } // layer
    } // GEMStation
  } // GEMRegion

  return layer_vector;
}


const reco::Track* GEMEfficiencyAnalyzer::getTrack(const reco::Muon& muon) {
  const reco::Track* track = nullptr;

  if (is_cosmics_) {
    if (muon.outerTrack().isNonnull())
      track = muon.outerTrack().get();

  } else {
    // beams, i.e. pp or heavy ions
    if (use_global_muon_ and muon.globalTrack().isNonnull())
      track = muon.globalTrack().get();

    else if ((not use_global_muon_) and muon.outerTrack().isNonnull())
      track = muon.outerTrack().get();
  }

  return track;
}


std::pair<TrajectoryStateOnSurface, DetId>
GEMEfficiencyAnalyzer::getStartingState(
    const reco::TransientTrack& transient_track,
    const GEMLayerData& layer,
    const edm::ESHandle<GlobalTrackingGeometry>& geometry) {

  TrajectoryStateOnSurface starting_state;
  DetId starting_id;

  if (use_global_muon_) {
    std::tie(starting_state, starting_id) = findStartingState(transient_track, layer, geometry);

  } else {
    // if outer track
    const reco::Track& track = transient_track.track();
    const bool is_insideout = isInsideOut(track);

    const DetId inner_id{(is_insideout ? track.outerDetId() : track.innerDetId())};
    if ((inner_id.det() == DetId::Detector::Muon) and (inner_id.subdetId() == MuonSubdetId::GEM)) {
      std::tie(starting_state, starting_id) = findStartingState(transient_track, layer, geometry);

    } else {
      starting_id = std::move(inner_id);
      if (is_insideout)
        starting_state = std::move(transient_track.outermostMeasurementState());
      else
        starting_state = std::move(transient_track.innermostMeasurementState());
    }
  }
 
  return std::make_pair(starting_state, starting_id);
}


std::pair<TrajectoryStateOnSurface, DetId>
GEMEfficiencyAnalyzer::findStartingState(
    const reco::TransientTrack& transient_track,
    const GEMLayerData& layer,
    const edm::ESHandle<GlobalTrackingGeometry>& geometry) {

  GlobalPoint starting_point;
  DetId starting_id;
  float min_distance = 1e12;
  bool found = false;
    
  for (auto rechit = transient_track.recHitsBegin(); rechit != transient_track.recHitsEnd(); rechit++) {
    const DetId&& det_id = (*rechit)->geographicalId();
    const int subdet_id = det_id.subdetId();

    // FIXME temporary
    if ((det_id.det() == DetId::Detector::Muon) and (subdet_id == MuonSubdetId::GEM))
      continue;

    const GeomDet* det = geometry->idToDet(det_id);
    const GlobalPoint&& global_position = det->toGlobal((*rechit)->localPosition());
    const float distance = std::abs(layer.surface->localZclamped(global_position));
    if (distance < min_distance) {
      found = true;
      min_distance = distance;
      starting_point = global_position;
      starting_id = det_id;
    }
  }

  TrajectoryStateOnSurface starting_state;
  if (found) {
    starting_state = std::move(transient_track.stateOnSurface(starting_point));
  }
  return std::make_pair(starting_state, starting_id);
}


bool GEMEfficiencyAnalyzer::skipLayer(const reco::Track* track,
                                      const GEMLayerData& layer) {
  const bool is_same_region = track->eta() * layer.region > 0;

  bool skip = false;
  if (is_cosmics_) {
    float p2_in = track->innerMomentum().mag2();
    float p2_out = track->outerMomentum().mag2();
    if (isInsideOut(*track))
      std::swap(p2_in, p2_out);
    const bool is_outgoing = p2_in > p2_out;

    skip = (is_outgoing xor is_same_region);

  } else {
    // beam scenario
    skip = not is_same_region;

  }

  return skip;
}


bool GEMEfficiencyAnalyzer::checkBounds(const GlobalPoint& global_point, const Plane& plane) {
  const LocalPoint&& local_point = plane.toLocal(global_point);
  const LocalPoint local_point_2d(local_point.x(), local_point.y(), 0.0f);
  return plane.bounds().inside(local_point_2d);
}


const GEMEtaPartition* GEMEfficiencyAnalyzer::findEtaPartition(
    const GlobalPoint& global_point,
    const std::vector<const GEMChamber*>& chamber_vector) {

  const GEMEtaPartition* bound = nullptr;
  for (const GEMChamber* chamber : chamber_vector) {
    if (not checkBounds(global_point, chamber->surface())) 
      continue;

    for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
      if (checkBounds(global_point, eta_partition->surface())) {
        bound = eta_partition;
        break;
      }
    } // GEMEtaPartition
  } // GEMChamber

  return bound;
}


const GEMRecHit* GEMEfficiencyAnalyzer::findMatchedHit(
    const GlobalPoint& dest_global_pos,
    const GEMRecHitCollection::range& range,
    const GEMEtaPartition* eta_partition) {

  const float dest_r = dest_global_pos.perp();
  const float dest_phi = dest_global_pos.barePhi();

  const GEMRecHit* closest_hit = nullptr;
  float min_rdphi{rdphi_cut_};

  for (auto hit = range.first; hit != range.second; ++hit) {
    const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit->localPosition());
    const float dphi = TVector2::Phi_mpi_pi(dest_phi - hit_global_pos.barePhi());
    const float rdphi = std::fabs(dest_r * dphi);
    if (rdphi <= min_rdphi) {
      min_rdphi = rdphi;
      closest_hit = &(*hit);
    }
  }

  return closest_hit;
}

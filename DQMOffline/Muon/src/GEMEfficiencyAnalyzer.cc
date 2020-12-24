#include "DQMOffline/Muon/interface/GEMEfficiencyAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

#include "TVector2.h"


GEMEfficiencyAnalyzer::GEMEfficiencyAnalyzer(const edm::ParameterSet& pset) : GEMOfflineDQMBase(pset) {
  rechit_token_ = consumes<GEMRecHitCollection>(pset.getParameter<edm::InputTag>("recHitTag"));
  muon_token_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"));

  auto muon_service_parameter = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  muon_service_ = new MuonServiceProxy(muon_service_parameter, consumesCollector());

  use_global_muon_ = pset.getUntrackedParameter<bool>("useGlobalMuon");

  rdphi_cut_ = static_cast<float>(pset.getParameter<double>("RDeltaPhiCut"));

  pt_binning_ = pset.getUntrackedParameter<std::vector<double> >("ptBinning");
  eta_nbins_ = pset.getUntrackedParameter<int>("etaNbins");
  eta_low_ = pset.getUntrackedParameter<double>("etaLow");
  eta_up_ = pset.getUntrackedParameter<double>("etaUp");

  folder_ = pset.getUntrackedParameter<std::string>("folder");

  // FIXME looks ugly
  title_ = (use_global_muon_ ? "Global Muon" : "Standalone Muon");
  matched_title_ = title_ + TString::Format(" (R#times|#phi_{Muon} - #phi_{Hit}| < %.1f)", rdphi_cut_);
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
  desc.add<double>("RDeltaPhiCut", 5.0); // TODO need to be tuned
  desc.addUntracked<std::vector<double> >("ptBinning", {20. ,30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 200.});
  desc.addUntracked<int>("etaNbins", 7);
  desc.addUntracked<double>("etaLow", 1.5);
  desc.addUntracked<double>("etaUp", 2.2);

  desc.addUntracked<bool>("useGlobalMuon", true); // FIXME looks ugly..
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
    const int region_number = region->region();
    const char* region_sign = region_number > 0 ? "+" : "-";

    for (const GEMStation* station : region->stations()) {
      const int station_number = station->station();

      const MEMapKey1 key1{region_number, station_number};
      const auto&& station_name_suffix = TString::Format("_ge%s%d1", region_sign, station_number);
      const auto&& station_title_suffix = TString::Format(" : GE%s%d/1", region_sign, station_number);
      bookDetectorOccupancy(ibooker, station, key1, station_name_suffix, station_title_suffix);
      bookOccupancy(ibooker, key1, station_name_suffix, station_title_suffix);

      const int num_etas = getNumEtaPartitions(station);
      for (int ieta = 1; ieta <= num_etas; ieta++) {
        const MEMapKey2 key2{region_number, station_number, ieta};
        const TString&& ieta_name_suffix = station_name_suffix + Form("_ieta%d", ieta);
        const TString&& ieta_title_suffix = station_title_suffix + Form(", i#eta = %d", ieta);
        bookResolution(ibooker, key2, ieta_name_suffix, ieta_title_suffix);
      }  // ieta
    }  // station
  }    // region
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
                                          const MEMapKey1& key,
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

void GEMEfficiencyAnalyzer::bookResolution(DQMStore::IBooker& ibooker,
                                           const MEMapKey2& key,
                                           const TString& name_suffix,
                                           const TString& title_suffix) {
  ibooker.setCurrentFolder(folder_ + "/Resolution");
  BookingHelper helper(ibooker, name_suffix, title_suffix);

  me_residual_rdphi_[key] = helper.book1D("residual_rdphi", title_, 50, -5.0, 5.0, "Residual R#Delta#phi [cm]");
  me_residual_y_[key] = helper.book1D("residual_y", title_, 60, -12.0, 12.0, "Residual in Local Y [cm]");

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

  const auto&& surface_vector = buildSurfaces(gem);

  for (const reco::Muon& muon : *muon_view) {
    const reco::Track* track = nullptr;

    if (use_global_muon_ and muon.globalTrack().isNonnull()) {
      track = muon.globalTrack().get();

    } else if ((not use_global_muon_) and muon.outerTrack().isNonnull()) {
      track = muon.outerTrack().get();
    }

    if (track == nullptr) {
      edm::LogError(log_category_) << "failed to get muon track" << std::endl;
      continue;
    }

    const reco::TransientTrack&& transient_track = transient_track_builder->build(track);
    if (not transient_track.isValid()) {
      edm::LogInfo(log_category_) << "failed to build TransientTrack" << std::endl;
      continue;
    }

    for (auto [surface, chamber_vector] : surface_vector) {
      // Skip propagation inn the opposite direction.
      if (muon.eta() * (*surface).eta() < 0)
        continue;

      const TrajectoryStateOnSurface&& tsos =
          propagator->propagate(transient_track.outermostMeasurementState(), *surface);
      if (not tsos.isValid()) {
        continue;
      }

      const LocalPoint&& tsos_local_pos = tsos.localPosition();
      const LocalPoint tsos_local_pos_2d(tsos_local_pos.x(), tsos_local_pos.y(), 0.0f);
      if (not (*surface).bounds().inside(tsos_local_pos_2d)) {
        continue;
      }

      const GlobalPoint&& tsos_global_pos = tsos.globalPosition();
      const GEMEtaPartition* eta_partition = findEtaPartition(tsos_global_pos, chamber_vector);

      const GEMDetId&& gem_id = eta_partition->id();

      const MEMapKey1 key1{gem_id.region(), gem_id.station()};
      const MEMapKey2 key2{gem_id.region(), gem_id.station(), gem_id.roll()};

      const int chamber_bin = getDetOccXBin(gem_id, gem);

      const float muon_abs_eta = std::fabs(muon.eta());
      const float muon_phi_degree = toDegree(muon.phi());

      fillME(me_detector_, key1, chamber_bin, gem_id.roll());
      fillME(me_muon_pt_, key1, muon.pt());
      fillME(me_muon_eta_, key1, muon_abs_eta);
      fillME(me_muon_phi_, key1, muon_phi_degree);

      const GEMRecHit* matched_hit = findMatchedHit(
          tsos_global_pos, rechit_collection->get(gem_id), eta_partition);

      if (matched_hit == nullptr) {
        continue;
      }

      fillME(me_detector_matched_, key1, chamber_bin, gem_id.roll());
      fillME(me_muon_pt_matched_, key1, muon.pt());
      fillME(me_muon_eta_matched_, key1, muon_abs_eta);
      fillME(me_muon_phi_matched_, key1, muon_phi_degree);

      const LocalPoint&& hit_local_pos = matched_hit->localPosition();
      const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit_local_pos);

      // FIXME got the following error
      // const float dphi = reco::deltaPhi(tsos_global_pos.phi(), hit_global_pos.phi());
      // /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_11_2_0_pre10/src/DataFormats/Math/interface/deltaPhi.h:19:17:
      // the type 'const Geom::Phi<float, Geom::MinusPiToPi>' of 'constexpr' variable 'o2pi' is not literal
      const float dphi = TVector2::Phi_mpi_pi(tsos_global_pos.barePhi() - hit_global_pos.barePhi());
      const float residual_rdphi = tsos_global_pos.perp() * dphi;
      const float residual_y = tsos_local_pos.y() - hit_local_pos.y();

      const LocalError&& tsos_local_err = tsos.localError().positionError();
      const LocalError&& hit_local_err = matched_hit->localPositionError();

      const BoundPlane& bound_plane = eta_partition->surface();
      const GlobalError& tsos_global_err = ErrorFrameTransformer().transform(tsos_local_err, bound_plane);
      const GlobalError& hit_global_err = ErrorFrameTransformer().transform(hit_local_err, bound_plane);

      const float global_phi_err = std::hypot(tsos_global_err.phierr(tsos_global_pos), hit_global_err.phierr(hit_global_pos));

      const float pull_rdphi = residual_rdphi / global_phi_err;
      const float pull_y = residual_y / std::sqrt(tsos_local_err.yy() + hit_local_err.yy());

      fillME(me_residual_rdphi_, key2, residual_rdphi);
      fillME(me_residual_y_, key2, residual_y);

      fillME(me_pull_rdphi_, key2, pull_rdphi);
      fillME(me_pull_y_, key2, pull_y);
    }  // GEMChamber
  }    // Muon
}

const GEMRecHit* GEMEfficiencyAnalyzer::findMatchedHit(
    const GlobalPoint& tsos_global_pos,
    const GEMRecHitCollection::range& range,
    const GEMEtaPartition* eta_partition) {

  const float tsos_r = tsos_global_pos.perp();
  const float tsos_phi = tsos_global_pos.phi();

  const GEMRecHit* closest_hit = nullptr;
  float min_rdphi{rdphi_cut_};

  for (auto hit = range.first; hit != range.second; ++hit) {
    const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit->localPosition());
    const float dphi = tsos_phi - hit_global_pos.phi();
    const float rdphi = std::fabs(tsos_r * dphi);
    if (rdphi <= min_rdphi) {
      min_rdphi = rdphi;
      closest_hit = &(*hit);
    }
  }

  return closest_hit;
}


std::vector<std::pair<ReferenceCountingPointer<Disk>, std::vector<const GEMChamber*> > >
GEMEfficiencyAnalyzer::buildSurfaces(const edm::ESHandle<GEMGeometry>& gem) {
  std::vector<std::pair<ReferenceCountingPointer<Disk>, std::vector<const GEMChamber*> > > disk_vector;

  for (const GEMRegion* region : gem->regions()) {
    for (const GEMStation* station : region->stations()) {

      // layer - chambers
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
        const GEMChamber* chamber = chamber_vector.front();

        Surface::PositionType position(0.f, 0.f, chamber->position().z());
        Surface::RotationType rotation;

        auto [rmin, rmax] = chamber->surface().rSpan();
        auto [zmin, zmax] = chamber->surface().zSpan();
        auto bounds = new SimpleDiskBounds(rmin, rmax, zmin, zmax);

        ReferenceCountingPointer<Disk> disk(new Disk(position, rotation, bounds));

        disk_vector.emplace_back(disk, chamber_vector);

      } // layer
    } // GEMStation
  } // GEMRegion

  return disk_vector;
}


const GEMEtaPartition* GEMEfficiencyAnalyzer::findEtaPartition(
    const GlobalPoint& global_point,
    const std::vector<const GEMChamber*>& chamber_vector) {

  const GEMEtaPartition* bound = nullptr;
  for (const GEMChamber* chamber : chamber_vector) {
    if (not checkBounds(global_point, chamber)) 
      continue;

    for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
      if (checkBounds(global_point, eta_partition)) {
        bound = eta_partition;
        break;
      }
    } // GEMEtaPartition
  } // GEMChamber

  return bound;
}


bool GEMEfficiencyAnalyzer::checkBounds(const GlobalPoint& global_point, const GeomDet* det) {
  const LocalPoint&& local_point = det->toLocal(global_point);
  const LocalPoint local_point_2d(local_point.x(), local_point.y(), 0.0f);
  return det->surface().bounds().inside(local_point_2d);
}

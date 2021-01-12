#include "DQMOffline/Muon/interface/GEMEfficiencyAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
// #include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Validation/MuonGEMHits/interface/GEMValidationUtils.h"
#include "Validation/MuonHits/interface/MuonHitHelper.h"


#include "TVector2.h"


GEMEfficiencyAnalyzer::GEMEfficiencyAnalyzer(const edm::ParameterSet& pset) : GEMOfflineDQMBase(pset) {
  rechit_token_ = consumes<GEMRecHitCollection>(pset.getParameter<edm::InputTag>("recHitTag"));
  muon_token_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"));

  const edm::ParameterSet&& muon_service_parameter = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  muon_service_ = new MuonServiceProxy(muon_service_parameter, consumesCollector());

  is_cosmics_ = pset.getUntrackedParameter<bool>("isCosmics");
  use_global_muon_ = pset.getUntrackedParameter<bool>("useGlobalMuon");
  use_fiducial_cut_ = pset.getUntrackedParameter<bool>("useFiducialCut");
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
  desc.addUntracked<bool>("useFiducialCut", false);
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

  ibooker.setCurrentFolder(folder_ + "/Debug");
  debug_me_residual_rdphi_ = ibooker.book1D("residual_rdphi", "", 100, 0.0f, 10.0f);

  bookEfficiencyMomentum(ibooker, gem);
  bookEfficiencyChamber(ibooker, gem);
  bookEfficiencyEtaPartition(ibooker, gem);
  bookResolution(ibooker, gem);
  debugBookEfficiencyStrip(ibooker, gem);
}


void GEMEfficiencyAnalyzer::bookEfficiencyMomentum(
    DQMStore::IBooker& ibooker,
    const edm::ESHandle<GEMGeometry>& gem) {
  ibooker.setCurrentFolder(folder_ + "/Efficiency");

  const std::string pt_x_title = "Muon p_{T} [GeV]";
  const std::string eta_x_title = "Muon |#eta|";
  const std::string phi_x_title = "Muon #phi [degree]";

  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();

    const GEMDetId&& key = getReStKey(region_id, station_id);
    // FIXME waiting for yeckang:gemValidationTools PR
    // const auto&& name_suffix = GEMUtils::getSuffixName(region_id, station_id);
    // const auto&& title_suffix = GEMUtils::getSuffixTitle(region_id, station_id);
    const TString&& name_suffix = TString::Format("_GE%+.2d", region_id * (station_id * 10 + 1));
    const TString&& title_suffix = TString::Format(" : GE%+.2d", region_id * (station_id * 10 + 1));

    const TString&& title = title_ + title_suffix;
    const TString&& matched_title = title_ + title_suffix;

    TH1F* h_muon_pt = new TH1F("muon_pt" + name_suffix, title, pt_binning_.size() -1, &pt_binning_[0]);
    me_muon_pt_[key] = ibooker.book1D(h_muon_pt->GetName(), h_muon_pt);

    TH1F* h_muon_pt_matched = new TH1F("muon_pt_matched" + name_suffix, matched_title, pt_binning_.size() -1, &pt_binning_[0]);
    me_muon_pt_matched_[key] = ibooker.book1D(h_muon_pt_matched->GetName(), h_muon_pt_matched);

    me_muon_eta_[key] = ibooker.book1D("muon_eta" + name_suffix, title, eta_nbins_, eta_low_, eta_up_);
    me_muon_eta_matched_[key] = ibooker.book1D("muon_eta_matched" + name_suffix, matched_title, eta_nbins_, eta_low_, eta_up_);

    me_muon_phi_[key] = ibooker.book1D("muon_phi" + name_suffix, title, 108, -5., 355.);
    me_muon_phi_matched_[key] = ibooker.book1D("muon_phi_matched" + name_suffix, matched_title, 108, -5., 355.);

    me_muon_pt_[key]->setAxisTitle(pt_x_title);
    me_muon_pt_matched_[key]->setAxisTitle(pt_x_title);

    me_muon_eta_[key]->setAxisTitle(eta_x_title);
    me_muon_eta_matched_[key]->setAxisTitle(eta_x_title);

    me_muon_phi_[key]->setAxisTitle(phi_x_title);
    me_muon_phi_matched_[key]->setAxisTitle(phi_x_title);

  }  // station
}


void GEMEfficiencyAnalyzer::bookEfficiencyChamber(
    DQMStore::IBooker& ibooker,
    const edm::ESHandle<GEMGeometry>& gem) {
  ibooker.setCurrentFolder(folder_ + "/Efficiency");

  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();

    const std::vector<const GEMSuperChamber*>&& superchambers = station->superChambers();
    if (not checkRefs(superchambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
      return;
    }

    // Number of GEMChambers per each layer
    // FIXME
    const int num_chambers = superchambers.size();
    for (const GEMChamber* chamber : superchambers[0]->chambers()) {
      const int layer_id = chamber->id().layer();

      // FIXME waiting for yeckang:gemValidationTools PR
      // const auto&& name_suffix = GEMUtils::getSuffixName(region_id, station_id, layer_id);
      // const auto&& title_suffix = GEMUtils::getSuffixTitle(region_id, station_id, layer_id);
      const TString&& name_suffix = TString::Format("_GE%+.2d_L%d", region_id * (station_id * 10 + 1), layer_id);
      const TString&& title_suffix = TString::Format(" : GE%+.2d Layer %d", region_id * (station_id * 10 + 1), layer_id);
      const GEMDetId&& key = getReStLaKey(chamber->id());

      me_chamber_[key] = ibooker.book1D("muon_chamber" + name_suffix, title_ + title_suffix, num_chambers, 0.5, num_chambers + 0.5);
      me_chamber_matched_[key] = ibooker.book1D("muon_chamber_matched" + name_suffix, matched_title_ + title_suffix, num_chambers, 0.5, num_chambers + 0.5);

      me_chamber_[key]->setAxisTitle("Chamber");
      me_chamber_matched_[key]->setAxisTitle("Chamber");
    } // layer
  } // station
}


void GEMEfficiencyAnalyzer::bookEfficiencyEtaPartition(
    DQMStore::IBooker& ibooker,
    const edm::ESHandle<GEMGeometry>& gem) {
  ibooker.setCurrentFolder(folder_ + "/Efficiency");

  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();

    const GEMDetId&& key = getReStKey(region_id, station_id);
    // FIXME waiting for yeckang:gemValidationTools PR
    // const auto&& name_suffix = GEMUtils::getSuffixName(region_id, station_id);
    // const auto&& title_suffix = GEMUtils::getSuffixTitle(region_id, station_id);
    const TString&& name_suffix = TString::Format("_GE%+.2d", region_id * (station_id * 10 + 1));
    const TString&& title_suffix = TString::Format(" : GE%+.2d", region_id * (station_id * 10 + 1));

    const std::vector<const GEMSuperChamber*>&& superchambers = station->superChambers();
    if (not checkRefs(superchambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
      return;
    }

    // the number of GEMChambers per GEMStation
    const int num_ch = superchambers.size() * superchambers.front()->nChambers();
    // the number of eta partitions per GEMChamber
    const int num_etas = getNumEtaPartitions(station);

    me_detector_[key] = ibooker.book2D("detector" + name_suffix, title_ + title_suffix, num_ch, 0.5, num_ch + 0.5, num_etas, 0.5, num_etas + 0.5);
    me_detector_matched_[key] =
        ibooker.book2D("detector_matched" + name_suffix, matched_title_ + title_suffix, num_ch, 0.5, num_ch + 0.5, num_etas, 0.5, num_etas + 0.5);

    setDetLabelsEta(me_detector_[key], station);
    setDetLabelsEta(me_detector_matched_[key], station);

  }  // station
}



void GEMEfficiencyAnalyzer::bookResolution(
    DQMStore::IBooker & ibooker,
    const edm::ESHandle<GEMGeometry>& gem) {

  // const std::string rdphi_x_title = "R_{muon}#times(#phi_{muon}-#phi_{hit}) [cm #times rad.]";
  // const std::string residual_y_x_title = "Residual in Local Y [cm]";
  // const std::string pull_phi_x_title = "Pull in Global #phi";
  // const std::string pull_y_x_title = "Pull in Local Y";

  ibooker.setCurrentFolder(folder_ + "/Resolution");
  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();

    // FIXME
    const std::vector<const GEMSuperChamber*>&& superchambers = station->superChambers();
    if (not checkRefs(superchambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
      return;
    }

    const std::vector<const GEMChamber*>& chambers = superchambers[0]->chambers();
    if (not checkRefs(chambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMChamber ptrs" << std::endl;
      return;
    }

    for (const GEMEtaPartition* eta_partition : chambers[0]->etaPartitions()) {
      const int roll_id = eta_partition->id().roll();
      const GEMDetId&& key = getReStEtKey(eta_partition->id());
      // FIXME waiting for yeckang:gemValidationTools PR
      // const TString&& name_suffix = GEMUtils::getSuffixName(region_id, station_id) + TString::Format("_R%d", roll_id);
      // const TString&& title_suffix = GEMUtils::getSuffixTitle(region_id, station_id) + TString::Format(" Roll %d", roll_id);
      const TString&& name_suffix = TString::Format("_GE%+.2d_R%d", region_id * (station_id * 10 + 1), roll_id);
      const TString&& title = matched_title_ + TString::Format(" : GE%+.2d Roll %d", region_id * (station_id * 10 + 1), roll_id);

      me_residual_rdphi_[key] = ibooker.book1D("residual_rdphi" + name_suffix, title, 50, -rdphi_cut_, rdphi_cut_);
      me_residual_rdphi_[key]->setAxisTitle("R_{muon}#times(#phi_{muon}-#phi_{hit}) [cm]");

      me_residual_y_[key] = ibooker.book1D("residual_y" + name_suffix, title, 60, -12.0, 12.0);
      me_residual_y_[key]->setAxisTitle("Residual in Local Y [cm]");

      me_pull_phi_[key] = ibooker.book1D("pull_phi" + name_suffix, title, 60, -3.0, 3.0);
      me_pull_phi_[key]->setAxisTitle("Pull in Global #phi");

      me_pull_y_[key] = ibooker.book1D("pull_y" + name_suffix, title, 60, -3.0, 3.0);
      me_pull_y_[key]->setAxisTitle("Pull in Local Y");

    }  // ieta
  }  // station
}


void GEMEfficiencyAnalyzer::debugBookEfficiencyStrip(
    DQMStore::IBooker & ibooker,
    const edm::ESHandle<GEMGeometry>& gem) {

  ibooker.setCurrentFolder(folder_ + "/Efficiency");
  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();

    // FIXME
    const std::vector<const GEMSuperChamber*>&& superchambers = station->superChambers();
    if (not checkRefs(superchambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
      return;
    }

    const std::vector<const GEMChamber*>& chambers = superchambers[0]->chambers();
    if (not checkRefs(chambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMChamber ptrs" << std::endl;
      return;
    }

    for (const GEMEtaPartition* eta_partition : chambers[0]->etaPartitions()) {
      const int roll_id = eta_partition->id().roll();
      const GEMDetId&& key = getReStEtKey(eta_partition->id());
      const TString&& name_suffix = TString::Format("_GE%+.2d_R%d", region_id * (station_id * 10 + 1), roll_id);
      const TString&& title_suffix = TString::Format(" : GE%+.2d Roll %d", region_id * (station_id * 10 + 1), roll_id);

      const int nstrips = eta_partition->nstrips();

      me_strip_[key] = ibooker.book1D("strip" + name_suffix, title_ + title_suffix, nstrips, 0, nstrips);
      me_strip_matched_[key] = ibooker.book1D("strip_matched" + name_suffix, matched_title_ + title_suffix, nstrips, 0, nstrips);

      me_strip_[key]->setAxisTitle("strip");
      me_strip_matched_[key]->setAxisTitle("strip");

    }  // ieta
  }  // station
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

  const std::vector<GEMLayerData>&& layer_vector = buildGEMLayers(gem);

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

    for (const GEMLayerData& layer : layer_vector) {
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
      const TrajectoryStateOnSurface&& dest_state = propagator->propagate(start_state, *(layer.surface));
      if (not dest_state.isValid()) {
        const char* msg = Form("failed to propagate a muon from DetId=(%d, %d) to GE%+.2d Layer %d",
                        static_cast<int>(start_id.det()),
                        start_id.subdetId(),
                        layer.region * (layer.station * 10 + 1),
                        layer.layer);
        edm::LogInfo(log_category_) << msg << std::endl;
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

      const LocalPoint&& dest_local_pos = eta_partition->toLocal(dest_global_pos);
      const float strip = eta_partition->strip(eta_partition->toLocal(dest_global_pos));

      // make it tunable
      if (use_fiducial_cut_ and ((strip < 20) or (strip > 364))) {

        edm::LogInfo(log_category_) << "fiducial" << std::endl;
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

      fillME(me_strip_, rse_key, strip);

      // const GEMRecHit* matched_hit = findMatchedHit(
      //     dest_global_pos, rechit_collection->get(gem_id), eta_partition);
      const auto&& [matched_hit, debug_rdphi] = findMatchedHit(
          dest_global_pos, rechit_collection->get(gem_id), eta_partition);

      if (matched_hit == nullptr) {
        continue;
      }

      debug_me_residual_rdphi_->Fill(debug_rdphi);

      fillME(me_detector_matched_, rs_key, chamber_bin, gem_id.roll());
      fillME(me_muon_pt_matched_, rs_key, muon.pt());
      fillME(me_muon_eta_matched_, rs_key, muon_abs_eta);
      fillME(me_muon_phi_matched_, rs_key, muon_phi_degree);
      fillME(me_chamber_matched_, rsl_key, gem_id.chamber());

      fillME(me_strip_matched_, rse_key, strip);

      const LocalPoint&& hit_local_pos = matched_hit->localPosition();
      const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit_local_pos);

      const float dphi = TVector2::Phi_mpi_pi(dest_local_pos.barePhi() - hit_local_pos.barePhi());
      // const float dphi = dest_local_pos.phi() - hit_local_pos.phi();
      const float residual_rdphi = dest_global_pos.perp() * dphi;
      const float residual_y = dest_local_pos.y() - hit_local_pos.y();

      const LocalError&& dest_local_err = dest_state.localError().positionError();
      const LocalError&& hit_local_err = matched_hit->localPositionError();

      const BoundPlane& bound_plane = eta_partition->surface();
      const GlobalError& dest_global_err = ErrorFrameTransformer().transform(dest_local_err, bound_plane);
      const GlobalError& hit_global_err = ErrorFrameTransformer().transform(hit_local_err, bound_plane);

      const float global_phi_err = std::sqrt(dest_global_err.phierr(dest_global_pos) + hit_global_err.phierr(hit_global_pos));

      // FIXME
      const float pull_phi = dphi / global_phi_err;
      const float pull_y = residual_y / std::sqrt(dest_local_err.yy() + hit_local_err.yy());

      fillME(me_residual_rdphi_, rse_key, residual_rdphi);
      fillME(me_residual_y_, rse_key, residual_y);

      fillME(me_pull_phi_, rse_key, pull_phi);
      fillME(me_pull_y_, rse_key, pull_y);
    }  // layer
  }    // Muon
}


std::vector<GEMEfficiencyAnalyzer::GEMLayerData>
GEMEfficiencyAnalyzer::buildGEMLayers(const edm::ESHandle<GEMGeometry>& gem) {
  std::vector<GEMLayerData> layer_vector;

  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
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

      // zmin -= layer_z;
      // zmax -= layer_z;

      // the bounds from min and max R and Z in the local coordinates.
      // but disk seems to ignore z
      SimpleDiskBounds* bounds = new SimpleDiskBounds(rmin, rmax, zmin - layer_z, zmax - layer_z);
      const Disk::DiskPointer&& layer = Disk::build(position, rotation, bounds);

      layer_vector.emplace_back(layer, chamber_vector, region_id, station_id, layer_id);

    } // layer
  } // GEMStation

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
    // if ((inner_id.det() == DetId::Detector::Muon) and (inner_id.subdetId() == MuonSubdetId::GEM)) {
    if (MuonHitHelper::isGEM(inner_id.rawId())) {
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
    // const int subdet_id = det_id.subdetId();

    // FIXME temporary
    // if ((det_id.det() == DetId::Detector::Muon) and (subdet_id == MuonSubdetId::GEM))
    if (MuonHitHelper::isGEM(det_id.rawId()))
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


/*
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
*/


std::pair<const GEMRecHit*, float> GEMEfficiencyAnalyzer::findMatchedHit(
    const GlobalPoint& dest_global_pos,
    const GEMRecHitCollection::range& range,
    const GEMEtaPartition* eta_partition) {

  const float dest_r = dest_global_pos.perp();

  const LocalPoint&& dest_local_pos = eta_partition->toLocal(dest_global_pos);
  const float dest_local_phi = dest_local_pos.barePhi();

  const GEMRecHit* closest_hit = nullptr;
  float min_rdphi = 1e6;

  for (auto hit = range.first; hit != range.second; ++hit) {
    const float dphi = TVector2::Phi_mpi_pi(dest_local_phi - hit->localPosition().barePhi());
    const float rdphi = std::fabs(dest_r * dphi);
    if (rdphi <= min_rdphi) {
      min_rdphi = rdphi;
      closest_hit = &(*hit);
    }
  }

  return std::make_pair(closest_hit, min_rdphi);
}

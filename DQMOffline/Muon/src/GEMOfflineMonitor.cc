#include "DQMOffline/Muon/interface/GEMOfflineMonitor.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Validation/MuonGEMHits/interface/GEMValidationUtils.h"


GEMOfflineMonitor::GEMOfflineMonitor(const edm::ParameterSet& pset) : GEMOfflineDQMBase(pset) {
  digi_token_ = consumes<GEMDigiCollection>(pset.getParameter<edm::InputTag>("digiTag"));
  rechit_token_ = consumes<GEMRecHitCollection>(pset.getParameter<edm::InputTag>("recHitTag"));
  lumi_scalers_token_ = consumes<LumiScalersCollection>(pset.getParameter<edm::InputTag>("lumiScalersTag"));
  vertex_token_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertexTag"));
  do_digi_occupancy_ = pset.getUntrackedParameter<bool>("doDigiOccupancy");
  do_hit_occupancy_ = pset.getUntrackedParameter<bool>("doHitOccupancy");
  do_hit_rate_ = pset.getUntrackedParameter<bool>("doHitRate");

  uninitialized_area_ = true;
}

GEMOfflineMonitor::~GEMOfflineMonitor() {}

void GEMOfflineMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("digiTag", edm::InputTag("muonGEMDigis"));
  desc.add<edm::InputTag>("recHitTag", edm::InputTag("gemRecHits"));
  // TODO check tag
  desc.add<edm::InputTag>("lumiScalersTag", edm::InputTag("scalersRawToDigi"));
  // TODO check tag
  desc.add<edm::InputTag>("vertexTag", edm::InputTag("offlinePrimaryVertices"));
  desc.addUntracked<std::string>("logCategory", "GEMOfflineMonitor");
  desc.addUntracked<bool>("doDigiOccupancy", true);
  desc.addUntracked<bool>("doHitOccupancy", true);
  desc.addUntracked<bool>("doHitRate", true);
  descriptions.add("gemOfflineMonitorDefault", desc);
}


void GEMOfflineMonitor::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& run, edm::EventSetup const& setup) {
  edm::ESHandle<GEMGeometry> gem;
  setup.get<MuonGeometryRecord>().get(gem);
  if (not gem.isValid()) {
    edm::LogError(log_category_) << "GEMGeometry is invalid" << std::endl;
    return;
  }

  if (do_digi_occupancy_)
    bookDigiOccupancy(ibooker, gem);

  if (do_hit_occupancy_)
    bookHitOccupancy(ibooker, gem);

  if (do_hit_rate_)
    bookHitRate(ibooker, gem);
}


void GEMOfflineMonitor::bookDigiOccupancy(DQMStore::IBooker& ibooker,
                                         const edm::ESHandle<GEMGeometry>& gem) {
  ibooker.setCurrentFolder("GEM/GEMOfflineMonitor/DigiOccupancy");

  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();
    const GEMDetId&& key = getReStKey(region_id, station_id);
    // FIXME waiting for yeckang:gemValidationTools PR
    // const auto&& rs_name = GEMUtils::getSuffixName(region_id, station_id);
    // const auto&& rs_title = GEMUtils::getSuffixTitle(region_id, station_id);
    const auto&& name_suffix = TString::Format("_GE%+.2d", region_id * (station_id * 10 + 1));
    const auto&& title_suffix = TString::Format(" : GE%+.2d", region_id * (station_id * 10 + 1));

    const auto&& superchambers = station->superChambers();
    if (not checkRefs(superchambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
      return;
    }

    // per station
    const int num_superchambers = superchambers.size();
    const int num_chambers = num_superchambers * superchambers.front()->nChambers();
    // the numer of VFATs per GEMEtaPartition
    const int max_vfat = getMaxVFAT(station->station());
    // the number of eta partitions per GEMChamber
    const int num_etas = getNumEtaPartitions(station);
    // the number of VFATs per GEMChamber
    const int num_vfat = num_etas * max_vfat;

    me_digi_det_[key] = 
        ibooker.book2D("digi_det" + name_suffix, "Digi Occupancy" + title_suffix, 
                       num_chambers, 0.5, num_chambers + 0.5, num_vfat, 0.5, num_vfat + 0.5);
    setDetLabelsVFAT(me_digi_det_[key], station);
  }  // station
}

void GEMOfflineMonitor::bookHitOccupancy(DQMStore::IBooker& ibooker,
                                         const edm::ESHandle<GEMGeometry>& gem) {
  ibooker.setCurrentFolder("GEM/GEMOfflineMonitor/HitOccupancy");

  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();

    const GEMDetId&& key = getReStKey(region_id, station_id);
    // FIXME waiting for yeckang:gemValidationTools PR
    // const auto&& rs_name = GEMUtils::getSuffixName(region_id, station_id);
    // const auto&& rs_title = GEMUtils::getSuffixTitle(region_id, station_id);
    const auto&& name_suffix = TString::Format("_GE%+.2d", region_id * (station_id * 10 + 1));
    const auto&& title_suffix = TString::Format(" : GE%+.2d", region_id * (station_id * 10 + 1));

    const auto&& superchambers = station->superChambers();
    if (not checkRefs(superchambers)) {
      edm::LogError(log_category_) << "failed to get a valid vector of GEMSuperChamber ptrs" << std::endl;
      return;
    }

    // per station
    const int num_superchambers = superchambers.size();
    const int num_chambers = num_superchambers * superchambers.front()->nChambers();
    // the number of eta partitions per GEMChamber
    const int num_etas = getNumEtaPartitions(station);

    me_hit_det_[key] =
        ibooker.book2D("hit_det" + name_suffix, "Hit Occupancy" + title_suffix, num_chambers, 0.5, num_chambers + 0.5, num_etas, 0.5, num_etas + 0.5);
    setDetLabelsEta(me_hit_det_[key], station);
  }  // station
}


void GEMOfflineMonitor::bookHitRate(DQMStore::IBooker& ibooker,
                                    const edm::ESHandle<GEMGeometry>& gem) {
  ibooker.setCurrentFolder("GEM/GEMOfflineMonitor/HitRate/Source");

  me_hit_rate_num_events_ = ibooker.bookInt("num_events");

  for (const GEMStation* station : gem->stations()) {
    // TODO if GE11

    // FIXME dirty...
    std::map<std::string, bool> found = {{"Odd", false}, {"Even", false}};
    for (const GEMSuperChamber* super_chamber : station->superChambers()) {
      if (found["Odd"] and found["Even"])
        break;
      const std::string chamber_parity = super_chamber->id().chamber() % 2 == 0 ? "Even" : "Odd";
      if (found[chamber_parity])
        continue;
      found[chamber_parity] = true;

      for (const GEMChamber* chamber : super_chamber->chambers()) {
        for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
          const GEMDetId&& key = getKey(eta_partition->id());

          // FIXME waiting for yeckang:gemValidationTools PR
          // const auto&& name_suffix = GEMUtils::getSuffixName(region_id, station_id, layer_id, roll_id) + "_" + chamber_parity_name;
          // const auto&& title_suffix = GEMUtils::getSuffixTitle(region_id, station_id, layer_id, roll_id) + " " + chamber_parity_name;
          const TString&& name_suffix = TString::Format(
              "_GE%+.2d_L%d_R%d_%s",
              key.region() * (key.station() * 10 + 1),
              key.layer(), key.roll(), chamber_parity.c_str());
          const TString&& title_suffix = TString::Format(
              " : GE%+.2d Layer %d Roll %d %s",
              key.region() * (key.station() * 10 + 1),
              key.layer(), key.roll(), chamber_parity.c_str());

          const TString&& source_name = "source" + name_suffix;
          const TString&& source_title = "Hit Multiplicity per Vertex Multiplicity for Background Hit Rate" + title_suffix;

          me_hit_rate_source_[key] = ibooker.book1D(source_name, source_title, 201, -0.5, 200.5);
          me_hit_rate_source_[key]->setAxisTitle("Number of Vertices", 1);
          me_hit_rate_source_[key]->setAxisTitle("Number of Hits", 2);

          const TString&& area_name = "area" + name_suffix;
          me_hit_rate_area_[key] = ibooker.bookFloat("area" + name_suffix);
        } // GEMEtaPartition (roll)
      } // GEMChamber (layer)
    } // GEMSuperChamber (chamber)
  }  // station
}


void GEMOfflineMonitor::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<GEMDigiCollection> digi_collection;
  if (do_digi_occupancy_) {
    event.getByToken(digi_token_, digi_collection);
    if (not digi_collection.isValid()) {
      edm::LogError(log_category_) << "GEMDigiCollection is invalid!" << std::endl;
      return;
    }
  }

  edm::Handle<GEMRecHitCollection> rechit_collection;
  if (do_hit_occupancy_ or do_hit_rate_) {
    event.getByToken(rechit_token_, rechit_collection);
    if (not rechit_collection.isValid()) {
      edm::LogError(log_category_) << "GEMRecHitCollection is invalid" << std::endl;
      return;
    }
  }

  edm::Handle<LumiScalersCollection> lumi_scalers_collection;
  event.getByToken(lumi_scalers_token_, lumi_scalers_collection);
  if (not lumi_scalers_collection.isValid()) {
    edm::LogError(log_category_) << "LumiScalersCollection is invalid" << std::endl;
  }

  edm::Handle<reco::VertexCollection> vertex_collection;
  if (do_hit_rate_) {
    event.getByToken(vertex_token_, vertex_collection);
    if (not vertex_collection.isValid()) {
      edm::LogError(log_category_) << "VertexCollection is invalid" << std::endl;
    }
  }

  edm::ESHandle<GEMGeometry> gem;
  setup.get<MuonGeometryRecord>().get(gem);
  if (not gem.isValid()) {
    edm::LogError(log_category_) << "GEMGeometry is invalid" << std::endl;
    return;
  }

  if (do_digi_occupancy_)
    doDigiOccupancy(gem, digi_collection);

  if (do_hit_occupancy_)
    doHitOccupancy(gem, rechit_collection);

  if (do_hit_rate_)
    doHitRate(gem, rechit_collection, vertex_collection);
}


void GEMOfflineMonitor::doDigiOccupancy(
    const edm::ESHandle<GEMGeometry>& gem,
    const edm::Handle<GEMDigiCollection>& digi_collection) {
  for (auto range_iter = digi_collection->begin(); range_iter != digi_collection->end(); range_iter++) {
    const GEMDetId& gem_id = (*range_iter).first;
    const GEMDigiCollection::Range& range = (*range_iter).second;

    const GEMDetId&& rs_key = getReStKey(gem_id);
    for (auto digi = range.first; digi != range.second; ++digi) {
      const int chamber_bin = getDetOccXBin(gem_id, gem);
      const int vfat_number = getVFATNumberByStrip(gem_id.station(), gem_id.roll(), digi->strip());

      fillME(me_digi_det_, rs_key, chamber_bin, vfat_number);
    } // digi
  } // range
}


void GEMOfflineMonitor::doHitOccupancy(
    const edm::ESHandle<GEMGeometry>& gem,
    const edm::Handle<GEMRecHitCollection>& rechit_collection) {
  for (auto hit = rechit_collection->begin(); hit != rechit_collection->end(); hit++) {
    const GEMDetId&& gem_id = hit->gemId();
    const GEMDetId&& rs_key = getReStKey(gem_id);
    const int chamber_bin = getDetOccXBin(gem_id, gem);
    fillME(me_hit_det_, rs_key, chamber_bin, gem_id.roll());
  }
}


void GEMOfflineMonitor::doHitRate(
    const edm::ESHandle<GEMGeometry>& gem,
    const edm::Handle<GEMRecHitCollection>& rechit_collection,
    const edm::Handle<reco::VertexCollection>& vertex_collection) {

  int num_vtx = 0;
  for (const reco::Vertex& vertex : *vertex_collection) {
    if (not vertex.isValid()) {
      continue;
    }
    num_vtx++;
  }

  for (auto hit = rechit_collection->begin(); hit != rechit_collection->end(); hit++) {
    const GEMDetId&& key = getKey(hit->gemId());
    fillME(me_hit_rate_source_, key, num_vtx);
  }

  // TODO initialization...?
  me_hit_rate_num_events_->Fill(me_hit_rate_num_events_->getIntValue() + 1);

  if (uninitialized_area_) {
    uninitialized_area_ = false;

    for (const GEMStation* station : gem->stations()) {
      // 0: even, 1: odd
      std::map<int, bool> found = {{0, false}, {1, false}};
      for (const GEMSuperChamber* super_chamber : station->superChambers()) {
        if (found[0] and found[1]) {
          break;
        }

        const int chamber_pairty = super_chamber->id().chamber() % 2;
        if (found[chamber_pairty]) {
          continue;
        }
        found[chamber_pairty] = true;

        for (const GEMChamber* chamber : super_chamber->chambers()) {
          for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
            const GEMDetId&& key = getKey(eta_partition->id());

            const StripTopology& strip_topology = eta_partition->specificTopology();
            const float area = strip_topology.stripLength() * strip_topology.pitch() * eta_partition->nstrips();

            fillME(me_hit_rate_area_, key, area);
          } // GEMEtaPartition
        } // GEMChamber
      } // GEMSuperChamber
    }  // GEMStation
  } // if (uninitialized_area_)
}

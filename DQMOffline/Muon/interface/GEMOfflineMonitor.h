#ifndef DQMOffline_Muon_GEMOfflineMonitor_h
#define DQMOffline_Muon_GEMOfflineMonitor_h

#include "DQMOffline/Muon/interface/GEMOfflineDQMBase.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


class GEMOfflineMonitor : public GEMOfflineDQMBase {
public:
  explicit GEMOfflineMonitor(const edm::ParameterSet &);
  ~GEMOfflineMonitor() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

protected:
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event &event, const edm::EventSetup &eventSetup) override;

private:
  void bookDigiOccupancy(DQMStore::IBooker &, const edm::ESHandle<GEMGeometry>&);
  void bookHitOccupancy(DQMStore::IBooker &, const edm::ESHandle<GEMGeometry>&);
  void bookHitRate(DQMStore::IBooker &, const edm::ESHandle<GEMGeometry>&);

  void doDigiOccupancy(const edm::ESHandle<GEMGeometry>&,
              const edm::Handle<GEMDigiCollection>&);
  void doHitOccupancy(const edm::ESHandle<GEMGeometry>&,
                const edm::Handle<GEMRecHitCollection>&);
  void doHitRate(const edm::ESHandle<GEMGeometry>&,
                 const edm::Handle<GEMRecHitCollection>&,
                 const edm::Handle<reco::VertexCollection>&);

  edm::EDGetTokenT<GEMDigiCollection> digi_token_;
  edm::EDGetTokenT<GEMRecHitCollection> rechit_token_;
  edm::EDGetTokenT<LumiScalersCollection> lumi_scalers_token_;
  edm::EDGetTokenT<reco::VertexCollection> vertex_token_;

  std::string log_category_;

  bool do_digi_occupancy_;
  bool do_hit_occupancy_;
  bool do_hit_rate_;

  bool uninitialized_area_;

  MEMap me_digi_det_; // TH2F, region-station
  MEMap me_hit_det_; // TH2F, region-station

  MEMap me_hit_rate_source_; // TH1F, (region, station, chamber parity, ieta)
  MEMap me_hit_rate_area_;
  MonitorElement* me_hit_rate_num_events_; // Int
};

#endif  // DQMOffline_Muon_GEMOfflineMonitor_h

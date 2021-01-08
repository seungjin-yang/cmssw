#ifndef DQMOffline_Muon_GEMOfflineHarvester_h
#define DQMOffline_Muon_GEMOfflineHarvester_h

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include <vector>
#include <string>

class GEMOfflineHarvester : public DQMEDHarvester {
public:
  GEMOfflineHarvester(const edm::ParameterSet&);
  ~GEMOfflineHarvester() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &);
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

private:
  void doHitRate(DQMStore::IBooker&, DQMStore::IGetter&);


  std::string folder_;
  std::string log_category_;
  double timing_window_;
};



#endif  // DQMOffline_Muon_GEMOfflineHarvester_h

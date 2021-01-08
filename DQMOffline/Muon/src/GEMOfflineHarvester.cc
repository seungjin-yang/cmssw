#include "DQMOffline/Muon/interface/GEMOfflineHarvester.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

GEMOfflineHarvester::GEMOfflineHarvester(const edm::ParameterSet& pset) {
  folder_ = pset.getUntrackedParameter<std::string>("folder");
  log_category_ = pset.getUntrackedParameter<std::string>("logCategory");
  timing_window_ = pset.getParameter<double>("timingWindow");
}

GEMOfflineHarvester::~GEMOfflineHarvester() {}


void GEMOfflineHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>("folder", "GEM/GEMOfflineMonitor");
  desc.addUntracked<std::string>("logCategory", "GEMOfflineHarvester");
  desc.add<double>("timingWindow", 1E-7); // 4*25 ns
  descriptions.add("gemOfflineHarvesterDefault", desc);
}

void GEMOfflineHarvester::doHitRate(DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter) {
  const std::string hit_rate_folder = folder_ + "/HitRate/";
  const std::string source_folder = hit_rate_folder + "/Source/";
  igetter.setCurrentFolder(source_folder);

  const std::string source_prefix = "source_";
  const std::string area_prefix = "area_";
  const std::string num_events_name = "num_events";

  MonitorElement* me_num_events = igetter.get(source_folder + num_events_name);
  if (me_num_events == nullptr) {
    edm::LogError(log_category_) << "failed to get " << source_folder + num_events_name << std::endl;
    return;
  }
  const int64_t num_events = me_num_events->getIntValue();
  const double data_taking_time = num_events * timing_window_;

  // map<key, pair<source, area> >
  std::map<std::string, std::pair<const MonitorElement*, const MonitorElement*> > me_pairs;
  for (const std::string& name : igetter.getMEs()) {
    if (name == num_events_name) continue;

    const std::string fullpath = source_folder + name;
    const MonitorElement* me = igetter.get(fullpath);
    if (me == nullptr) {
      edm::LogError(log_category_) << "failed to get " << fullpath << std::endl;
      continue;
    }

    const bool is_source = name.find(source_prefix) != std::string::npos;
    const bool is_area = name.find(area_prefix) != std::string::npos;

    std::string key = name;
    if (is_source)
      key.erase(key.find(source_prefix), source_prefix.length());
    else if (is_area)
      key.erase(key.find(area_prefix), area_prefix.length());
    else {
      edm::LogError(log_category_) << "unexpected" << std::endl;
      continue;
    }

    if (me_pairs.find(key) == me_pairs.end()) {
      me_pairs[key] = {nullptr, nullptr};
    }

    if (is_source)
      me_pairs[key].first = me;
    if (is_area)
      me_pairs[key].second = me;
  } // igetter.getMEs()

  ibooker.setCurrentFolder(hit_rate_folder);
  for (auto&& [key, value] : me_pairs) {
    const auto& [me_source, me_area] = value;
    if (me_source == nullptr) {
      edm::LogError(log_category_) << "source is missing" << std::endl;
      continue;
    }

    if (me_area == nullptr) {
      edm::LogError(log_category_) << "area is missing" << std::endl;
      continue;
    }

    TH1F* h_source = me_source->getTH1F();
    if (h_source == nullptr) {
      edm::LogError(log_category_) << "failed to get TH1F from " << key << std::endl;
      continue;
    }

    // weired but it returns double
    const double area = me_area->getFloatValue();
    if (area <= 0) {
      auto msg = Form("got invalid area=%f", area);
      edm::LogError(log_category_) << msg << std::endl;
      continue;
    }

    // TODO replacement?
    const std::string name = "hit_rate_" + key;
    // TODO
    const std::string title = "Hit Rate";

    TH1F* h_hit_rate = dynamic_cast<TH1F*>(h_source->Clone(name.c_str()));
    h_hit_rate->SetTitle(title.c_str());
    h_hit_rate->GetYaxis()->SetTitle("Hit Rate [Hz/cm^{2}]");

    // TODO check scale
    const double scale = 1.0f / (data_taking_time * area);
    h_hit_rate->Scale(scale);

    ibooker.book1D(name, h_hit_rate);
  }
}



void GEMOfflineHarvester::dqmEndJob(DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter) {
  doHitRate(ibooker, igetter);
}

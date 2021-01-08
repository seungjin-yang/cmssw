#ifndef DQMOffline_Muon_GEMOfflineDQMBase_h
#define DQMOffline_Muon_GEMOfflineDQMBase_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "CondFormats/GEMObjects/interface/GEMeMap.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"


class GEMOfflineDQMBase : public DQMEDAnalyzer {
public:
  using MEMap = std::map<GEMDetId, dqm::impl::MonitorElement*>;

  explicit GEMOfflineDQMBase(const edm::ParameterSet&);
 
  inline int getVFATNumber(const int, const int, const int);
  inline int getVFATNumberByStrip(const int, const int, const int);
  inline int getMaxVFAT(const int);
  inline int getDetOccXBin(const int, const int, const int);

  // Re: region / St: station, La: layer, Ch: chamber parity, Et: eta partition
  inline GEMDetId getReStKey(const int, const int);
  inline GEMDetId getReStKey(const GEMDetId&);
  inline GEMDetId getReStLaKey(const GEMDetId&);
  inline GEMDetId getReStEtKey(const GEMDetId&);
  inline GEMDetId getKey(const GEMDetId&); // == getReStLaChEtKey

  int getDetOccXBin(const GEMDetId&, const edm::ESHandle<GEMGeometry>&);
  void setDetLabelsVFAT(MonitorElement*, const GEMStation*);
  void setDetLabelsEta(MonitorElement*, const GEMStation*);
  // the number of eta partitions per GEMChamber
  int getNumEtaPartitions(const GEMStation*);

  // FIXME template<typename T>
  void fillME(MEMap& me_map, const GEMDetId& key, const float x);
  void fillME(MEMap& me_map, const GEMDetId& key, const float x, const float y);

  template <typename T>
  inline bool checkRefs(const std::vector<T*>&);

  inline float toDegree(float);

  std::string log_category_;

  class BookingHelper {
  public:
    BookingHelper(DQMStore::IBooker& ibooker, const TString& name_suffix, const TString& title_suffix)
        : ibooker_(&ibooker), name_suffix_(name_suffix), title_suffix_(title_suffix) {}

    ~BookingHelper() {}

    MonitorElement* book1D(TString name,
                           TString title,
                           int nbinsx,
                           double xlow,
                           double xup,
                           TString x_title = "",
                           TString y_title = "Entries") {
      name += name_suffix_;
      title += title_suffix_ + ";" + x_title + ";" + y_title;
      return ibooker_->book1D(name, title, nbinsx, xlow, xup);
    }

    MonitorElement* book1D(TString name,
                           TString title,
                           std::vector<double>& x_binning,
                           TString x_title = "",
                           TString y_title = "Entries") {
      name += name_suffix_;
      title += title_suffix_ + ";" + x_title + ";" + y_title;
      TH1F* h_obj = new TH1F(name, title, x_binning.size() - 1, &x_binning[0]);
      return ibooker_->book1D(name, h_obj);
    }

    MonitorElement* book2D(TString name,
                           TString title,
                           int nbinsx,
                           double xlow,
                           double xup,
                           int nbinsy,
                           double ylow,
                           double yup,
                           TString x_title = "",
                           TString y_title = "") {
      name += name_suffix_;
      title += title_suffix_ + ";" + x_title + ";" + y_title;
      return ibooker_->book2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
    }

  private:
    DQMStore::IBooker* ibooker_;
    const TString name_suffix_;
    const TString title_suffix_;
  };  // BookingHelper
};

inline int GEMOfflineDQMBase::getMaxVFAT(const int station) {
  if (GEMSubDetId::station(station) == GEMSubDetId::Station::GE0)
    return GEMeMap::maxVFatGE0_;
  else if (GEMSubDetId::station(station) == GEMSubDetId::Station::GE11)
    return GEMeMap::maxVFatGE11_;
  else if (GEMSubDetId::station(station) == GEMSubDetId::Station::GE21)
    return GEMeMap::maxVFatGE21_;
  else
    return -1;
}

inline int GEMOfflineDQMBase::getVFATNumber(const int station, const int ieta, const int vfat_phi) {
  const int max_vfat = getMaxVFAT(station);
  return max_vfat * (ieta - 1) + vfat_phi;
}

inline int GEMOfflineDQMBase::getVFATNumberByStrip(const int station, const int ieta, const int strip) {
  const int vfat_phi = (strip % GEMeMap::maxChan_) ? strip / GEMeMap::maxChan_ + 1 : strip / GEMeMap::maxChan_;
  return getVFATNumber(station, ieta, vfat_phi);
}

inline int GEMOfflineDQMBase::getDetOccXBin(const int chamber, const int layer, const int n_chambers) {
  return n_chambers * (chamber - 1) + layer;
}

inline float GEMOfflineDQMBase::toDegree(float radian) {
  float degree = radian / M_PI * 180.f;
  if (degree < -5.f)
    degree += 360.f;
  return degree;
}

template <typename T>
inline bool GEMOfflineDQMBase::checkRefs(const std::vector<T*>& refs) {
  if (refs.empty())
    return false;
  if (refs.front() == nullptr)
    return false;
  return true;
}

inline GEMDetId GEMOfflineDQMBase::getReStKey(const int region, const int station) {
  // region, ring, station, layer, chamber, roll
  return GEMDetId{region, 1, station, 0, 0, 0};
}

inline GEMDetId GEMOfflineDQMBase::getReStKey(const GEMDetId& id) {
  return getReStKey(id.region(), id.station());
}

inline GEMDetId GEMOfflineDQMBase::getReStLaKey(const GEMDetId& id) {
  return GEMDetId{id.region(), 1, id.station(), id.layer(), 0, 0};
}

inline GEMDetId GEMOfflineDQMBase::getReStEtKey(const GEMDetId& id) {
  return GEMDetId(id.region(), 1, id.station(), 0, 0, id.roll());
}

inline GEMDetId GEMOfflineDQMBase::getKey(const GEMDetId& id) {
  return GEMDetId{id.region(), 1, id.station(), id.layer(), id.chamber() % 2, id.roll()};
}

#endif  // DQMOffline_Muon_GEMOfflineDQMBase_h

#ifndef THINSLICESAMPLE_hh
#define THINSLICESAMPLE_hh

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <map>
#include <sstream>
#include <stdexcept>

#include "fhiclcpp/ParameterSet.h"

namespace protoana {
std::string PreciseToString(const double val, const int n = 2);


class ThinSliceSample {
 public:
  ThinSliceSample(std::string name, int flux_type,
                  const std::vector<fhicl::ParameterSet> & selections,
                  const std::vector<double> & incident_bins,
                  bool is_signal = false, std::pair<double, double> range = {0., 0.});

  ~ThinSliceSample(){};

  void SetFactor(double f) {fFactor = f;};

  const std::map<int, TH1 *> & GetSelectionHists() const {
    return fSelectionHists;
  };

  TH1 * GetSelectionHist(int id) {
    return fSelectionHists.at(id);
  };

  TH1D & GetIncidentHist() {
    return fIncidentHist;
  };

  TH1D & GetRebinnedIncidentHist() {
    return fIncidentHistRebinned;
  };

  TH1 * GetRebinnedSelectionHist(int id) {
    return fSelectionHistsRebinned.at(id);
  };

  const std::string & GetName() const {
    return fSampleName;
  };

  const int & GetFluxType() const {
    return fFluxType;
  };

  const double & GetNominalFlux() const {
    return fNominalFlux;
  };

  void AddFlux(double val = 1.) {
    fNominalFlux += val;
  };

  void FillIncidentHist(const std::vector<double> & vals) {
    for (size_t i = 0; i < vals.size(); ++i) {
      fIncidentHist.Fill(vals.at(i));
    }
  };

  void FillSelectionHist(int id, double val) {
    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      fSelectionHists.at(id)->Fill(val);
    }
  };

  template <size_t N> void FillSelectionHist(int id, const double (& vals)[N]) {
    if (N < 1 || N > 3) {
      std::string message = "Error: trying to fill hists with too many values";
      message += std::to_string(N);
      throw std::runtime_error(message);
    }

    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      if (N == 1) {
        fSelectionHists.at(id)->Fill(vals[0]);
      }
      else if (N == 2) {
        ((TH2D*)fSelectionHists.at(id))->Fill(vals[0], vals[1]);
      }
      else if (N == 3) {
        ((TH3D*)fSelectionHists.at(id))->Fill(vals[0], vals[1], vals[2]);
      }
    }
  }

  void ScaleHists(double val) {
    fIncidentHist.Scale(val);
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      it->second->Scale(val);
    }
  };

  void SetDataMCScale(double val) {
    fDataMCScale = val;
    ScaleHists(fDataMCScale);
  };

  void SetFactorAndScale(double val) {
    ResetFactor();
    fFactor = val;
    ScaleHists(val);
  };

  void ResetFactor() {
    ScaleHists(1./fFactor);
    fFactor = 1.;
  };

  bool CheckIsSignal() {return fIsSignal;};
  bool CheckInSignalRange(double val) {return ((fRange.first < val) &&
                                               (val <= fRange.second));};
  const std::pair<double, double> & GetRange() const {return fRange;};

  void RefillRebinnedHists();
  void MakeRebinnedHists();

 private:
  double fFactor = 1.;
  std::string fSampleName;
  int fFluxType;
  double fNominalFlux = 0.;
  double fDataMCScale = 1.;
  bool fIsSignal;
  std::pair<double, double> fRange;

  void Rebin1D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin2D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin3D(TH1 * sel_hist, TH1 * rebinned);
  std::map<int, TH1 *> fSelectionHists;
  TH1D fIncidentHist;
  std::map<int, TH1 *> fSelectionHistsRebinned;
  TH1D fIncidentHistRebinned;
  bool fMadeRebinned = false;

};

}
#endif
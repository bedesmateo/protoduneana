#ifndef PDSPTHINSLICEFITTER_hh
#define PDSPTHINSLICEFITTER_hh

#include <string>
#include <vector>
#include <map>

#include "TTree.h"
#include "TFile.h"
#include "THStack.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TRandom3.h"

#include "ThinSliceSample.h"
#include "ThinSliceDataSet.h"
#include "ThinSliceDriver.h"

namespace protoana {

class PDSPThinSliceFitter {
 public:
  PDSPThinSliceFitter(std::string fcl_file, std::string output_file);
  void BuildMCSamples();
  void SaveMCSamples();
  void GetNominalFluxes();
  void BuildDataHists();
  void InitializeMCSamples();
  void CompareDataMC(bool post_fit = false);
  void ScaleMCToData();
  void RunFitAndSave();
  ~PDSPThinSliceFitter();

 private:
  void Configure(std::string fcl_file);
  void DefineFitFunction();
  void MakeMinimizer();
  void ParameterScans();

  ThinSliceDriver * fThinSliceDriver;
  std::map<int, std::vector<ThinSliceSample>> fSamples;
  ThinSliceDataSet fDataSet;
  std::map<int, bool> fIsSignalSample;
  TFile fMCFile;
  TTree * fMCTree;
  TFile fDataFile;
  TTree * fDataTree;
  TFile fOutputFile;
  ROOT::Math::Functor fFitFunction;
  std::unique_ptr<ROOT::Math::Minimizer> fMinimizer;


  std::map<int, TH1D> fSelectedDataHists;
  std::map<int, TH1D> fRebinnedSelectedDataHists;
  TH1D fIncidentDataHist;
  TH1D fRebinnedIncidentDataHist;

  THStack * fNominalIncidentMCStack;
  THStack * fPostFitIncidentMCStack;
  std::map<int, THStack *> fNominalSelectedMCStacks;
  std::map<int, THStack *> fPostFitSelectedMCStacks;

  std::map<int, double> fNominalFluxes;
  std::map<int, std::vector<double>> fFluxesBySample;
  double fDataFlux;
  double fMCDataScale = 1.;

  std::map<int, std::vector<double>> fSignalParameters;
  std::map<int, std::vector<std::string>> fSignalParameterNames;
  size_t fTotalSignalParameters;

  std::map<int, double> fFluxParameters;
  std::map<int, std::string> fFluxParameterNames;
  size_t fTotalFluxParameters;

  TRandom3 fRNG;

  //Configurable members
  std::string fMCFileName;
  std::string fDataFileName;
  std::string fTreeName;
  std::vector<fhicl::ParameterSet> fSelectionSets;
  std::vector<fhicl::ParameterSet> fSampleSets;
  std::map<int, std::string> fFluxTypes;
  int fMaxCalls;
  unsigned int fNScanSteps;
  double fTolerance, fLowerLimit, fUpperLimit;
  bool fReducedIncidentChi2;
  std::vector<std::pair<int, int>> fPlotStyle;
  bool fPlotRebinned;
  bool fRandomStart;
  std::string fDriverName;
  std::string fAnalysis;
  fhicl::ParameterSet fAnalysisOptions;
  
  std::vector<double> fIncidentRecoBins;
  //////////////////////////
};

}
#endif
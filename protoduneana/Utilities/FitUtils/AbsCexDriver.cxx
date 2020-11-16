#include "AbsCexDriver.h"
#include "TCanvas.h"
#include "THStack.h"

#include "ThinSliceDriverFactory.h"
DECLARE_THINSLICEDRIVER_FACTORY_NS(protoana::AbsCexDriver, protoana, AbsCexDriver)

protoana::AbsCexDriver::~AbsCexDriver() {}

protoana::AbsCexDriver::AbsCexDriver(
    const fhicl::ParameterSet & extra_options)
    : ThinSliceDriver(extra_options) {}

void protoana::AbsCexDriver::BuildMCSamples(
    TTree * tree,
    std::map<int, std::vector<ThinSliceSample>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<double>> & fluxes_by_sample) {
  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                            &true_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                            &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                            &reco_beam_incidentEnergies);

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (samples.find(sample_ID) == samples.end())
      continue;

    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID);
    bool is_signal = signal_sample_checks.at(sample_ID);

    if (!is_signal) {
      ThinSliceSample & sample = samples_vec.at(0);
      int flux_type = sample.GetFluxType();
      nominal_fluxes[flux_type] += 1.;
      fluxes_by_sample[sample_ID][0] += 1.;
      sample.AddFlux();
      if (reco_beam_incidentEnergies->size()) {
        sample.FillIncidentHist(*reco_beam_incidentEnergies);
        double energy[1] = {reco_beam_interactingEnergy};
        sample.FillSelectionHist(selection_ID, energy);
      }
    }
    else {
      //Iterate through the true bins and find the correct one
      for (size_t j = 0; j < samples_vec.size(); ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(true_beam_interactingEnergy)) {
          int flux_type = sample.GetFluxType();
          nominal_fluxes[flux_type] += 1.;
          fluxes_by_sample[sample_ID][j] += 1.;
          sample.AddFlux();
          if (reco_beam_incidentEnergies->size()) {
            sample.FillIncidentHist(*reco_beam_incidentEnergies);
            sample.FillSelectionHist(selection_ID, reco_beam_interactingEnergy);
          }
          break;
        }
      }
    }
  }
}

void protoana::AbsCexDriver::BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set) {
  int selection_ID; 
  double reco_beam_interactingEnergy;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                        &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                        &reco_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        selected_hists[selection_ID]->Fill(reco_beam_interactingEnergy);
      }
    }
  }
}

std::pair<double, size_t> protoana::AbsCexDriver::CalculateChi2(
    std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set) {

  double chi2 = 0.;
  size_t nPoints = 0;

  std::map<int, TH1 *> & selected_data_hists = data_set.GetSelectionHists();
  for (auto it = selected_data_hists.begin();
       it != selected_data_hists.end(); ++it) {
    TH1D * data_hist = (TH1D*)it->second;
    int selection_ID = it->first;
    for (int i = 1; i <= data_hist->GetNbinsX(); ++i) {
      double data_val = data_hist->GetBinContent(i);
      double data_err = data_hist->GetBinError(i);

      double mc_val = 0.;
      //Go through all the samples and get the values from mc
      for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
        std::vector<ThinSliceSample> & samples_vec = it2->second;
        for (size_t j = 0; j < samples_vec.size(); ++j) {
          ThinSliceSample & sample = samples_vec[j];
          mc_val += sample.GetSelectionHist(selection_ID)->GetBinContent(i);
        }
      }
      chi2 += (std::pow((data_val - mc_val), 2) / std::pow(data_err, 2));
      ++nPoints;
    }
  }

  //Go through the incident hist. Calculate it as its own point.
  //Then add reduced incident chi2 to full chi2 
  if (!fExtraOptions.get<bool>("SkipIncidentChi2")) {
    double incident_chi2 = 0.;
    size_t incident_points = 0;
    TH1D & incident_data_hist = data_set.GetIncidentHist();
    for (int i = 1; i <= incident_data_hist.GetNbinsX(); ++i) {
      double data_val = incident_data_hist.GetBinContent(i);
      double data_err = incident_data_hist.GetBinError(i);

      double mc_val = 0.;
      //Go through all the samples and get the values from mc
      for (auto it = samples.begin(); it != samples.end(); ++it) {
        std::vector<ThinSliceSample> & samples_vec = it->second;
        for (size_t j = 0; j < samples_vec.size(); ++j) {
          ThinSliceSample & sample = samples_vec[j];
          mc_val += sample.GetIncidentHist().GetBinContent(i);
        }
      }

      incident_chi2 += (std::pow((data_val - mc_val), 2) / std::pow(data_err, 2));
      ++incident_points;
    }

    if (fExtraOptions.get<bool>("ReducedChi2")) {
      //Add reduced incident chi2 to full chi2, increment total nPoints
      chi2 += incident_chi2 / incident_points;
      ++nPoints;
    }
    else {
      chi2 += incident_chi2;
      nPoints += incident_points;
    }
  }

  return {chi2, nPoints};
}

void protoana::AbsCexDriver::CompareSelections(
    std::map<int, std::vector<ThinSliceSample>> & samples,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    bool plot_rebinned,
    bool post_fit) {

  output_file.cd();
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    TH1D * data_hist = (TH1D*)it->second;
    data_hist->SetLineColor(kBlack);
    data_hist->SetMarkerColor(kBlack);
    data_hist->SetMarkerStyle(20);   

    std::string canvas_name = "c";
    canvas_name += (post_fit ? "PostFit" : "Nominal") +
                   data_set.GetSelectionName(selection_ID);
    TCanvas cSelection(canvas_name.c_str(), "");
    cSelection.SetTicks();

    //Build the mc stack
    std::string stack_name = (post_fit ? "PostFit" : "Nominal") +
                             data_set.GetSelectionName(selection_ID) +
                             "Stack";
    THStack mc_stack(stack_name.c_str(), "");
    size_t iColor = 0;
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        TH1D * sel_hist
            = (TH1D*)(plot_rebinned ?
                      it2->second.at(i).GetRebinnedSelectionHist(selection_ID) :
                      it2->second.at(i).GetSelectionHist(selection_ID));

        std::pair<int, int> color_fill = GetColorAndStyle(iColor, plot_style);
        sel_hist->SetFillColor(color_fill.first);
        sel_hist->SetFillStyle(color_fill.second);
        sel_hist->SetLineColor(kBlack);
        mc_stack.Add(sel_hist);
        ++iColor;
      }
    }
    mc_stack.Write();


    mc_stack.Draw("hist");
    std::string title = data_set.GetSelectionName(selection_ID) +
                        ";Reconstructed KE (MeV)";
    mc_stack.SetTitle(title.c_str());
    mc_stack.GetHistogram()->SetTitleSize(.04, "X");
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");

    cSelection.Write();

    //Get the full incident hist from stack
    TList * l = (TList*)mc_stack.GetHists();
    TH1D * hMC = (TH1D*)l->At(0)->Clone();
    for (int i = 1; i < l->GetSize(); ++i) {
      hMC->Add((TH1D*)l->At(i));
    }

    std::string ratio_name = data_set.GetSelectionName(selection_ID) + "Ratio" +
                             (post_fit ? "PostFit" : "Nominal");
    TH1D * hRatio
        = (TH1D*)data_hist->Clone(ratio_name.c_str());
    hRatio->Divide(hMC);
    hRatio->Write(); 
  }
}
////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEelectronWireAna
// Plugin Type: analyzer (art v3_03_01)
// File:        ProtoDUNEelectronWireAna_module.cc
// 
// from cetlib version v3_08_00. 
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimChannel.h" // Ewerton
#include "larreco/Calorimetry/CalorimetryAlg.h" // Ewerton
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "dune/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "dune/Calib/LifetimeCalib.h" 
#include "dune/CalibServices/LifetimeCalibService.h" 
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"

#include "TFile.h"
#include "TProfile.h"
#include "TH1.h"
#include "TTree.h"
#include <stdio.h>
#include <stdlib.h> 
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

const int MAXHits = 6000;//max number of hits

class ProtoDUNEelectronWireAna;


class ProtoDUNEelectronWireAna : public art::EDAnalyzer {
public:
  explicit ProtoDUNEelectronWireAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtoDUNEelectronWireAna(ProtoDUNEelectronWireAna const&) = delete;
  ProtoDUNEelectronWireAna(ProtoDUNEelectronWireAna&&) = delete;
  ProtoDUNEelectronWireAna& operator=(ProtoDUNEelectronWireAna const&) = delete;
  ProtoDUNEelectronWireAna& operator=(ProtoDUNEelectronWireAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;
  void beginJob() override;

private:

  void Initialize();
  void FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle);

  // Declare member data here.
  protoana::ProtoDUNEDataUtils dataUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEBeamlineUtils beamlineUtil;
 
  geo::GeometryCore const * fGeometry = &*(art::ServiceHandle<geo::Geometry>());
  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  const art::InputTag fBeamModuleLabel;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fShowerCaloTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fWireTag;
  std::string fSimChannelTag;
  std::string fRawDigitTag;
  calo::CalorimetryAlg fCalorimetryAlg;
  int fVerbose;

  TTree *fTree;

  // Tree variables
  int fRun;
  int fSubRun;
  int fevent;

  // Beam track
  std::vector<int> fcerenkovStatus;
  double fbeamtrackMomentum;
  std::vector<double> fbeamtrackDir;

  // Reconstructed tracks/showers
  int fprimaryIstrack;
  int fprimaryIsshower;
  int fprimaryNHits;
  std::vector<double> fprimaryStartDirection;
  int    fprimaryShower_nHits; //collection only
  std::vector<int>    fprimaryShower_hit_w;
  std::vector<double> fprimaryShower_hit_q;
  std::vector<double> fprimaryShower_hit_t; 
  std::vector<double> fprimaryShower_hit_X;
  std::vector<int>    fprimaryShower_hit_ch;
  std::vector<double> fprimaryShower_hit_Y;
  std::vector<double> fprimaryShower_hit_Z;
  std::vector<double> fprimaryShower_hit_nElectrons;

  double fLifetime; //us
  double fdet_frequency; // Frequency in MHz

  // Ewerton (piece of code to select track-like beam particles)
  Int_t tree_Track_index= 0;
  Int_t tree_Shower_index= 0;

//  int fNSimChannels =0;
//  int fNTrackIDEs=0;

  std::vector<int>    fprimary_Shower_wire_w;
  std::vector<double> fprimary_Shower_wire_ch; //charge 
  std::vector<double> fprimary_Shower_wire_ch_corr;
  std::vector<double> fprimary_Shower_wire_X;  
  std::vector<double> fprimary_Shower_wire_Z;  
  std::vector<double> fprimary_Shower_wire_Y; 
//  std::vector<std::vector<std::pair<int,double>>> fprimary_Shower_wire_time_charge;
/*
  std::vector<double> fprimary_Shower_simchannel_electrons_t1t2window; // number of ionization electrons on simchannel in window (t1,t2)
  std::vector<std::vector<double>> fprimary_Shower_simchannel_simide_energyFrac;   // fraction of hit energy from the particle on channel 
  std::vector<std::vector<double>> fprimary_Shower_simchannel_simide_energy;       // energy from the particle on channel
  std::vector<std::vector<double>> fprimary_Shower_simchannel_simide_numElectrons;  // number of electrons from particle on channel 
*/
  std::vector<double>  fprimary_Shower_MCwire_E;  //energy in MeV
  std::vector<int>  fprimary_Shower_MCwire_w;
  std::vector<double> fprimaryStartPosition;

  TH1D *rawwf[3][800];
  TH1D *deconwf[3][800];
  double charge[3][800];
  double deconchg[3][800];
};

void ProtoDUNEelectronWireAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("WireAna","Recob::Wire");
  fTree->Branch("run",                        &fRun,"run/I");
  fTree->Branch("subrun",                     &fSubRun,"subrun/I");
  fTree->Branch("event",                      &fevent,"event/I");

  fTree->Branch("cerenkovStatus",             &fcerenkovStatus);
  fTree->Branch("beamtrackMomentum",          &fbeamtrackMomentum,      "beamtrackMomentum/D");
  fTree->Branch("beamtrackDir",               &fbeamtrackDir);


  fTree->Branch("primaryIstrack",             &fprimaryIstrack,         "primaryIstrack/I");
  fTree->Branch("primaryIsshower",            &fprimaryIsshower,        "primaryIsshower/I");
  fTree->Branch("primaryNHits",               &fprimaryNHits,           "primaryNHits/I");
  fTree->Branch("primaryStartPosition",       &fprimaryStartPosition);
  fTree->Branch("primaryStartDirection",      &fprimaryStartDirection);

  fTree->Branch("primaryShower_nHits",        &fprimaryShower_nHits,    "primaryShower_nHits/I");
  fTree->Branch("primaryShower_hit_q",        &fprimaryShower_hit_q);
  fTree->Branch("primaryShower_hit_w",        &fprimaryShower_hit_w);
  fTree->Branch("primaryShower_hit_t",        &fprimaryShower_hit_t);
  fTree->Branch("primaryShower_hit_X",        &fprimaryShower_hit_X);
  fTree->Branch("primaryShower_hit_Y",        &fprimaryShower_hit_Y);
  fTree->Branch("primaryShower_hit_Z",        &fprimaryShower_hit_Z);
  fTree->Branch("primaryShower_hit_nElectrons",&fprimaryShower_hit_nElectrons);
  fTree->Branch("primaryShower_hit_ch",       &fprimaryShower_hit_ch);

  // Ewerton (added from CalibrationdEdXPDSP_module.cc)
  fTree->Branch("lifetime",  &fLifetime, "lifetime/D");

  fTree->Branch("primary_Shower_wire_w",           &fprimary_Shower_wire_w);
  fTree->Branch("primary_Shower_wire_ch",          &fprimary_Shower_wire_ch);
  fTree->Branch("primary_Shower_wire_ch_corr",     &fprimary_Shower_wire_ch_corr);
  fTree->Branch("primary_Shower_wire_X",           &fprimary_Shower_wire_X);
  fTree->Branch("primary_Shower_wire_Z",           &fprimary_Shower_wire_Z);
  fTree->Branch("primary_Shower_wire_Y",           &fprimary_Shower_wire_Y);
//  fTree->Branch("primary_Shower_wire_time_charge", &fprimary_Shower_wire_time_charge);
/*
  fTree->Branch("NSimChannels",               &fNSimChannels,           "NSimChannels/I");
  fTree->Branch("NTrackIDEs",                 &fNTrackIDEs,             "NTrackIDEs/I");
  fTree->Branch("primary_Shower_simchannel_electrons_t1t2window", &fprimary_Shower_simchannel_electrons_t1t2window);
  fTree->Branch("primary_Shower_simchannel_simide_energyFrac",    &fprimary_Shower_simchannel_simide_energyFrac);
  fTree->Branch("primary_Shower_simchannel_simide_energy",        &fprimary_Shower_simchannel_simide_energy);
  fTree->Branch("primary_Shower_simchannel_simide_numElectrons",  &fprimary_Shower_simchannel_simide_numElectrons);
*/
  fTree->Branch("primary_Shower_MCwire_w",                        &fprimary_Shower_MCwire_w);
  fTree->Branch("primary_Shower_MCwire_E",                        &fprimary_Shower_MCwire_E);

 
  // Ewerton (piece of code to select track-like beam particles)
  fTree->Branch("Track_index",  &tree_Track_index);
  fTree->Branch("Shower_index", &tree_Shower_index);

  //fTree->Branch("charge", charge, "charge[3][800]/D");
  //fTree->Branch("deconchg", deconchg, "deconchg[3][800]/D");
/*
   for (int i = 0; i<3; ++i){
     for (int j = 0; j<800; ++j){
       rawwf[i][j] = tfs->make<TH1D>(Form("rawwf_%d_%d",i,j),Form("Plane %d Wire %d",i,j),6000,0,6000);
       deconwf[i][j] = tfs->make<TH1D>(Form("deconwf_%d_%d",i,j),Form("Plane %d Wire %d",i,j),6000,0,6000);
     }
   }
*/
}


ProtoDUNEelectronWireAna::ProtoDUNEelectronWireAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  dataUtil(p.get<fhicl::ParameterSet>("DataUtils")),
  beamlineUtil(p.get<fhicl::ParameterSet>("BeamLineUtils")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fShowerCaloTag(p.get<std::string>("ShowerCalorimetryTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fWireTag(p.get<std::string>("WireTag")),
  fSimChannelTag(p.get<std::string>("SimChannelTag")),
  fRawDigitTag(p.get<std::string>("RawDigitTag")),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg")),
  fVerbose(p.get<int>("Verbose"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ProtoDUNEelectronWireAna::analyze(art::Event const& evt)
{

  Initialize(); 
  fRun = evt.run();
  fSubRun = evt.subRun();
  fevent  = evt.id().event(); 
   
  std::cout<<"          INITIAL CHECK"<<std::endl;

  //======================================
  // Begin ProtoDUNEelectronAnaTree_module
  //======================================
  bool beamTriggerEvent = false;
  if(!evt.isRealData())
    {
    // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);

    // Also get the reconstructed beam information in the MC - TO DO
    const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
//    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    if(geantGoodParticle != 0x0)
      {
      beamTriggerEvent = true;
      fbeamtrackMomentum = geantGoodParticle->P();
      }

    if( abs(geantGoodParticle->PdgCode()) != 11 ) return;
    }//MC

  //			000
  std::cout<<"		FIRST CHECK"<<std::endl;
  //			000  

  else
    { //data
    //                    000
    std::cout<<"          SECOND CHECK"<<std::endl;
    //                    000
    // For data we can see if this event comes from a beam trigger
    beamTriggerEvent = dataUtil.IsBeamTrigger(evt);
    if( !beamTriggerEvent ) return;
    //                    000
    std::cout<<"          THIRD CHECK"<<std::endl;
    //                    000
    art::Handle< std::vector<beam::ProtoDUNEBeamEvent> > pdbeamHandle;
    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    if(evt.getByLabel(fBeamModuleLabel, pdbeamHandle))
      art::fill_ptr_vector(beaminfo, pdbeamHandle);
    else
      {
      std::cout<<"No beam information from "<<fBeamModuleLabel<<std::endl;
      } 
    if(beaminfo.size())
      {
      if( beamlineUtil.IsGoodBeamlineTrigger( evt ) )
        {  
        // Get Cerenkov
        fcerenkovStatus.push_back(beaminfo[0]->GetCKov0Status());
        fcerenkovStatus.push_back(beaminfo[0]->GetCKov1Status());
        auto & tracks = beaminfo[0]->GetBeamTracks();
        if(!tracks.empty())
          {
          fbeamtrackDir.push_back(tracks[0].StartDirection().X());
          fbeamtrackDir.push_back(tracks[0].StartDirection().Y());
          fbeamtrackDir.push_back(tracks[0].StartDirection().Z());
          }
        auto & beammom = beaminfo[0]->GetRecoBeamMomenta();
        if(!beammom.empty()) fbeamtrackMomentum = beammom[0];
        } //good beam trigger
      }
    }//for data
  //                    000
  std::cout<<"          THIRD CHECK"<<std::endl;
  //                    000
  //check for reco pandora stuff
  art::Handle<std::vector<recob::PFParticle>> recoParticleHandle;
  if( !evt.getByLabel(fPFParticleTag,recoParticleHandle) ) return;

  // Get all of the PFParticles, by default from the "pandora" product
  auto recoParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  // We'd like to find the beam particle. Pandora tries to do this for us, so let's use the PFParticle utility 
  // to look for it. Pandora reconstructs slices containing one (or sometimes more) primary PFParticles. These
  // are tagged as either beam or cosmic for ProtoDUNE. This function automatically considers only those
  // PFParticles considered as primary
  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
  for(const recob::PFParticle* particle : pfParticles){
    //			OOO
    std::cout<< particle  <<std::endl:
    std::cout<<"pfParticles = "<< pfParticles <<std::endl:
    //			000
    FillPrimaryPFParticle(evt, particle);
    
    // For now only consider the first primary track. Need a proper treatment if more than one primary particles are found.
    break;
  }

///*
  //------------------------------------------------------------
  // Ewerton (added from CalibrationdEdXPDSP_module.cc)
  //------------------------------------------------------------
  // Electron lifetime from database calibration service provider
  art::ServiceHandle<calib::LifetimeCalibService> lifetimecalibHandler;
  calib::LifetimeCalibService & lifetimecalibService = *lifetimecalibHandler; 
  calib::LifetimeCalib *lifetimecalib = lifetimecalibService.provider();
  
  if (evt.isRealData()) {
    fLifetime = lifetimecalib->GetLifetime()*1000.0; // [ms]*1000.0 -> [us]
    std::cout << "\nuse lifetime from database(run " << fRun << ")   " << fLifetime << "\n\n";
  } 
  else {
    fLifetime = detprop->ElectronLifetime(); // [us] 
    std::cout << "\nuse lifetime from DetectorProperties(run " << fRun << ")   " << fLifetime << "\n\n";
  }
  //------------------------------------------------------------
  // End Ewerton (added from CalibrationdEdXPDSP_module.cc)
  //------------------------------------------------------------

  //------------------------------------------------------------
  // End Ewerton (piece of code to select track-like beam particles)
  //------------------------------------------------------------
  Int_t Shower_index=0;
  Int_t Track_index=0;
  if( recoParticles->size() > 0 ){
    const art::FindManyP<recob::Track> findTracks(recoParticles,evt,fTrackerTag);
    const std::vector<art::Ptr<recob::Track>> pfpTracks = findTracks.at(recoParticles->at(0).Self());
    if(pfpTracks.size() != 0){Track_index=1;}
    const art::FindManyP<recob::Shower> findShowers(recoParticles,evt,fShowerTag);
    const std::vector<art::Ptr<recob::Shower>> pfpShowers = findShowers.at(recoParticles->at(0).Self());
    if(pfpShowers.size() != 0){Shower_index=1;}
  }
  if(Track_index==1 && Shower_index==0){std::cout << "Track-like   \n\n";}
  if(Track_index==0 && Shower_index==1){std::cout << "Shower-like  \n\n";}
  if(Track_index==1 && Shower_index==1){std::cout << "Ambiguity    \n\n";}
  if(Track_index==0 && Shower_index==0){std::cout << "Unclassified \n\n";}
  tree_Track_index=Track_index;
  tree_Shower_index=Shower_index;
  //------------------------------------------------------------
  // End Ewerton (piece of code to select track-like beam particles)
  //------------------------------------------------------------
  //======================================
  // End ProtoDUNEelectronAnaTree_module
  //======================================

  //check for reco pandora stuff
//  art::Handle<std::vector<recob::PFParticle>> recoParticleHandle;
//  if( !evt.getByLabel(fPFParticleTag,recoParticleHandle) ) return;

  // sim channel info (Ewerton)
  art::Handle<std::vector<sim::SimChannel>> simChannelHandle;
  if( !evt.isRealData() ) {if( !evt.getByLabel(fSimChannelTag,simChannelHandle) ) return;}


  //std::vector<int> hit_w;
  std::map<int,std::vector<double>> hit_w_and_t1;
  std::map<int,std::vector<double>> hit_w_and_t2;
  std::map<int,std::vector<double>> hit_w_and_y;
  std::map<int,std::vector<double>> hit_w_and_x;
//  std::vector<const recob::PFParticle*> pfParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
//std::cout << "\n event " << fevent << ": [PFParticlesFromBeamSlice] pfParticles.size()= " << pfParticles.size() << "\n\n";
  const simb::MCParticle* mcparticle = NULL;
  bool doAna = false;
  for(const recob::PFParticle* particle : pfParticles){
     const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
     if( thisShower == 0x0 ) continue;
     if( !evt.isRealData() ){
       mcparticle = truthUtil.GetMCParticleFromRecoShower(*thisShower, evt, "pandoraShower");
       if( abs(mcparticle->PdgCode()) != 11 ) return;
     }
     fprimaryStartPosition.push_back(thisShower->ShowerStart().X());
     fprimaryStartPosition.push_back(thisShower->ShowerStart().Y());
     fprimaryStartPosition.push_back(thisShower->ShowerStart().Z());
     doAna = true; 
     const std::vector<const recob::Hit*> sh_hits = showerUtil.GetRecoShowerHits(*thisShower, evt, fShowerTag);
     art::FindManyP<recob::Wire> wFromHits(sh_hits,evt,"hitpdune");
     art::FindManyP<recob::SpacePoint> spFromShowerHits(sh_hits,evt,fPFParticleTag);
     for( size_t j=0; j<sh_hits.size() && j<MAXHits; ++j){
       if( sh_hits[j]->WireID().Plane != 2 ) continue;
       std::vector<art::Ptr<recob::Wire>> wires = wFromHits.at(j);
       //hit_w.push_back(wires[0]->Channel());
//std::cout << "\n event " << fevent << ": shower hit " << j << " | wire " << wires[0]->Channel() << ": PeakTimeMinusRMS=" << sh_hits[j]->PeakTimeMinusRMS(5.0);
//std::cout << "\n event " << fevent << ": shower hit " << j << " | wire " << wires[0]->Channel() << ": PeakTimePlusRMS=" << sh_hits[j]->PeakTimePlusRMS(5.0);
       hit_w_and_t1[wires[0]->Channel()].push_back(sh_hits[j]->PeakTimeMinusRMS(5.0));
       hit_w_and_t2[wires[0]->Channel()].push_back(sh_hits[j]->PeakTimePlusRMS(5.0));
//std::cout << "\n sh_hits[j]->WireID().TPC=" << sh_hits[j]->WireID().TPC;
       hit_w_and_x[wires[0]->Channel()].push_back(detprop->ConvertTicksToX(sh_hits[j]->PeakTime(),sh_hits[j]->WireID().Plane,sh_hits[j]->WireID().TPC,0));
       std::vector<art::Ptr<recob::SpacePoint>> sp = spFromShowerHits.at(j); 
       if(!sp.empty()){
          hit_w_and_y[wires[0]->Channel()].push_back(sp[0]->XYZ()[1]);
       }
       else hit_w_and_y[wires[0]->Channel()].push_back(fprimaryStartPosition[1]); //use vtx if no sp 
     }
     break;
  }
 
  //one wire can have various hits so lets remove duplicate wires
  //std::sort(hit_w.begin(),hit_w.end());
  //hit_w.erase(std::unique(hit_w.begin(),hit_w.end()), hit_w.end()); 
  

  if( doAna ){
    auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >(fWireTag);
    //auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >("wclsdatasp:gauss"); //new reco 
    auto w1 = hit_w_and_t1.begin();
    auto w2 = hit_w_and_t2.begin();
    auto x  = hit_w_and_x.begin();
    auto y  = hit_w_and_y.begin(); 

    while( w1 != hit_w_and_t1.end()){
      int it_w = w1->first;
      int n_hits = w1->second.size();
// std::cout<<it_w<<" "<<n_hits<<" "<<w1->second[n_hits-1]<<" "<<w2->second[n_hits-1]<<std::endl;
      double t1 =  w1->second[0]; //first hit
      double t2 =  w2->second[n_hits-1]; //last hit
//std::cout << "\n t1=" << t1;
//std::cout << "\n t2=" << t2;
//std::cout << "\n deltat=" << t2-t1;
      double x_w, y_w;
      if( x->second.size() < 1 ){           //<<<<<<---- Ewerton: should it be >1?
        x_w = (x->second[0]-x->second[n_hits-1])/2.0;
        y_w = (y->second[0]-y->second[n_hits-1])/2.0;
      }
      else {
        x_w = x->second[0];
        y_w = y->second[0];
      }
      for(auto & wire : * wires){
         int channel_no = wire.Channel();
         int plane = fGeometry->View(wire.Channel()); 
         if( plane != 2 ) continue;
         std::vector< geo::WireID > wireID= fGeometry->ChannelToWire(channel_no);
         const geo::WireGeo* pwire = fGeometry->WirePtr(wireID.at(0)); //for collection plane there is one wire per channel
         TVector3 xyzWire = pwire->GetCenter<TVector3>();
         if( it_w == channel_no ){  
           double charge =0.0;
           double charge_corr =0.0;
//           std::vector<std::pair<int,double>> tick_charge;
           for(size_t i = 0; i < wire.Signal().size(); ++i){
              if( i > t1 && i < t2 ) {
                charge += wire.Signal()[i];  //<<<<<<---- Ewerton: comparing index to ticks? 
                double tick_in_us = i*500; // tick in ns
                tick_in_us /= 1000; // tick in us
                charge_corr += wire.Signal()[i]*exp(tick_in_us/fLifetime);
//                tick_charge.push_back(std::make_pair(i,wire.Signal()[i]));
              }
           }
           fprimary_Shower_wire_Z.push_back(xyzWire.Z()); 
           fprimary_Shower_wire_Y.push_back(y_w); 
           fprimary_Shower_wire_X.push_back(x_w); 
           fprimary_Shower_wire_ch.push_back(charge);
           fprimary_Shower_wire_ch_corr.push_back(charge_corr);
           fprimary_Shower_wire_w.push_back(channel_no);
//           fprimary_Shower_wire_time_charge.push_back(tick_charge);
           // --- Ewerton ---
           double this_charge = 0;
           //int _wire = wireID[0].Wire;
           for(size_t i = 0; i < wire.Signal().size(); ++i){
              if( i > t1 && i < t2 ) {
                this_charge += wire.Signal()[i];  //<<<<<<---- Ewerton: comparing index to ticks? 
                //deconwf[plane][_wire]->SetBinContent(i, wire.Signal()[i]);
              }
           }   
/*
           if( !evt.isRealData() ){
             unsigned int isimch=0;
             fNSimChannels=0;
             for(auto simchannel : *simChannelHandle) {
               // Get all ionization electrons in interval t1<tick<t2
               double total_electrons=0; // number of electrons
               if((int)simchannel.Channel()==it_w) {
                 for(double tick=t1; tick<=t2; tick++) { // <<--- tick steps of 1?
                   total_electrons += simchannel.Charge(tick);
                 }
                 fprimary_Shower_simchannel_electrons_t1t2window.push_back(total_electrons);
//                 fprimary_Shower_simchannel_electrons_t1t2window[isimch] = total_electrons;
                 // test
                 std::vector<sim::TrackIDE> trackIDEs = simchannel.TrackIDEs(t1,t2);
                 std::vector<double> ide_energyFrac;
                 std::vector<double> ide_energy;
                 std::vector<double> ide_numElectrons;
                 unsigned int itrkide=0;
                 fNTrackIDEs=0;
                 for(auto ide : trackIDEs) {
                  ide_energyFrac.push_back(ide.energyFrac);
                  ide_energy.push_back(ide.energy);
                  ide_numElectrons.push_back(ide.numElectrons);
                   itrkide++;
                   fNTrackIDEs++;
                 }
                 fprimary_Shower_simchannel_simide_energyFrac.push_back(ide_energyFrac);
                 fprimary_Shower_simchannel_simide_energy.push_back(ide_energy);
                 fprimary_Shower_simchannel_simide_numElectrons.push_back(ide_numElectrons);

                 break;
               }
//std::cout << "\n total_electrons=" << total_electrons << "\n\n";
               isimch++;
               fNSimChannels++;
             }
           }
           // --- End Ewerton ---    
*/       
           break;
         }
      }// all recob::wire
      w1 ++; w2 ++;
      x ++; y ++;
    }//wire from shower 

//std::cout << "\n fprimary_Shower_wire_time_charge.size()=" << fprimary_Shower_wire_time_charge.size() << "\n\n";

    //look at MC info if available
    if(!evt.isRealData()){ 
      auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
      const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
      art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
      int trackid = mcparticle->TrackId();
      if( mcparticle->TrackId() == geantGoodParticle->TrackId() ){
        if( doAna ){   //we have a reco shower
          if(evt.getByLabel("largeant", simchannelHandle)){
            for(auto const& simchannel : (*simchannelHandle)){
               if(fGeometry->View(simchannel.Channel()) != 2) continue;
               auto const& alltimeslices = simchannel.TDCIDEMap();
               double EperCh=0;
               for(auto const& tslice : alltimeslices){
	          auto const& simide = tslice.second;
	          // Loop over energy deposits
	          for(auto const& eDep : simide){
	             if(eDep.trackID == trackid || eDep.trackID == -trackid){
                       EperCh += eDep.energy;
                     }
	          }
               } 
               if( EperCh== 0 ) continue;//save only channeles with ID
               fprimary_Shower_MCwire_w.push_back(simchannel.Channel());
               fprimary_Shower_MCwire_E.push_back(EperCh);
            }
          }
        }//do analysis of MC
      }//is MC same as beam? 
    }//is MC?
  }//shower


//================TEST===================
/*
  for (int i = 0; i<3; ++i){
    for (int j = 0; j<800; ++j){
      charge[i][j] = 0;
      deconchg[i][j] = 0;
    }
  }
*/
  // Implementation of required member function here.
  art::Handle< std::vector<raw::RawDigit> > rawListHandle;
  std::vector<art::Ptr<raw::RawDigit> > rawlist;
  if (evt.getByLabel(fRawDigitTag, rawListHandle)) //tpcrawdecoder:daq
    art::fill_ptr_vector(rawlist, rawListHandle);

  art::Handle < std::vector < recob::Wire > > wireListHandle;
  std::vector < art::Ptr < recob::Wire > > wires;
  if (evt.getByLabel(fWireTag, wireListHandle)) {//wclsdatanfsp:gauss
    art::fill_ptr_vector(wires, wireListHandle);
  }
/*
  auto const* geo = lar::providerFrom<geo::Geometry>();

 for (unsigned int ich = 0; ich < rawlist.size(); ++ich){
   const auto & digitVec = rawlist[ich];
   const auto & wireid = geo->ChannelToWire(digitVec->Channel());
   if (wireid[0].TPC != 1) continue;
   //if (wireid[0].Plane != 2) continue;
   int wire = wireid[0].Wire;
   int plane = wireid[0].Plane;
   double this_charge = 0;
   std::vector<short> rawadc(6000);
   raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
   for (size_t itck = 0; itck < rawadc.size(); ++itck){
     this_charge += rawadc[itck] - digitVec->GetPedestal();
     rawwf[plane][wire]->SetBinContent(itck, rawadc[itck] - digitVec->GetPedestal());
//std::cout << "\n rawadc[itck] - digitVec->GetPedestal()=" << rawadc[itck] - digitVec->GetPedestal();
   }
   charge[plane][wire] = this_charge;
 }
*/
/*
 for (unsigned int ich = 0; ich < wires.size(); ++ich){
   const auto & wireid = geo->ChannelToWire(wires[ich]->Channel());
   if (wireid[0].TPC != 1) continue;
   //if (wireid[0].Plane != 2) continue;
   int wire = wireid[0].Wire;
   int plane = wireid[0].Plane;
   double this_charge = 0;
   const auto & signal = wires[ich]->Signal();

   for (size_t itck = 0; itck < signal.size(); ++itck){
     this_charge += signal[itck];
     deconwf[plane][wire]->SetBinContent(itck, signal[itck]);
//std::cout << "\n signal[itck]=" << signal[itck];
   }
   deconchg[plane][wire] = this_charge;
 }
*/
//================END TEST===================


//ElectronsToADC= 6.8906513e-3;

//  fTree->Fill();
  if(beamTriggerEvent) fTree->Fill();
}

// -----------------------------------------------------------------------------
void ProtoDUNEelectronWireAna::FillPrimaryPFParticle(art::Event const & evt, const recob::PFParticle* particle){

  // NHits associated with this pfParticle
  fprimaryNHits = (pfpUtil.GetPFParticleHits(*particle,evt,fPFParticleTag)).size();

  // "particle" is the pointer to our beam particle. The recob::Track or recob::Shower object
  // of this particle might be more helpful. These return null pointers if not track-like / shower-like
  const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*particle, evt,fPFParticleTag,fTrackerTag);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);

//  const simb::MCParticle* mcparticle = NULL; 
  if(thisTrack != 0x0){
    fprimaryIstrack                    = 1;
    fprimaryIsshower                   = 0;
    fprimaryStartDirection.push_back(thisTrack->Trajectory().StartDirection().X());
    fprimaryStartDirection.push_back(thisTrack->Trajectory().StartDirection().Y());
    fprimaryStartDirection.push_back(thisTrack->Trajectory().StartDirection().Z());

  } // end is track
  else if(thisShower != 0x0){
    fprimaryIstrack                     = 0;
    fprimaryIsshower                    = 1;
    fprimaryStartDirection.push_back(thisShower->Direction().X());
    fprimaryStartDirection.push_back(thisShower->Direction().Y());
    fprimaryStartDirection.push_back(thisShower->Direction().Z());
   
    const std::vector<const recob::Hit*> sh_hits = showerUtil.GetRecoShowerHits(*thisShower, evt, fShowerTag);
    art::FindManyP<recob::SpacePoint> spFromShowerHits(sh_hits,evt,fPFParticleTag);

    int idx =0;
 //   fprimaryShowerCharge =0.0;
    for( size_t j=0; j<sh_hits.size() && j<MAXHits; ++j){
       if( sh_hits[j]->WireID().Plane != 2 ) continue;
       art::FindManyP<recob::Wire> wFromHits(sh_hits,evt,"hitpdune");
       std::vector<art::Ptr<recob::Wire>> wires = wFromHits.at(j);
       const geo::WireGeo* pwire = fGeometry->WirePtr(sh_hits[j]->WireID());
       TVector3 xyzWire = pwire->GetCenter<TVector3>();
//       double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]; 
//       fprimaryShowerCharge += sh_hits[j]->Integral();
       fprimaryShower_hit_w.push_back(sh_hits[j]->WireID().Wire);
       fprimaryShower_hit_t.push_back(sh_hits[j]->PeakTime());
       fprimaryShower_hit_q.push_back(sh_hits[j]->Integral());        
       fprimaryShower_hit_X.push_back(detprop->ConvertTicksToX(sh_hits[j]->PeakTime(),sh_hits[j]->WireID().Plane,sh_hits[j]->WireID().TPC,0));
       fprimaryShower_hit_Z.push_back(xyzWire.Z());  
       fprimaryShower_hit_nElectrons.push_back(fCalorimetryAlg.ElectronsFromADCArea(sh_hits[j]->Integral(), sh_hits[j]->WireID().Plane));  
       fprimaryShower_hit_ch.push_back(wires[0]->Channel()); //only one channel per wire at collection plane
       std::vector<art::Ptr<recob::SpacePoint>> sp = spFromShowerHits.at(j); 

       if(!sp.empty() ){
         //fprimaryShower_hit_X[idx]= sp[0]->XYZ()[0];
         fprimaryShower_hit_Y.push_back(sp[0]->XYZ()[1]);
         //fprimaryShower_hit_Z[idx]= sp[0]->XYZ()[2];
       }
    
       idx ++;
    }
    fprimaryShower_nHits = idx; //only collection hits
    // Calorimetry only colleciton plane
    std::vector<anab::Calorimetry> calovector = showerUtil.GetRecoShowerCalorimetry(*thisShower, evt, fShowerTag, fShowerCaloTag);
    if(calovector.size() != 3 && fVerbose > 0)
      std::cerr << "WARNING::Calorimetry vector size for primary is = " << calovector.size() << std::endl;

  } // end is shower

}

// -----------------------------------------------------------------------------

void ProtoDUNEelectronWireAna::Initialize(){
  fRun =-999;
  fSubRun =-999;
  fevent =-999;

  fbeamtrackDir.clear();
  fprimaryStartDirection.clear();

  fprimaryShower_nHits =-999;
  fprimaryShower_hit_w.clear();
  fprimaryShower_hit_q.clear();
  fprimaryShower_hit_t.clear();
  fprimaryShower_hit_X.clear();
  fprimaryShower_hit_Y.clear();
  fprimaryShower_hit_Z.clear();
  fprimaryShower_hit_nElectrons.clear();
  fprimaryShower_hit_ch.clear();

  fcerenkovStatus.clear();

  fdet_frequency = -999.0;

  fbeamtrackMomentum = -999.0;
  fprimaryIstrack = -999;
  fprimaryIsshower = -999;

  fprimaryNHits =-999;

  fprimary_Shower_wire_w.clear();
  fprimary_Shower_wire_Y.clear();
  fprimary_Shower_wire_X.clear();
  fprimary_Shower_wire_Z.clear();
  fprimary_Shower_wire_ch.clear();
  fprimary_Shower_wire_ch_corr.clear();
//  fprimary_Shower_wire_time_charge.clear();
/*
  fprimary_Shower_simchannel_electrons_t1t2window.clear();
  fprimary_Shower_simchannel_simide_energyFrac.clear();
  fprimary_Shower_simchannel_simide_energy.clear();
  fprimary_Shower_simchannel_simide_numElectrons.clear();
*/
  fprimary_Shower_MCwire_w.clear();
  fprimary_Shower_MCwire_E.clear();
  fprimaryStartPosition.clear();

//  fNSimChannels =-999;
//  fNTrackIDEs =-999;


}
DEFINE_ART_MODULE(ProtoDUNEelectronWireAna)

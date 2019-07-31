/**
 * @file PixSimAnaTree_module.cc
 * @brief Analyzer module for * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

// Framework includes 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes 
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

// ROOT includes 
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"

// pixsim includes
#include "pixsim/Geometry/DetectorGeometryService.h"

namespace pixsim
{

class PixSimAnaTree : public art::EDAnalyzer 
{
public:
  explicit PixSimAnaTree(fhicl::ParameterSet const & p);
  virtual ~PixSimAnaTree();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void reconfigure(fhicl::ParameterSet const & p);

private:

  // Function used to reset all the variables  
  void ResetVars();
  
  // Storing information into TTree
  TTree* fTree;
   
  // Storing Run Information 
  int    run;			    ///< Run Number
  int    subrun;			///< SubRun Number
  int    event;			  ///< Event Number

  // Storing Geant4 MC Truth Information 
  int    geant_list_size;		///< Number of Geant4 particles tracked

  // Neutrino information
  std::vector<int> nu_pdg;
  std::vector<int> nu_ndau;
  std::vector<int> nu_ccnc;
  std::vector<int> nu_mode;
  std::vector<int> nu_inttype;
  std::vector<double> nu_StartPointx; 
  std::vector<double> nu_StartPointy;
  std::vector<double> nu_StartPointz;
  std::vector<double> nu_EndPointx; 
  std::vector<double> nu_EndPointy;
  std::vector<double> nu_EndPointz;
  std::vector<double> nu_Vertexx;
  std::vector<double> nu_Vertexy;
  std::vector<double> nu_Vertexz;
  std::vector<double> nu_StartPx;
  std::vector<double> nu_StartPy;
  std::vector<double> nu_StartPz;
  std::vector<double> nu_EndPx;
  std::vector<double> nu_EndPy;
  std::vector<double> nu_EndPz;

  std::vector<double> ides_x;
  std::vector<double> ides_y;
  std::vector<double> ides_z;
  std::vector<double> ides_energy;
  std::vector<double> ides_tid;
  std::vector<double> ides_numElectrons;

  // Other particle information
  std::vector<int>    TrackId;
  std::vector<int>    PDG;
  std::vector<double> StartEnergy;
  std::vector<double> EndEnergy;
  std::vector<double> StartPx;
  std::vector<double> StartPy;
  std::vector<double> StartPz;
  std::vector<double> EndPx;
  std::vector<double> EndPy;
  std::vector<double> EndPz;
  std::vector<double> StartPointx; 
  std::vector<double> StartPointy;
  std::vector<double> StartPointz;
  std::vector<double> EndPointx; 
  std::vector<double> EndPointy;
  std::vector<double> EndPointz;
  std::vector<int>    Mother;
  std::vector<std::vector<int>> Daughters;

  std::vector<int> isPrimary;	    
  std::vector<std::string> Process;

  std::string fTreeName;
  std::string fG4ModuleLabel;
  std::string fGenieModuleLabel;
  std::string fMCTrackModuleLabel;
};

//--------------------------------------------------------------------
PixSimAnaTree::PixSimAnaTree(fhicl::ParameterSet const & pset) 
 : EDAnalyzer{pset}
{
  this->reconfigure(pset);
}

//--------------------------------------------------------------------
PixSimAnaTree::~PixSimAnaTree()
{}

//--------------------------------------------------------------------
void PixSimAnaTree::reconfigure(fhicl::ParameterSet const & pset)
{
  fTreeName            = pset.get< std::string >("TreeName", "anatree");
  fG4ModuleLabel       = pset.get< std::string >("G4ModuleLabel", "largeant");
  fGenieModuleLabel    = pset.get< std::string >("GenieModuleLabel", "generator");
  fMCTrackModuleLabel = pset.get< std::string >("MCTrackModuleLabel", "mcreco");
  return;
}

//--------------------------------------------------------------------
void PixSimAnaTree::analyze(art::Event const & evt)
{
  // Reset variables
  ResetVars();

  // Detector properties service 
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // MC particle list
  art::Handle< std::vector<simb::MCParticle> > plistHandle;
  evt.getByLabel(fG4ModuleLabel, plistHandle);
  auto plist = *plistHandle;

  // MC truth 
  art::Handle< std::vector<simb::MCTruth> > tlistHandle;
  evt.getByLabel(fGenieModuleLabel, tlistHandle);
  auto tlist = *tlistHandle;

  // G truth 
  art::Handle< std::vector<simb::GTruth> > glistHandle;
  evt.getByLabel(fGenieModuleLabel, glistHandle);
  auto glist = *glistHandle;

  // Sim channels
  art::Handle< std::vector<sim::SimChannel> > SimListHandle;
  std::vector<art::Ptr<sim::SimChannel> > Simlist;
  if(evt.getByLabel(fG4ModuleLabel, SimListHandle))
       { art::fill_ptr_vector(Simlist, SimListHandle); }



  run    = evt.run();
  subrun = evt.subRun();
  event  = evt.id().event();
   
  std::cout<<std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout << "Run = "         << run 
            << ", SubRun = "    << subrun 
            << ", Evt = "       << event    << std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout<<std::endl;
   
  // Setting a boolian to only output MC info if this is MC-info 
  bool isdata = false;
  if (evt.isRealData()) {isdata = true;}
  else isdata = false;
   
  //
  // Filling Sim channels 
  //
  if (!isdata)
  {
    // Loop over channels
    for (int iCh = 0; iCh < (int)Simlist.size(); iCh++)
    {
      const auto& TDCIDEs = Simlist.at(iCh)->TDCIDEMap(); 
      for (const auto& TDCinfo : TDCIDEs)
      {
        for (const auto& ide : TDCinfo.second)
        {
          ides_tid.push_back(ide.trackID);
          ides_numElectrons.push_back(ide.numElectrons);
          ides_energy.push_back(ide.energy); 
          ides_x.push_back(ide.x);
          ides_y.push_back(ide.y);
          ides_z.push_back(ide.z);
        }
      }
    } 

    // Fill neutrino information
    if (tlist.size() > 1) throw cet::exception("PixSimAnaTree") << "MCTruth list greater than 1\n";
    if (glist.size() > 1) throw cet::exception("PixSimAnaTree") << "GTruth list greater than 1\n";
    for (size_t iT = 0; iT < tlist.size(); iT++)
    {
      simb::MCNeutrino nu   = tlist[iT].GetNeutrino();
      simb::MCParticle mcnu = nu.Nu(); 
    
      nu_pdg.push_back(mcnu.PdgCode());
      nu_ndau.push_back(mcnu.NumberDaughters());
      nu_ccnc.push_back(nu.CCNC());
      nu_mode.push_back(nu.Mode());
      nu_inttype.push_back(nu.InteractionType());
      nu_StartPointx.push_back(mcnu.Position().Vect().X());
      nu_StartPointy.push_back(mcnu.Position().Vect().Y());
      nu_StartPointz.push_back(mcnu.Position().Vect().Z());
      nu_EndPointx  .push_back(mcnu.EndPosition().Vect().X());
      nu_EndPointy  .push_back(mcnu.EndPosition().Vect().Y());
      nu_EndPointz  .push_back(mcnu.EndPosition().Vect().Z());
      nu_StartPx    .push_back(mcnu.Momentum().Vect().X());
      nu_StartPy    .push_back(mcnu.Momentum().Vect().Y());
      nu_StartPz    .push_back(mcnu.Momentum().Vect().Z());
      nu_EndPx      .push_back(mcnu.EndMomentum().Vect().X());
      nu_EndPy      .push_back(mcnu.EndMomentum().Vect().Y());
      nu_EndPz      .push_back(mcnu.EndMomentum().Vect().Z());
    }//<-- End loop over truth information

    for (size_t iT = 0; iT < tlist.size(); iT++)
    {
      simb::GTruth truth = glist[iT];
      nu_Vertexx.push_back(truth.fVertex.Vect().X());
      nu_Vertexy.push_back(truth.fVertex.Vect().Y());
      nu_Vertexz.push_back(truth.fVertex.Vect().Z());
    }
    
    int geant_particle(0);
    
    // Determine the number of primary particles from geant 
    for( unsigned int i = 0; i < plist.size(); ++i )
	  {
	    geant_particle++;
	  }//<---End i loop

    geant_list_size = geant_particle;
     
    // Looping over all the Geant4 particles 
    for(unsigned int i = 0; i < plist.size(); ++i)
    {
      // Check if primary
	    if(plist[i].Process()=="primary"){isPrimary.push_back(1);}
	    else                             {isPrimary.push_back(0);}
	   
      // Get the process
      Process.push_back(plist[i].Process());
    
	    // Saving other info 
      PDG.push_back(plist[i].PdgCode());
      Mother.push_back(plist[i].Mother());
      TrackId.push_back(plist[i].TrackId());
      StartEnergy.push_back(plist[i].E());
      EndEnergy.push_back(plist[i].EndE());

	    // Saving the start and end Px, Py, Pz info 
	    StartPx.push_back(plist[i].Px());
	    StartPy.push_back(plist[i].Py());
	    StartPz.push_back(plist[i].Pz());
	    EndPx.push_back(plist[i].EndPx());
	    EndPy.push_back(plist[i].EndPy());
	    EndPz.push_back(plist[i].EndPz());
	   
	    // Saving the Start and End Point for this particle 
	    StartPointx.push_back(plist[i].Vx());
	    StartPointy.push_back(plist[i].Vx());
	    StartPointz.push_back(plist[i].Vx());
	    EndPointx.push_back(plist[i].EndPosition()[0]);
	    EndPointy.push_back(plist[i].EndPosition()[1]);
	    EndPointz.push_back(plist[i].EndPosition()[2]);
 	   
	    // Saving the number of Daughters for this particle 
      std::vector<int> tempDtr;
      for(unsigned int iDtr = 0; iDtr < plist.size(); ++iDtr)
      {
        if (plist[iDtr].Mother() == plist[i].TrackId())
        {
          tempDtr.push_back(plist[iDtr].TrackId());
        }
      }
      Daughters.push_back(tempDtr);
    }//<--End loop on geant particles   
  }//<---End checking if this is MC

  //
  // Filling MC track info
  //
  if (!isdata)
  {
    art::Handle<std::vector<sim::MCTrack>> mctrackh;
    evt.getByLabel(fMCTrackModuleLabel, mctrackh);

  
    // Check to make sure the MCShower Handle is valid 
    if (mctrackh.isValid())
    {
      //  Looping over all MC Showers (using uBooNEAna method) 
      for (std::vector<sim::MCTrack>::const_iterator imctrk = mctrkh->begin(); imctrk != mctrkh->end(); ++imctrk)
      {

        const sim::MCTrack &mctrk = *imctrk;

        mctrk_origin.push_back(mctrk.Origin());
        mctrk_pdg.push_back(mctrk.PdgCode());
        mctrk_TrackId.push_back(mctrk.TrackID());
        mctrk_startX.push_back(mctrk.Start().X());
        mctrk_startY.push_back(mctrk.Start().Y());
        mctrk_startZ.push_back(mctrk.Start().Z());
        mctrk_endX.push_back(mctrk.End().X());
        mctrk_endY.push_back(mctrk.End().Y());
        mctrk_endZ.push_back(mctrk.End().Z());

        mctrk_Motherpdg.push_back(mctrk.MotherPdgCode());
        mctrk_MotherTrkId .push_back(mctrk.MotherTrackID());
        mctrk_MotherstartX.push_back(mctrk.MotherStart().X());
        mctrk_MotherstartY.push_back(mctrk.MotherStart().Y());
        mctrk_MotherstartZ.push_back(mctrk.MotherStart().Z());
        mctrk_MotherendX.push_back(mctrk.MotherEnd().X());
        mctrk_MotherendY.push_back(mctrk.MotherEnd().Y());
        mctrk_MotherendZ.push_back(mctrk.MotherEnd().Z());
        mctrk_Ancestorpdg.push_back(mctrk.AncestorPdgCode());
        mctrk_AncestorTrkId .push_back(mctrk.AncestorTrackID());
        mctrk_AncestorstartX.push_back(mctrk.AncestorStart().X());
        mctrk_AncestorstartY.push_back(mctrk.AncestorStart().Y());
        mctrk_AncestorstartZ.push_back(mctrk.AncestorStart().Z());
        mctrk_AncestorendX.push_back(mctrk.AncestorEnd().X());
        mctrk_AncestorendY.push_back(mctrk.AncestorEnd().Y());
        mctrk_AncestorendZ.push_back(mctrk.AncestorEnd().Z());
      } //<---End imctrk iterator loop

    } //<---Only going in if the handle is valid
  }   //<---End only looking at MC information

  fTree->Fill();
}//<---End analyze()

//--------------------------------------------------------------------
void PixSimAnaTree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str(),  fTreeName.c_str());
  fTree->Branch("run",                         &run,"run/I");
  fTree->Branch("subrun",                      &subrun,"subrun/I");
  fTree->Branch("event",                       &event,"event/I");

  fTree->Branch("geant_list_size",             &geant_list_size,"geant_list_size/I");
  fTree->Branch("PDG",                         &PDG);
  fTree->Branch("StartEnergy",                 &StartEnergy);
  fTree->Branch("EndEnergy",                   &EndEnergy);
  fTree->Branch("StartPx",                     &StartPx);
  fTree->Branch("StartPy",                     &StartPy);
  fTree->Branch("StartPz",                     &StartPz);
  fTree->Branch("EndPx",                       &EndPx);
  fTree->Branch("EndPy",                       &EndPy);
  fTree->Branch("EndPz",                       &EndPz);
  fTree->Branch("StartPointx",                 &StartPointx);
  fTree->Branch("StartPointy",                 &StartPointy);
  fTree->Branch("StartPointz",                 &StartPointz);
  fTree->Branch("EndPointx",                   &EndPointx);
  fTree->Branch("EndPointy",                   &EndPointy);
  fTree->Branch("EndPointz",                   &EndPointz);
  fTree->Branch("Process",                     &Process);
  fTree->Branch("Mother",                      &Mother);
  fTree->Branch("TrackId",                     &TrackId);
  fTree->Branch("isPrimary",                   &isPrimary); 

  fTree->Branch("nu_pdg",     &nu_pdg        );
  fTree->Branch("nu_ndau",    &nu_ndau       );
  fTree->Branch("nu_ccnc",    &nu_ccnc       );
  fTree->Branch("nu_mode",    &nu_mode       );
  fTree->Branch("nu_inttype", &nu_inttype    );

  fTree->Branch("nu_StartPointx",  &nu_StartPointx);
  fTree->Branch("nu_StartPointy",  &nu_StartPointy);
  fTree->Branch("nu_StartPointz",  &nu_StartPointz);
  fTree->Branch("nu_EndPointx",    &nu_EndPointx);
  fTree->Branch("nu_EndPointy",    &nu_EndPointy);
  fTree->Branch("nu_EndPointz",    &nu_EndPointz);
  fTree->Branch("nu_Vertexx",      &nu_Vertexx);
  fTree->Branch("nu_Vertexy",      &nu_Vertexy);
  fTree->Branch("nu_Vertexz",      &nu_Vertexz);
  fTree->Branch("nu_StartPx",      &nu_StartPx);
  fTree->Branch("nu_StartPy",      &nu_StartPy);
  fTree->Branch("nu_StartPz",      &nu_StartPz);
  fTree->Branch("nu_EndPx",        &nu_EndPx);
  fTree->Branch("nu_EndPy",        &nu_EndPy);
  fTree->Branch("nu_EndPz",        &nu_EndPz);

  fTree->Branch("ides_x", &ides_x);
  fTree->Branch("ides_y", &ides_y);
  fTree->Branch("ides_z", &ides_z);
  fTree->Branch("ides_energy", &ides_energy);
  fTree->Branch("ides_tid", &ides_tid);
  fTree->Branch("ides_numElectrons", &ides_numElectrons);
}

//--------------------------------------------------------------------
void PixSimAnaTree::ResetVars()
{
  TrackId.clear();
  PDG.clear();         
  StartEnergy.clear(); 
  EndEnergy.clear();   
  StartPx.clear();     
  StartPy.clear();     
  StartPz.clear();     
  EndPx.clear();       
  EndPy.clear();       
  EndPz.clear();       
  StartPointx.clear(); 
  StartPointy.clear(); 
  StartPointz.clear(); 
  EndPointx.clear();   
  EndPointy.clear();   
  EndPointz.clear();   
  Mother.clear();      
  Daughters.clear();   
  isPrimary.clear();	  
  Process.clear();     
  
  ides_x.clear();
  ides_y.clear();
  ides_z.clear();
  ides_energy.clear();
  ides_tid.clear();
  ides_numElectrons.clear();

  nu_pdg.clear();
  nu_ndau.clear();
  nu_ccnc.clear();
  nu_mode.clear();
  nu_inttype.clear();
  nu_StartPointx.clear(); 
  nu_StartPointy.clear();
  nu_StartPointz.clear();
  nu_EndPointx.clear(); 
  nu_EndPointy.clear();
  nu_EndPointz.clear();
  nu_Vertexx.clear();
  nu_Vertexy.clear();
  nu_Vertexz.clear();
  nu_StartPx.clear();
  nu_StartPy.clear();
  nu_StartPz.clear();
  nu_EndPx.clear();
  nu_EndPy.clear();
  nu_EndPz.clear();

  run = -99999;
  subrun = -99999;
  event = -99999;

  geant_list_size=-999;
}

}

DEFINE_ART_MODULE(pixsim::PixSimAnaTree)

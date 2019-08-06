////////////////////////////////////////////////////////////////////////
// Class:       EmptyEventFilter
// Plugin Type: filter (art v3_02_06)
// File:        EmptyEventFilter_module.cc
//
// Generated at Wed Jul 31 20:27:15 2019 by Hunter Sullivan using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/Simulation/SimChannel.h"

#include <memory>

namespace pixsim
{

class EmptyEventFilter;


class EmptyEventFilter : public art::EDFilter {
public:
  explicit EmptyEventFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EmptyEventFilter(EmptyEventFilter const&) = delete;
  EmptyEventFilter(EmptyEventFilter&&) = delete;
  EmptyEventFilter& operator=(EmptyEventFilter const&) = delete;
  EmptyEventFilter& operator=(EmptyEventFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


EmptyEventFilter::EmptyEventFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
}

bool EmptyEventFilter::filter(art::Event& evt)
{
  // Sim channels
  art::Handle< std::vector<sim::SimChannel> > SimListHandle;
  std::vector<art::Ptr<sim::SimChannel> > Simlist;
  if(evt.getByLabel("largeant", SimListHandle))
       { art::fill_ptr_vector(Simlist, SimListHandle); }

  // Setting a boolian to only output MC info if this is MC-info 
  bool isdata = false;
  if (evt.isRealData()) {isdata = true;}
  else isdata = false;
   
  //
  // Filling Sim channels 
  //
  std::vector<int> ides_tid;
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
        }
      }
    }
    if (ides_tid.size() > 1000) return true;
    else return false;
  }
  else return true; 
}

void EmptyEventFilter::beginJob()
{
  // Implementation of optional member function here.
}

void EmptyEventFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(EmptyEventFilter)
}

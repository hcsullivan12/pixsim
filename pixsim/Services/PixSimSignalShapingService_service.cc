////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceT1034_service.cc
/// \author H. Greenlee 
////////////////////////////////////////////////////////////////////////

#include "pixsim/Services/PixSimSignalShapingService.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/LArFFT.h"
#include "TFile.h"
#include <fstream>

//----------------------------------------------------------------------
// Constructor.
util::PixSimSignalShapingService::PixSimSignalShapingService(const fhicl::ParameterSet& pset,
                                                           art::ActivityRegistry& /* reg */)
: fInit(false)
{
  reconfigure(pset);
}


//----------------------------------------------------------------------
// Destructor.
util::PixSimSignalShapingService::~PixSimSignalShapingService()
{
}


//----------------------------------------------------------------------
// Reconfigure method.
void util::PixSimSignalShapingService::reconfigure(const fhicl::ParameterSet& pset)
{
  // Reset initialization flag.
  fInit = false;
  
  // Reset kernels.
  fColSignalShaping.Reset();
} 

//--------------------------------------------------------------------
const util::SignalShaping&
util::PixSimSignalShapingService::SignalShaping() const
{
  if(!fInit) init();

  // We always return collection
  return fColSignalShaping;
}

//----------------------------------------------------------------------
void util::PixSimSignalShapingService::init()
{
  if (fInit) return;
  fInit = true;

  SetElectResponse();

  // Configure convolution kernels.
  fColSignalShaping.AddResponseFunction(fElectResponse);
  fColSignalShaping.save_response();
  fColSignalShaping.set_normflag(false);
}

//----------------------------------------------------------------------
double util::PixSimSignalShapingService::ElectronsToCurrent(double const& el) const
{
  auto const *detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>{}->provider();
  // Convert and scale to pA
  // Sampling rate in ns
  // qConv = 1.6e-19 * 1e12 * 1e9
  double qConv = 160.;
  return el*qConv/detprop->SamplingRate();
}

//----------------------------------------------------------------------
void util::PixSimSignalShapingService::SetElectResponse()
{
  art::ServiceHandle<util::LArFFT> fft;
  int nticks = fft->FFTSize();
  fElectResponse.resize(nticks, 0.);

  double max(0);

  // Let's take our response function to be gaussian
  TF1 g("electgaus", "gaus");
  g.SetParameters(1.0, 0.0, 100.0);
  float time = -3*g.GetParameter(2);

  for (size_t i = 0; i < fElectResponse.size(); i++)
  {
    time += i;
    fElectResponse[i] = g.Eval(time);
  }
}

namespace util {

  DEFINE_ART_SERVICE(PixSimSignalShapingService)

}

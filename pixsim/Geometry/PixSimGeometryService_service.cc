//////////////////////////////////////////////////////////////////////
/// \file  PixSimGeometryService_service.cxx
/// \brief Service for pixel geometry information.
///
/// \author  hsulliva@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "PixSimGeometryService.h"

namespace pixgeo 
{

PixSimGeometryService::PixSimGeometryService
  (fhicl::ParameterSet const& pset, art::ActivityRegistry&)
  {
    fProvider = std::make_unique<pixgeo::PixSimGeometry>(pset);
  }


void PixSimGeometryService::reconfigure(fhicl::ParameterSet const& pset)
{

}
}
DEFINE_ART_SERVICE_INTERFACE_IMPL(pixgeo::PixSimGeometryService, geo::DetectorGeometryService)


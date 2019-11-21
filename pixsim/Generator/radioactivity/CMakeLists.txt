art_make(LIB_LIBRARIES nusimdata_SimulationBase
                       cetlib
                       cetlib_except
                       ${FHICLCPP}
                       ${CLHEP}
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${ROOT_CORE}
                       ${ROOT_RIO}
                       ${ROOT_HIST}
                       ${ROOT_PHYSICS}
          MODULE_LIBRARIES larsim_PhotonPropagation
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           larsim_Simulation
                           lardataobj_Simulation
                           larcoreobj_SummaryData
                           pixsim_Geometry
                           nurandom_RandomUtils_NuRandomService_service
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Provenance
                           ${MF_MESSAGELOGGER}
                           ${CLHEP}
                           ${ROOT_EG}
                           ${ROOT_GEOM}
        )

install_headers()
install_fhicl()
install_source()

#install headers
install(FILES iga.hh
              bsplinebasis.hh
              bsplinepatch.hh
              vtkfile.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga)

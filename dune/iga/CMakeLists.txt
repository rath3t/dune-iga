#install headers
file(
        GLOB_RECURSE headers
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        CONFIGURE_DEPENDS *.hh
)
install(FILES ${headers}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga)

add_subdirectory(test)
cmake_minimum_required( VERSION 2.6 )
find_package( brainvisa-cmake REQUIRED )
BRAINVISA_PROJECT()

BRAINVISA_FIND_PACKAGE( aims-free REQUIRED )

include( "${AIMS-FREE_USE_FILE}" )

add_subdirectory( src )

find_package( Doxygen )
IF ( DOXYGEN_FOUND )
  set( _doxygenInput "${CMAKE_BINARY_DIR}/include/cortical_surface"  )
  set(aims_version "${aims-free_VERSION_MAJOR}.${aims-free_VERSION_MINOR}")
  set( DOXYFILE_TAGFILES "${CMAKE_BINARY_DIR}/share/doc/aimsalgo-${aims_version}/doxygen/aimsalgo.tag=../../aimsalgo-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/aimsdata-${aims_version}/doxygen/aimsdata.tag=../../aimsdata-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/graph-${aims_version}/doxygen/graph.tag=../../graph-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/cartobase-${aims_version}/doxygen/cartobase.tag=../../cartobase-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/cartodata-${aims_version}/doxygen/cartodata.tag=../../cartodata-${aims_version}/doxygen")
  BRAINVISA_GENERATE_DOXYGEN_DOC( _doxygenInput
                                  COMPONENT "cortical_surface" )
ENDIF ( DOXYGEN_FOUND )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()


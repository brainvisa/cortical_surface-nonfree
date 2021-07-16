cmake_minimum_required( VERSION 2.6 )
find_package( brainvisa-cmake REQUIRED )
BRAINVISA_PROJECT()

BRAINVISA_FIND_PACKAGE( aims-free REQUIRED )

include( "${AIMS-FREE_USE_FILE}" )

math( EXPR result "${aims-free_VERSION_MINOR}+1" )
BRAINVISA_DEPENDENCY( RUN DEPENDS "aims-free" RUN ">= ${aims-free_VERSION};<< ${aims-free_VERSION_MAJOR}.${result}" )

# don't warn on obsolete classes within this library: we still provide obsolete
# features so we have to build them without warnings
add_definitions( "-DAIMSDATA_CLASS_NO_DEPREC_WARNING=1" )

add_subdirectory( src )

find_package( Doxygen )
IF ( DOXYGEN_FOUND )
  set( _doxygenInput ${HEADERS} )
  set(aims_version "${aims-free_VERSION_MAJOR}.${aims-free_VERSION_MINOR}")
  set( DOXYFILE_TAGFILES "${CMAKE_BINARY_DIR}/share/doc/aimsalgo-${aims_version}/doxygen/aimsalgo.tag=../../aimsalgo-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/aimsdata-${aims_version}/doxygen/aimsdata.tag=../../aimsdata-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/graph-${aims_version}/doxygen/graph.tag=../../graph-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/cartobase-${aims_version}/doxygen/cartobase.tag=../../cartobase-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/cartodata-${aims_version}/doxygen/cartodata.tag=../../cartodata-${aims_version}/doxygen")
  BRAINVISA_GENERATE_DOXYGEN_DOC( _doxygenInput
                                  INPUT_PREFIX "${CMAKE_BINARY_DIR}/include/cortical_surface" 
                                  COMPONENT "cortical_surface" )
  add_dependencies(cortical_surface-doxygen aimsalgo-doxygen graph-doxygen cartobase-doxygen cartodata-doxygen )
ENDIF ( DOXYGEN_FOUND )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()


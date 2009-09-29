cmake_minimum_required( VERSION 2.6 )
find_package( brainvisa-cmake REQUIRED )
BRAINVISA_PROJECT()

BRAINVISA_FIND_PACKAGE( aims-free REQUIRED )

include( "${AIMS-FREE_USE_FILE}" )

add_subdirectory( src )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()


TEMPLATE = lib
TARGET = cortical_surface${BUILDMODEEXT}

#!include ../../config-cpp-library

HEADERS = \
  mesh/isoLine.h \
  mesh/median.h \
  mesh/pointDistance.h \
  mesh/linkPath.h \
  mesh/meshToMeshResample.h \
  surfacereferential/corticalConstraints.h \
  surfacereferential/corticalReferential.h \
  surfacereferential/corticalTools.h \
  surfacereferential/shortestPath.h \
  surfacereferential/sulcusCleaner.h \
  surfacereferential/sulcusCorticalSnake.h \
  surfacereferential/sulcusCorticalSnake_energy.h \


SOURCES = \
  mesh/isoLine.cc \
  mesh/median.cc \
  mesh/pointDistance.cc \
  surfacereferential/corticalConstraints.cc \
  surfacereferential/corticalReferential.cc \
  surfacereferential/corticalTools.cc \
  surfacereferential/sulcusCleaner.cc \
  surfacereferential/sulcusCorticalSnake.cc \
  surfacereferential/sulcusCorticalSnake_energy.cc \


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
  surfacereferential/conformalMapping.h \
  structuralanalysis/iograph.h \
  structuralanalysis/meshdistance.h \
  structuralanalysis/sites.h \
  structuralanalysis/cliques.h \
  structuralanalysis/minimization.h \
  structuralanalysis/icm.h \
  structuralanalysis/anneal.h \
  structuralanalysis/old_anneal.h \
  structuralanalysis/cluster.h \
  structuralanalysis/validation.h \


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
  surfacereferential/conformalMapping.cc \
  structuralanalysis/iograph.cc \
  structuralanalysis/meshdistance.cc \
  structuralanalysis/sites.cc \
  structuralanalysis/cliques.cc \
  structuralanalysis/minimization.cc \
  structuralanalysis/icm.cc \
  structuralanalysis/anneal.cc \
  structuralanalysis/old_anneal.cc \
  structuralanalysis/cluster.cc \
  structuralanalysis/validation.cc

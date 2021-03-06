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
  surfacereferential/sulcalLinesGeodesic.h \
  surfacereferential/gyri/basictex_operations.h \
  surfacereferential/gyri/compconn_operations.h \
  surfacereferential/gyri/constraints_operations.h \
  surfacereferential/gyri/gyri_operations.h \
  surfacereferential/gyri/intersec_operations.h \
  surfacereferential/gyri/mesh_operations.h \
  surfacereferential/gyri/misctex_operations.h \
  surfacereferential/gyri/model_operations.h \
  surfacereferential/gyri/vector_operations.h \
  surfacereferential/gyri/vertices_operations.h \
  surfacereferential/gyri/verif_operations.h \
  surfacereferential/gyri/gyri_parameterization.h \
  structuralanalysis/iograph.h \
  structuralanalysis/sites.h \
  structuralanalysis/cliques.h \
  structuralanalysis/minimization.h \
  structuralanalysis/icm.h \
  structuralanalysis/anneal.h \
  structuralanalysis/cluster.h \
  structuralanalysis/validation.h \
  structuralanalysis/representation.h \
  structuralanalysis/blobs.h \
  structuralanalysis/subjectdata.h \
  structuralanalysis/region.h \
  structuralanalysis/texturetoblobs.h \ 
  structuralanalysis/meshdistance.h \
  structuralanalysis/primalsketch_operations.h 


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
  surfacereferential/sulcalLinesGeodesic.cc \
  surfacereferential/gyri/basictex_operations.cc \
  surfacereferential/gyri/compconn_operations.cc \
  surfacereferential/gyri/constraints_operations.cc \
  surfacereferential/gyri/gyri_operations.cc \
  surfacereferential/gyri/intersec_operations.cc \
  surfacereferential/gyri/mesh_operations.cc \
  surfacereferential/gyri/misctex_operations.cc \
  surfacereferential/gyri/model_operations.cc \
  surfacereferential/gyri/vector_operations.cc \
  surfacereferential/gyri/verif_operations.cc \
  surfacereferential/gyri/vertices_operations.cc \
  surfacereferential/gyri/gyri_parameterization.cc \
  structuralanalysis/iograph.cc \
  structuralanalysis/cliques.cc \
  structuralanalysis/minimization.cc \
  structuralanalysis/icm.cc \
  structuralanalysis/anneal.cc \
  structuralanalysis/cluster.cc \
  structuralanalysis/validation.cc \
  structuralanalysis/representation.cc \
  structuralanalysis/blobs.cc \
  structuralanalysis/subjectdata.cc \
  structuralanalysis/region.cc \
  structuralanalysis/texturetoblobs.cc \
  structuralanalysis/meshdistance.cc \
  structuralanalysis/primalsketch_operations.cc 
  

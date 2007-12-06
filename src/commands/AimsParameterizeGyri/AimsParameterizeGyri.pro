TEMPLATE = app
TARGET = AimsParameterizeGyri

#!include ../../../config-cpp-command

SOURCES =       \
      main.cc  \
	   vector_operations.cc  \
      basictex_operations.cc  \
      model_operations.cc  \
      vertices_operations.cc  \
      constraints_operations.cc  \
      gyri_operations.cc  \
      intersec_operations.cc  \
      compconn_operations.cc  \
      mesh_operations.cc  \
      verif_operations.cc  \
      misctex_operations.cc  \
      gyri_parameterization.cc
HEADERS =   \
      vector_operations.h  \
      basictex_operations.h  \
      model_operations.h  \
      vertices_operations.h  \
      constraints_operations.h  \
      gyri_operations.h  \
      intersec_operations.h  \
      compconn_operations.h  \
      mesh_operations.h  \
      verif_operations.h  \
      misctex_operations.h  \
      gyri_parameterization.h



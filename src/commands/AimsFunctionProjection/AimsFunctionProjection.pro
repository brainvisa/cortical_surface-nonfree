TEMPLATE = app
TARGET = AimsFunctionProjection

#!include ../../../config-cpp-command

SOURCES =  createKernels.cc     \
           projectOntoMesh.cc    \
      main.cc
      
HEADERS = createKernels.h  \
         projectOntoMesh.h 


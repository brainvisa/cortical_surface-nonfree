TEMPLATE = app
TARGET = AimsCheck2DCoordinates

#!include ../../../config-cpp-command

SOURCES = \
          main.cc \
          mesh_operations.cc \
          vector_operations.cc 

INCLUDES = \
          mesh_operations.h \
          vector_operations.h 

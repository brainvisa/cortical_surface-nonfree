TEMPLATE = app
TARGET = surfStructuralAnalysis

#!include ../../../config-cpp-command


SOURCES = iograph.cc \
          meshdistance.cc \
          sites.cc \
          cliques.cc \
          minimization.cc \
          icm.cc \
          anneal.cc \
          old_anneal.cc \
          cluster.cc \
          main.cc

HEADERS = iograph.h \
          meshdistance.h \
          sites.h \
          cliques.h \
          minimization.h \
          icm.h \
          anneal.h \
          old_anneal.h \
          cluster.h
          
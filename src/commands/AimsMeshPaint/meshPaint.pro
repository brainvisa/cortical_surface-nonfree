TEMPLATE = app
TARGET = AimsMeshPaint

#!include ../../../config-cpp-command

QT = core \
    gui \
    opengl

HEADERS = glwidget.h \
    meshpaint.h \
    vector.h \
    trackball.h
          
SOURCES = glwidget.cc \
    meshpaint.cc \
    trackball.cc \
    main.cc

INCLUDEPATH += /usr/lib/qt4/mkspecs/linux-g++ \
/usr/lib/qt4/include/QtCore \
/usr/lib/qt4/include/QtGui \
/usr/lib/qt4/include/QtOpenGL \
/usr/lib/qt4/include

LIBS += -lQtOpenGL -lQtGui -lGLU -lGL -lanatomist
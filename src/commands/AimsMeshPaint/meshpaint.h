#ifndef MESHPAINT_H
#define MESHPAINT_H

#include <QMainWindow>
#include <QtGui>
#include "glwidget.h"
#include <cstdlib>
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <string.h>
#include <aims/mesh/surfaceOperation.h>
#include <anatomist/application/globalConfig.h>
#include <anatomist/application/settings.h>
#include <QMainWindow>

using namespace anatomist;
using namespace aims;
using namespace carto;
using namespace std;

class MeshPaint : public QMainWindow
{
    Q_OBJECT

public:
    MeshPaint(AimsSurfaceTriangle mesh,TimeTexture<float> tex,string colorMap);
    ~MeshPaint();

private slots:
    void paintBrush();
    void colorPicker();
    void trackball();
    void keyPressEvent( QKeyEvent* _event );
private:
    void createActions();
    void createToolBars();

    AimsSurfaceTriangle _mesh;
    TimeTexture<float> _tex;
    string _colorMap;

    GLWidget *glWidget;
    QToolBar *paintToolBar;
    QAction *colorPickerAction;
    QAction *paintBrushAction;
    QAction *trackballAction;
};

#endif // MESHPAINT_H

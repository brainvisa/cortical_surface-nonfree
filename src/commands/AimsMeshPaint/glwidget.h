#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QBrush>
#include <QFont>
#include <QImage>
#include <QPen>
#include <QGLWidget>
#include <QTime>
#include <QWheelEvent>
#include <iostream>
#include <QtGui>
#include <QDrag>
#include <QtOpenGL>
#include <stdlib.h>
#include <math.h>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/def/path.h>
#include <aims/def/general.h>
#include <cartobase/stream/fileutil.h>
#include <cartobase/stream/directory.h>
#include <cartobase/config/paths.h>
#include <cartobase/config/version.h>
#include <aims/config/aimsdata_config.h>
#include <cmath>
#include <anatomist/surface/texture.h>
#include <aims/utility/converter_texture.h>
#include <float.h>
#include "trackball.h"


using namespace carto;
using namespace aims;
using namespace std;
using namespace anatomist;

QT_BEGIN_NAMESPACE
class QPaintEvent;
class QWidget;
class QDragEnterEvent;
class QDropEvent;
QT_END_NAMESPACE

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget (QWidget *parent) : QGLWidget(QGLFormat(QGL::SampleBuffers), parent){};
    ~GLWidget (){};
};

template<typename T>
class myGLWidget : public GLWidget
{
public:
    myGLWidget (QWidget *parent, string adressTexIn,string adressMeshIn,string adressTexOut,string colorMap,string dataType);
    ~myGLWidget ();

    int getMode () const { return _mode; }
    float getZoom () const { return _zoom; }
    float getTranslate () const { return _trans; }

    void changeMode (int mode);

    void unitize (AimsSurfaceTriangle as, Point3df *meshCenter, float *meshScale);
    int buildDisplayList (AimsSurfaceTriangle as,int mode);
    int computeNearestVertexFromPolygonPoint (Point3df position, int poly, AimsSurfaceTriangle as);
    void computeIndexColor (AimsSurfaceTriangle as);
    void copyBackBuffer2Texture (void);

    GLuint loadColorMap (const char * filename);
    void drawColorMap (void);
    void keyPressEvent (QKeyEvent *event);

    void setZoom(float z);
    void setTranslate(float t);
//
protected:
    void initializeGL ();
    void paintGL ();
    void resizeGL (int width, int height);
    void mousePressEvent (QMouseEvent *event);
    void mouseMoveEvent (QMouseEvent *event);
    void mouseReleaseEvent (QMouseEvent *event);
    void wheelEvent (QWheelEvent *event);
    void dragEnterEvent (QDragEnterEvent *event);
    void dragMoveEvent (QDragMoveEvent *event);
    void dropEvent (QDropEvent *event);
    QPointF pixelPosToViewPos (const QPointF& p);
    void drawScenetoBackBuffer (void);

    GLuint checkIDpolygonPicked (int x, int y);
    Point3df check3DpointPicked (int x, int y);
    void projectionPerspective (void);
    void projectionOrtho (void);
    void trackBallTransformation(void);
    void drawPrimitivePicked (void);

private:
//    void drawInfos (QPainter *painter, string t);
    void setupViewport (int width, int height);
//
    int _mode;
    float _zoom ;
    float _trans;
    bool _resized;

    GLUquadricObj *_qobj_cursor;

    int _frames;
    QTime _time;

    string _dataType;
    AimsSurfaceTriangle _mesh;
    TimeTexture<T> _tex;

    string _adressTexIn;
    string _adressMeshIn;
    string _adressTexOut;
    string _colorMap;

    TrackBall _trackBall;

    GLuint _listMeshRender;
    GLuint _listMeshPicking;
    GLuint _listMeshParcelation;
    GLuint _listMeshSmooth;

    Point3df _meshCenter;
    float _meshScale;

    std::vector<int> _indexTexture;
    GLubyte *backBufferTexture;

    unsigned char* dataColorMap;

    Point3df _point3Dpicked;
    Point3df _vertexNearestpicked;
    GLuint _indexPolygon;
    GLuint _indexVertex;
    float _textureValue;
    bool _wireframe;
    bool _parcelation;

    GLuint _IDcolorMap;
};

#endif

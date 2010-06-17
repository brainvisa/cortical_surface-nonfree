#include <iostream>
#include <QtGui>
#include <QDrag>
#include <QtOpenGL>
#include <stdlib.h>
#include <math.h>
#include "glwidget.h"

/* enums */
enum
{
  X, Y, Z, W
};

GLWidget::GLWidget(QWidget *parent, AimsSurfaceTriangle as,TimeTexture<float> tex, string colorMap) :
QGLWidget(QGLFormat(QGL::SampleBuffers), parent), _mesh(as), _tex(tex), _colorMap(colorMap)
{
  _zoom = -2.0;
  _trans = 0.0;
  _mode = 1;
  _listMeshRender = 0;
  _listMeshPicking = 0;
  _indexVertex = 0;
  _indexPolygon = 0;

  _wireframe = false;

  //_trackBall = TrackBall(0.05f, gfx::Vector3f::vector(0, 1, 0), TrackBall::Sphere);
  _trackBall = TrackBall(0.0f, gfx::Vector3f::vector(0, 1, 0),TrackBall::Plane);
  setAutoFillBackground(false);
  setMinimumSize(800, 600);
  setWindowTitle(tr("painting on the mesh"));
  setAutoBufferSwap(false);
  setAcceptDrops(true);
}

GLWidget::~GLWidget()
{
}

void GLWidget::unitize(AimsSurfaceTriangle as, Point3df *meshCenter, float *meshScale)
{
  GLfloat maxx, minx, maxy, miny, maxz, minz;
  GLfloat w, h, d;

  const vector<Point3df> & vert = as.vertex();
  /* get the max/mins */
  maxx = minx = vert[0][X];
  maxy = miny = vert[0][Y];
  maxz = minz = vert[0][Z];

  for (int j = 0; j < (int) vert.size(); ++j)
  {
    if (maxx < vert[j][X])
      maxx = vert[j][X];
    if (minx > vert[j][X])
      minx = vert[j][X];
    if (maxy < vert[j][Y])
      maxy = vert[j][Y];
    if (miny > vert[j][Y])
      miny = vert[j][Y];
    if (maxz < vert[j][Z])
      maxz = vert[j][Z];
    if (minz > vert[j][Z])
      minz = vert[j][Z];
  }
  /* calculate mesh width, height, and depth */
  w = fabs(maxx) + fabs(minx);
  h = fabs(maxy) + fabs(miny);
  d = fabs(maxz) + fabs(minz);
  /* calculate center of the mesh */
  (*meshCenter)[X] = (maxx + minx) / 2.0;
  (*meshCenter)[Y] = (maxy + miny) / 2.0;
  (*meshCenter)[Z] = (maxz + minz) / 2.0;
  /* calculate unitizing scale factor */
  *meshScale = 2.0 / max(max(w, h), d);
  /* translate around center then scale */
}

int GLWidget::buildDisplayList(AimsSurfaceTriangle as, int mode)
{
  int j;
  GLuint list = 0;
  list = glGenLists(1);
  glNewList(list, GL_COMPILE);

  vector<Point3df> & vert = as.vertex();
  const vector<Point3df> & norm = as.normal();

  vector<AimsVector<uint, 3> > & tri = as.polygon();

  rc_ptr<Texture1d>	tex( new Texture1d );
  Converter<TimeTexture<float>, Texture1d>	c;
  c.convert( _tex, *tex );
  ATexture	*ao = new ATexture;
  ao->setTexture( tex );
  ao->normalize();

  const float* t = ao->textureCoords();

  for (j = 0; j < (int) vert.size(); ++j)
    {
      vert[j][X] -= (_meshCenter)[X];
      vert[j][Y] -= (_meshCenter)[Y];
      vert[j][Z] -= (_meshCenter)[Z];
      vert[j][X] *= (_meshScale);
      vert[j][Y] *= (_meshScale);
      vert[j][Z] *= (_meshScale);
    }

  glBegin( GL_TRIANGLES);
  for (j = 0; j < (int) tri.size(); ++j)
  {
	if (mode)
	  glColor3ub(_indexTexture[3 * j],
			  _indexTexture[3 * j + 1],
				  _indexTexture[3 * j + 2]);
	else
      glColor3ub(255, 255, 255);

	glTexCoord2d(t[tri[j][0]],0);
    glNormal3f(norm[tri[j][0]][0], norm[tri[j][0]][1], norm[tri[j][0]][2]);
    glVertex3f(vert[tri[j][0]][0], vert[tri[j][0]][1], vert[tri[j][0]][2]);

    glTexCoord2d(t[tri[j][1]],0);
    glNormal3f(norm[tri[j][1]][0], norm[tri[j][1]][1], norm[tri[j][1]][2]);
    glVertex3f(vert[tri[j][1]][0], vert[tri[j][1]][1], vert[tri[j][1]][2]);

    glTexCoord2d(t[tri[j][2]],0);
    glNormal3f(norm[tri[j][2]][0], norm[tri[j][2]][1], norm[tri[j][2]][2]);
    glVertex3f(vert[tri[j][2]][0], vert[tri[j][2]][1], vert[tri[j][2]][2]);
  }

  glEnd();
  glEndList();

  cout << "ID list created = " << list << endl;

//  _qobj_cursor = gluNewQuadric();
//  gluQuadricDrawStyle(_qobj_cursor, GLU_FILL);
//  gluQuadricNormals(_qobj_cursor, GLU_SMOOTH);

  return list;
}

void GLWidget::setZoom(float z)
{
  _zoom = z;
  updateGL();
}

void GLWidget::setTranslate(float t) {
	_trans = t;
	updateGL();
}

void GLWidget::changeMode(int mode) {
	_mode = mode;
}

void GLWidget::initializeGL()
{
  glEnable( GL_LIGHTING);
  glEnable( GL_LIGHT0);
  glEnable( GL_DEPTH_TEST);
  glShadeModel( GL_SMOOTH);
  static GLfloat lightPosition[4] = { 0.5, 5.0, 7.0, 1.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
  glClearColor(1, 1, 1, 1);

  unitize(_mesh, &_meshCenter, &_meshScale);
  cout << "mesh scale = " <<  _meshScale << "\n";
  cout << "mesh center = (" <<  _meshCenter[0] << "," << _meshCenter[1] << "," << _meshCenter[2] << ")\n";

  computeIndexColor(_mesh);

  loadColorMap(_colorMap.c_str());

  cout << "nb triangle = " << _mesh.polygon().size() << endl;
  cout << "nb vertex = " << _mesh.vertex().size() << endl;

  _listMeshRender = buildDisplayList(_mesh,0);
  _listMeshPicking = buildDisplayList(_mesh,1);

//  for(std::size_t i=0;i<_indexTexture.size();++i)
//	  std::cout << _indexTexture[i] << " " << std::endl;

}

void GLWidget::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasUrls())
  {
    if (event->mimeData()->urls()[0].toString().startsWith("file:///"))
    {
      // Check extension
      QString path = event->mimeData()->urls()[0].toString();
      QString ext = path.mid(path.lastIndexOf(".") + 1);
      //std::cout << ext.toStdString() << ext.compare("gii") << std::endl;
      if (ext.compare("gii") == 0)
      {
        event->acceptProposedAction();
      }
    }
  }
}

void GLWidget::dragMoveEvent(QDragMoveEvent *event)
{
  if (event->mimeData()->hasText())
  {
    std::cout << "drag move" << std::endl;
  }
  event->accept();
}

void GLWidget::dropEvent(QDropEvent *event)
{
  QString fName;

  if (event->mimeData()->hasUrls())
  {
    if (event->mimeData()->urls()[0].toString().startsWith("file:///"))
    {
      // Check extension
      QString path = event->mimeData()->urls()[0].toString();
      QString ext = path.mid(path.lastIndexOf(".") + 1);
      if (ext.compare("gii") == 0)
      {
        fName = event->mimeData()->urls()[0].toLocalFile();
        std::cout << fName.toStdString() << std::endl;
      }
    }
  }
}

QPointF GLWidget::pixelPosToViewPos(const QPointF& p)
{
  return QPointF(2.0 * float(p.x()) / width() - 1.0, 1.0 - 2.0 * float(p.y())/ height());
}

Point3df GLWidget::check3DpointPicked (int x, int y)
{
  Point3df p;

  drawScenetoBackBuffer();
  glReadBuffer(GL_BACK);
  glFinish();

  projectionPerspective();
  trackBallTransformation();

  GLdouble modelview[4 * 4];
  GLdouble projection[4 * 4];
  GLint viewport[4];
  GLdouble meshx = 0.0;
  GLdouble meshy = 0.0;
  GLdouble meshz = 0.0;
  GLfloat wz = 0.0;
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);

  glReadPixels(x, viewport[3] - y, 1, 1,
		GL_DEPTH_COMPONENT, GL_FLOAT, &wz);

  gluUnProject((GLdouble) x,(GLdouble)(viewport[3] - y), (GLdouble) wz, modelview,
		projection, viewport, &meshx, &meshy, &meshz);

  p[0] = meshx;
  p[1] = meshy;
  p[2] = meshz;

//  _point3Dpicked[0] = _meshCenter[0] + (float)meshx/_meshScale;
//  _point3Dpicked[1] = _meshCenter[1] + (float)meshy/_meshScale;
//  _point3Dpicked[2] = _meshCenter[2] + (float)meshz/_meshScale;

  return p;
}

GLuint GLWidget::checkIDpolygonPicked (int x, int y)
{
  GLuint r,g,b;
  r = backBufferTexture[3 * (height() - y) * width() + 3 * x];
  g = backBufferTexture[3 * (height() - y) * width() + 3 * x +1];
  b = backBufferTexture[3 * (height() - y) * width() + 3 * x +2];

  return (GLuint)(b + 256 * g + 256 * 256 * r);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
  if (event->buttons()==Qt::LeftButton & _mode == 1)
  {
	_trackBall.start();
    _trackBall.push(pixelPosToViewPos(event->pos()),gfx::Quaternionf::identity());
    event->accept();
    updateGL();
  }

  if (event->buttons()==Qt::LeftButton & _mode == 2)
  {
	int indexPolygon = checkIDpolygonPicked (event->x(),event->y());
    //cout << "ID polygon : " << indexPolygon << endl;
    _point3Dpicked = check3DpointPicked(event->x(),event->y());
    _trackBall.stop();
  }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
  if (event->isAccepted())
    return;

  if (event->buttons() & Qt::LeftButton & _mode == 1)
  {
    _trackBall.release(pixelPosToViewPos(event->pos()),gfx::Quaternionf::identity());
    event->accept();
    updateGL();
  }

  if (event->buttons() & Qt::LeftButton & _mode == 2)
  {
    event->accept();
    updateGL();
  }
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
  if (event->buttons()==Qt::LeftButton && _mode == 1)
  {
    _trackBall.move(pixelPosToViewPos(event->pos()),gfx::Quaternionf::identity());
    event->accept();
    updateGL();
  }

  if (_mode == 2)
  {
	_trackBall.stop();
	_indexPolygon = checkIDpolygonPicked (event->x(),event->y());
    //cout << "ID polygon : " << _indexPolygon << endl;
	_point3Dpicked = check3DpointPicked(event->x(),event->y());
	//cout << "3D : " << _point3Dpicked[0] << " " << _point3Dpicked[1] << " " << _point3Dpicked[2] << " " << endl;

	Point3df p;
	p[0] = _meshCenter[0] + (float)_point3Dpicked[0]/_meshScale;
	p[1] = _meshCenter[1] + (float)_point3Dpicked[1]/_meshScale;
	p[2] = _meshCenter[2] + (float)_point3Dpicked[2]/_meshScale;

	_indexVertex = computeNearestVertexFromPolygonPoint( p, _indexPolygon, _mesh);

	//cout << "3D coord vertex value = " << _vertexNearestpicked[X] << " " << _vertexNearestpicked[Y] << " " << _vertexNearestpicked[Z] << "\n" ;

	updateGL();
  }
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
  switch (event->key())
  {
    case Qt::Key_Plus:
      break;
    case Qt::Key_Minus:
      break;
    case Qt::Key_Left:
      break;
    case Qt::Key_Right:
      break;
    case Qt::Key_Down:
      break;
    case Qt::Key_Up:
      break;
    case Qt::Key_W :
      _wireframe = !_wireframe;
      break;

    default:
      QWidget::keyPressEvent(event);
  }

  updateGL();
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
  int numDegrees = event->delta() / 8;
  float numSteps = (float) (numDegrees / 300.);

  if (_mode == 1)
  {
	_trackBall.push(pixelPosToViewPos(event->pos()),gfx::Quaternionf::identity());

	if (event->orientation() == Qt::Horizontal)
    {
      setTranslate(_trans + numSteps);
    }
    else
    {
      setZoom(_zoom + numSteps);
    }

  event->accept();
  }
}

void GLWidget::drawColorMap (void)
{
  glPushAttrib( GL_ALL_ATTRIB_BITS );

  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glShadeModel(GL_FLAT);

  glDisable(GL_DEPTH_TEST);

  glBegin( GL_QUADS );
      glColor4f(1,1,1,0.75);
      glVertex2d(5,5);
      glVertex2d(200,5);
      glVertex2d(200,105);
      glVertex2d(5,105);
  glEnd();

  glEnable( GL_TEXTURE_2D );
  glBindTexture( GL_TEXTURE_2D,_IDcolorMap );

  glBegin( GL_QUADS );
    glColor4f(1,1,1,0.75);
    glTexCoord2d(0,0);
    glVertex2d(10,10);
    glTexCoord2d(0,0);
    glVertex2d(20,10);
    glTexCoord2d(1,0);
    glVertex2d(20,100);
    glTexCoord2d(1,0);
    glVertex2d(10,100);
  glEnd();

  glPopAttrib();
}

void GLWidget::drawPrimitivePicked (void)
{
  glPushAttrib( GL_ALL_ATTRIB_BITS );

  glDisable( GL_TEXTURE_2D );
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);

  GLUquadricObj *quadric;
  quadric = gluNewQuadric();
  gluQuadricDrawStyle(quadric, GLU_FILL );

  glColor3d(1,0,0);
  glTranslatef(_point3Dpicked[0],_point3Dpicked[1],_point3Dpicked[2]);
  gluSphere( quadric , 0.001 , 36 , 18 );
  glTranslatef(-_point3Dpicked[0],-_point3Dpicked[1],-_point3Dpicked[2]);

  glColor3d(1,1,1);
  glTranslatef(_vertexNearestpicked[0],_vertexNearestpicked[1],_vertexNearestpicked[2]);
  gluSphere( quadric , 0.001 , 36 , 18 );
  glTranslatef(-_vertexNearestpicked[0],-_vertexNearestpicked[1],-_vertexNearestpicked[2]);

  gluDeleteQuadric(quadric);
//printf("angle %f ",(float)(180.*acos(model->facetnorms[3*triangle->findex + 2]))/M_PI);
  	 //glRotatef(acos(model->facetnorms[3*triangle->findex + 2]),1,0,0);
  	 //glutSolidSphere(0.001, 15, 15);

//  	 glMultMatrixf(matrice_rotation);
//
//  	 glColor3d(1,1,0);
//  	 gluCylinder(qobj_cone,0,0.002,0.004,20,1);
//  	 glColor3d(0,1,1);
//  	 gluCylinder(qobj_cone,0,0.001,0.01,20,1);
//
//
//  glPushMatrix();
//
//
//  glPopMatrix();

  glPopAttrib();
}

void GLWidget::projectionOrtho()
{
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, width(), 0, height(), -1,1);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void GLWidget::projectionPerspective()
{
  glMatrixMode( GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(50, (float) (width() / height()), 0.001, 5);
  glMatrixMode( GL_MODELVIEW);
  glLoadIdentity();
}

void GLWidget::trackBallTransformation(void)
{
  glTranslatef(_trans, 0.0, _zoom);
  gfx::Matrix4x4f m;
  _trackBall.rotation().matrix(m);
  glMultMatrixf(m.bits());
}

void GLWidget::paintGL()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  projectionPerspective();

  glEnable( GL_DEPTH_TEST);
  glBindTexture( GL_TEXTURE_2D,_IDcolorMap );


  glPushMatrix();
    trackBallTransformation();

    glEnable( GL_TEXTURE_2D );
    glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);
	if (_listMeshRender != 0)
	      glCallList(_listMeshRender);
	glDisable(GL_POLYGON_OFFSET_FILL);

	glDisable( GL_TEXTURE_2D );
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);

	if (_listMeshRender != 0 && _wireframe)
	      glCallList(_listMeshRender);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);


  glPopMatrix();

  glPushMatrix();
    trackBallTransformation();
    drawPrimitivePicked();
  glPopMatrix();

  glDisable(GL_DEPTH_TEST);
  glDisable( GL_TEXTURE_2D );

  glPushMatrix();
    projectionOrtho();
    drawColorMap();
  glPopMatrix();

  QPainter painter;
  painter.begin(this);
  QString framesPerSecond;
  framesPerSecond.setNum(_frames / (_time.elapsed() / 1000.0), 'f', 2);
  painter.setPen(Qt::black);
  painter.drawText(30, height()-20, framesPerSecond + " fps");

  QString indexPolygon;
  indexPolygon.setNum((int)_indexPolygon,10);
  painter.drawText(30, height()-80, "ID polygon = " + indexPolygon);

  QString indexVertex;
  indexVertex.setNum((int)_indexVertex,10);
  painter.drawText(30, height()-60, "ID vertex = " + indexVertex);

  QString textureValue;
  if (_indexVertex < _mesh.vertex().size())
    textureValue.setNum(_tex[0].item((int)_indexVertex), 'f', 4);
  painter.drawText(30, height()-40, "texture value = " + textureValue);

  painter.end();
  swapBuffers();

  if (!(_frames % 100))
  {
    _time.start();
    _frames = 0;
  }

  _frames++;
}

void GLWidget::resizeGL(int width, int height)
{
  setupViewport(width, height);
}

void GLWidget::setupViewport(int width, int height)
{
  glViewport(0, 0, width, height);
}

int GLWidget::computeNearestVertexFromPolygonPoint( Point3df position, int poly, AimsSurfaceTriangle as)
{
  int index_nearest_vertex,index_min= 0;
  Point3df			pt[3];
  uint				v[3];

  const vector<Point3df>	& vert  =  as.vertex();
  vector< AimsVector<uint,3> >	& tri = as.polygon();

  if (poly < (int)tri.size() && poly > 0)
  {
    v[0] = tri[poly][0];
    v[1] = tri[poly][1];
    v[2] = tri[poly][2];

    pt[0] = vert[v[0]];
    pt[1] = vert[v[1]];
    pt[2] = vert[v[2]];

  //cout << "nb poly= " << tri.size() << endl;

  //compute the nearest polygon vertex
  float min,dist_min = FLT_MAX;

  for (int i = 0 ; i < 3 ; i++)
	{
	  min = (float) sqrt (
	  (position[0]-pt[i][0])*(position[0]-pt[i][0]) +
	  (position[1]-pt[i][1])*(position[1]-pt[i][1]) +
	  (position[2]-pt[i][2])*(position[2]-pt[i][2]) );

	  if (min < dist_min)
	  {
		  dist_min = min;
		  index_min = i;
	  }
	}
  }

  index_nearest_vertex = v[index_min];

  _vertexNearestpicked = pt[index_min];

  _vertexNearestpicked[X] -= (_meshCenter)[X];
  _vertexNearestpicked[Y] -= (_meshCenter)[Y];
  _vertexNearestpicked[Z] -= (_meshCenter)[Z];
  _vertexNearestpicked[X] *= (_meshScale);
  _vertexNearestpicked[Y] *= (_meshScale);
  _vertexNearestpicked[Z] *= (_meshScale);

  //cout << "3D coord vertex value = " << _vertexNearestpicked[X] << " " << _vertexNearestpicked[Y] << " " << _vertexNearestpicked[Z] << "\n" ;

  return index_nearest_vertex;
}

void GLWidget::computeIndexColor (AimsSurfaceTriangle as)
{
  int i;
  int r, g, b;

  _indexTexture.clear();

  for (i = 0; i < (int)(as.polygon().size()); i++)
  {
    r = (i / 256) / 256;
    g = (i / 256) % 256;
    b = i % 256;
    _indexTexture.push_back(r);
    _indexTexture.push_back(g);
    _indexTexture.push_back(b);
  }
}

void GLWidget::drawScenetoBackBuffer (void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  projectionPerspective();
  trackBallTransformation();

  glPushAttrib( GL_ALL_ATTRIB_BITS );

  glEnable( GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  glDisable( GL_TEXTURE_2D );

  if (_listMeshPicking != 0)
    glCallList(_listMeshPicking);

  glPopAttrib();
}

void GLWidget::copyBackBuffer2Texture (void)
{
  drawScenetoBackBuffer();
  glReadBuffer(GL_BACK);

  if (backBufferTexture != NULL)
    free(backBufferTexture);

  backBufferTexture = (GLubyte*) malloc((width()*height()) * 3* sizeof(GLubyte));
  glReadPixels(0, 0, width(), height(), GL_RGB,GL_UNSIGNED_BYTE, backBufferTexture);
  glFinish();
}

GLuint GLWidget::loadColorMap( const char * filename)
{

  int i;
  int width, height;
  unsigned char* data;
  unsigned char* data2;

  FILE * file;
  char textureName[256];

  const Path	& p = Path::singleton();
  char		s = FileUtil::separator();

  string path = p.globalShared() + s + "aims-"
			+ carto::cartobaseShortVersion() + s + "Rgb" + s + filename;

  cout << "colorMap : " << path << endl;

  file = fopen( path.c_str(), "rb" );

  if ( file == NULL ) return 0;

  // read texture data
  fscanf(file,"%s\n%d\nRed\n",&textureName,&width);

  height = 1;
  data = (unsigned char *) malloc( width * height * 3 * sizeof(unsigned char)  );
  data2 = (unsigned char *) malloc( width * height * 3 * sizeof(unsigned char)  );

  for (i = 0 ; i < width ; i++)
    fscanf(file,"%d ",&data[i]);

  fscanf(file,"\nGreen\n");

  for (i = 0 ; i < width ; i++)
    fscanf(file,"%d ",&data[256 + i]);

  fscanf(file,"\nBlue\n");
  for (i = 0 ; i < width ; i++)
    fscanf(file,"%d ",&data[256 + 256 + i]);

  for (i = 0 ; i < width ; i++)
  {
    data2[3*i] = data[i];
    data2[3*i+1] = data[256 + i ];
    data2[3*i+2] = data[512 + i ];
  }

  fclose( file );

  // allocate a texture name
  glGenTextures( 1, &_IDcolorMap );

  // select our current texture
  glBindTexture( GL_TEXTURE_2D,_IDcolorMap);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  // select modulate to mix texture with color for shading
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
  glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, width, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, data2);

  return 1;
}

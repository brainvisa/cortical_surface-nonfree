#ifndef MESHPAINT_H
#define MESHPAINT_H

#include <QMainWindow>
#include <QtGui>
#include "glwidget.h"
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <aims/mesh/surfaceOperation.h>
#include <anatomist/application/globalConfig.h>
#include <anatomist/application/settings.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>

using namespace anatomist;
using namespace aims;
using namespace carto;
using namespace std;

class MeshPaint : public QMainWindow
{
  Q_OBJECT

public:
  MeshPaint();
  ~MeshPaint();

protected :
  void resizeEvent(QResizeEvent *event ) {
	  event->accept();
	  trackballAction->setChecked(true);
	  paintBrushAction->setChecked(false);
	  colorPickerAction->setChecked(false);
  }

  virtual void changeMode(int mode){};
  virtual void saveTexture(void){};

private slots:

  void trackball()
  {
    //std::cout << "trackball" << std::endl;
    paintBrushAction->setChecked(false);
  	colorPickerAction->setChecked(false);
  	changeMode(1);
  }

  void colorPicker()
  {
    //std::cout << "color picker" << std::endl;
    paintBrushAction->setChecked(false);
    trackballAction->setChecked(false);
    changeMode(2);
  }

  void paintBrush()
  {
    //std::cout << "paint brush" << std::endl;
    colorPickerAction->setChecked(false);
    trackballAction->setChecked(false);
    changeMode(3);
  }

  void save()
    {
    saveTexture();
    }

private :
  void createActions();
  void createToolBars();
  QToolBar *paintToolBar;
  QAction *colorPickerAction;
  QAction *paintBrushAction;
  QAction *trackballAction;
  QAction *saveAction;
};

template<typename T>
class myMeshPaint : public MeshPaint
{
public:
    myMeshPaint(string adressTexIn,string adressMeshIn,string adressTexOut,string colorMap, string dataType);
    ~myMeshPaint();

    void changeMode(int mode);
    void saveTexture(void);
    void keyPressEvent( QKeyEvent* event );

private :
    string _adressTexIn;
    string _adressMeshIn;
    string _adressTexOut;
    string _colorMap;
    string _dataType;

    myGLWidget<T> *glWidget;
};

#endif // MESHPAINT_H

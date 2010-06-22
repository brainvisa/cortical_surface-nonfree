#include "meshpaint.h"

template<typename T>
myMeshPaint<T>::myMeshPaint(string adressTexIn,string adressMeshIn,string adressTexOut,string colorMap, string dataType)
	:_adressTexIn(adressTexIn),_adressMeshIn(adressMeshIn),_adressTexOut(adressTexOut),_colorMap(colorMap),_dataType(dataType)
{
  QRect r = geometry();
  r.moveCenter(QApplication::desktop()->availableGeometry().center());
  setGeometry(r);

  glWidget = new myGLWidget<T> (this,adressTexIn,adressMeshIn,adressTexOut,colorMap,dataType);

  //glWidget->setFocusPolicy(Qt::StrongFocus);

  setCentralWidget(glWidget);
}

template<typename T>
myMeshPaint<T>::~myMeshPaint()
{
}

MeshPaint::MeshPaint()
{
  createActions();
  createToolBars();
  show();
}

MeshPaint::~MeshPaint()
{
}

template<typename T>
void myMeshPaint<T>::changeMode(int mode)
{
  cout << "mode = " << mode << endl;
  glWidget->changeMode(mode);
}

template<typename T>
void myMeshPaint<T>::keyPressEvent( QKeyEvent* event )
{
  glWidget->keyPressEvent( event );
}

void MeshPaint::createActions()
{
  string iconname = Settings::globalPath() + "/icons/meshPaint/colorpicker.png";

  colorPickerAction = new QAction(QIcon(iconname.c_str()), tr("&ColorPicker"), this);
  colorPickerAction->setShortcut(tr("c"));
  colorPickerAction->setStatusTip(tr("ColorPicker"));
  colorPickerAction->setCheckable(true);
  connect(colorPickerAction, SIGNAL(triggered()), this, SLOT(colorPicker()));

  iconname = Settings::globalPath() + "/icons/meshPaint/paintbrush.png";

  paintBrushAction = new QAction(QIcon(iconname.c_str()), tr("&PaintBrush"), this);
  paintBrushAction->setShortcut(tr("b"));
  paintBrushAction->setStatusTip(tr("PaintBrush"));
  paintBrushAction->setCheckable(true);
  connect(paintBrushAction, SIGNAL(triggered()), this, SLOT(paintBrush()));

  iconname = Settings::globalPath() + "/icons/meshPaint/zoom.png";

  trackballAction = new QAction(QIcon(iconname.c_str()), tr("&trackballAction"), this);
  trackballAction->setShortcut(tr("s"));
  trackballAction->setStatusTip(tr("trackball"));
  trackballAction->setCheckable(true);
  trackballAction->setChecked(true);
  connect(trackballAction, SIGNAL(triggered()), this, SLOT(trackball()));
}

void MeshPaint::createToolBars()
{
  paintToolBar = addToolBar(tr("PaintToolBar"));
  paintToolBar->setIconSize(QSize(32, 32));
  paintToolBar->addAction(trackballAction);
  paintToolBar->addAction(colorPickerAction);
  paintToolBar->addAction(paintBrushAction);
}

template class myMeshPaint<float>;
template class myMeshPaint<short>;


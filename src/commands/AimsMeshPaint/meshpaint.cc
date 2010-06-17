#include "meshpaint.h"

MeshPaint::MeshPaint(AimsSurfaceTriangle mesh,TimeTexture<float> tex, string colorMap)
	:_mesh(mesh),_tex(tex),_colorMap(colorMap)
{
  createActions();
  createToolBars();

  glWidget = new GLWidget(this,_mesh,_tex,_colorMap);

  setCentralWidget(glWidget);
}

MeshPaint::~MeshPaint()
{
}

void MeshPaint::trackball()
{
  std::cout << "trackball" << std::endl;
  colorPickerAction->setChecked(false);
  paintBrushAction->setChecked(false);

  glWidget->changeMode(1);
}

void MeshPaint::colorPicker()
{
  std::cout << "color picker" << std::endl;
  paintBrushAction->setChecked(false);
  trackballAction->setChecked(false);

  glWidget->changeMode(2);
  glWidget->copyBackBuffer2Texture();
}

void MeshPaint::paintBrush()
{
  std::cout << "paint brush" << std::endl;
  colorPickerAction->setChecked(false);
  trackballAction->setChecked(false);

  glWidget->changeMode(3);
}

void MeshPaint::keyPressEvent( QKeyEvent* _event )
{
  glWidget->keyPressEvent( _event );
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

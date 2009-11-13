/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 *
 *  Just my own little binary for various purposes
 */

#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <aims/utility/utility_g.h>
#include <aims/mesh/mesh_g.h>
#include <cortical_surface/mesh/isoLine.h>
#include <aims/io/io_g.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileMesh1, fileMesh2;
  uint point;

  AimsApplication    app( argc, argv, "From a point on a mesh find the closest point on another mesh" );
  app.addOption( fileMesh1, "-m1", "source mesh" );
  app.alias( "--mesh1", "-m1" );
  app.addOption( fileMesh2, "-m2", "target mesh" );
  app.alias( "--mesh2", "-m2" );
  app.addOption( point, "-p", "point index on source mesh" );
  app.alias( "--point", "-p" );
  app.initialize();
  
  cout << "reading source triangulation   : " << flush;
  AimsSurfaceTriangle surface1;
  Reader<AimsSurfaceTriangle> surf1R( fileMesh1 );
  surf1R >> surface1;
  cout << "done" << endl;
  cout << "reading target triangulation   : " << flush;
  AimsSurfaceTriangle surface2;
  Reader<AimsSurfaceTriangle> surf2R( fileMesh2 );
  surf2R >> surface2;
  cout << "done" << endl;

  Point3df point1, pointTemp;
  std::vector<Point3df> vert1=surface1.vertex();
  std::vector<Point3df> vert2=surface2.vertex();

  point1=vert1[point];
  uint ind2;
  uint nv1=vert1.size();
  uint nv2=vert2.size();

  float distMin=10000.0, dist;

  cout << "Looking for the target point" <<endl;
  for (uint i=0; i<nv2; i++)
  {
    pointTemp=vert2[i];
    dist=(pointTemp-point1).norm();
    if (dist < distMin)
    {
      distMin=dist;
      ind2=i;
    }
  }
  cout << "Generating textrues" <<endl;
  TimeTexture<short> tex1(1,nv1), tex2(1,nv2);
  for (uint i=0; i<nv1; i++)
  {
    if (i==point) tex1[0].item(i)=100;
    else tex1[0].item(i)=0;
  }
  for (uint i=0; i<nv2; i++)
  {
    if (i==ind2) tex2[0].item(i)=100;
    else tex2[0].item(i)=0;
  }
  std::cout << "Point index on target surface: " << ind2 << std::endl; 
  Writer<TimeTexture<short> > texW1("texSource.tex");
  Writer<TimeTexture<short> > texW2("textarget.tex");
  texW1.write(tex1);
  texW2.write(tex2);
  return EXIT_SUCCESS;
}


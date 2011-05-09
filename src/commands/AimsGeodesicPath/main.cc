/*
 *  Copyright (C) 1997-2005 CEA
 *
 *  This software and supporting documentation were developed by
 *
 *      Laboratoire LSIS, equipe LXAO,
 *      Marseille, France
 *
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 */
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/process.h>
#include <string.h>
#include <fstream>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/geodesicpath/geodesicPath.h>

using namespace aims;
using namespace carto;
using namespace std;

int main( int argc, const char** argv )
{
  try
  {
  string meshFileIn;string FileOut="./out";
  int method = 0;
  int strain = 3;
  unsigned source,target;

  AimsApplication app( argc, argv, "Compute the shortest path between two vertex" );

  app.addOption( meshFileIn, "-i", "mesh" );
  app.alias( "--input", "-i" );

  app.addOption( source, "-s", "index of source vertex" );
  app.alias( "--source", "-s" );

  app.addOption( target, "-t", "index of target vertex" );
  app.alias( "--target", "-t" );

  app.addOption( FileOut, "-o", "output file without extension file specified (.tex or .mesh)" );
  app.alias( "--output", "-o" );

  app.addOption( method, "-c", "constraintType:\n\"0\" -> no constraint\n"
      "\"1\" -> constrained sulci\n\"2\" -> constrained gyri\n\"3\" -> Exact geodesic path",true);
  app.alias( "--constraint", "-c" );

  app.addOption( strain, "-st", "strain parameter (3 by default)",true );
  app.alias( "--strain", "-st" );

  app.initialize();

  // read triangulation
  cout << "reading triangulation   : " << flush;
  AimsSurfaceTriangle surface;
  Reader<AimsSurfaceTriangle> triR( meshFileIn );
  triR >> surface;
  cout << "done" << endl;

  GeodesicPath sp(surface,method,strain);

  vector<int> pathIndex;

  cout << "set of index vertex" << endl;
  pathIndex = sp.shortestPath_1_1_ind(source,target);
  for (int i = 0; i < pathIndex.size(); i++)
  cout << pathIndex[i] << " ";
  cout << "done" << endl;

  cout << "set of 3D coord vertex" << endl;
  vector<Point3df> path3D;
  path3D = sp.shortestPath_1_1_xyz(source,target);
  for (int i = 0; i < path3D.size(); i++)
  cout << path3D[i][0] << " " << path3D[i][1] << " "<< path3D[i][2] << "\n";
  cout << "done" << endl;

  if (method==3)
  {
    AimsSurfaceTriangle meshOut, *tmpMeshOut;
    tmpMeshOut = new AimsSurfaceTriangle;
    int i;
    for (i = 0; i < path3D.size() - 1; i++)
    {
      tmpMeshOut = SurfaceGenerator::sphere(path3D[i], 0.25 ,20 );
      SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
      delete tmpMeshOut;

      tmpMeshOut = SurfaceGenerator::cylinder( path3D[i],path3D[i+1], 0.2, 0.2, 12, false, true );
      SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
      delete tmpMeshOut;
    }

    tmpMeshOut = SurfaceGenerator::sphere(path3D[i], 0.2 ,10 );
    SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
    delete tmpMeshOut;

    FileOut = FileOut + ".mesh";
    Writer<AimsSurfaceTriangle> wm(FileOut);
    wm.write(meshOut);
  }
  else
  {
    TimeTexture<float> texOut(1, surface.vertex().size() );

    sp.shortestPath_1_1_tex(source,target,180,texOut);

    FileOut = FileOut + ".tex";
    Writer<TimeTexture<float> > texW(FileOut);
    texW << texOut;

    double len = sp.shortestPath_1_1_len (source,target);
    cout << "path length = " << len << endl;

    vector<unsigned> targetList;
    targetList.push_back(40);
    targetList.push_back(400);
    targetList.push_back(340);
    targetList.push_back(403);
    sp.shortestPath_1_N_ind(source,targetList,&target,&len);
    cout << "best target = " << target << " length = " << len << endl;

  }

  return( 0 );
  }

  catch( user_interruption & )
  {
  }
  catch( exception & e )
  {
  cerr << e.what() << endl;
  }
  return 1;
}

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

using namespace aims;
using namespace carto;
using namespace std;

#include <aims/geodesicpath/geodesicPath.h>

int main( int argc, const char** argv )
{
  try
  {
    string meshFileIn;
    string FileOut= "./out";
    int method = 0;
    int strain = 3;
    unsigned source,target;

    AimsApplication    app( argc, argv, "Compute the shortest path between two vertex" );

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

    app.addOption( strain, "-st", "strain parameter (5 by default)",true );
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

    pathIndex = sp.shortestPathIndiceVextex(source,target);

    for (int i = 0; i < pathIndex.size(); i++)
      cout << pathIndex[i] << " " ;

    vector<Point3df> path3D;

    path3D = sp.shortestPathCoordVextex(source,target);

    for (int i = 0; i < path3D.size(); i++)
      cout << path3D[i][0] << " " << path3D[i][1] << " "<< path3D[i][2] << "\n";

    TimeTexture<float> texOut(1, surface.vertex().size() );

    sp.shortestPath2Texture(source,target,180, texOut);

    Writer<TimeTexture<float> > texW(FileOut);
    texW << texOut;

    float len = sp.shortestPathLength (source,target);
    cout << "length = " << len << endl;

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

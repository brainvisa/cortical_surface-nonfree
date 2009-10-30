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
  string fileOut, fileMesh, fileGyri;
  float gyrus;

  AimsApplication    app( argc, argv, "Extract a particular gyrus from a mesh and its gyri texture" );
  app.addOption( fileMesh, "-m", "input mesh" );
  app.alias( "--mesh", "-m" );
  app.addOption( fileGyri, "-t", "gyri texture" );
  app.alias( "--tex", "-t" );
  app.addOption( gyrus, "-g", "gyrus index" );
  app.alias( "--gyrus", "-g" );
  app.addOption( fileOut, "-o", "output gyrus mesh" );
  app.alias( "--out", "-o" );
  app.initialize();
  
  cout << "reading triangulation   : " << flush;
  AimsSurfaceTriangle surface;
  Reader<AimsSurfaceTriangle> triR( fileMesh );
  triR >> surface;
  cout << "done" << endl;

  cout << "reading gyri texture   : " << flush;
  Reader<TimeTexture<float> > texGyriR( fileGyri );
  TimeTexture<float> texG;
  texGyriR >> texG ;
  cout << "done " << endl;

  std::vector< Point3df > vertG;
  std::vector< AimsVector< uint, 3 > >  polyG;

  std::vector<Point3df> vert=surface.vertex();
  uint nv=vert.size();
  std::vector< AimsVector< uint, 3 > > poly=surface.polygon();
  uint np=poly.size();
  uint i;
  uint count=0;
  std::map<uint, uint> trans;
  
  cout << "Extracting patch" << endl;
  for (i=0; i<nv; i++)
  {
    if (fabs(texG[0].item(i) - gyrus)<0.1)
    {
      vertG.push_back(vert[i]);
      trans[i]=count;
      count++;
    }
  }
  for (i=0; i<np; i++)
  {
    if ( (fabs(texG[0].item(poly[i][0]) - gyrus)< 0.1)
      && (fabs(texG[0].item(poly[i][1]) - gyrus)< 0.1)
      && (fabs(texG[0].item(poly[i][2]) - gyrus)< 0.1) )
    {
      polyG.push_back( AimsVector< uint, 3 >(trans[poly[i][0]], trans[poly[i][1]], trans[poly[i][2]] ) );
    }
  }
  
  AimsSurfaceTriangle surfaceG;
  surfaceG.vertex()=vertG;
  surfaceG.polygon()=polyG;
   
  // recomputing normals
  cout << "Recomputing normals" << endl;
  surface.updateNormals();
  
  cout << "saving gyrus patch" << flush;
  Writer<AimsSurfaceTriangle> triW( fileOut );
  triW << surfaceG;
  cout << "done" << endl;

  return EXIT_SUCCESS;
}


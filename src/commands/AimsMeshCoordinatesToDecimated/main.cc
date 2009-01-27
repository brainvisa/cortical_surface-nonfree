/*
 *  Copyright (C) 1997-2005 CEA
 *
 *  This software and supporting documentation were developed by
 *
 *      Laboratoire LSIS, �quipe LXAO,
 *      Marseille, France
 *
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 */

#include <aims/mesh/geometric.h>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>
#include <aims/math/math_g.h>
#include <aims/vector/vector.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/curv.h>
#include <aims/distancemap/meshdistance_d.h>
#include <aims/distancemap/meshmorphomat_d.h>
#include <aims/connectivity/meshcc_d.h>
#include <aims/morphology/morphology_g.h>
#include <cortical_surface/surfacereferential/shortestPath.h>
#include <cortical_surface/mesh/linkPath.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/getopt/getopt2.h>
#include <aims/connectivity/meshcc.h>
#include <iostream>
#include <iomanip>

using namespace aims;
using namespace carto;
using namespace std;

int main( int argc, const char** argv )
{
  try
    {
      string  meshfile, meshDeci, texIn_x, texIn_y, texOut_x, texOut_y;

      AimsApplication    app( argc, argv, "Transfer 2D coordinates of a mesh to a decimated version of this mesh" );
      app.addOption( meshfile, "-i", "input mesh" );
      app.alias( "--input", "-i" );
      app.addOption(meshDeci, "-d", "decimatedMesh" );
      app.alias( "--deci", "-d");
      app.addOption( texIn_x, "-ix", "input x coordinate" );
      app.alias( "--inputx", "-ix" );
      app.addOption( texIn_y, "-iy", "input y coordinate" );
      app.alias( "--inputy", "-iy" );
       app.addOption( texOut_x, "-ox", "output x coordinate" );
      app.alias( "--outputx", "-ox" );
      app.addOption( texOut_y, "-oy", "output y coordinate" );
      app.alias( "--outputy", "-oy" );
      app.initialize();

      //
      // read triangulation

      cout << "reading triangulations   : " << endl;
      AimsSurfaceTriangle surfaceIn;
      Reader<AimsSurfaceTriangle> triR( meshfile );
      triR >> surfaceIn;
      AimsSurfaceTriangle surfaceDeci;
      Reader<AimsSurfaceTriangle> deciR( meshDeci );
      deciR >> surfaceDeci;
      cout << "done" << endl;
      
      uint ni, no;
      ni=surfaceIn.vertex().size();
      no=surfaceDeci.vertex().size();
      
      std::cout << "Mesh : " << ni << " vertices, and decimated : " << no << " vertices" << std::endl;
      
      TimeTexture<float> xout(1, no);
      TimeTexture<float> yout(1, no);
      
      TimeTexture<float> xin;
      TimeTexture<float> yin;
      Reader<TimeTexture<float> > xr(texIn_x);
      xr >> xin;
      Reader<TimeTexture<float> > yr(texIn_y);
      yr >> yin;
      
      uint i, j, neigh;
      double d, dmin;
      Point3df p1, p2;
      
      std::cout << "Computing new values" << std::endl;
      for (i=0; i<no; i++)
      {
          dmin=100000.0;
          p1=surfaceDeci.vertex()[i];
          for (j=0; j<ni; j++)
          {
               p2=surfaceIn.vertex()[j];
               d=sqrt( ( p1[0]-p2[0] )*( p1[0]-p2[0] ) + ( p1[1]-p2[1] )*( p1[1]-p2[1] ) + ( p1[2]-p2[2] )*( p1[2]-p2[2] ) );
               if (d<dmin)
               {
                    dmin=d; neigh=j;
               }
          }
          xout[0].item(i)=xin[0].item(neigh);
          yout[0].item(i)=yin[0].item(neigh);
      }
      std::cout << "OK, writing textures" << std::endl;
      Writer<TimeTexture<float> > xw(texOut_x);
      xw.write(xout);
      Writer<TimeTexture<float> > yw(texOut_y);
      yw.write(yout);
      std::cout << "Done" << std::endl;
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

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
  string fileOut, fileMesh, fileLat, fileLon;

  AimsApplication    app( argc, argv, "Check that a 2D surface-based coordinate system is OK (unicity)" );
  try
  {
     app.addOption( fileMesh, "-m", "input mesh" );
     app.alias( "--mesh", "-m" );
     app.addOption( fileLat, "-x", "latitude texture" );
     app.alias( "--xcoord", "-x" );
     app.addOption( fileLon, "-y", "longitude texture" );
     app.alias( "--ycoord", "-y" );
     app.addOption( fileOut, "-o", "output sign texture" );
     app.alias( "--out", "-o" );
     
     app.initialize();
     
     cout << "reading triangulation   : " << flush;
     AimsSurfaceTriangle surface;
     Reader<AimsSurfaceTriangle> triR( fileMesh );
     triR >> surface;
     cout << "done" << endl;
     
     cout << "reading textures   : " << flush;
     Reader<TimeTexture<float> > texLonR( fileLon );
     TimeTexture<float> texLon;
     texLonR >> texLon ;
     Reader<TimeTexture<float> > texLatR( fileLat );
     TimeTexture<float> texLat;
     texLatR >> texLat ;
     cout << "done " << endl;
     
     
     cout << "computing neighbours  " << endl;
     vector<set<uint> >  neigh = SurfaceManip::surfaceNeighbours(surface);
     
     uint i, nVert=surface.vertex().size();
     TimeTexture<float> texOut(1, nVert);
     vector<Point3df>  vert=surface.vertex();

     cout << "computing normals" << endl;
     surface.updateNormals();
     vector<Point3df> norm=surface.normal();
     cout << "OK. Starting check" << endl;

     for (i=0; i< nVert; i++)
     {
          set<uint> v=neigh[i];
          set<uint>::iterator vIt=v.begin();
          uint j, jlat, jlon;
          float glat, glon, latM=-1000.0, lonM=-1000.0;
          
          for ( ; vIt!=v.end(); ++vIt)
          {
               j=*vIt;
               glat=texLat[0].item(j)-texLat[0].item(i);
               glon=texLon[0].item(j)-texLon[0].item(i);
               if (glat>latM) {latM=glat; jlat=j;}
               if (glon>lonM) {lonM=glon; jlon=j;}
          }
          Point3df gradLat, gradLon, vp, n=norm[i];
          float sign;
          gradLat=vert[jlat]-vert[i]; gradLat=gradLat/gradLat.norm();
          gradLon=vert[jlon]-vert[i]; gradLon=gradLon/gradLon.norm();
          n=n/n.norm();
          vp=vectProduct(gradLat, gradLon); 
          sign=vp.dot(n);
//           if (sign>0)
//                texOut[0].item(i)=1;
//           else if (sign<0)
//                texOut[0].item(i)=-1;
//           else
//                texOut[0].item(i)=0;
          texOut[0].item(i)=sign;
     }
     Writer<TimeTexture<float> > texOutW( fileOut );
     texOutW << texOut ;
     return EXIT_SUCCESS;
  }
  catch( user_interruption & )
  {
  }
  catch( ... )
  {
    throw;
  }
  return EXIT_FAILURE;
}


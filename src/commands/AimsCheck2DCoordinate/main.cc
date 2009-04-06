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
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>
#include <aims/distancemap/meshmorphomat_d.h>
#include <aims/connectivity/meshcc.h>
#include <aims/connectivity/meshcc_d.h>

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
     TimeTexture<short> texOut(1, nVert);
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
               Point3df dist=vert[j]-vert[i];
               glat=(texLat[0].item(j)-texLat[0].item(i))/dist.norm();
               glon=(texLon[0].item(j)-texLon[0].item(i))/dist.norm();
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
          if (sign>=0)
               texOut[0].item(i)=0;
          else if (sign<0)
               texOut[0].item(i)=1;
          else
               texOut[0].item(i)=0;
/*          texOut[0].item(i)=sign;*/
     }

     cout << "Postprocessing check" << endl;
     for (i=0; i<nVert; i++)
     {
          set<uint> v=neigh[i];
          set<uint>::iterator vIt=v.begin();
          uint j, count=0;
          if (texOut[0].item(i) == 1)
          {
               for ( ; vIt!=v.end(); ++vIt)
               {
                    j=*vIt;
                    if (texOut[0].item(i) == texOut[0].item(j))
                         count++;
               }
               if (count<=1)
               {
                    texOut[0].item(i) = 0;
               }
          }
     }
     cout << "OK. Writing texture texOut.tex" << endl;
     Writer<TimeTexture<short> > texOutW( "texOut.tex" );
     texOutW << texOut ;
     
     cout << "Closing result" << endl;

     TimeTexture<short> dil, closed;
     float sizeClosing=6.0;
     
     dil[0] = MeshDilation<short>( surface[0], texOut[0], short(0), -1 , sizeClosing, true);
     cout << "OK. Writing texture dil.tex" << endl;
     Writer<TimeTexture<short> > dilW( "dil.tex" );
     dilW << dil ;
     closed[0] = MeshErosion<short>( surface[0], dil[0], short(0), -1 , sizeClosing, true);
     cout << "OK. Writing texture closed.tex" << endl;
     Writer<TimeTexture<short> > closedW( "closed.tex" );
     closedW << closed ;
     
     cout << "Extracting connected components" << endl;
     
     TimeTexture<short> compoOut;
     compoOut[0] = AimsMeshLabelConnectedComponent<short>( surface[0], closed[0], 0, 0);

     std::map<short, std::vector<uint> > compos;
     short val;
     for (i=0; i<nVert; i++)
     {
          val=compoOut[0].item(i);
          if (val>0)
          {
               if (compos.find(val) == compos.end())
               {
                    compos[val]=std::vector<uint>();
                    compos[val].push_back(i);
               }
               else
               {
                    compos[val].push_back(i);
               }

          }
     }

     cout << "Found " << compos.size() << " components";
     cout << "Processing them" << endl;

     TimeTexture<short> texResult(1,nVert);
     for (i=0; i<nVert; i++) texResult[0].item(i)=0;
     
     uint j,k, l, sc;
     float lat, lon, latMax, latMin, lonMax, lonMin;
     for (i=1; i<=compos.size(); i++)
     {
          cout << "Composante " << i << endl;
          lonMax=latMax=0.0; lonMin=latMin=360.0;
          std::vector<uint> compo=compos[i];
          sc=compo.size();
          cout << "\t size " << sc << endl;
          
          for (j=0; j<sc; j++)
          {
               k=compo[j];
               lat=texLat[0].item(k);
               lon=texLon[0].item(k);
               if (lat<latMin) latMin=lat;
               if (lat>latMax) latMax=lat;
               if (lon<lonMin) lonMin=lon;
               if (lon>lonMax) lonMax=lon;
          }
          cout << "bounding box : (" << latMin << ", " << lonMin << ") - (" << latMax << ", " << lonMax << ")" << endl;
          for (l=0; l<nVert; l++)
          {
               lat=texLat[0].item(l);
               lon=texLon[0].item(l);
               if ( (lat<=latMax) && (lat>=latMin) && (lon<=lonMax) && (lon>=lonMin))
                    texResult[0].item(l)=i;
          }
     }

     cout << "OK. Writing texture " << fileOut << endl;
     Writer<TimeTexture<short> > compoOutW( fileOut );
     compoOutW << texResult ;
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


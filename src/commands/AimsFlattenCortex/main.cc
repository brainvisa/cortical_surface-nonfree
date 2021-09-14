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
#include <aims/data/data.h>

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileOut, fileMesh, fileLat, fileLon, fileTex, fileIma;
  
  AimsApplication    app( argc, argv, "Create a flat representation of a cortical mesh if it has been parameterised" );
  app.addOption( fileMesh, "-m", "input mesh" );
  app.alias( "--mesh", "-m" );
  app.addOption( fileLat, "-x", "latitude texture" );
  app.alias( "--xcoord", "-x" );
  app.addOption( fileLon, "-y", "longitude texture" );
  app.alias( "--ycoord", "-y" );
  app.addOption (fileTex, "-t", "texture");
  app.alias( "--tex", "-t");
  app.addOption( fileOut, "-om", "output mesh" );
  app.alias( "--outM", "-om" );
  app.addOption( fileIma, "-oi", "output ima" );
  app.alias( "--outI", "-oi" );

  try
  {
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

    Reader<TimeTexture<float> > texImaR(fileTex);
    TimeTexture<float> texIma;
    texImaR >> texIma;
    cout << "done " << endl;

  std::vector<Point3df> & vert=surface.vertex();
  int nv=vert.size();
  int i;
  float x,y;
  float lon, lat;

  cout << "Computing new coordinates" << endl;
  for (i=0; i<nv; ++i)
  {
      lon=texLon[0].item(i);
      lat=texLat[0].item(i);
      if (lon<=180) x=180.0 -lon;
      else if (lon >180) x=540.0-lon;
  //      x=lon;
      y=180.0-lat;
      vert[i][0]=x;
      vert[i][1]=y;
      vert[i][2]=texIma[0].item(i)*10.0;
  }

  cout << "Cleaning triangles" << endl;
  vector< AimsVector<uint,3> > poly = surface.polygon();
  vector< AimsVector<uint,3> > newPoly;
  int np=poly.size(), nbNew=0;;
  std::set<uint> removed;
  for (i=0; i<np; ++i)
  {
      uint v1, v2, v3;
      v1=poly[i][0];
      v2=poly[i][1];
      v3=poly[i][2];

      if ( ( (fabs(texLon[0].item(v1))>10) && (fabs(texLon[0].item(v2))>10) &&  (fabs(texLon[0].item(v3))>10) &&
              (  ((texLon[0].item(v1)<=180) && (texLon[0].item(v2) >180))      // this test has to be improved
            || ((texLon[0].item(v2)<=180) && (texLon[0].item(v1) >180))  // this is a quick try
            || ((texLon[0].item(v1)<=180) && (texLon[0].item(v3) >180))
            || ((texLon[0].item(v3)<=180) && (texLon[0].item(v1) >180))
            || ((texLon[0].item(v2)<=180) && (texLon[0].item(v3) >180))
            || ((texLon[0].item(v3)<=180) && (texLon[0].item(v2) >180)) ) )
            || (fabs(texLat[0].item(v1))<0.001) || (fabs(texLat[0].item(v2))<0.001) || (fabs(texLat[0].item(v3))<0.001)
            || (fabs(texLat[0].item(v1) - 180.0)<0.001) || (fabs(texLat[0].item(v2) - 180.0)<0.001) || (fabs(texLat[0].item(v3) - 180.0)<0.001) )
      {}
      else
      {
            nbNew++;
            AimsVector<uint,3> triangle(v1,v2,v3);
            newPoly.push_back(triangle);
      }
  }

  cout << nbNew << " triangles have been kept (total = " << np << ")" << endl;

  surface.polygon()=newPoly;


  // recomputing normals
  cout << "Recomputing normals" << endl;
  surface.updateNormals();
  /*
  std::vector< Point3df > & norm=surface.normal();
  std::vector< Point3df >::iterator itNorm=norm.begin();
  for ( ; itNorm!=norm.end(); ++itNorm)
  {
      (*itNorm)[0]= -(*itNorm)[0];
      (*itNorm)[1]= -(*itNorm)[1];
      (*itNorm)[2]= -(*itNorm)[2];
  }*/

  cout << "Computing image from texture" << endl;
  int index1, index2, index3;
  float x1, x2, x3, y1, y2, y3;
  float fx, fy;
  int flag=0;
  int dx=80, dy=30;
  int ix, iy;
  AimsData<float>  ima(360-(2*dx), 180-(2*dy), 1);
  for (iy=dy; iy< 180-dy; ++iy)
  {
      cout << iy-dy << "/" << (180-2*dy) << endl;
      for (ix=dx; ix<360-dx; ++ix)
      {
            fx=(float) ix;
            fy=(float) iy;
            vector< AimsVector<uint,3> >::iterator itPoly=newPoly.begin();
            flag=0;
            for ( ; (itPoly != newPoly.end()) && (flag==0); ++itPoly)
            {
                AimsVector<uint, 3> triangle=*itPoly;
                index1=triangle[0];
                index2=triangle[1];
                index3=triangle[2];

                x1=vert[index1][0]; y1=vert[index1][1];
                x2=vert[index2][0]; y2=vert[index2][1];
                x3=vert[index3][0]; y3=vert[index3][1];

  //                x1=texLon[0].item(index1);y1=texLat[0].item(index1);
  //                x2=texLon[0].item(index2);y2=texLat[0].item(index2);
  //                x3=texLon[0].item(index3);y3=texLat[0].item(index3);

                float v1=((x2-x1)*(fy-y1) - (y2-y1)*(fx-x1))*((fx-x1)*(y3-y1) - (fy-y1)*(x3-x1));
                float v2=((x1-x3)*(fy-y3) - (y1-y3)*(fx-x3))*((fx-x3)*(y2-y3) - (fy-y3)*(x2-x3));
                if ((v1>=0) && (v2>=0))
                      flag=1;
            }
            if (flag==0)
            {
                cout << "Triangle not found for x=" << ix << " and y=" << iy << endl;
                exit(1);
            }
            else
            {
                double t1, t2, t3;
                double sum;
                double l1, l2, l3;
                if (((fabs(x2-x1)<0.0001) && (fabs(x3-x1)<0.0001))
                      || ((fabs(y2-y1)<0.0001) && (fabs(y3-y1)<0.0001)))
                // three points of the triangle have the same x or y-coordinate
                // weird but it can happen on the main coordinate axis
                // (bottom ridge of a sulcus for instance)
                {
                      t1=t2=t3=1.0/3.0;
                }
                else
                {
                      l1=sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
                      l2=sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1));
                      l3=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

                      t1=fabs((x2-x3)*(fy-y3) - (y2-y3)*(fx-x3))/l1;
                      t2=fabs((x3-x1)*(fy-y1) - (y3-y1)*(fx-x1))/l2;
                      t3=fabs((x1-x2)*(fy-y2) - (y1-y2)*(fx-x2))/l3;
                      sum=t1+t2+t3; t1/=sum; t2/=sum; t3/=sum;
                }
                ima(ix-dx, iy-dy, 0)=(float) t1*texIma[0].item(index1) + t2*texIma[0].item(index2) + t3*texIma[0].item(index3);
            }
      }
  }

  cout << "saving triangulation    : " << flush;
  Writer<AimsSurfaceTriangle> triW( fileOut );
  triW << surface;
  Writer<AimsData<float> > imaW( fileIma );
  imaW << ima;
  cout << "done" << endl;
  }
  catch( user_interruption & )
  {
  }
  catch( exception & e )
  {
    cerr << argv[ 0 ] << ": " << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


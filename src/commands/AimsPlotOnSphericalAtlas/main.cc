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
#include <aims/mesh/surface.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace carto;
using namespace std;


int main( int argc, const char* argv[] )
{
  Reader<AimsSurfaceTriangle> triR;
  Reader<TimeTexture<float> > lonR;
  Reader<TimeTexture<float> > latR;
  Writer<AimsSurfaceTriangle> sphW;

  string fileIn;

  AimsApplication app( argc, argv,
                       "coordinate file (ASCII) contains one line per point, "
                       "each line contains only longitude and latitude "
                       "separated by a space:\n"
                       "\n"
                       " u1 v1\n"
                       " u2 v2\n"
                       " u3 v3\n"
                       " ..." );
  app.addOption( fileIn, "-i", "coordinates file" );
  app.addOption( triR, "-m", "atlas mesh" );
  app.addOption( lonR, "-u", "longitude texture" );
  app.addOption( latR, "-v", "latitude texture" );
  app.addOption( sphW, "-o", "output mesh" );
  app.alias( "--input", "-i" );
  app.alias( "--longitude", "-u" );
  app.alias( "--latitude", "-v" );
  app.alias( "--output", "-o" );

  try
  {
    app.initialize();

    //
    // read triangulation
    //
    cout << "reading atlas mesh   : " << flush;
    AimsSurfaceTriangle atlas;
    triR >> atlas;
    cout << "done" << endl;

    cout << "reading atlas coordinates   : " << flush;
    TimeTexture<float> longi, lati;
    lonR >> longi;
    latR >> lati;
    cout << "done" << endl;

    cout << "Reading points file   : " << flush;
    FILE *pointList;
    pointList=fopen(fileIn.c_str(), "r");
    int n=0;
    float u,v;
    std::vector<std::pair<float,float> > points;

    while (fscanf(pointList, "%f %f\n", &u, &v) != EOF)
    {
        if ((u<0) || (u>360))
        {
              cerr << "Point with longitude outside [0; 360]" << endl;
              return(EXIT_FAILURE);
        }
        if ((v<0) || (v>180))
        {
              cerr << "Point with latitude outside [0; 180]" << endl;
              return(EXIT_FAILURE);
        }
        points.push_back(std::pair<float, float>(u,v));
        n++;
    }
    fclose(pointList);

    cout << "found " << n << " points" << endl;

    cout << "Generating mesh    : " << flush;
    AimsSurfaceTriangle pointsM;
    std::vector<Point3df> vertices=atlas.vertex();
    int nv=vertices.size();
    int i,j, plot;
    float s,t;
    SurfaceManip merging;
    double dmin, d;

    for (i=0; i<n; i++)
    {
        u=points[i].first;
        v=points[i].second;
        dmin=1000.0;
        // finding the closet point on the mesh
        for (j=0; j<nv; j++)
        {
              s=longi.item(j);
              t=lati.item(j);
              d=((u-s)*(u-s) + (v-t)*(v-t));
              if (d<dmin)
              {
                  dmin=d;
                  plot=j;
              }
        }
        //generating little sphere
          AimsSurfaceTriangle	*msh;
          msh = SurfaceGenerator::sphere(vertices[plot], 1.5, 40);
          merging.meshMerge<3, Void>(pointsM, *msh);
          delete msh;
    }
    cout << "done" << endl;

    sphW.write(pointsM);

    return EXIT_SUCCESS;
  }
  catch( user_interruption & )
  {
  }
  catch( exception & e )
  {
    cerr << e.what() << endl;
  }
  return EXIT_FAILURE;
}


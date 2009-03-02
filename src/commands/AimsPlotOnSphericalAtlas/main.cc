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
#include <aims/getopt/getopt.h>
#include <aims/utility/utility_g.h>
#include <aims/mesh/mesh_g.h>
#include <aims/io/io_g.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>

using namespace aims;
using namespace std;

BEGIN_USAGE( usage )
  "--------------------------------------------------------------------------",
  "AimsPlotOnSphericalAtlas  -i[nput]  <coordinates file>                    ",
  "                          -m        <atlas mesh>                          ",
  "                          -u        <longitude texture>                   ",
  "                          -v        <latitude texture>                    ",
  "                          -o        <output mesh>                         ",
  "                         [-h[elp]]                                        ",
  "                                                                          ",
  " coordinate file (ASCII) contains one line per point, each line           ",
  " contains only longitude and latitude separated by a space :              ",
  "                                                                          ",
  " u1 v1                                                                    ",
  " u2 v2                                                                    ",
  " u3 v3                                                                    ",
  " ...                                                                      ",
  "--------------------------------------------------------------------------",
END_USAGE


void Usage( void )
{
  AimsUsage( usage );
}


int main( int argc, char* argv[] )
{
  char *fileIn = NULL, *fileMesh = NULL, *fileLon = NULL, *fileLat = NULL, *fileOut=NULL;

  AimsOption opt[] = {
  { 'h',"help"      ,AIMS_OPT_FLAG  ,( void* )Usage           ,AIMS_OPT_CALLFUNC,0},
  { 'i',"input"     ,AIMS_OPT_STRING,&fileIn         ,0                ,1},
  { 'm',"mesh"      ,AIMS_OPT_STRING,&fileMesh, 0, 1},
  { 'u',"longitude" ,AIMS_OPT_STRING,&fileLon, 0, 1},
  { 'v',"latitude"  ,AIMS_OPT_STRING,&fileLat, 0, 1},
  { 'o',"output"    ,AIMS_OPT_STRING,&fileOut, 0, 1}, 
  };


  AimsParseOptions( &argc, argv, opt, usage );

  //
  // read triangulation
  //
  cout << "reading atlas mesh   : " << flush;
  AimsSurfaceTriangle atlas;
  Reader<AimsSurfaceTriangle> triR( fileMesh );
  triR >> atlas;
  cout << "done" << endl;
  
  cout << "reading atlas coordinates   : " << flush;  
  TimeTexture<float> longi, lati;
  Reader<TimeTexture<float> > lonR( fileLon);
  lonR >> longi;
  Reader<TimeTexture<float> > latR( fileLat);
  latR >> lati;
  cout << "done" << endl;

  cout << "Reading points file   : " << flush;
  FILE *pointList;
  pointList=fopen(fileIn, "r");
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
   
  Writer<AimsSurfaceTriangle> sphW(fileOut);
  sphW.write(pointsM);

  return EXIT_SUCCESS;
}


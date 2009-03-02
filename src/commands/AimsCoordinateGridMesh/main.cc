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
  string fileOut, fileMesh, fileLat, coord, fileLon;
  float diam=0.25;

  AimsApplication    app( argc, argv, "generate iso-parameter lines from a mesh and textures" );
  try
  {
    app.addOption( fileMesh, "-m", "input mesh" );
    app.alias( "--mesh", "-m" );
    app.addOption( fileLat, "-x", "latitude texture (or -x for sulci)" );
    app.alias( "--xcoord", "-x" );
    app.addOption( fileLon, "-y", "longitude texture (or -y for sulci)" );
    app.alias( "--ycoord", "-y" );
    app.addOption( coord, "-c", "coordinates (r=regular, c=constraints, s=sillon)" );
    app.alias( "--coord", "-c" );
    app.addOption( fileOut, "-o", "output mesh" );
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

//   cout << "Parsing coordinate list from " << fileCoord << endl;
//   std::list<int> latitude;
//   std::list<int> longitude;
//   char flag[3];
//   int value;
//   FILE *coord;
//   if ((coord=fopen(fileCoord, "r")) == NULL)
//   {
//      cerr << "Cannot open file " << fileCoord << endl;
//      exit(EXIT_FAILURE);
//   }
// 
//   while (!feof(coord))
//   {
//      fscanf(coord, "%s %i\n", flag, &value);
//      if (strcmp("lon", flag)==0)
//      {
//           cout << "inserting " << flag << " with value " << value ; fflush(stdout);
//           longitude.push_back(value);
//           cout << "...OK" << endl;
//      }
//      else if (strcmp("lat", flag)==0)
//      {
//           cout << "inserting " << flag << " with value " << value ; fflush(stdout);
//           latitude.push_back(value);
//           cout << "...OK" << endl;
//      }
//      else
//      {
//           cerr << "Flag unknown in coordinates file" << endl;
//           exit(EXIT_FAILURE);
//      }
//   }
// fclose(coord);


  std::list<short> latitude;
  std::list<short> longitude;
  cout << "Building list of meridian and parallels" << endl;
  // if (strcmp("r", coord)==0)
  if (coord=="r")
  {
          cout << "Chosing regular grid coordinates" << endl;
          for (int x=0; x<=180; x+=20)
            latitude.push_back((short)x);
      latitude.push_back(30); latitude.push_back(150);
      for (int y=0; y<360; y+=5)
              longitude.push_back((short)y);
  }
  else  if (coord=="s")
  {
          cout << "Chosing regular grid coordinates for sulci" << endl;
          for (int x=0; x<200; x+=20)
            latitude.push_back((short)x);
      for (int y=0; y<100; y+=5)
              longitude.push_back((short)y);
  }
  //  else if (strcmp("c", coord)==0)
  else if (coord=="c")
  {
      cout << "Chosing constraint coordinates" << endl;
      latitude.push_back(0); latitude.push_back(30); latitude.push_back(55); latitude.push_back(81); latitude.push_back(92); latitude.push_back(106); latitude.push_back(150); latitude.push_back(180); 
      longitude.push_back(0); longitude.push_back(16); longitude.push_back(40); longitude.push_back(61); longitude.push_back(281); longitude.push_back(297); longitude.push_back(339); 
  }
  else if (coord=="d")
  {
          cout << "Chosing regular grid coordinates" << endl;
          for (int x=0; x<= 180; x+=5)
                  latitude.push_back((short)x);
          latitude.push_back(30); latitude.push_back(150);
          for (int y=0; y<360; y+=5)
                  longitude.push_back((short)y);
  }

    // GENERER LES TUBES ICI

  cout << "starting tube generation" << endl;
  std::list<short>::iterator coordIt;
  AimsSurfaceTriangle grille;
/*  AimsSegments grille;*/
  IsoLine mer(surface, texLat);
  mer.radius1=diam;
  mer.radius2=diam;
  for (coordIt=latitude.begin(); coordIt!=latitude.end(); ++coordIt)
  {
      AimsSurfaceTriangle meridien;
/*      AimsSegments meridien;*/
      cout << "val = " <<  *coordIt << endl;
      meridien=mer.makeTubes((short)*coordIt);
/*      meridien=mer.makeLine();*/
      SurfaceManip::meshMerge(grille, meridien);
  }
  IsoLine par(surface, texLon);
  par.radius1=diam;
  par.radius2=diam;
  for (coordIt=longitude.begin(); coordIt!=longitude.end(); ++coordIt)
  {
/*      AimsSegments parallele;*/
      AimsSurfaceTriangle parallele;
      cout << "val = " <<  *coordIt << endl;
      parallele=par.makeTubes((short)*coordIt);
/*      parallele=par.makeLine();*/
      SurfaceManip::meshMerge(grille, parallele);
  }

  cout << "saving triangulation    : " << flush;
  Writer<AimsSurfaceTriangle> triW( fileOut );
/*  Writer<AimsSegments> triW( fileOut );*/
  triW << grille;
  cout << "done" << endl;

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


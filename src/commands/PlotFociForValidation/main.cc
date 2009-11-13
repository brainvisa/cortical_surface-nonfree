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
#include <cortical_surface/mesh/pointDistance.h>
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
     string  meshfile, fociFile, fileLat="", fileLon="", fileOut;
     uint mode=0;

     AimsApplication    app( argc, argv, "Plot a set of foci on a surface, from node indexes or coordinates, and compute stats on their dispersion in the 2D space" );
     app.addOption( meshfile, "-i", "input mesh" );
     app.alias( "--input", "-i" );
     app.addOption( fileOut, "-o", "output texture");
     app.alias( "--outut", "-o" );
     app.addOption(fociFile, "-f", "foci file" );
     app.alias( "--foci", "-f");
     app.addOption( mode, "-m", "mode : 0=from indexes, 1=from coordinates", 0);
     app.alias( "--mode", "-m" );
     app.addOption( fileLat, "-x", "latitude texture", "" );
     app.alias( "--xcoord", "-x" );
     app.addOption( fileLon, "-y", "longitude texture", "" );
     app.alias( "--ycoord", "-y" );
     app.initialize();

      //
      // read triangulation

     AimsSurfaceTriangle surfaceIn;
     Reader<AimsSurfaceTriangle> triR( meshfile );
     triR >> surfaceIn;

     TimeTexture< float > lat, lon;
     if (mode==1)
     {
      if ((fileLat!="") && (fileLon!="") && (mode==1))
      {    
        Reader< TimeTexture<float> > latW(fileLat);
        latW >> lat;
        Reader< TimeTexture<float> > lonW(fileLon);
        lonW >> lon;
      }
     }
     else if (mode != 0) 
     {
        std::cerr << "Mode incompatible with coordinate files" << std::endl;
        exit(EXIT_FAILURE);
     }

     std::vector< std::pair<float, float> > listeCoord;
     std::vector< uint > listeIndex;
     std::vector<Point3df> vert=surfaceIn.vertex();

     uint nVert=vert.size();
     TimeTexture< short > texOut(1, nVert);
     for (uint i=0; i<nVert; i++)
          texOut[0].item(i)=0;
     float x, y;
     float xv, yv;
     float dista, distMin=10000.0;
     uint a_index;
     FILE *foci=fopen(fociFile.c_str(), "r");
     uint n_points;
     if (mode==0)
     {
        uint a;
        fscanf(foci, "%u\n", &a);
        listeIndex.push_back(a);
        while (!feof(foci))
        {
            fscanf(foci, "%u\n", &a);
            listeIndex.push_back(a);
        }
        n_points=listeIndex.size();
     }
     else if (mode==1)
     {
        fscanf(foci, "%f %f\n", &x, &y);
        listeCoord.push_back(std::pair<float, float>(x,y));
        while (!feof(foci))
        {
              fscanf(foci, "%f %f\n", &x, &y);
              listeCoord.push_back(std::pair<float, float>(x,y));
        }
        n_points=listeCoord.size();
  
        for ( uint i=0; i<n_points; i++)
        {
          distMin=10000.0;
          x=listeCoord[i].first;
          y=listeCoord[i].second;
/*          std::cout << "(" << x << "' " << y << ") -> ";*/
          for ( uint j=0; j<nVert; j++)
          {
               xv=lat[0].item(j);
               yv=lon[0].item(j);
               dista=sqrt( (x-xv)*(x-xv) + (y-yv)*(y-yv) );
               if (dista<distMin)
               {
                    distMin=dista;
                    a_index=j;
               }
          }
          listeIndex.push_back(a_index);
/*          std::cout << a_index << std::endl;*/
        }
     }
     else 
     {
        std::cerr << "Wrong mode" << std::endl;
        exit(EXIT_FAILURE);
     }
/*     std::cout << "Computing distances" << std::endl;*/
     std::vector<uint>::iterator listeIndexIt=listeIndex.begin(), listIndexIt2;

     MeshPointDistance dist(surfaceIn);
     std::map<uint, float> mapDist;
//      std::cout << "Computing distances" << std::endl;
     float dAv=0;
     float d;

     for (uint i=0; i<n_points; i++)
          for (uint j=i+1; j<n_points; j++)
          {
/*               cout << i << " - " << j << endl;*/
               d=dist.compute(listeIndex[i], listeIndex[j]);
/*               std::cout << "\t" << i << " -> " << j << " : " << d << std::endl;*/
//           mapDist.insert(std::pair<uint, float>(*listeIndexIt, d));
               dAv+=d;
          }

/*     std::cout << "Main loop over" << endl;*/
//      for ( listeIndexIt=listeIndex.begin() ; listeIndexIt!=listeIndex.end(); ++listeIndexIt)
//      {
//           cout << (*listeIndexIt) << endl;
//           d=dist.compute((*listeIndexIt), a_index);
//           mapDist.insert(std::pair<uint, float>(*listeIndexIt, d));
//           dAv+=d;
//      }


     float sum=(n_points*(n_points-1))/2.0;
/*     cout << "Sum=" << sum << endl;*/
     dAv/=(float) sum;
/*     std::map<uint, float>::iterator dIt=mapDist.begin();*/
//      for ( ; dIt!=mapDist.end(); ++dIt)
//      {
//           std::cout << dIt->first << " : " << dIt->second << std::endl;
//      }
/*     std::cout << "Number of distances : " << sum << std::endl;*/
     std::cout << "Average distance :" << dAv << std::endl;

     for (uint i=0; i<nVert; i++)
      texOut[0].item(i)=0;
     for (uint i=0; i<n_points; i++)
      texOut[0].item(listeIndex[i])=100;
          
/*     cout << "OK. Writing texture " << endl;*/
     Writer<TimeTexture<short> > texOutW( fileOut );
     texOutW << texOut ;
/*     std::cout << "Done" << std::endl;*/
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

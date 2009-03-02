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

uint closer(AimsSurfaceTriangle surface, Point3df p);

int main( int argc, const char** argv )
{
     string fileOut, fileMesh, fileX, fileY;

     AimsApplication    app( argc, argv, "Compute normalized depth profile from " );
     try
     {
          app.addOption( fileMesh, "-m", "input mesh" );
          app.alias( "--mesh", "-m" );
          app.addOption( fileX, "-x", "x coordinate texture" );
          app.alias( "--xcoord", "-x" );
          app.addOption( fileY, "-y", "y coordinate texture" );
          app.alias( "--ycoord", "-y" );
          app.addOption( fileOut, "-o", "output text file" );
          app.alias( "--out", "-o" );
          
          app.initialize();
          
          cout << "reading triangulation   : " << flush;
          AimsSurfaceTriangle surface;
          Reader<AimsSurfaceTriangle> triR( fileMesh );
          triR >> surface;
          cout << "done" << endl;
          
          cout << "reading textures   : " << flush;
          Reader<TimeTexture<float> > texYR( fileY );
          TimeTexture<float> texY;
          texYR >> texY ;
          Reader<TimeTexture<float> > texXR( fileX );
          TimeTexture<float> texX;
          texXR >> texX ;
          cout << "done " << endl;
          
          // GENERER LES SEGMENTS ICI ET CALCULER LES PROFONDEURS

          cout << "Computing iso-coordinates" << endl;
          std::list<short>::iterator coordIt;
/*          std::map<uint, AimsSegments> isoseg;*/
          std::map<uint, float> depth;
          IsoLine par(surface, texY);
          std::cout << "Processing values" << std::endl;
          uint i;
          for (i=1; i<=100; i++)
          {
               AimsSegments parallele;
               cout <<  i << std::endl;
               parallele=par.makeLine( i );
/*               isoseg[i]=parallele;*/
               float d1=0.0, d2=0.0;
               std::vector<Point3df>  vert=parallele.vertex();
               std::vector<AimsVector<uint,2> > poly=parallele.polygon();
               std::vector<AimsVector<uint,2> >::iterator polyIt=poly.begin();
               for ( ; polyIt!=poly.end(); ++polyIt)
               {
                    uint v1=(*polyIt)[0];
                    uint v2=(*polyIt)[1];
                    Point3df p1=vert[v1];
                    Point3df p2=vert[v2];
                    uint i1=closer(surface, p1);
                    uint i2=closer(surface, p2);
                    float x1=texX.item(i1), x2=texX.item(i2);
                    float d=norm(p1-p2);
                    if (x1<100.0)
                    {
                         if (x2<=100.0) d1+=d;
                         else {d1+=d; d2+=d;}
                    }
                    else
                    {
                         if (x2>=100.0) d2+=d;
                         else {d1+=d; d2+=d;}
                    }
               }
               std::cout << d1 << "\t" << d2 << std::endl;
               if (d1<=d2) depth.insert(std::pair<uint, float>(i, d1));
               else depth.insert(std::pair<uint, float>(i, d2));
          }

          cout << "done" << endl;
          
          ofstream fout(fileOut.c_str());
          for (i=0; i<=100; i++)
          {
               fout << i << "\t" << depth[i] << std::endl;
          }
          fout.close();

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

//---------------------------------------------------------------------------------

uint closer(AimsSurfaceTriangle surface, Point3df p)
{
     uint imin, i;
     float dmin=1000.0, d;
     std::vector<Point3df>  vert=surface.vertex();
     for (i=0; i<vert.size(); i++)
     {
          d=norm(p-vert[i]);
          if (d<dmin)
          {
               dmin=d; imin=i;
          }
     }
     return(imin);
}



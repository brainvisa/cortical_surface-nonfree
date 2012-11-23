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
     string fileOut, fileMesh, fileX, fileY, fileC;

     AimsApplication    app( argc, argv, "Compute depth and profile curve from SC " );
     try
     {
          app.addOption( fileMesh, "-m", "input mesh" );
          app.alias( "--mesh", "-m" );
          app.addOption( fileX, "-x", "x coordinate texture" );
          app.alias( "--xcoord", "-x" );
          app.addOption( fileY, "-y", "y coordinate texture" );
          app.alias( "--ycoord", "-y" );
          app.addOption( fileC, "-d", "dist to plan file" );
          app.alias( "--dist", "-d");
          app.addOption( fileOut, "-o", "output text file" );
          app.alias( "--out", "-o" );

          app.initialize();
          
          cout << "reading triangulation " << fileMesh << ": " << flush;
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
          Reader<TimeTexture<float> > texCR( fileC );
          TimeTexture<float> texC;
          texCR >> texC ;


          cout << "done " << endl;
          
          // GENERER LES SEGMENTS ICI ET CALCULER LES PROFONDEURS

          cout << "Computing iso-coordinates" << endl;
          std::list<short>::iterator coordIt;
/*          std::map<uint, AimsSegments> isoseg;*/
          std::map<uint, float> depth;
          std::map<uint, float> curvM, curvM2, sig, smooth;
          IsoLine par(surface, texY);
          std::cout << "Processing values" << std::endl;
          uint i;
          for (i=1; i<=100; i++)
          {
               AimsSegments parallele;
               // cout <<  i << std::endl;
               parallele=par.makeLine( i );
/*               isoseg[i]=parallele;*/
               float d1=0.0, d2=0.0;
               std::vector<Point3df>  vert=parallele.vertex();
               std::vector<AimsVector<uint,2> > poly=parallele.polygon();
               std::vector<AimsVector<uint,2> >::iterator polyIt=poly.begin();
               float curv=0.0, c;
               float curv2=0.0, c2;

               int ncurv2=0;
               int ncurv=0;
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
                    c=texC.item(i1);
                    if (x1<100.0)
                    {
                         if (x2<=100.0) d1+=d;
                         else {d1+=d; d2+=d;}
                         if ((x1>=0) && (x1<=100))
                         {
                        	 curv2+=c;
                             ncurv2++;
                         }
                    }
                    else
                    {
                         if (x2>=100.0) d2+=d;
                         else {d1+=d; d2+=d;}
                         if ((x1>=100) && (x1<=200))
                         {
                             curv+=c;
                             ncurv++;
                         }
                    }
               }
               if (ncurv!=0)
            	   curv=curv/float(ncurv);
               else curv=0.0;
               if (ncurv2!=0)
                   curv2=curv2/float(ncurv2);
               else curv2=0.0;
               //std::cout << d1 << "\t" << d2 << std::endl;
               if (d1<=d2) depth.insert(std::pair<uint, float>(i, d1));
               else depth.insert(std::pair<uint, float>(i, d2));
               curvM.insert(std::pair<uint, float>(i, curv));
               curvM2.insert(std::pair<uint, float>(i, curv2));
          }

          cout << "done" << endl;
          
          ofstream fout(fileOut.c_str());
          float minC=100.0;
          uint imin=0, imax;
          for (i=0; i<=100; i++)
          {
              fout << i << "\t" << depth[i] << "\t" << (curvM[i]+curvM2[i])/2.0 << std::endl;
          }
          fout.close();

          // smoothing before looking for extrema

          cout << "Smoothing profile before detecting landmarks" << endl;
       	  for (i=0; i<=100; i++)
       	  {
       		  smooth.insert(std::pair<uint, float>(i, curvM[i]));
       		  sig.insert(std::pair<uint, float>(i, curvM[i]));
       	  }

       	  float lapl;
          for (uint t=0; t<24; t++)
           {
        	  for (i=0; i<=100; i++)
        		  sig[i]=smooth[i];

         	  for (i=0; i<=100; i++)
         	  {
         		  if (i==0)
         			  lapl=sig[1]-sig[0];
         		  else if (i==100)
         			  lapl=sig[99]-sig[100];
         		  else
         			  lapl=sig[i-1]-2*sig[i]+sig[i+1];
         		  smooth[i]=sig[i]+0.25*lapl*0.5;
         	  }
           }

          cout << "looking for extrema" << endl;

          for (i=0; i<=50; i++)
          {
        	  if ((smooth[i] < minC) && (smooth[i] < smooth[i+1]) && (smooth[i] < smooth[i-1]))
        	  {
        		  imin=i;
        		  minC=smooth[i];
        	  }
          }

          imax=imin;
          for (i=imin+1; i<100; i++)
          {
        	  if ((smooth[i] > smooth[i+1]) && (smooth[i] > smooth[i-1]))
        	  {
        		  imax=i;
        		  break;
        	  }
          }
          cout << "Landmark L1: " << imin << endl;
          cout << "Landmark L2: " << imax << endl;
          cout << "\n\n\n" << endl;
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



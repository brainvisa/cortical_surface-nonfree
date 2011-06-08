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

#include <cstdlib>
#include <aims/mesh/geometric.h>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>
#include <aims/math/math_g.h>
#include <aims/vector/vector.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/curv.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>
#include <aims/connectivity/meshcc.h>
#include <aims/morphology/morphology_g.h>
#include <cortical_surface/surfacereferential/shortestPath.h>
#include <cortical_surface/mesh/linkPath.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/getopt/getopt2.h>
#include <aims/connectivity/meshcc.h>
#include <iostream>
#include <iomanip>

#include <aims/geodesicpath/geodesicPath.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace aims::meshdistance;

#define TOP2BOTTOM 0
#define BACK2FRONT 1
#define RIDGE_TOP 100
#define RIDGE_BOT 200
#define POLE_N 500
#define POLE_S 600

// SEE THE COMMENTS IN FRENCH AT THE END FOR DETAILS OF THE ALGORITHM.
// BUILDING OF THE CONSTRAINTS (POLES, TOP AND BOTTOM RIDGES) IS QUITE 
// A HEAVY TEDIOUS, BUT ROBUST, MACHINERY.


// IT HAS NOW CHANGED TO USE SHORTEST PATH CONSTRAINED BY GEOMETRY
// A LOT OF WORK BY ARNAUD LETROTER BEHIND IT
// HOPEFULLY LESS MESSY AND MORE ROBUST
// A LOT OF COMMENTS TO TAKE OUT

int main( int argc, const char** argv )
{
  try
    {
      string  meshfile, texFile_x, texFile_y, bottomfile, hullfile;
      string  method = "boix";
      float deltaT=0.05;
      float stop=0.000003;
      float offset=0.0;
      short dil=1;

/*      float tCurv=0.95;*/
      int orientation;

      AimsApplication    app( argc, argv, "Parameterize the mesh of a sulcus" );
      app.addOption( meshfile, "-i", "input mesh" );
      app.alias( "--input", "-i" );
      app.addOption( bottomfile, "-b", "sulcus bottom point image" );
      app.alias( "--bottom", "-b" );
      app.addOption( hullfile, "-t", "sulcus top point (junction with brain hull) image" );
      app.alias( "--top", "-t" );
      app.addOption( orientation, "-o", "sulcus main orientation : top->bottom = 0, back->front = 1" );
      app.alias( "--orientation", "-o" );
      app.addOption( texFile_x, "-ox", "1st output texture (depth)" );
      app.alias( "--outputx", "-ox" );
      app.addOption( texFile_y, "-oy", "2nd output texture ('along' the sulcus)" );
      app.alias( "--outputy", "-oy" );
      app.addOption( method, "-m", "curvature computation method (fem/boix/barycenter, default=boix)", "boix" );
      app.alias( "--method", "-m" );
      app.addOption( deltaT, "-d", "diffusion iteration step (default=0.05)", 0.05 );
      app.alias( "--deltaT", "-d" );
      app.addOption( stop, "-s", "laplacian variation stopping criterion (default=0.000003)", 0.000003 );
      app.alias( "--stop", "-s" );
      app.addOption( dil, "-di", "dilation of the ridges (1=yes (default), 0=no)", 1);
      app.alias("--dilation", "-di");
      app.addOption( offset, "-mo", "Morphological offset between dilation and erosion of ridges for extrem cases (default=0, otherwise should be 1.0)", 0.0);
      app.initialize();

      if ( (orientation != TOP2BOTTOM) && (orientation != BACK2FRONT) )
      {
           cerr << "-o option must be set to 0 or 1" << endl;
           return(1);
      }

      //
      // read triangulation
      //

//      std::cerr << "Reading all files" << std::endl;

      cout << "reading triangulation   : " << flush;
      AimsSurfaceTriangle surface;
      Reader<AimsSurfaceTriangle> triR( meshfile );
      triR >> surface;
      cout << "done" << endl;

      cout << "reading bottom image  : " << bottomfile << flush;
      AimsData<short> bottom;
      Reader<AimsData<short> > bottomR(bottomfile );
      bottomR >> bottom;
      cout << "done" << endl;
      cout << "reading hull image  : " << flush;
      AimsData<short> hull;
      Reader<AimsData<short> > hullR(hullfile);
      hullR >> hull;
      cout << "done" << endl;

      int x,y,z,sx,sy,sz;
      float dx, dy, dz;

      sx=bottom.dimX(); sy=bottom.dimY(); sz=bottom.dimZ();
      dx=bottom.sizeX(); dy=bottom.sizeY(); dz=bottom.sizeZ();

      AimsData<short> bottomDil(sx, sy, sz, 1, 1), hullDil(sx, sy, sz, 1, 1);
      bottomDil.setSizeXYZT(dx, dy, dz); hullDil.setSizeXYZT(dx, dy, dz);
      std::cout << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << endl;
      std::cout << "sx=" << sx << ", sy=" << sy << ", sz=" << sz << endl;

      uint i, ns=surface.vertex().size();
      TimeTexture<short> texBot(1,ns), texHull(1,ns), texPole(1,ns),
                         topDilation(1,ns), botDilation(1,ns),
                         topClosing(1,ns), botClosing(1,ns),
                         poleDilation(1, ns),
                         topLine(1,ns), botLine(1,ns);

      GeodesicPath spGeo(surface,3,30);
      GeodesicPath spGyri(surface,2,30);
      uint pS, pN;

//-------------------------------------------------------------------
      // TEST : LOOKING FOR EXTREMITIES IN A NEW WAY.

 /*     float min=1000.0, max=-1000.0;
      std::vector<uint> ext1, ext2;
      if (orientation==TOP2BOTTOM)
      {
    	  float z;
    	  for (i=0; i<ns; i++)
    	  {
    		  z=surface.vertex()[i][2];
    		  if (z<min) min=z;
    		  if (z>max) max=z;
    	  }
    	  cerr << "Found min=" << min << ", and max=" << max << endl;
    	  float t=(max-min)/5.0;
    	  for (i=0; i<ns; i++)
    	  {
    		  z=surface.vertex()[i][2];
    	      if (z<(min+t)) ext1.push_back(i);
    	      if (z>(max-t)) ext2.push_back(i);
    	  }
      }
      else if (orientation == BACK2FRONT)
      {
       	  float y;
       	  for (i=0; i<ns; i++)
       	  {
       		  y=surface.vertex()[i][1];
       		  if (y<min) min=y;
       		  if (y>max) max=y;
       	  }
       	  cerr << "Found min=" << min << ", and max=" << max << endl;
       	  float t=(max-min)/5.0;
       	  for (i=0; i<ns; i++)
       	  {
       		  y=surface.vertex()[i][1];
        	  if (y<(min+t)) ext1.push_back(i);
        	  if (y>(max-t)) ext2.push_back(i);
       	  }
      }
      TimeTexture<short> testExt(1,ns);
      for (i=0; i< ns; i++)
    	  testExt[0].item(i)=0;

      uint j, i1, i2;
      GeodesicPath sp(surface,2,5);
      float l, lmax=0.0;

      std::vector<uint>::iterator extIt1=ext1.begin(), extIt2=ext2.begin();
      cerr << "sizes: " << ext1.size() << " and " << ext2.size() << endl;
      int nb=ext1.size() * ext2.size();
      int cnt=0;

      for (; extIt1!=ext1.end(); ++extIt1)
    	  for (extIt2=ext2.begin() ; extIt2!=ext2.end(); ++extIt2)
    	  {
    		  cnt++;
    		  cerr << cnt << "/" << nb << endl;
    		  i=*extIt1;
    		  j=*extIt2;
    		  l=sp.shortestPathLength(i,j);
    		  if (l>lmax)
    		  {
    			  lmax=l;
    		      i1=i;
    		      i2=j;
    		  }
    	  }
      cerr << "Found " << i1 << ", " << i2 << ", and length " << lmax << endl;

      vector<int> listIndexVertexPathSP;

      listIndexVertexPathSP = sp.shortestPathIndiceVextex(i1,i2);
      for (i = 0; i < listIndexVertexPathSP.size(); i++)
    	  testExt[0].item(listIndexVertexPathSP[i]) = 100;


  //    for (; extIt!=ext1.end(); ++extIt)
  //  	  testExt[0].item(*extIt)=10;

//      for (extIt=ext2.begin(); extIt!=ext2.end(); ++extIt)
  //        	  testExt[0].item(*extIt)=20;

      Writer<TimeTexture<short> > wext("/Users/olivier/Desktop/lmax.tex");
      wext << testExt;
      return( 0 );*/
      /*
      cerr << "Looking for extremities (experimental) " << endl;
      uint j, i1, i2;
      GeodesicPath sp(surface,2,5);
      TimeTexture<short> testExt(1,ns);
      float l, lmax=0.0;
      for (i=0; i<ns; i++)
    	  for (j=i+1; j<ns; j++)
    	  {
    		  l=sp.shortestPathLength(i,j);
    		  if (l>lmax)
    		  {
    			  lmax=l;
    			  i1=i;
    			  i2=j;
    		  }
    	  }
      cerr << "Found " << i1 << ", " << i2 << ", and length " << lmax << endl;
      cerr << "Writing texture /Users/olivier/Dekstop/lmax.tex" << endl;

      vector<int> listIndexVertexPathSP;

      listIndexVertexPathSP = sp.shortestPathIndiceVextex(i1,i2);
      for (i=0; i< ns; i++)
    	  testExt[0].item(i)=0;

      for (unsigned i = 0; i < listIndexVertexPathSP.size(); i++)
          testExt[0].item(listIndexVertexPathSP[i]) = 100;
      Writer<TimeTexture<short> > wext("/Users/olivier/Dekstop/lmax.tex");
      wext << testExt;

      cerr << "Done" << endl;
      */
//--------------------------------------------------------------------

      //bottomDil=AimsMorphoDilation(bottom, 2.0);
      for (int z=0; z<sz; z++)
     	  for (int y=0; y<sy; y++)
     		  for (int x=0; x<sx; x++)
     		  {
					  bottomDil(x,y,z)=0;
//					  bottomDil(x,y,z)=bottom(x,y,z);
					  hullDil(x,y,z)=0;
//					  hullDil(x,y,z)=hull(x,y,z);
     		  }

      for (int z=1; z<sz-1; z++)
    	  for (int y=1; y<sy-1; y++)
    		  for (int x=1; x<sx-1; x++)
    		  {
    			  if (bottom(x,y,z)!=0)
    			  {
    				  bottomDil(x,y,z)=1;
    				  if (dil==1)
    				  {
    					  bottomDil(x,y,z+1)=1;
    					  bottomDil(x,y,z-1)=1;
    					  bottomDil(x,y-1,z)=1;
    					  bottomDil(x,y-1,z-1)=1;
    					  bottomDil(x,y-1,z+1)=1;
    					  bottomDil(x,y+1,z)=1;
    					  bottomDil(x,y+1,z-1)=1;
    					  bottomDil(x,y+1,z+1)=1;
    					  bottomDil(x-1,y,z)=1;
    					  bottomDil(x-1,y,z-1)=1;
    					  bottomDil(x-1,y,z+1)=1;
    					  bottomDil(x-1,y-1,z)=1;
    					  bottomDil(x-1,y-1,z-1)=1;
    					  bottomDil(x-1,y-1,z+1)=1;
    					  bottomDil(x-1,y+1,z)=1;
    					  bottomDil(x-1,y+1,z-1)=1;
    					  bottomDil(x-1,y+1,z+1)=1;
    					  bottomDil(x+1,y,z)=1;
    					  bottomDil(x+1,y,z-1)=1;
    					  bottomDil(x+1,y,z+1)=1;
    					  bottomDil(x+1,y-1,z)=1;
    					  bottomDil(x+1,y-1,z-1)=1;
    					  bottomDil(x+1,y-1,z+1)=1;
    					  bottomDil(x+1,y+1,z)=1;
    					  bottomDil(x+1,y+1,z-1)=1;
    					  bottomDil(x+1,y+1,z+1)=1;
    				  }
    			  }
       			  if (hull(x,y,z)!=0)
       			  {
					  hullDil(x,y,z)=1;
					  hullDil(x,y,z+1)=1;
					  hullDil(x,y,z-1)=1;
					  hullDil(x,y-1,z)=1;
					  hullDil(x,y-1,z-1)=1;
					  hullDil(x,y-1,z+1)=1;
					  hullDil(x,y+1,z)=1;
					  hullDil(x,y+1,z-1)=1;
					  hullDil(x,y+1,z+1)=1;
					  hullDil(x-1,y,z)=1;
					  hullDil(x-1,y,z-1)=1;
					  hullDil(x-1,y,z+1)=1;
					  hullDil(x-1,y-1,z)=1;
					  hullDil(x-1,y-1,z-1)=1;
					  hullDil(x-1,y-1,z+1)=1;
					  hullDil(x-1,y+1,z)=1;
					  hullDil(x-1,y+1,z-1)=1;
					  hullDil(x-1,y+1,z+1)=1;
					  hullDil(x+1,y,z)=1;
					  hullDil(x+1,y,z-1)=1;
					  hullDil(x+1,y,z+1)=1;
					  hullDil(x+1,y-1,z)=1;
					  hullDil(x+1,y-1,z-1)=1;
					  hullDil(x+1,y-1,z+1)=1;
					  hullDil(x+1,y+1,z)=1;
					  hullDil(x+1,y+1,z-1)=1;
					  hullDil(x+1,y+1,z+1)=1;
				  }
    		  }

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // computing intersection of mesh with bottom and hull image for selection of "ridges"
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//      std::cerr << "Ridges 1" << std::endl;

     cerr << "Detecting top and bottom ridges (new style)" << endl;
      
      for (i=0; i<ns; i++)
      {
          texPole[0].item(i)=short(0);
          texBot[0].item(i)=short(0);
          texHull[0].item(i)=short(0);
          Point3df vert=surface.vertex()[i];
          x=int(floor(vert[0]/dx));
          y=int(floor(vert[1]/dy));
          z=int(floor(vert[2]/dz));
          if (  (hullDil(x,y,z)!=0) || (hullDil(x,y,z+1)!=0) || (hullDil(x,y+1,z)!=0) || (hullDil(x,y+1,z+1)!=0)
             || (hullDil(x+1,y,z)!=0) || (hullDil(x+1,y,z+1)!=0) || (hullDil(x+1,y+1,z)!=0) || (hullDil(x+1,y+1,z+1)!=0) )
          {
              texHull[0].item(i)=short(RIDGE_TOP);
          }
          if (  (bottomDil(x,y,z)!=0) || (bottomDil(x,y,z+1)!=0) || (bottomDil(x,y+1,z)!=0) || (bottomDil(x,y+1,z+1)!=0)
             || (bottomDil(x+1,y,z)!=0) || (bottomDil(x+1,y,z+1)!=0) || (bottomDil(x+1,y+1,z)!=0) || (bottomDil(x+1,y+1,z+1)!=0) )
          {
              texBot[0].item(i)=short(RIDGE_BOT);
          }
      }

           // adding a second (different) pass for robustness

//     for (i=0; i<ns; i++)
//     {
//          texHull[0].item(i)=0;
//          texBot[0].item(i)=0;
//     }
     float fx, fy, fz;
     int count=0;
     for (z=0; z<sz; z++)
          for (y=0; y<sy; y++)
               for (x=0; x<sx; x++)
               {
                    if (hullDil(x,y,z) != 0)
                    {
                         float dist, distMin=10000.0;
                         uint imin=0;
                         for (i=0; i<ns; i++)
                         {
                              Point3df vert=surface.vertex()[i];
                              fx=vert[0]/dx;
                              fy=vert[1]/dy;
                              fz=vert[2]/dz;
                              dist=sqrt((fx-x)*(fx-x) + (fy-y)*(fy-y) + (fz-z)*(fz-z));
                              if (dist<distMin)
                              {
                                   distMin=dist; imin=i;
                              }
                         }
                         texHull[0].item(imin)=short(RIDGE_TOP);
                    }
                    if (bottomDil(x,y,z) != 0)
                    {
                         float dist, distMin=10000.0;
                         uint imin=0;
                         for (i=0; i<ns; i++)
                         {
                              Point3df vert=surface.vertex()[i];
                              fx=vert[0]/dx;
                              fy=vert[1]/dy;
                              fz=vert[2]/dz;
                              dist=sqrt((fx-x)*(fx-x) + (fy-y)*(fy-y) + (fz-z)*(fz-z));
                              if (dist<distMin)
                              {
                                   distMin=dist; imin=i;
                              }
                         }
                         texBot[0].item(imin)=short(RIDGE_BOT);
                         count ++;
                    }
               }

     cout << "Done... moving on to postprocessing of ridges" << endl;
     cout << "Done " << count << " texBot points" << endl;

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // dilation of ridges
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//      std::cerr << "Ridges 2" << std::endl;


      topDilation[0]=MeshDilation<short>( surface[0], texHull[0], short(0), -1, 6.0, true);
      botDilation[0]=MeshDilation<short>( surface[0], texBot[0], short(0), -1, 6.0, true);

      topClosing[0]=MeshErosion<short>( surface[0], topDilation[0], short(0), -1, 6.0-offset, true);// was 6.0
      botClosing[0]=MeshErosion<short>( surface[0], botDilation[0], short(0), -1, 6.0-offset, true);


      // New Stuff: direct pole computation:

      cerr << "Detecting extremities" << endl;
      TimeTexture<short> extremities(1,ns);
      std::vector<uint> extrV, extr_temp;
      for (i=0; i<ns; i++)
      {
    	  if ((topClosing[0].item(i) != 0) && (botClosing[0].item(i) != 0))
    	  {
    		  extremities[0].item(i)=0;
    		  extr_temp.push_back(i);
    	  }
    	  else
    		  extremities[0].item(i)=0;
      }

      std::vector<uint>::iterator ve1, ve2;

      cerr << "Size subset before clean: " << extr_temp.size() << endl;
      float div=5.0;
      float min=1000.0, max=-1000.0;
      if (orientation==TOP2BOTTOM)
      {
     	  float z;
     	  for (i=0; i<ns; i++)
     	  {
     		  z=surface.vertex()[i][2];
     		  if (z<min) min=z;
     		  if (z>max) max=z;
     	  }
     	  cerr << "Found min=" << min << ", and max=" << max << endl;
     	  float t=(max-min)/div;
     	  cerr << "Setting t at " << t << endl;

     	  for (ve1=extr_temp.begin(); ve1!=extr_temp.end(); ++ve1)
     	  {
     		 z=surface.vertex()[*ve1][2];
     		 if ((z<(min+t)) || (z>(max-t)))
     		 {
     			extrV.push_back(*ve1);
     			extremities[0].item(*ve1)=100;
     		 }
     	  }
       }
       else if (orientation == BACK2FRONT)
       {
        	  float y;
        	  for (i=0; i<ns; i++)
        	  {
        		  y=surface.vertex()[i][1];
        		  if (y<min) min=y;
        		  if (y>max) max=y;
        	  }
        	  cerr << "Found min=" << min << ", and max=" << max << endl;
        	  float t=(max-min)/div;
        	  for (ve1=extr_temp.begin(); ve1!=extr_temp.end(); ++ve1)
         	  {
         		 y=surface.vertex()[*ve1][1];
         		 if ((y<(min+t)) || (y>(max-t)))
         		 {
         			extremities[0].item(*ve1)=100;
         			extrV.push_back(*ve1);
         		 }
         	  }
       }

      cerr << "Size subset after clean: " << extrV.size() << endl;

      int nb=0; //extrV.size() * extrV.size();
      int cnt=0;
      double l, lmax=0.0;
      uint j, i1, i2;
      nb=extrV.size();

      float Tsplit=(max-min)/3.0;
      cerr << "Splitting in two (split=" << Tsplit << ")" << endl;
      std::vector<uint> extrV1, extrV2;
      ve1=(extrV.begin());
      i1=*ve1;
      extrV1.push_back(i1);
      ++ve1;
      //float xi=surface.vertex()[i1][0];
      //float yi=surface.vertex()[i1][1];
      int indo;
      if (orientation==TOP2BOTTOM)
    	  indo=2;
      else if (orientation == BACK2FRONT)
    	  indo=1;
      float zi=surface.vertex()[i1][indo];

      for (; ve1 != extrV.end(); ++ve1)
      {
          //float xj=surface.vertex()[*ve1][0];
          //float yj=surface.vertex()[*ve1][1];
          float zj=surface.vertex()[*ve1][indo];
          //if (sqrt((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi)) < Tsplit)
          if (fabs(zi-zj) < Tsplit)
        	  extrV1.push_back(*ve1);
          else
        	  extrV2.push_back(*ve1);
      }

      cerr << "Split in " << extrV1.size() << " and " << extrV2.size() << endl;

      for (ve1=extrV1.begin(); ve1!=extrV1.end(); ++ve1)
      {
  		  //cnt++;
    	  //if ((cnt%100)==0)
    		  //cerr << cnt << "/" << extrV1.size() << endl;
    	  i=*ve1;
    	  spGeo.longestPath_1_N_ind(i, extrV2, &j, &l, 0);
    	  if (l>=lmax)
    	  {
    		  lmax=l;
    	      i1=i;
    	      i2=j;
    	  }

      }

      extremities[0].item(i1)=200;
      extremities[0].item(i2)=200;

      if (orientation==TOP2BOTTOM)
      {
           if ((surface.vertex()[i1])[2] > (surface.vertex()[i2])[2])
           {
                pN=i2;
                pS=i1;
           }
           else
           {
                pN=i1;
                pS=i2;
           }
      }
      else if (orientation == BACK2FRONT)
      {
           if ((surface.vertex()[i1])[1] > (surface.vertex()[i2])[1])
           {
                pN=i2;
                pS=i1;
           }
           else
           {
                pN=i1;
                pS=i2;
           }
      }
      vector<int> listIndexVertexPathSP;

      cerr << "found index i1=" << i1 << ", i2=" << i2 << ", and lmax=" << lmax << endl;

      cerr << "OK (index " << i1 << " and " << i2 << ")" << endl;
      cerr << "Computing shortest gyri path" << endl;

      listIndexVertexPathSP = spGyri.shortestPath_1_1_ind(i1,i2);
      l=spGyri.shortestPath_1_1_len(i1,i2);
      for (i = 0; i < listIndexVertexPathSP.size(); i++)
    	  extremities[0].item(listIndexVertexPathSP[i]) = 300;

      cerr << "Length=" << listIndexVertexPathSP.size() << endl;



     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // skeletization/thinnning of ridges
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      vector< list<unsigned> > neighbourso(ns);
      neighbourso = AimsMeshOrderNode(surface[0]);
      topLine[0]=MeshSkeletization<short> ( surface[0], topClosing[0], short(RIDGE_TOP), short(0), neighbourso );
      botLine[0]=MeshSkeletization<short> ( surface[0], botClosing[0], short(RIDGE_BOT), short(0), neighbourso );


      TimeTexture<short> topRidge(1,ns), botRidge(1,ns);

      //listIndexVertexPathSP = spGyri.shortestPath_1_1_1_ind(pS,candT, pN);



      listIndexVertexPathSP = spGyri.shortestPath_1_1_ind(pS, pN, topClosing);
      cerr << "first SP done" << endl;

      for (i = 0; i < listIndexVertexPathSP.size(); i++)
      {
    	  topRidge[0].item(listIndexVertexPathSP[i]) = RIDGE_TOP;
    	  extremities[0].item(listIndexVertexPathSP[i]) = 500;
      }

      cerr << "first loop done" << endl;

      //listIndexVertexPathSP = spGyri.shortestPath_1_1_1_ind(pS,candB, pN);
      listIndexVertexPathSP = spGyri.shortestPath_1_1_ind(pS, pN, botClosing);

      cerr << "second SP done" << endl;

      for (i = 0; i < listIndexVertexPathSP.size(); i++)
      {
    	  botRidge[0].item(listIndexVertexPathSP[i]) = RIDGE_BOT;
    	  extremities[0].item(listIndexVertexPathSP[i]) = 600;

      }
      cerr << "Ridges detected" << endl;


     // Here the poles are removed from the ridges 
     topRidge[0].item(pS)=0;
     topRidge[0].item(pN)=0;
     botRidge[0].item(pS)=0;
     botRidge[0].item(pN)=0;
  
     cout << "Generating constraints" << endl;

     TimeTexture<float> coord_x(1,ns), coord_y(1,ns);
     TimeTexture<float> init(1,ns), dist(1,ns);
     TimeTexture<short> poleCC(1,ns), topCC, botCC;


     cout << "Building constraints sets" << endl;

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // Building the two sets of constraints
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      std::set<unsigned int> pole1, pole2, topR, botR;
      std::set<unsigned int>::iterator p1It, p2It;
      float p1v=0.0, p2v=100.0;
      for (i=0; i<ns; i++)
      {
           if (i==pN)
           {
                dist[0].item(i)=0;
                coord_x[0].item(i)=-100.0;
                coord_y[0].item(i)=p1v;
                pole1.insert(i);
           }
           else if (i==pS)
           {
                dist[0].item(i)=10;
                coord_x[0].item(i)=-100.0;
                coord_y[0].item(i)=p2v;
                pole2.insert(i);
           }
           else if (topRidge[0].item(i)>0)
           {
                dist[0].item(i)=0;
                coord_x[0].item(i)=0.0;
                coord_y[0].item(i)=50.0;
                topR.insert(i);
           }
          else if (botRidge[0].item(i)>0)
          {
                dist[0].item(i)=0;
                coord_x[0].item(i)=100.0;
                coord_y[0].item(i)=50.0;
                botR.insert(i);
          }
          else
          {
                dist[0].item(i)=0;
                coord_x[0].item(i)=50.0;
                coord_y[0].item(i)=50.0;
          }
     }


      cerr << "Diffusion of coordinates" << endl;


      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Diffusion
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      // distance map for y initialisation

      init[0]=meshdistance::MeshDistance<float>(surface[0], dist[0], true);
      float distM=-10;
      for (i=0; i<ns; i++)
           if (init[0].item(i)>distM)
                distM=init[0].item(i);
      for (i=0; i<ns; i++)
           if (fabs(coord_y[0].item(i) - 50)<1)
                coord_y[0].item(i)=100.0 - (100.0*init[0].item(i)/distM);

      // Getting the list of weights for laplacian computation

      std::map<unsigned, std::set< std::pair<unsigned,float> > > weights;
      weights=AimsMeshWeightFiniteElementLaplacian( surface[0], 0.9);

      // propagating the y coordinate, pole to pole

      TimeTexture<float> laplacian(1, ns), laplacian360(1,ns);

      cerr << "Diffusing y coordinate" << endl;
      int iteration=0;
      std::vector<float> laplM;
      float lMax=0, lMaxP=0;
      int flagOut=0;
      float Lmean=0, LmeanP=-10000;
      while (flagOut<1)
      {
           lMax=-10000.0;
           Lmean=0;
           laplacian[0]=AimsMeshLaplacian( coord_y[0], weights);
           for (i=0; i<ns; i++)
           {
		    coord_y[0].item(i)=coord_y[0].item(i)+deltaT*laplacian[0].item(i);
                if ((coord_y[0].item(i) != 0) && (coord_y[0].item(i) != 100) && (laplacian[0].item(i) > lMax))
                     lMax=laplacian[0].item(i);
                Lmean+=laplacian[0].item(i);

           }
           Lmean=Lmean/float(ns);

           p1It=pole1.begin(), p2It=pole2.begin();
           for ( ; p1It!=pole1.end() ; ++p1It)
                coord_y[0].item(*p1It)=p1v;
           for ( ; p2It!=pole2.end() ; ++p2It)
                coord_y[0].item(*p2It)=p2v;

           iteration++;

           if ((iteration%200) == 0)
           {
        	   cout << fabs(Lmean-LmeanP) << " | " << flush;
                if (fabs(Lmean - LmeanP)<=stop)
                {
                     cout << "reached at iteration " << iteration << endl;
                     flagOut++;
                }
                lMaxP=lMax;
                LmeanP=Lmean;
           }
     }

      cout << "\nDiffusion of y coordinate stopped after " << iteration << " iterations" << endl;

      // for propagation in the other direction, poles have to be virtually 'removed'
     // in the same time we identify the two sides of the origin meridian, that are concerned
     // by the periodicity of the coordinates
     cout << "Preparing data for diffusion of x coordinate" << endl;
      TimeTexture<short> sides(1,ns), split(1,ns);
      for (i=0; i<ns; i++)
      {
          sides[0].item(i)=0;
          split[0].item(i)=0;
      }
      std::map<unsigned, std::set< std::pair<unsigned,float> > >::iterator itW=weights.begin();
      std::map<unsigned, std::set< std::pair<unsigned, float> > > newWeights;
      for ( ; itW!=weights.end(); ++itW)
      {
           i=(*itW).first;
           std::set< std::pair<unsigned,float> > setW;
           std::set< std::pair<unsigned,float> >::iterator itN=weights[i].begin();
           for ( ; itN!=weights[i].end(); ++itN)
           {
               unsigned neigh=(*itN).first;
               float poids=(*itN).second;
               if ((neigh == pS) || (neigh == pN))
                    poids=0.0;
               setW.insert(std::pair<unsigned,float>(neigh, poids));
               if ((topRidge[0].item(i)>0) && (topRidge[0].item(neigh) == 0) && (neigh != pS)  && (neigh != pN))
               {
                    sides[0].item(neigh)=10;
               }
           }
           newWeights[i]=setW;
      }

      // splitting 'sides' in two separate components
      

//      cout << "Writing Sides" << endl;
//      Writer<TimeTexture<short> > sideW("sides.tex");
//      sideW.write(sides);
//      cout << "Done" << endl;
     cout << "splitting ... " << flush;
     split[0]=AimsMeshLabelConnectedComponent<short>( surface[0], sides[0], 9, 0);
     cout << "splitting done" << endl;


     cout << "Getting nb connected components : " << flush;
     unsigned nbSides=AimsMeshLabelNbConnectedComponent<short>( surface[0], sides[0], 10 );
     cout << "found " << nbSides << " components" << endl;
     
     if (nbSides < 2)
     {
          cerr << "Finding " << nbSides << " sides to the principal meridian (should be 2). Exiting ..." << endl;
          exit(1);
     }
     else if (nbSides > 2)
     {
          std::vector<int> conn;
          for (i=0; i<nbSides; i++)
               conn.push_back(0);
          cerr << "Finding " << nbSides << " sides to the principal meridian. Keeping the two biggest ones" << endl;
          for (i=0; i<ns; i++)
               if (split[0].item(i) > 0)
                    conn[split[0].item(i) - 1]++;
          int label1, label2, n1=0, n2=0;
          for (i=0; i<nbSides; i++)
          {
               if (conn[i] >= n2)
               {
                    if (conn[i] >= n1)
                    {
                         n1=conn[i];
                         label1=i+1;
                    }
                    else
                    {
                         n2=conn[i];
                         label2=i+1;
                    }
               }
          }
          cerr << "keeping labels " << label1 << " and " << label2 << endl;
          for (i=0; i<ns; i++)
          {
               if ((split[0].item(i) != label1) && (split[0].item(i) != label2))
                    sides[0].item(i)=0;
          }
          split[0]=AimsMeshLabelConnectedComponent<short>( surface[0], sides[0], 9, 0);
     }

     //propagating the x coordinate
     cerr << "Diffusing x coordinate" << endl;
     iteration=0;
     lMax=0; lMaxP=0;
     Lmean=0, LmeanP=-10000;
     flagOut=0;
     std::vector<Point3df> vertV=surface[0].vertex();

     FILE *laplF;
//      laplF=fopen("laplacienX.txt", "w");

     float x1=0.0, x2=0.0, y1=0.0, y2=0.0, z1=0.0, z2=0.0;
     short l1, l2;
     set<uint> set1, set2;

          // here we use the orientation to decide which side of the sulci
          // will have x between 0 and 100
     int n1=0, n2=0;
     for (uint i=0; i<ns ; i++)
     {
          Point3df vert=surface.vertex()[i];
          if (split[0].item(i)==1)
          {
               x1+=vert[0];
               y1+=vert[1];
               z1+=vert[2];
               n1++;
               set1.insert(i);
          }
          else if (split[0].item(i)==2)
          {
               x2+=vert[0];
               y2+=vert[1];
               z2+=vert[2];
               n2++;
               set2.insert(i);
          }
     }
     x1/=float(n1); y1/=float(n1); z1/=float(n1);
     x2/=float(n2); y2/=float(n2); z2/=float(n2);
     if (orientation==TOP2BOTTOM)
     {
          if (y1<y2)
          {
               l1=1; l2=2;
          }
          else
          {
               l1=2; l2=1;
          }
     }
     else if (orientation == BACK2FRONT)
     {
          if (z1>z2)
          {
               l1=1; l2=2;
          }
          else
          {
               l1=2; l2=1;
          }
     }
     set<uint>::iterator it1=set1.begin();
     set<uint>::iterator it2=set2.begin();
     for ( ; it1!=set1.end(); ++it1)
               split[0].item(*it1)=l1;
     for ( ; it2!=set2.end(); ++it2)
               split[0].item(*it2)=l2;


//     Writer<TimeTexture<short> >  verifsplit("/Users/olivier/Desktop/split.tex");
//     verifsplit.write( split);

     while (flagOut<1)
     {
          lMax=-10000.0;
          Lmean=0;
          // computing the 'periodic' laplacian
          std::map<unsigned, std::set< std::pair<unsigned,float> > >::const_iterator il,el;
          int cpt=0;
          std::set< std::pair<unsigned,float> >::iterator ip,ep;
          unsigned node, voisin;
          Point3df temp, temp_current;
          float L,weight;

          for (il=newWeights.begin(), el=newWeights.end(); il!=el; ++il, cpt++)
          {
               node = il->first;
               L = 0;
               temp=vertV[cpt];

               for ( ip = (il->second).begin(), ep = (il->second).end(); ip != ep; ++ip    )
               {
                    temp_current=vertV[(ip->first)];
                    voisin = ip->first;
                    weight = ip->second;

                    // periodicity management around the meridian
                    //Side chnge with left or right hemisphere -- not managed yet
                    if(( split[0].item(node)==1) && (topR.find(voisin)!=topR.end()) )
                    //1 for left hemisphere, 2 for right hemisphere
                         L += weight * (200.0 - coord_x[0].item(node) );
                    else if( (split[0].item(node)==2) && (topR.find(voisin)!=topR.end()) )
                         //2 for left hemisphere, 1 for righthemisphere
                         L += weight * (-coord_x[0].item(node));
                    //for all other nodes
                    else
                         L += weight * (- coord_x[0].item(node) + coord_x[0].item(voisin));
               }
               laplacian360[0].item(node) = L;
               Lmean+=L;
          }
          Lmean=Lmean/float(newWeights.size());

          for (i=0; i<ns; i++)
           {
                coord_x[0].item(i)=coord_x[0].item(i)+deltaT*laplacian360[0].item(i);
                if ((coord_x[0].item(i) != 0) && (coord_x[0].item(i) != 100) && (laplacian360[0].item(i) > lMax))
                     lMax=laplacian360[0].item(i);
           }

           std::set<unsigned>::iterator trIt=topR.begin(), brIt=botR.begin();
           for ( ; trIt!=topR.end() ; ++trIt)
                coord_x[0].item(*trIt)=0;
           for ( ; brIt!=botR.end() ; ++brIt)
                coord_x[0].item(*brIt)=100;
/*           fprintf(laplF, "%.6f\n", Lmean);*/
           iteration++;
           if ((iteration%200) == 0)
           {
                cout << fabs(Lmean-LmeanP) << " | " << flush;
                if (fabs(Lmean-LmeanP)<=stop)
                {
                     cout << "reached at iteration " << iteration << endl;
                     flagOut++;
                }
                lMaxP=lMax;
                LmeanP=Lmean;
           }
     }

//     Texture1d tempBug(1,ns);
//     for (i=0; i<ns; i++)
//    	 tempBug[0].item(i)=0;
//     std::set<unsigned>::iterator trIt=topR.begin(), brIt=botR.begin();
//     for ( ; trIt!=topR.end() ; ++trIt)
//    	 tempBug[0].item(*trIt)=100;
//     for ( ; brIt!=botR.end() ; ++brIt)
//         tempBug[0].item(*brIt)=200;
//     Writer<Texture1d>  tbW( "/Users/olivier/Desktop/tempBug.tex" );
//     tbW.write( tempBug );

     /*     fclose(laplF);*/
      cout << "\nDiffusion of x coordinate stopped after " << iteration << " iterations" << endl;

      cerr << "writing texture " << texFile_x << " : " << flush;
      Writer<Texture1d>  texxW( texFile_x );
      texxW.write( coord_x);
      cout << "done " << endl;
      cerr << "writing texture " << texFile_y << " : " << flush;
      Writer<Texture1d>  texyW( texFile_y );
      texyW.write( coord_y);
      cerr << "done " << endl;
      return( 0 );
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

/* QUELQUES DETAILS SUR L'USINE A GAZ

Top et Bottom ridges :
----------------------

- intersection buckets/maillage
- => textures texBot et texHull à valeurs 0 ou RIDGE_TOP/RIDGE_BOT

- dilatation des deux textures => topDilation/botDilation
- fermeture des deux textures => topClosing/botClosing

- squelettisation => topLine/botLine

- enlève les branches avec un plus court chemin :
     - sélection start/end points avec le "plus long plus court chemin"
     => ibn/ibs et itn/its
     - plus court chemin itn->its ou ibn->ibs en restant sur la valeur RIDGE_TOP ou RIDGE_BOT

=> textures topRidge et botRIdge


poles : 
-------

- 2 points calculés comme le milieu de [itn, ibn] et [its, ibs] => on prend les 2 points du maillage les plus près de ca => poleN et poleS
- texture de poles à valeurs POLE_N et POLE_S sur les deux points : texPole
- dilatation de texPole => poleDilation et deux listes de points, nord et sud. 
- le couple de points de pole et sud les plus éloignés sont les poles définitifs et poleDIlation ne contient plus qu'eux. 


finalisation du méridien 0/360 :
------------------------------

- On relie le topRidge aux poles
- On ferme le résultat
- Pour chacun on refait suqlette puis le plus court chemin d'un pole à l'autre.



*/

/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 *
 *  checking and fixing duplicate coordinates in the 2D surface referential
 */

#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <aims/utility/utility_g.h>
#include <aims/mesh/mesh_g.h>
#include <cortical_surface/mesh/isoLine.h>
#include <cortical_surface/mesh/pointDistance.h>
#include <aims/io/io_g.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/mesh/geometric.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>
#include <aims/distancemap/meshmorphomat_d.h>
#include <aims/connectivity/meshcc.h>
#include <aims/connectivity/meshcc_d.h>
#include <aims/math/random.h>
#include "gyriModel.h"

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileOut, fileMesh, gyriIn;
  float weight=1.0;
  int opti=1;
  int smooth=2; // was 80

  AimsApplication    app( argc, argv, "Build a regularized version of gyri from the 2D coordinate system. TEXTURE MUST BE IN FLOAT !!!" );
  try
  {
     app.addOption( fileMesh, "-m", "input mesh" );
     app.alias( "--mesh", "-m" );
     app.addOption( gyriIn, "-i", "input gyri texture" );
     app.alias( "--in", "-i" );
     app.addOption( fileOut, "-o", "output gyri texture" );
     app.alias( "--out", "-o" );
     app.addOption( weight, "-w", "data-driven term weight", 0.1);
     app.alias( "--weight", "-w" );
     app.addOption( opti, "-a", "annealing (0=no, ICM; 1=yes", 1);
     app.alias( "--anneal", "-a");
     app.addOption(smooth, "-s", "Smoothing of probability maps (default=80)", 50);
     app.alias( "--smooth", "-s");
     app.initialize();

     std::cout << "Data-driven potential weight : " << weight << std::endl;
     cout << "reading triangulation   : " << flush;
     AimsSurfaceTriangle surface;
     Reader<AimsSurfaceTriangle> triR( fileMesh );
     triR >> surface;
     cout << "done" << endl;
     
     cout << "reading gyri   : " << flush;
     Reader<TimeTexture<float> > texGyriR( gyriIn);
     TimeTexture<float> texGyriF;
     texGyriR >> texGyriF;
     cout << "done " << endl;
     
     float T0=100.0;
     float kT=0.999;

     uint size=surface.vertex().size();
     
     // debug
     cout << "size vertices: " << size << endl; 

     // Somebody tell me why gyri are in float textures !!!
     TimeTexture<int> texGyri(1, size);
     for (uint i=0; i<size; i++)
     {
    	 texGyri[0].item(i)=(int) floor(texGyriF[0].item(i)+0.2);
     }


     cout << "computing neighbours  " << endl;
     vector<set<uint> >  neigh = SurfaceManip::surfaceNeighbours(surface);
     
     TimeTexture<short> texOut(1, size);

     //float sig=0.3;

     std::cout << "Attempt to compute mean curvatures by averaging normal curvatures" << std::endl;

     surface.updateNormals();
     TimeTexture<float> normalMean(1, size);
     TimeTexture<float> normalMax(1, size);

     for (uint i=0; i<size; i++)
     {
    	 float nk=0.0;
    	 std::set<uint> voisins=neigh[i];
    	 std::set<uint>::iterator voisIt=voisins.begin();
    	 float nx, ny, nz, vx, vy, vz;
    	 float maxk=-100.0;
    	 for ( ; voisIt != voisins.end(); voisIt++)
    	 {
    		 nx=((surface.normal())[i])[0];
    		 ny=((surface.normal())[i])[1];
    		 nz=((surface.normal())[i])[2];
    		 vx=((surface.vertex())[*voisIt])[0] - ((surface.vertex())[i])[0];
    		 vy=((surface.vertex())[*voisIt])[1] - ((surface.vertex())[i])[1];
    		 vz=((surface.vertex())[*voisIt])[2] - ((surface.vertex())[i])[2];
    		 float k=(float) 2*((nx*vx)+(ny*vy)+(nz*vz))/float(vx*vx + vy*vy + vz*vz);
    		 if (k>maxk) maxk=k;
    		 nk+=k;
    	 }
    	 nk /= (float)voisins.size();
    	 normalMean[0].item(i)=nk;
    	 normalMax[0].item(i)=maxk;
//    	 if (maxk>0)
//    		 normalMax[0].item(i)=exp(-(maxk*maxk)/(2*sig*sig));
//    	 else normalMax[0].item(i)=1.0;
     }
     cout << "OK. Writing textures normalMeanCurvatureMap and normalMaxCurvatureMap " << endl;
     Writer<TimeTexture<float> > normW( "normalMeanCurvatureMap" );
     normW << normalMean ;
     Writer<TimeTexture<float> > maxW( "normalMaxCurvatureMap" );
     maxW << normalMax ;

 	GyriRegularization regul(surface, texGyri, normalMax, weight, smooth);
 	if (opti==1) regul.runAnnealing(T0, kT);
 	else regul.runICM();

// 	regul.debugCliques();
//
// 	regul.runICMdebug(14644);

	std::cout << "Writing texture " << fileOut << std::endl;
 	regul.writeGyri(fileOut);


//	for (uint i=0; i<size; i++)
//		bandeCh[0].item(i)=labelMap[i].size();
//
//     cout << "OK. Writing texture bandeChanges" << endl;
//     x
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


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

  AimsApplication    app( argc, argv, "Build a regularized version of gyri from the 2D coordinate system" );
  try
  {
     app.addOption( fileMesh, "-m", "input mesh" );
     app.alias( "--mesh", "-m" );
     app.addOption( gyriIn, "-i", "input gyri texture" );
     app.alias( "--in", "-i" );
     app.addOption( fileOut, "-o", "output gyri texture" );
     app.alias( "--out", "-o" );
     app.addOption( weight, "-w", "data-driven term weight", 1.0);
     app.alias( "--weight", "-w" );
     app.addOption( opti, "-a", "annealing (0=no, ICM; 1=yes", 1);
     app.alias( "--anneal", "-a"),
     
     app.initialize();
     

     std:cout << "Testing the uniform random draw" << std::endl;

     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;
     std::cout << UniformRandom() << std::endl;



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
     
     uint size=surface.vertex().size();

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

 	GyriRegularization regul(surface, texGyri, normalMax, weight);
 	if (opti==1) regul.runAnnealing(2000.0, 0.99);
 	else regul.runICM();


//     for (uint i=0; i<size; i++)
//    	 texOut[0].item(i)=0;
//
//			TimeTexture<int> gyri(1,size), bandeCh(1, size);
//			std::map<uint, std::set<int> > labelMap;
//
//			std::map< std::string, std::set<float> > map_global;
//
//			const char *adr_cor= adr_corl.c_str();
//
//			FILE *corres; //correspondance between contraints names and real values
//
//			if ((corres=fopen(adr_cor, "r")) == NULL)
//			{
//				cerr << "Cannot open file " << adr_cor << endl;
//				exit(EXIT_FAILURE);
//			}
//
//			std::string arg1;
//
//			//int val_contraint;
//			int val_projection;
//
//			while (!feof(corres))
//			{
//				char c[40];
//				char a[5];
//				fscanf(corres, "%s %s %i\n", a, c, &val_projection);
//				arg1=a;
//				map_global[a].insert(val_projection);
//			}
//
//			float u=0;
//			float v=0;
//
//			std::set<float>::const_iterator it_lon, it_lat, it_lon2, it_lat2;
//
//			for(uint i=0;i<size;i++)
//			{
//				gyri[0].item(i)=0;
//				int cpt=2;
//
//				u=texLat[0].item(i);
//				v=texLon[0].item(i);
//
//				it_lat=map_global["lat"].begin();
//				for(; it_lat!=map_global["lat"].end(); ++it_lat)
//				{
//					it_lat2=it_lat;
//					it_lat2++;
//					if( it_lat2!=map_global["lat"].end() )
//					{
//						it_lon=map_global["lon"].begin();
//						for(; it_lon!=map_global["lon"].end(); ++it_lon)
//						{
//							it_lon2=it_lon;
//							it_lon2++;
//							if( it_lon2!=map_global["lon"].end() )
//							{
//
//
//								if( u>(*it_lat) && u<=(*it_lat2) && v>(*it_lon) && v<=(*it_lon2) )
//								{
//									gyri[0].item(i)=cpt;
//								}
//								cpt++;
//							}
//						}
//					}
//				}
//
//				//Pole cingulaire
//
//				if(u<=30)
//					gyri[0].item(i)=1;
//
//			}
//
//			// deuxième passe pour augmenter la liste des labels possibles.
//			for(uint i=0;i<size;i++)
//			{
//				int cpt=2;
//				labelMap[i]=std::set<int>();
//				labelMap[i].insert(gyri[0].item(i));
//				u=texLat[0].item(i);
//				v=texLon[0].item(i);
//
//				it_lat=map_global["lat"].begin();
//				for(; it_lat!=map_global["lat"].end(); ++it_lat)
//				{
//					it_lat2=it_lat;
//					it_lat2++;
//					if( it_lat2!=map_global["lat"].end() )
//					{
//						it_lon=map_global["lon"].begin();
//						for(; it_lon!=map_global["lon"].end(); ++it_lon)
//						{
//							it_lon2=it_lon;
//							it_lon2++;
//							if( it_lon2!=map_global["lon"].end() )
//							{
//								if ( ( (u>((*it_lat)-bande)) && (u<=(*it_lat)) && v>((*it_lon)-bande) && v<=((*it_lon2)+bande) )
//								||   ((u>(*it_lat2)) && (u<=((*it_lat2)+bande)) && v>((*it_lon)-bande) && v<=((*it_lon2)+bande))
//								||   ((u>((*it_lat)-bande)) && (u<=((*it_lat2)+bande)) && (v>((*it_lon)-bande)) && v<=(*it_lon))
//								||   ((u>((*it_lat)-bande)) && (u<=((*it_lat2)+bande)) && (v>(*it_lon2)) && (v<=((*it_lon2)+bande))) )
//								{
//									labelMap[i].insert(cpt);
//								}
//								cpt++;
//							}
//						}
//					}
//				}
//				if(u<=30)
//				{
//					labelMap[i]=std::set<int>();
//					labelMap[i].insert(1);
//				}
//			}
//
//
////--------------------------------------------------------


//	for (uint i=0; i<size; i++)
//		bandeCh[0].item(i)=labelMap[i].size();
//
//     cout << "OK. Writing texture bandeChanges" << endl;
//     Writer<TimeTexture<int> > texResultW( "bandeChanges" );
//     texResultW << bandeCh ;
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


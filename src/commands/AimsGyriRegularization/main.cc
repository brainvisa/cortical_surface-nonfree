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
#include "gyriModel.h"

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileOut, fileMesh, fileLat, fileLon, adr_corl;
  float bande;

  AimsApplication    app( argc, argv, "Build a regularized version of gyri from the 2D coordinate system" );
  try
  {
     app.addOption( fileMesh, "-m", "input mesh" );
     app.alias( "--mesh", "-m" );
     app.addOption( fileLat, "-x", "latitude texture" );
     app.alias( "--xcoord", "-x" );
     app.addOption( fileLon, "-y", "longitude texture" );
     app.alias( "--ycoord", "-y" );
  	 app.addOption( adr_corl, "-a", "input Correspondance File");
	 app.alias( "--inCor", "-a" );
	 app.addOption( bande, "-b", "bandwidth");
	 app.alias( "--band", "-b"),
     app.addOption( fileOut, "-o", "output gyri texture" );
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

     cout << "Bandwith: " << bande << endl;
     
     cout << "computing neighbours  " << endl;
     vector<set<uint> >  neigh = SurfaceManip::surfaceNeighbours(surface);
     
     uint size=surface.vertex().size();
     TimeTexture<short> texOut(1, size);


     TimeTexture<float>	curvM, curvG, curvB, k1(1, size);

     CurvatureFactory CF;
     Curvature *cM = CF.createCurvature(surface,"boix");
     Curvature *cG = CF.createCurvature(surface,"boixgaussian");
     Curvature *cB = CF.createCurvature(surface,"barycenter");

     float sig=0.2;

     curvB[0] = cB->doIt();
     cB->regularize(curvB[0], 0.98);
     for (uint i=0; i<size; i++)
     {
    	 if (curvB[0].item(i) >=0)
    		 k1[0].item(i)=1/sqrt(2*M_PI*sig*sig);
    	 else
    		 k1[0].item(i)=(1.0/sqrt(2*M_PI*sig*sig))*exp(-(curvB[0].item(i)*curvB[0].item(i))/(2*sig*sig));
     }
/*     cout << "processing curvatures" << flush;
     curvM[0] = cM->doIt();
     curvG[0] = cG->doIt();
     cM->regularize(curvM[0], 0.95);
     cG->regularize(curvG[0], 0.95);
     std::cout << "Computing max curvature" << std::endl;
     for (uint i=0; i<size; i++)
     {
    	 float H=curvM[0].item(i);
    	 float G=curvG[0].item(i);
    	 float det=(H*H)-G;
    	 if (det>=0)
    	 {
    		 float kmax=H+sqrt(det);
    		 float kmin=H-sqrt(det);
    		 float keff;
    		 if (fabs(kmax)>=fabs(kmin)) keff=kmax;
    		 else keff=kmin;
    		 k1[0].item(i)=keff;
    //		 if (keff>=0)
   // 			 k1[0].item(i)=0;
    //		 else k1[0].item(i)=fabs(keff);
    	 }
    	 else
    		 k1[0].item(i)=0.0;
     } */


     cout << "OK. Writing textures curvatureMap and curvatureValue " << endl;
     Writer<TimeTexture<float> > texCurvW( "curvatureMap" );
     texCurvW << k1 ;
     Writer<TimeTexture<float> > texCurvBW( "curvatureValue" );
	 texCurvBW << curvB ;

     for (uint i=0; i<size; i++)
    	 texOut[0].item(i)=0;

			TimeTexture<int> gyri(1,size), bandeCh(1, size);
			std::map<uint, std::set<int> > labelMap;

			std::map< std::string, std::set<float> > map_global;

			const char *adr_cor= adr_corl.c_str();

			FILE *corres; //correspondance between contraints names and real values

			if ((corres=fopen(adr_cor, "r")) == NULL)
			{
				cerr << "Cannot open file " << adr_cor << endl;
				exit(EXIT_FAILURE);
			}

			std::string arg1;

			int val_contraint;
			int val_projection;

			while (!feof(corres))
			{
				char c[40];
				char a[5];
				fscanf(corres, "%s %s %i\n", a, c, &val_projection);
				arg1=a;
				map_global[a].insert(val_projection);
			}

			float u=0;
			float v=0;

			std::set<float>::const_iterator it_lon, it_lat, it_lon2, it_lat2;

			for(uint i=0;i<size;i++)
			{
				gyri[0].item(i)=0;
				int cpt=2;

				u=texLat[0].item(i);
				v=texLon[0].item(i);

				it_lat=map_global["lat"].begin();
				for(; it_lat!=map_global["lat"].end(); ++it_lat)
				{
					it_lat2=it_lat;
					it_lat2++;
					if( it_lat2!=map_global["lat"].end() )
					{
						it_lon=map_global["lon"].begin();
						for(; it_lon!=map_global["lon"].end(); ++it_lon)
						{
							it_lon2=it_lon;
							it_lon2++;
							if( it_lon2!=map_global["lon"].end() )
							{
								
								
								if( u>(*it_lat) && u<=(*it_lat2) && v>(*it_lon) && v<=(*it_lon2) )
								{
									gyri[0].item(i)=cpt;
								}
								cpt++;
							}
						}
					}
				}

				//Pole cingulaire

				if(u<=30)
					gyri[0].item(i)=1;

			}
			
			// deuxième passe pour augmenter la liste des labels possibles.
			for(uint i=0;i<size;i++)
			{
				int cpt=2;
				labelMap[i]=std::set<int>();
				labelMap[i].insert(gyri[0].item(i));
				u=texLat[0].item(i);
				v=texLon[0].item(i);

				it_lat=map_global["lat"].begin();
				for(; it_lat!=map_global["lat"].end(); ++it_lat)
				{
					it_lat2=it_lat;
					it_lat2++;
					if( it_lat2!=map_global["lat"].end() )
					{
						it_lon=map_global["lon"].begin();
						for(; it_lon!=map_global["lon"].end(); ++it_lon)
						{
							it_lon2=it_lon;
							it_lon2++;
							if( it_lon2!=map_global["lon"].end() )
							{
								if ( ( (u>((*it_lat)-bande)) && (u<=(*it_lat)) && v>((*it_lon)-bande) && v<=((*it_lon2)+bande) )
								||   ((u>(*it_lat2)) && (u<=((*it_lat2)+bande)) && v>((*it_lon)-bande) && v<=((*it_lon2)+bande))
								||   ((u>((*it_lat)-bande)) && (u<=((*it_lat2)+bande)) && (v>((*it_lon)-bande)) && v<=(*it_lon))
								||   ((u>((*it_lat)-bande)) && (u<=((*it_lat2)+bande)) && (v>(*it_lon2)) && (v<=((*it_lon2)+bande))) )
								{
									labelMap[i].insert(cpt);
								}
								cpt++;
							}
						}
					}
				}
				if(u<=30)
				{
					labelMap[i]=std::set<int>();
					labelMap[i].insert(1);
				}
			}
			

//--------------------------------------------------------

	GyriRegularization regul(surface, labelMap, gyri, k1);
	regul.runICM();

	for (uint i=0; i<size; i++)
		bandeCh[0].item(i)=labelMap[i].size();

		
     cout << "OK. Writing texture bandeChanges" << endl;
     Writer<TimeTexture<int> > texResultW( "bandeChanges" );
     texResultW << bandeCh ;
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


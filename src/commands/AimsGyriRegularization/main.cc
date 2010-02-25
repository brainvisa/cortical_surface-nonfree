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
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>
#include <aims/distancemap/meshmorphomat_d.h>
#include <aims/connectivity/meshcc.h>
#include <aims/connectivity/meshcc_d.h>
#include "mesh_operations.h"

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileOut, fileMesh, fileLat, fileLon, adr_corl;

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
     
     cout << "computing neighbours  " << endl;
     vector<set<uint> >  neigh = SurfaceManip::surfaceNeighbours(surface);
     
     uint i, size=surface.vertex().size();
     TimeTexture<short> texOut(1, size);

     for (uint i=0; i<size; i++)
    	 texOut[0].item(i)=0;

     MeshPointNeighborhoodFromDistance testVoisinage(surface);

     cout << "doing node 5000  " << endl;

     set<uint> listeNoeuds=testVoisinage.compute( 5000, 3.0);
     set<uint>::iterator listIt=listeNoeuds.begin();
     for ( ; listIt!=listeNoeuds.end(); listIt++)
     {
    	 texOut[0].item(*listIt)=100;
     }

     cout << "doing node 10000  " << endl;

     listeNoeuds=testVoisinage.compute( 10000, 3.0);
	 listIt=listeNoeuds.begin();
	 for ( ; listIt!=listeNoeuds.end(); listIt++)
	 {
		 texOut[0].item(*listIt)=200;
	 }
     cout << "doing node 15000  " << endl;

	 listeNoeuds=testVoisinage.compute( 15000, 3.0);
	 listIt=listeNoeuds.begin();
	 for ( ; listIt!=listeNoeuds.end(); listIt++)
	 {
		 texOut[0].item(*listIt)=300;
	 }
	 cout << "OK. Writing texture test " << endl;
	 Writer<TimeTexture<short> > texTestW( "testOlive.3.0" );
	 texTestW << texOut ;
	 return EXIT_SUCCESS;


//-------------------------------------------------------

			TimeTexture<int> gyri(1,size);
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
								if ( ( (u>(*it_lat)-5.0) && (u<=(*it_lat)) && v>(*it_lon) && v<=(*it_lon2) )
								||   ((u>(*it_lat2)) && (u<=(*it_lat2)+5.0) && v>(*it_lon) && v<=(*it_lon2))
								||   ((u>(*it_lat)) && (u<=(*it_lat2)) && (v>(*it_lon)-5.0) && v<=(*it_lon))
								||   ((u>(*it_lat)) && (u<=(*it_lat2)) && (v>(*it_lon2)) && (v<=(*it_lon)+5.0)) )
								{
									labelMap[i].insert(cpt);
								}
								cpt++;
							}
						}
					}
				}
			}
			

//--------------------------------------------------------

	for (uint i=0; i<size; i++)
		gyri[0].item(i)=labelMap[i].size();	
		
		
     cout << "OK. Writing texture " << fileOut << endl;
     Writer<TimeTexture<int> > texResultW( fileOut );
     texResultW << gyri ;
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


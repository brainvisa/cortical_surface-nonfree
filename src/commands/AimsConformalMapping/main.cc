/*
 *  Copyright (C) 1997-2005 CEA
 *
 *  This software and supporting documentation were developed by
 *
 *      Laboratoire LSIS, �uipe LXAO,
 *      Marseille, France
 *
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 */
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <string.h>
#include <aims/def/path.h>
#include <cartobase/stream/fileutil.h>

#include <cortical_surface/surfacereferential/conformalMapping.h>

using namespace aims;
using namespace carto;
using namespace std;


int main( int argc, const char** argv )
{
	try
	{
		AimsSurfaceTriangle meshIn;
                AimsSurfaceTriangle meshOut;
                
                std::string addMeshIn, addMeshOut, addMeshTuette="nofile";
                float dt=0.1, de=0.001;

		AimsApplication    app( argc, argv, "SPherical conformalpParameterization of a mesh" );
		app.addOption( addMeshIn, "-i", "input mesh" );
		app.alias( "--input", "-i" );
		app.addOption( dt, "-dt", "time step (recommended : 0.1, nothing above that for stability)", 0.1 );
		app.alias( "--deltat", "-dt" );
		app.addOption( de, "-de", "energy variation threshold for convergence (recommended : 0.001)", 0.001 );
		app.alias( "--deltae", "-de" );
		app.addOption( addMeshOut, "-o", "output (spherical) mesh" );
		app.alias( "--output", "-o" );
                app.addOption( addMeshTuette, "-t", "Tuette map precomputed mesh", "nofile" );
		app.alias( "--tuette", "-t" );
		app.initialize();
		
                ConformalMapping map(addMeshIn, dt, de);
		
                std::cout << "Getting Conformal mapping" << std::endl;
                meshOut=map.GetConformalMapping(addMeshTuette);
                std::cout << "OK. Writing it" << std::endl;
                std::cout << "OK. Writing coordinates" << std::endl;
                map.WriteConformalCoordinates();
                std::cout << "OK" << std::endl;
		Writer<AimsSurfaceTriangle> meshW(addMeshOut);
		meshW.write(meshOut);
		
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

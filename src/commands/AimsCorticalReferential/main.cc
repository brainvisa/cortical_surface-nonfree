#include <cstdlib>
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <string.h>
#include <aims/io/aimsGraphW.h>
#include <aims/def/path.h>
#include <cartobase/stream/fileutil.h>

#include <cortical_surface/mesh/isoLine.h>
#include <cortical_surface/surfacereferential/corticalReferential.h>


int main(int argc, const char **argv)
{
    try
    {
    	std::string adress_mesh;
    	std::string adress_texMer;
    	std::string adress_texPar;
    	std::string adress_texCall;
    	std::string adress_texPoles;
    	std::string x;
    	std::string y;
    	float criter;
    	float dt;
    	float t;
    	int c;
    	int tBeta;
    	int choice;
    	bool doInsulaParameterization = true;
    	
    	AimsSurfaceTriangle mesh;
    	aims::AimsApplication app ( argc, argv, "Diffusion of mesh-based coordinates with constraints" );
    	app.addOption ( adress_mesh, "-i", "input Surface Mesh");
    	app.addOption ( adress_texPar, "-p", "input Texture Parallel Constraints" );
    	app.addOption ( adress_texMer, "-m", "input Texture Meridian Constraints" );
    	app.addOption ( adress_texCall, "-l", "input Texture Corpus Callosum Pole" );
    	app.addOption ( adress_texPoles, "-r", "input Texture Poles" );
    	app.addOption ( criter, "-c", "Max difference stop criterium for the diffusion process ( default = 1e-6 )" );
    	app.addOption ( dt, "-d", "Iterative step ( default = 0.2 )");
    	app.addOption ( c, "-b", "Context (in order to swap 0/360 around central sulcus if necessary)" );
    	app.addOption ( choice, "-f", "process choice : 1 latitude only 2 longitude only 3 none 0 latitude and longitude" );
    	app.addOption ( t, "-t", "data-driven weight" );
    	app.addOption ( tBeta, "-a", "Beta Map (1=Yes, 0=No)" );
    	app.addOption ( x, "-x", "output Latitude Texture" );
    	app.addOption ( y, "-y", "output Longitude Texture" );
    	app.addOption ( doInsulaParameterization, "--insula", "do parameterize insula", true );
    	app.initialize();

        app.alias ( "--inMesh", "-i" );
        app.alias ( "--inTexPar", "-p" );
        app.alias ( "--inTexMer", "-m" );
        app.alias ( "--inTexCall", "-l" );
        app.alias ( "--inTexPoles", "-r" );
        app.alias ( "--criter", "-c" );
        app.alias ( "--dt", "-d" );
        app.alias ( "--bool", "-b" );
        app.alias ( "--choice", "-f" );
        app.alias ( "--attache", "-t" );
        app.alias ( "--tBeta", "-a" );
        app.alias ( "--latitude", "-x" );
        app.alias ( "--longitude", "-y" );

    	aims::CorticalReferential cr ( adress_mesh, 
    	                               adress_texPar, 
    	                               adress_texMer, 
    	                               adress_texCall, 
    	                               adress_texPoles, 
    	                               criter, 
    	                               dt,
    	                               c,
    	                               choice,
    	                               t,
    	                               tBeta,
    	                               x,
    	                               y,
    	                               doInsulaParameterization );
   
    	cr.process();
    	return 0;

    }
    catch( carto::user_interruption & )
    {
    }
    catch( std::exception & e )
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return( EXIT_SUCCESS );
}


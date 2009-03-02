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

using namespace aims;
using namespace carto;
using namespace std;


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
	AimsSurfaceTriangle mesh;
	AimsApplication     app( argc, argv, "Diffusion de coordonnees sur une surface contrainte");
	app.addOption( adress_mesh, "-i", "inputMesh");
	app.alias( "--inMesh", "-i" );
	app.addOption( adress_texPar, "-p", "input Texture Parallel Constraints");
	app.alias( "--inTexPar", "-p" );
	app.addOption( adress_texMer, "-m", "input Texture Meridian Constraints");
	app.alias( "--inTexMer", "-m" );
	app.addOption( adress_texCall, "-l", "input Texture Corpus Callosum Pole", "pole_calleux_plein.tex");
	app.alias( "--inTexCall", "-l" );
	app.addOption( adress_texPoles, "-r", "input Texture Poles");
	app.alias( "--inTexPoles", "-r" );
	app.addOption( criter, "-c", "Limit for the diffusion process ( default = 1e-6 )");
	app.alias( "--criter", "-c" );
	app.addOption( dt, "-d", "Iterative step ( default = 0.2 )");
	app.alias( "--dt", "-d" );
	app.addOption( c, "-b", "context");
 	app.alias( "--bool", "-b" );
	app.addOption( choice, "-f", "choice");
 	app.alias( "--choice", "-f" );
	app.addOption( t, "-t", "attache");
	app.alias( "--attache", "-t" );
	app.addOption( tBeta, "-a", "Beta Map (1=Yes, 0=No)");
	app.alias( "--tBeta", "-a" );
	app.addOption( x, "-x", "latitude texture output");
 	app.alias( "--latitude", "-x" );
	app.addOption( y, "-y", "longitude texture output");
 	app.alias( "--longitude", "-y" );
	app.initialize();

	CorticalReferential cr(adress_mesh, adress_texPar, adress_texMer, adress_texCall, adress_texPoles, criter, dt,c,choice,t,tBeta,x,y);
	cr.process();
	return 0;

    }
   catch( user_interruption & )
    {
    }
  catch( exception & e )
    {
      cerr << e.what() << endl;
      return EXIT_FAILURE;
    }

  return( EXIT_SUCCESS );
}


#include <cstdlib>
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <string.h>
#include <aims/def/path.h>
#include <cartobase/stream/fileutil.h>

#include <cortical_surface/surfacereferential/sulcusCleaner.h>

using namespace aims;
using namespace carto;
using namespace std;


int main(int argc, const char **argv)
{
  try
  {
	  std::string adr_mesh;
	  std::string adr_calleux;
	  std::string adr_poles;
	  std::string adr_mer;
	  std::string adr_par;
	  std::string adr_corres;
	  std::string adr_file;
	  std::string adr_long_cleaned;
	  std::string adr_lat_cleaned;
	  std::string side;
	  int context;
	  float contr;
	  float curvature;
	  float elasticity;
			  
	AimsSurfaceTriangle mesh;
	AimsApplication     app( argc, argv, "Cortical Constraints Cleaning (for cortical surface coordinate system)");
	app.addOption( adr_mesh, "-m", "input Mesh");
	app.alias( "--inMesh", "-m" );
	app.addOption( adr_calleux, "-t", "input Texture Cingular Pole");
	app.alias( "--inTexCingularPole", "-t" );
	app.addOption( adr_poles, "-p", "Poles Textures");
	app.alias( "--inTexPoles", "-p" );
	app.addOption( adr_mer, "-x", "input Texture Meridian Constraints");
	app.alias( "--inTexMer", "-x" );
	app.addOption( adr_par, "-y", "input Texture Parallel Constraints");
	app.alias( "--inTexPar", "-y" );
	app.addOption( adr_corres, "-f", "input Constraint Correspondance File");
	app.alias( "--inCorrFile", "-f" );
	app.addOption( adr_file, "-g", "input Projected constraint Correspondance File");
	app.alias( "--inProjFile", "-g" );
	app.addOption( adr_long_cleaned, "-a", "outut Texture Cleaned Meridian Constraints");
	app.alias( "--outMer", "-a" );
	app.addOption( adr_lat_cleaned, "-b", "output Texture Cleaned Parallel Constraints");
	app.alias( "--outPar", "-b" );
	app.addOption( contr, "-i", "constraint distance parameter");
	app.alias( "--constraint_distance", "-i" );
	app.addOption( curvature, "-j", "curvature parameter");
	app.alias( "--curvature", "-j" );
	app.addOption( elasticity, "-k", "elasticity parameter");
	app.alias( "--elasticity", "-k" );
	app.addOption( context, "-c", "context");
	app.alias( "--context", "-c" );
	app.addOption( side, "-s", "side (left or right)");
	app.alias( "--side", "-s" );
	app.initialize();

	SulcusCleaner sc( adr_mesh, adr_calleux, adr_poles, adr_mer, adr_par, adr_long_cleaned, adr_lat_cleaned, adr_corres, adr_file, context, side, contr, curvature, elasticity );
	sc.processConstraints();
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


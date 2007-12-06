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

#include <aims/surfacereferential/corticalTools.h>

using namespace aims;
using namespace carto;
using namespace std;


int main(int argc, const char **argv)
{
  try
    {
	std::string adress_long;
	std::string adress_lat;
	std::string adress_output;
	std::string adr_cor;
	std::string side;
	
	AimsApplication     app( argc, argv, "Creation de gyri utilisant le systeme de coordonnees surfacique AimsCorticalReferential");
	app.addOption( adress_long, "-x", "input Longitude Texture");
	app.alias( "--inLon", "-x" );
	app.addOption( adress_lat, "-y", "input Latitude Texture");
	app.alias( "--inLat", "-y" );
	app.addOption( adr_cor, "-a", "input Correspondance File");
	app.alias( "--inCor", "-a" );
	app.addOption( adress_output, "-o", "output texture");
	app.alias( "--outTex", "-o" );
	app.initialize();

	//InitValues i(adress_dir, adress_mesh);
	
	drawGyri(adress_long,adress_lat,adress_output,adr_cor,side);
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


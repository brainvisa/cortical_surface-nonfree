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

#include <cortical_surface/surfacereferential/sulcalLinesGeodesic.h>

using namespace aims;
using namespace carto;
using namespace std;

int main(int argc, const char **argv)
{
  try
  {
    std::string adrMesh;

    std::string adrGeodesicDepth;

    std::string adrRootsLon;
    std::string adrRootsLat;

    std::string adrLinesOut;
    std::string adrBassinsOut;

    std::string adrLonGeodesicOut;
    std::string adrLatGeodesicOut;

	  int strain = 3;

    AimsSurfaceTriangle mesh;
    AimsApplication     app( argc, argv, "Cortical Sulcal Lines (for cortical surface coordinate system)");

    app.addOption( adrMesh, "-m", "input Mesh");
    app.alias( "--inMesh", "-m" );

    app.addOption( adrGeodesicDepth, "-d", "input Texture Geodesic Depth",true);
    app.alias( "--inTexGeoDepth", "-d" );

    app.addOption( adrRootsLon, "-im", "input Texture Meridian Constraints");
    app.alias( "--inTexMer", "-im" );
    app.addOption( adrRootsLat, "-ip", "input Texture Parallel Constraints");
    app.alias( "--inTexPar", "-ip" );

    app.addOption( adrLinesOut, "-o", "outut Texture all sulcal lines");
    app.alias( "--outTexLines", "-o" );

    app.addOption( adrBassinsOut, "-b", "outut Texture all bassins regions");
    app.alias( "--outTexBassins", "-b" );

    app.addOption( adrLonGeodesicOut, "-om", "outut Texture sulcal lines meridian constraints (longitude)");
    app.alias( "--outTexMer", "-om" );
    app.addOption( adrLatGeodesicOut, "-op", "output Texture sulcal lines parallel constraints (latitude)");
    app.alias( "--outTexPar", "-op" );

    app.addOption( strain, "-st", "strain parameter (3 by default)",true );
    app.alias( "--strain", "-st" );

    app.initialize();

    SulcalLinesGeodesic slg( adrMesh, adrGeodesicDepth, adrRootsLon, adrRootsLat, adrLinesOut, adrBassinsOut, adrLonGeodesicOut, adrLatGeodesicOut, strain );
    slg.run();

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


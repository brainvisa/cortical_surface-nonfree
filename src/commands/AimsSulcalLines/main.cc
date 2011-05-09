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
    string adrMesh;

    string adrRootsLon;
    string adrRootsLat;

    int extremeties_method = 1;
    int constraint_type = 1;

    string adrCurv = "";
    string adrGeodesicDepth = "";

	  int strain = 3;
	  float proba = 0.4;

	  bool save = false;

    AimsSurfaceTriangle mesh;
    AimsApplication     app( argc, argv, "Cortical Sulcal Lines (for cortical surface coordinate system)");

    app.addOption( adrMesh, "-i", "input Mesh");
    app.alias( "--inMesh", "-i" );

    app.addOption( adrRootsLon, "-lon", "input Texture Longitude Constraints");
    app.alias( "--inTexLon", "-lon" );
    app.addOption( adrRootsLat, "-lat", "input Texture Latitude Constraints");
    app.alias( "--inTexLat", "-lat" );

    app.addOption( extremeties_method, "-m", "extraction of extremities method :\n1 : projection crop by basins (by default)\n2 : map of probability\n",true);
    app.alias( "--inMethod", "-m" );

    app.addOption( adrCurv, "-c", "input Texture Curvature (barycenter curvature by default)",true);
    app.alias( "--inTexCurv", "-c" );

    app.addOption( adrGeodesicDepth, "-d", "input Texture Geodesic Depth",true);
    app.alias( "--inTexGeoDepth", "-d" );

    app.addOption( constraint_type, "-t", "constraint type (shortest path) :\n1 : on curvature map (by default)\n2 : on depth map\n",true);
    app.alias( "--inConstraint", "-t" );

    app.addOption( strain, "-st", "strain parameter (3 by default)",true );
    app.alias( "--strain", "-st" );

    app.addOption( proba, "-p", "threshold of probability (0.4 by default)",true );
    app.alias( "--proba", "-p" );

    app.addOption( save, "-s", "save all textures", true );
    app.alias( "--save", "-s" );

//    app.addOption( adrBassinsDepthNorm, "-bn", "input Texture Bassins Depth Normalized",true);
//    app.alias( "--inTexBassinsDepthNorm", "-bn" );

//    app.addOption( adrLinesOut, "-o", "outut Texture all sulcal lines",true);
//    app.alias( "--outTexLines", "-o" );
//
//    app.addOption( adrBassinsOut, "-b", "outut Texture all bassins regions",true);
//    app.alias( "--outTexBassins", "-b" );
//
//    app.addOption( adrLonGeodesicOut, "-om", "outut Texture sulcal lines meridian constraints (longitude)",true);
//    app.alias( "--outTexMer", "-om" );
//    app.addOption( adrLatGeodesicOut, "-op", "output Texture sulcal lines parallel constraints (latitude)",true);
//    app.alias( "--outTexPar", "-op" );
//


    app.initialize();

    SulcalLinesGeodesic slg( adrMesh, adrCurv, adrGeodesicDepth, adrRootsLon, adrRootsLat, extremeties_method, constraint_type, strain, proba, save );
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


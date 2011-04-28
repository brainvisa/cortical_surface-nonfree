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
    int segmentation_method = 1;

    string adrCurv = "";
    string adrGeodesicDepth = "";

//    std::string adrLinesOut;
//    std::string adrLinesLonOut;
//    std::string adrLinesLatOut;
//
//    std::string adrProbaLatOut;
//    std::string adrProbaLonOut;

	  int strain = 3;

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

    app.addOption( segmentation_method, "-s", "segmentation method :\n1 : on curvature map (by default)\n2 : on depth map\n",true);
    app.alias( "--inSegmentation", "-s" );

    app.addOption( strain, "-st", "strain parameter (3 by default)",true );
    app.alias( "--strain", "-st" );


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

    SulcalLinesGeodesic slg( adrMesh, adrCurv, adrGeodesicDepth, adrRootsLon, adrRootsLat, extremeties_method, segmentation_method, strain );
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


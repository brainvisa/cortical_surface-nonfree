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
    string adrMesh = "";
    string adrRootsBottom = "";
    string adrRootsLon = "";
    string adrRootsLat = "";
    string adrLabelBasins= "";
    string adrLabelSulcalines= "";

    string adrSulcalines= "";

    string adrSaveFolder = "";

    int extremeties_method = 1;
    int constraint_type = 1;

    string adrCurv = "";
    string adrGeodesicDepth = "";

    int strain = 15;
    float proba = 0.4;
    vector<float> proba_list(0.4);

    bool save = false;

    float curv_threshold = 0.0;

    string side;
    float threshold_size_basin = 50.0;

    int constraintValue = 1;

    int max_extremities = 2;

    AimsSurfaceTriangle mesh;

    AimsApplication     app( argc, argv, "Cortical Sulcal Lines (for cortical surface coordinate system)");

    app.addOption( adrMesh, "-i", "input Mesh");
    app.alias( "--inMesh", "-i" );

    app.addOption( extremeties_method, "-m", "extraction of extremities method :\n1 : projection crop by basins \n2 : map of probability (embc11 method)\n3 : "
    		"	map of density (NeuroImage method, by default) \n 4 : mGPDM (miccai) \n 5 : map of probability (basin user defined) : ",true);
    app.alias( "--inMethod", "-m" );

    app.addOption( adrSulcalines, "-o", "output sulcal lines texture (.tex)",true );
    app.alias( "--inSulcalines", "-o" );

    app.addOption( adrCurv, "-c", "input Texture Curvature (barycenter curvature by default)",true);
    app.alias( "--inTexCurv", "-c" );

    app.addOption( adrGeodesicDepth, "-d", "input Texture Geodesic Depth",true);
    app.alias( "--inTexGeoDepth", "-d" );

    app.addOption( constraint_type, "-t", "constraint type (shortest path) :\n1 : on curvature map (by default)\n2 : on depth map\n",true);
    app.alias( "--inConstraint", "-t" );

    app.addOption( strain, "-st", "strain parameter (15 by default)",true );
    app.alias( "--strain", "-st" );

    app.addOption( adrRootsBottom, "-b", "sulcus bottom point volume",true );
    app.alias( "--bottom", "-b" );

    app.addOption( adrLabelBasins, "-lb", "input label of basins (.txt)",true );
    app.alias( "--inLabelBasins", "-lb" );

    app.addOption( adrLabelSulcalines, "-ls", "input file : label-constraint correspondances of Sulcalines (.txt)",true );
    app.alias( "--inLabelSulcalines", "-ls" );

    app.addOption( curv_threshold, "-ct", "curvature threshold for basins segmentation (0.0 by default)",true );
    app.alias( "--curvthresh", "-ct" );

    app.addOption( adrSaveFolder, "-s", "folder path for save texture", true );
    app.alias( "--save", "-s" );

    app.addOption( side, "-si", "side of hemisphere (left, right, both)",true );
    app.alias( "--side", "-si" );

    app.addOption( threshold_size_basin, "-sb", "threshold of basins size (50 by default)",true );
    app.alias( "--size_basins", "-sb" );

    app.addOption( constraintValue, "-cv", "constraint value :\n1 Basin (by default) :\n2 LatLon value",true );
    app.alias( "--constraintvalue", "-cv" );

    app.addOption(max_extremities , "-e", "max number of extremities (2 by default)",true );
    app.alias( "--max_extremities", "-e" );

    app.addOptionSeries( proba_list, "-p", "threshold of probability (0.4 by default)", false) ;

    app.addOption( adrRootsLon, "-lon", "input Texture Longitude Constraints",true);
    app.alias( "--inTexLon", "-lon" );
    app.addOption( adrRootsLat, "-lat", "input Texture Latitude Constraints",true);
    app.alias( "--inTexLat", "-lat" );


    app.initialize();

    SulcalLinesGeodesic slg( adrMesh, adrCurv, adrGeodesicDepth, adrRootsLon, adrRootsLat, adrRootsBottom, adrLabelBasins,adrLabelSulcalines, adrSulcalines,
        extremeties_method, constraint_type, strain, proba_list, adrSaveFolder, curv_threshold,side,threshold_size_basin,constraintValue,max_extremities);

    if (extremeties_method < 5 )
      slg.run();

    if (extremeties_method == 5)
      slg.probaMap();


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


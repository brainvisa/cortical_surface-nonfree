#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <string.h>

#include <cortical_surface/mesh/isoLine.h>
//#include <cortical_surface/surfacereferential/autoConstraints.h>

using namespace aims;
using namespace carto;
using namespace std;/*
BEGIN_USAGE(usage)
  "-------------------------------------------------------------------------",
  "AimsMeshIsoLine          -i[nput] <mesh file input>                      ",
  "                         -t[exture]    <texture input>		    ",
  "                         -v[alue]    <isoline value>		            ",
  "                         -o[utput] <output mesh file>]                   ",
  "-------------------------------------------------------------------------",
  " Generates a mesh representing an isodensity line on a textured mesh.     ",
END_USAGE


//
// Usage
//
void Usage( void )
{
  AimsUsage( usage );
}*/



int main(int argc, const char **argv) //int argc, const char **argv)
{
	//DECLARATIONS
	std::string adressTex="./";
	std::string adressMesh="./";
	std::string adressOutput="./";
	int value;
	AimsSurfaceTriangle meshResult;
	AimsSurfaceTriangle mesh;
	TimeTexture<float> texOriginal;

	AimsApplication     app( argc, argv, "Create an isoline mesh (tube) for a textured mesh");
        try
        {
	app.addOption( adressMesh, "-i", "input mesh");
	app.alias( "--inputMesh", "-i" );
	app.addOption( adressTex, "-t", "input texture (TimeTexture<float>)");
	app.alias( "--inputTex", "-t" );
	app.addOption( adressOutput, "-o", "output mesh");
	app.alias( "--output", "-o" );
	app.addOption( value, "-v", "value of the isoline");
	app.alias( "--value", "-v" );

	app.initialize();

	Reader < AimsSurfaceTriangle > rm(adressMesh);
	rm.read( mesh );
	Reader < TimeTexture<float> > rt(adressTex);
	rt.read( texOriginal );

	IsoLine mt(mesh, texOriginal);
	meshResult=mt.makeTubes(value);

	Writer<AimsSurfaceTriangle> wm(adressOutput);
	wm.write(meshResult);
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

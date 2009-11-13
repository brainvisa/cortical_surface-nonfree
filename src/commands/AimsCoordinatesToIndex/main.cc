#include <cstdlib>
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/math/eigen.h>
#include <string.h>

using namespace aims;
using namespace carto;
using namespace std;

int main(int argc, const char **argv) //int argc, const char **argv)
{
  //DECLARATIONS
  std::string adressMesh="./";
  std::string texLat="./";
  std::string texLon="./";
  float px, py;

  AimsApplication     app( argc, argv, "get the index of a point from its 2D surface coordinates");
  try
  {
    app.addOption( adressMesh, "-i", "input mesh");
    app.alias( "--inputMesh", "-i" );
    app.addOption( texLat, "-l", "latitude texture");
    app.alias( "--lat", "-l" );
    app.addOption( texLon, "-L", "longitude texture");
    app.alias( "--lon", "-L" );
    app.addOption( px, "-x", "latitude");
    app.addOption( py, "-y", "longitude");
    app.initialize();

   
    std::cout << "Reading mesh" << endl;
    AimsSurfaceTriangle mesh;
    Reader < AimsSurfaceTriangle > rm(adressMesh);
    rm.read( mesh );
    std::cout << "OK" << endl;

    std::cout << "Reading textures" << endl;
    TimeTexture< float > lat, lon;
    Reader < TimeTexture< float > > latW(texLat);
    Reader < TimeTexture< float > > lonW(texLon);
    latW.read( lat ); lonW.read( lon );
    std::cout << "OK" << endl;

    uint ns=mesh.vertex().size();
    float x, y, dist, distMin=10000.0;
    uint index;

    for (uint i=0; i<ns; i++)
    {
      x=lat[0].item(i);
      y=lon[0].item(i);
      dist=sqrt( (x-px)*(x-px) + (y-py)*(y-py) );
      if (dist<distMin)
      {
        distMin=dist;
        index=i;
      }
    }

    std::cout << "Point index is : " << index << std::endl;
    return(0);
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

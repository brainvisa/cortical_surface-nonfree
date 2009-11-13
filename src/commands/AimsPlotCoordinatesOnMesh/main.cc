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
  std::string adressTexX="./";
  std::string adressTexY="./";
  std::string outTex="./";
  float px, py;	

  int local=0;

  AimsApplication     app( argc, argv, "Plot a point defined by 2D coordinates by creating a texture");
  try
  {
    app.addOption( adressMesh, "-i", "input mesh");
    app.alias( "--inputMesh", "-i" );
    app.addOption( adressTexX, "-l", "latitude texture");
    app.alias( "--lat", "-l" );
    app.addOption( adressTexY, "-L", "longitude texture");
    app.alias( "--lon", "-L" );
    app.addOption( px, "-x", "point latitude");
    app.alias( "--plat", "-x");
    app.addOption( py, "-y", "point longitude");
    app.alias( "--plon", "-y");
    app.addOption( outTex, "-o", "output texture", 0);
    app.alias( "--outT", "-o");
    app.initialize();

    if ((local != 0) && (local !=1))
    {
          cerr << "-l option must be 0 or 1" << endl;
          exit(1);
    }

    std::cout << "Reading mesh" << endl;
    AimsSurfaceTriangle mesh;
    Reader < AimsSurfaceTriangle > rm(adressMesh);
    rm.read( mesh );

    std::cout << "Reading coordinate textures" << endl;
    TimeTexture< float > texLat, texLon;
    Reader< TimeTexture< float > > latw(adressTexX), lonw(adressTexY);
    latw.read(texLat); lonw.read(texLon);

    uint ns=mesh.vertex().size();
    float dist, distMin=10000.0; 
    float x, y;
    uint index;
    for (uint i=0; i<ns; i++)
    {
      x=texLat[0].item(i); 
      y=texLon[0].item(i); 
      dist=sqrt((px-x)*(px-x) + (py-y)*(py-y));
      if (dist<distMin)
      {
	distMin=dist; index=i;
      }
    }
    std::cout << "index is " << index << std::endl;
    TimeTexture<short> out(1,ns);
    for (uint i=0;  i<ns; i++)
    {
      if (i==index)
	out[0].item(i)=100;
      else out[0].item(i)=0;
    }

    // �criture des textures

    Writer< TimeTexture<short> > wt(outTex);
    wt.write( out );
    std::cout << "Done" << endl;
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

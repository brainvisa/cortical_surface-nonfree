#include <cstdlib>
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <cortical_surface/mesh/meshToMeshResample.h>
#include <string.h>

using namespace aims;
using namespace carto;
using namespace std;

int main(int argc, const char **argv) //int argc, const char **argv)
{
     //DECLARATIONS
     std::string adressTexIn="./";
     std::string adressMesh="./";
     std::string adressTexOut="./";
     std::string adressAtlas="./";
     std::string adressIx="./";
     std::string adressIy="./";
     std::string adressAx="./";
     std::string adressAy="./";
     float px=0;

     AimsApplication     app( argc, argv, "Projet a texture from one mesh onto an atlas using spherical parameterization of both. x-coordinate is the longitude (with a period)");

     try
     {
      app.addOption( adressMesh, "-i", "input mesh");
      app.alias( "--inputMesh", "-i" );
      app.addOption( adressTexIn, "-t", "input texture (float)");
      app.alias( "--inputTex", "-t" );
      app.addOption( adressTexOut, "-o", "output texture");
      app.alias( "--outputTex", "-o" );
      app.addOption( adressAtlas, "-a", "atlas mesh");
      app.alias( "--atlasMesh", "-a" );
      app.addOption( adressIx, "-ix", "mesh x-coordinate texture");
      app.addOption( adressIy, "-iy", "mesh y-coordinate texture");
      app.addOption( adressAx, "-ax", "atlas x-coordinate texture");
      app.addOption( adressAy, "-ay", "atlas y-coordinate texture");
      app.addOption( px, "-px", "x-coord period (none=0)", 0);
  
      app.initialize();
     }
     catch( user_interruption &)
     {
       return EXIT_FAILURE;
     }
     catch( ... )
     {
       throw;
     }


     std::cout << "Reading all mesh and textures" << endl;
     AimsSurfaceTriangle mesh, atlas;
     TimeTexture<float> texIn, texIx, texIy, texAx, texAy, texOut;
     Reader < AimsSurfaceTriangle > rm(adressMesh);
     rm.read( mesh );
     Reader < AimsSurfaceTriangle > ra(adressAtlas);
     ra.read( atlas );
     Reader < TimeTexture<float> > rt(adressTexIn);
     rt.read( texIn );
     Reader < TimeTexture<float> > rix(adressIx);
     rix.read( texIx );
     Reader < TimeTexture<float> > riy(adressIy);
     riy.read( texIy );
     Reader < TimeTexture<float> > rax(adressAx);
     rax.read( texAx );
     Reader < TimeTexture<float> > ray(adressAy);
     ray.read( texAy );

     cout << "Building interpolator and computing new texture" << endl;
     Mesh2mesh projection(mesh, atlas, texIx, texIy, texAx, texAy, px);
     texOut=projection.sendTextureToTarget(texIn);

    std::cout << "Writing interpolated texture" << endl;
    // ECRITURE DE LA TEXTURE

    Writer< TimeTexture<float> > wt(adressTexOut);
    wt.write( texOut );

    return(0);
}

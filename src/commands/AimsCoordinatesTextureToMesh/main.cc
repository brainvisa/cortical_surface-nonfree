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
     std::string adressTexZ="./";
     std::string adressOut="./";

     AimsApplication     app( argc, argv, "Creates a mesh from 3 coordinates textures (x,y and z)");
     try
     {
      app.addOption( adressMesh, "-i", "input mesh");
      app.alias( "--inputMesh", "-i" );
      app.addOption( adressTexX, "-x", "texture_x");
      app.alias( "--inputx", "-x" );
      app.addOption( adressTexY, "-y", "texture_y");
      app.alias( "--inputy", "-y" );
      app.addOption( adressTexZ, "-z", "texture_z");
      app.alias( "--inputz", "-z" );
      app.addOption( adressOut, "-o", "output mesh");
      app.initialize();
  
      std::cout << "Reading mesh" << endl;
      AimsSurfaceTriangle mesh;
      Reader < AimsSurfaceTriangle > rm(adressMesh);
      rm.read( mesh );
      std::cout << "Reading coordinate textures" << endl;
      TimeTexture<float> texx, texy, texz;
      Reader < TimeTexture<float> > rtx(adressTexX);
      rtx.read( texx );
      Reader < TimeTexture<float> > rty(adressTexY);
      rty.read( texy );
      Reader < TimeTexture<float> > rtz(adressTexZ);
      rtz.read( texz );
  
      uint ns=mesh.vertex().size();
      if ((texx[0].nItem() != ns) || (texy[0].nItem() != ns) ||(texz[0].nItem() != ns))
      {
            cerr << "One of the coordinate textures does not have the same number of items than the mesh" << endl;
            exit(EXIT_FAILURE);
      }
  
      std::cout << "found " << ns << " vertices" << endl;
      float x, y, z;
  
      // Creation des textures.
  
      std::cout << "generating new mesh" << endl;
      std::vector<Point3df> vert=mesh.vertex();
      for (uint i=0; i<ns; i++)
      {
            x=texx[0].item(i);
            y=texy[0].item(i);
            z=texz[0].item(i);
            (mesh.vertex()[i])[0]=x;(mesh.vertex()[i])[1]=y;(mesh.vertex()[i])[2]=z;
      }
      
      mesh.updateNormals();
  
      // �criture des textures
      cout << "Writing new mesh" << endl;
      Writer< AimsSurfaceTriangle > mo(adressOut);
      mo.write( mesh );
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

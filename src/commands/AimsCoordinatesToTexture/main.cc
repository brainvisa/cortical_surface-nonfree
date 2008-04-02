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
  std::string adressTexZ="./";
  int local=0;

  AimsApplication     app( argc, argv, "Create 3 textures of node coordinates");
  try
  {
    app.addOption( adressMesh, "-i", "input mesh");
    app.alias( "--inputMesh", "-i" );
    app.addOption( adressTexX, "-x", "texture_x");
    app.alias( "--outx", "-x" );
    app.addOption( adressTexY, "-y", "texture_y");
    app.alias( "--outy", "-y" );
    app.addOption( adressTexZ, "-z", "texture_z");
    app.alias( "--outz", "-z" );
    app.addOption( local, "-l", "local coordinates transo (yes=1, default=0)", 0);
    app.alias( "--local", "-l");
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

    uint ns=mesh.vertex().size();

    std::cout << "found " << ns << " vertices" << endl;
    TimeTexture<float> texx(1,ns), texy(1,ns), texz(1,ns);
    float x, y, z;

    // if necessary, transforming coordinates in local referential.

    if (local ==1)
    {
          float A=0.0, B=0.0, C=0.0, D=0.0, E=0.0, F=0.0;
          std::vector<Point3df> vert=mesh.vertex();
          float gx=0.0, gy=0.0, gz=0.0;
          for (uint i=0; i<ns; i++)
          {
              Point3df pt=vert[i];
              x=pt[0];y=pt[1];z=pt[2];
              gx+=x; gy+=y; gz+=z;
          }
          gx/=ns; gy/=ns; gz/=ns;
          // cette partie est une simple normalisation des coordonnées 
          // centrée autour du barycentre
          for (uint i=0; i<ns; i++)
          {
              Point3df pt=vert[i];
              x=pt[0];y=pt[1];z=pt[2];
              (mesh.vertex()[i])[0]=x-gx;(vert[i])[1]=y-gy;(vert[i])[2]=z-gz;
          }
          //--------------------------------
          
//           for (uint i=0; i<ns; i++)
//           {
//                Point3df pt=vert[i];
//                x=pt[0];y=pt[1];z=pt[2];
//                A+=(y-gy)*(y-gy) + (z-gz)*(z-gz);
//                B+=(x-gx)*(x-gx) + (z-gz)*(z-gz);
//                C+=(x-gx)*(x-gx) + (y-gy)*(y-gy);
//                D+=(y-gy)*(z-gz);
//                E+=(x-gx)*(z-gz);
//                F+=(x-gx)*(y-gy);
//           }
//           cout << "Matrice d'inertie :" << endl;
//           cout << A << "   " << -F << "   " << -E << endl;
//           cout << -F << "   " << B << "   " << -D << endl;
//           cout << -E << "   " << -D << "   " << C << endl;
//           cout << endl;
// 
//           AimsData< float > inertia(3,3);
//           inertia(0,0)=A; inertia(1,0)=inertia(0,1)=-F; inertia(2,0)=inertia(0,2)=-E;
//           inertia(1,1)=B; inertia(1,2)=inertia(2,1)=-D; inertia(2,2)=C;
//      
//           AimsEigen<float> diago;
//           AimsData<float> eigval;
//           eigval=diago.doit(inertia);
//           int ex=eigval.dimX(), ey=eigval.dimY();
//           int i,j;
//           cout << "Eigenvalues : " << endl;
//           for (j=0; j<ey; j++)
//           {
//                for (i=0; i<ex ; i++)
//                     cout << eigval(i,j) << "   ";
//                cout << endl;
//           }
//           cout << endl;
//           cout << "Eigenvectors" << endl;
//           cout << inertia(0,0) << "   " << inertia(1,0) << "   " << inertia(2,0) << endl;
//           cout << inertia(0,1) << "   " << inertia(1,1) << "   " << inertia(2,1) << endl;
//           cout << inertia(0,2) << "   " << inertia(1,2) << "   " << inertia(2,2) << endl;

    }


    // Creation des textures.

    std::cout << "generating textures" << endl;
    std::vector<Point3df> vert=mesh.vertex();
    for (uint i=0; i<ns; i++)
    {
          Point3df pt=vert[i];
          x=pt[0];y=pt[1];z=pt[2];
          texx[0].item(i)=x;
          texy[0].item(i)=y;
          texz[0].item(i)=z;
    }

    // écriture des textures

    Writer< TimeTexture<float> > wtx(adressTexX);
    Writer< TimeTexture<float> > wty(adressTexY);
    Writer< TimeTexture<float> > wtz(adressTexZ);
    wtx.write( texx ); wty.write( texy ); wtz.write( texz );
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

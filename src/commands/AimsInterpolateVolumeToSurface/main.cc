#include <cstdlib>
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/mesh/mesh_g.h>
#include <string.h>
#include <aims/def/path.h>
#include <cartobase/stream/fileutil.h>

using namespace aims;
using namespace carto;
using namespace std;


int main(int argc, const char **argv)
{
  try
    {
	std::string adress_mesh;
	std::string adress_texOut;
	std::string adress_vol;

	AimsApplication     app( argc, argv, "Projection d'un volume sur une surface par simple interpolation lin�aire des noeuds de la surface dans le volume ");
	app.addOption( adress_mesh, "-im", "input Mesh");
	app.alias( "--inMesh", "-im" );
	app.addOption( adress_vol, "-iv", "input Volume");
	app.alias( "--inVol", "-iv" );
	app.addOption( adress_texOut, "-o", "output Texture");
	app.alias( "--outTex", "-o" );

	app.initialize();

	std::cout << "Reading volume and mesh" << std::endl;

	AimsData<short> func;
	Reader<AimsData<short> > fR(adress_vol);
	fR.read(func);
	AimsSurfaceTriangle surface;
	Reader<AimsSurfaceTriangle > sR(adress_mesh);
	sR.read(surface);

	float sx, sy, sz;
	int nx, ny, nz;
	uint nv;

	std::cout << "Reading input" << std::endl;

	sx=func.sizeX(); sy=func.sizeY(); sz=func.sizeZ();
	nx=func.dimX(); ny=func.dimY(); nz=func.dimZ();

	std::vector< Point3df > vert=surface.vertex();
	nv=vert.size();

	TimeTexture<float> fTex(1, nv);

	std::cout << "Interpolating" << std::endl;

	for (uint i=0; i< nv; i++)
	{
		float xf=vert[i][0],yf=vert[i][1],zf=vert[i][2];
		int x=int(floor(xf/sx));
		int y=int(floor(yf/sy));
		int z=int(floor(zf/sz));

		float tx=(xf-float(x))/sx;
		float ty=(yf-float(y))/sy;
		float tz=(zf-float(z))/sz;

		float val=float( (1-tx)*(1-ty)*(1-tz)* func(x,y,z)
				+  tx*(1-ty)*(1-tz)    * func(x+1,y,z)
				+  (1-tx)*ty*(1-tz)    * func(x,y+1,z)
				+  tx*ty*(1-tz)        * func(x+1,y+1,z)
				+  (1-tx)*(1-ty)*tz    * func(x,y,z+1)
				+  tx*(1-ty)*tz        * func(x+1,y,z+1)
				+  (1-tx)*ty*tz        * func(x,y+1,z+1)
				+  tx*ty*tz            * func(x+1,y+1,z+1)
				);
		fTex[0].item(i)=val;
	}

	std::cout << "Writing texture" << std::endl;
	Writer<TimeTexture<float> > texW(adress_texOut);
	texW.write(fTex);
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


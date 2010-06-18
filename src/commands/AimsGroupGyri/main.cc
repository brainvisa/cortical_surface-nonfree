/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 *
 *  checking and fixing duplicate coordinates in the 2D surface referential
 */

#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <aims/utility/utility_g.h>
#include <aims/mesh/mesh_g.h>
#include <aims/io/io_g.h>
#include <aims/mesh/surfaceOperation.h>

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileMesh, fileModel, gyriIn, gyriOut;

  AimsApplication    app( argc, argv, "Group some of the original cortical gyri to get a new gyri model." );
  try
  {
     app.addOption( fileMesh, "-m", "input mesh" );
     app.alias( "--mesh", "-m" );
     app.addOption( gyriIn, "-i", "input gyri texture" );
     app.alias( "--in", "-i" );
     app.addOption( gyriOut, "-o", "output gyri texture" );
     app.alias( "--out", "-o" );
     app.addOption( fileModel, "-g", "grouping model file" );
     app.alias( "--grouping", "-g" );
     app.initialize();

     cout << "reading triangulation   : " << flush;
     AimsSurfaceTriangle surface;
     Reader<AimsSurfaceTriangle> triR( fileMesh );
     triR >> surface;
     cout << "done" << endl;
     
     cout << "reading gyri   : " << flush;
     Reader<TimeTexture<float> > texGyriR( gyriIn);
     TimeTexture<float> texGyri;
     texGyriR >> texGyri;
     cout << "done " << endl;

     uint size=surface.vertex().size();
     uint i;

     TimeTexture<float> texOut(1, size);

     int label, lmax=0;

     for (i=0; i<size; i++)
     {
    	 label = int(floor(texGyri[0].item(i)));
    	 if (label > lmax) lmax=label;
     }
     
     std::map<int, int> grouping;
     for (int j=0; j<=lmax; j++)
     {
    	 grouping[j]=j;
     }

     TimeTexture<short> texGrouped(1, size);
     
     ifstream fin(fileModel.c_str());

     int one, two;
	 fin >> one >> two;

     while (!fin.eof())
     {
    	 cout << one << " - " << two << endl;
    	 if (two>one) grouping[two]=one;
    	 else grouping[one]=two;
    	 fin >> one >> two;
     }

     for (int j=0; j<=lmax; j++)
    	 if (j != grouping[j])
    		 cout << j << " -> " << grouping[j] << endl;

     for (i=0; i<size; i++)
     {
    	 label = int(floor(texGyri[0].item(i)));
    	 texOut[0].item(i)=(float) grouping[label];
     }

     Writer<TimeTexture<float> > texW(gyriOut);
     texW << texOut;

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


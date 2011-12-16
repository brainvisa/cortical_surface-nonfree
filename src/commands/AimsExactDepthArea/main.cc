/*
 *  Copyright (C) 1997-2005 CEA
 *
 *  This software and supporting documentation were developed by
 *
 *      Laboratoire LSIS, equipe LXAO,
 *      Marseille, France
 *
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 */
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/process.h>
#include <string.h>

#include <fstream>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>

#include <aims/geodesicpath/geodesicPath.h>
#include <cortical_surface/surfacereferential/sulcalLinesGeodesic.h>

using namespace aims;
using namespace carto;
using namespace std;


double mydistance(Point3df v1, Point3df v2)
  {
    double dx = v2[0] - v1[0];
    double dy = v2[1] - v1[1];
    double dz = v2[2] - v1[2];

    return sqrt(dx*dx + dy*dy + dz*dz);
  };

double mylength(vector<Point3df>& path)
{
  double length = 0;
  if(!path.empty())
  {
    for(unsigned i=0; i<path.size()-1; ++i)
    {
      length += mydistance(path[i],path[i+1]);
    }
  }
  return length;
}

int main( int argc, const char** argv )
{
try{
string meshFileIn;
std::string adrTexContour;
std::string adrTexDepth;
int value,method;

AimsApplication    app( argc, argv, "Compute the Exact Geodesic Depth Map" );

app.addOption( meshFileIn, "-im", "mesh" );
app.alias( "--inputMesh", "-im" );

app.addOption( adrTexContour, "-it", "input texture (contour area)" );
app.alias( "--inputTexContour", "-it" );

app.addOption( adrTexDepth, "-ot", "output texture (depth map)" );
app.alias( "--inputTexDepth", "-ot" );

app.addOption( value, "-v", "texture value contour" );
app.alias( "--inputValue", "-v" );

app.addOption( method, "-c", "constraintType:\n\"0\" -> no constraint\n"
    "\"1\" -> constrained sulci\n\"2\" -> constrained gyri\n\"3\" -> Exact geodesic path"
    "\n\"4\" -> Euclidian distance\n\"5\" -> z value",true);
app.alias( "--constraint", "-c" );

app.initialize();

// read triangulation
cout << "reading triangulation   : " << flush;
AimsSurfaceTriangle surface;
Reader<AimsSurfaceTriangle> triR( meshFileIn );
triR >> surface;
cout << "done" << endl;

cout << "reading input texture   : " << flush;

TimeTexture<float> texSource;
Reader < TimeTexture<float> > rits(adrTexContour);
rits.read( texSource );

vector<unsigned> listIndexVertexSource,listIndexVertexTarget;

TimeTexture<short> texSourceShort(1,texSource[0].nItem());
for( uint i = 0; i < texSource[0].nItem(); i++)
{
  if (texSource[0].item(i)>= (float)value - 0.01 && texSource[0].item(i)<= (float)value + 0.01)
  {
  listIndexVertexTarget.push_back(i);
  cout << i << "-" ;
  }

  listIndexVertexSource.push_back(i);
}


GeodesicPath gp(surface,method,0);

vector<Point3df> path3D;

vector<Point3df> & vert = surface.vertex();
Point3df p1,p2;

TimeTexture<float> texOut(1, surface.vertex().size() );

unsigned source_vertex_index;
unsigned target_vertex_index;
//
double distance_temp,min_distance;
float max_dist;
max_dist = 0;

for( uint i = 0; i < texSource[0].nItem(); i++)
{
  source_vertex_index = listIndexVertexSource[i];
  min_distance = 100000.;

  if (method < 3)
    {
    gp.shortestPath_1_N_ind(source_vertex_index,listIndexVertexTarget,&target_vertex_index,&min_distance);
    cout << source_vertex_index << " " << target_vertex_index << " " << min_distance << endl ;
    }
  if (method == 3)
    {
      for( uint j = 0; j < listIndexVertexTarget.size(); j++)
      {
        path3D = gp.shortestPath_1_1_xyz(source_vertex_index,listIndexVertexTarget[j]);
        distance_temp = mylength(path3D);
        min_distance = std::min(min_distance,distance_temp);
      }

      cout << source_vertex_index << " " << min_distance << endl ;

//    for (int i = 0; i < path3D.size(); i++)
//      cout << path3D[i][0] << " " << path3D[i][1] << " "<< path3D[i][2] << "\n";

    }

  if (method == 4)
    {
    p1 = vert[i];

    for( uint j = 0; j < listIndexVertexTarget.size(); j++)
      {
        p2 = vert[listIndexVertexTarget[j]];
        distance_temp = mydistance(p1,p2);
        min_distance = std::min(min_distance,distance_temp);
      }

      cout << i << " " << min_distance << endl ;

//    for (int i = 0; i < path3D.size(); i++)
//      cout << path3D[i][0] << " " << path3D[i][1] << " "<< path3D[i][2] << "\n";


      texOut.item(i) = (float)min_distance;

    }

  if (method == 5)
    {
    texOut.item(i) = -vert[i][1];
    max_dist = max(max_dist,texOut.item(i));
    }

}

if (method == 5)
  {
  for( uint i = 0; i < texOut[0].nItem(); i++)
  {
    texOut.item(i) = (float)(texOut.item(i)/max_dist);
    texOut.item(i) = texOut.item(i)*texOut.item(i)*texOut.item(i);
  }
  }


Writer<TimeTexture<float> > texW(adrTexDepth);
texW << texOut;

return( 0 );
}

catch( user_interruption & )
{
}
catch( exception & e )
{
cerr << e.what() << endl;
}
return 1;
}


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
//#include <aims/io/writer.h>
#include <aims/io/process.h>
//#include <aims/io/finder.h>
#include <string.h>

//#include <aims/io/aimsGraphW.h>
//#include <aims/def/path.h>
//#include <cartobase/stream/fileutil.h>

#include <aims/geodesicpath/geodesic_mesh.h>
#include <aims/geodesicpath/geodesic_algorithm_dijkstra.h>
#include <aims/geodesicpath/geodesic_algorithm_subdivision.h>
#include <aims/geodesicpath/geodesic_algorithm_exact.h>

#include <fstream>

#include <aims/mesh/curv.h>
//#include <aims/data/data.h>
#include <aims/mesh/surfaceOperation.h>
//#include <aims/mesh/geometric.h>
//#include <aims/data/data_g.h>
//#include <aims/io/io_g.h>
//#include <aims/math/math_g.h>
//#include <aims/vector/vector.h>
//#include <aims/mesh/texture.h>
//#include <aims/distancemap/meshdistance_d.h>
//#include <aims/distancemap/distancemap_g.h>
//#include <aims/morphology/morphology_g.h>

#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace geodesic;

int main( int argc, const char** argv )
{
try{
string meshFileIn;
std::string adrTexSource;
std::string adrTexTarget;
int strain = 3;

int value;


AimsApplication app( argc, argv, "Compute the shortest path between two vertex" );

app.addOption( meshFileIn, "-im", "mesh" );
app.alias( "--inputMesh", "-im" );

app.addOption( adrTexSource, "-itS", "input texture source" );
app.alias( "--inputTexSource", "-itS" );

app.addOption( adrTexTarget, "-itT", "input texture target" );
app.alias( "--inputTexTarget", "-itT" );

app.addOption( value, "-v", "input texture value" );
app.alias( "--inputValue", "-v" );

app.initialize();

// read triangulation
cout << "reading triangulation   : " << flush;
AimsSurfaceTriangle surface;
Reader<AimsSurfaceTriangle> triR( meshFileIn );
triR >> surface;
cout << "done" << endl;

cout << "reading input texture   : " << flush;

TimeTexture<float> texSource;
Reader < TimeTexture<float> > rits(adrTexSource);
rits.read( texSource );

TimeTexture<float> texTarget;
Reader < TimeTexture<float> > ritt(adrTexTarget);
ritt.read( texTarget );

cout << "done" << endl;

// compute and copy curvature
/*
TimeTexture<float> texCurv;
cout << "compute texture curvature : ";
texCurv = TimeTexture<float>(1, surface.vertex().size());
texCurv = AimsMeshCurvature(surface[0]);
cout << "done" << endl;

float *f = (float*) malloc (texCurv[0].nItem() * sizeof(float));
for( uint i = 0; i < texCurv[0].nItem(); i++)
{
  f[i] = (float)(texCurv[0].item(i));
}
*/
// copy vertex and faces vector
std::vector<double> pointsSP;
std::vector<unsigned> facesSP;
vector<Point3df> & vert = surface.vertex();
vector<AimsVector<uint, 3> > & tri = surface.polygon();
pointsSP.resize(3*vert.size());
facesSP.resize(3*tri.size());

for (uint j = 0; j < (int) vert.size(); j++)
{
  pointsSP[3*j] = vert[j][0];
  pointsSP[3*j+1] = vert[j][1];
  pointsSP[3*j+2] = vert[j][2];
}
for (uint j = 0; j < (int) tri.size(); j++)
{
  facesSP[3*j] = tri[j][0];
  facesSP[3*j+1] = tri[j][1];
  facesSP[3*j+2] = tri[j][2];
}

// compute adjacence graph
geodesic::Mesh meshSP;
cout << "compute adjacence graph : ";

meshSP.initialize_mesh_data(pointsSP,facesSP, NULL ,0,0);

cout << "done" << endl;

std::vector<unsigned> listIndexVertexSource;

cout << "source points" << endl;
for( uint i = 0; i < texSource[0].nItem(); i++)
{
  if ( (texSource[0].item(i) > 0.3 && texSource[0].item(i)<=1))
    {
    listIndexVertexSource.push_back(i);
    texSource[0].item(i) = 1.0;
    cout << i << " " << texSource[0].item(i) << endl;
    }

  if (texSource[0].item(i)==value)
    {
    texSource[0].item(i) = (float)(texSource[0].item(i))/value;
    cout << i << " " << texSource[0].item(i) << endl;
    listIndexVertexSource.push_back(i);
    }

}

std::vector<unsigned> listIndexVertexTarget;
cout << "target points" << endl;

for( uint i = 0; i < texTarget[0].nItem(); i++)
{
  if ( (texTarget[0].item(i) > 0.3 && texTarget[0].item(i)<=1))
  {
  listIndexVertexTarget.push_back(i);
  texTarget[0].item(i) = 1.0;
  cout << i << " " << texTarget[0].item(i) << endl;
  }

  if (texTarget[0].item(i) == value)
  {
      listIndexVertexTarget.push_back(i);
      texTarget[0].item(i) = (float)(texTarget[0].item(i))/value;
      cout << i << " " << texTarget[0].item(i) << endl;
  }


}

//unsigned source = 5,target = 100;

cout << "compute mean and standard deviation distance :" << endl;

std::vector<geodesic::SurfacePoint> sources;

geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;
dijkstra_algorithm = new geodesic::GeodesicAlgorithmDijkstra(&meshSP);

unsigned source_vertex_index;

//std::vector<geodesic::SurfacePoint> SPath;
//SPath.clear();


double min_distance_temp;
double distance_temp;

std::vector<double> min_distance_source_target(listIndexVertexSource.size(), 0.0);

double max_distance_s_t = 0;
double moy_distance_s_t = 0;
double ecart_type_distance = 0;
int j_min;

cout << "\nsource -> target : " << endl;

for (unsigned i = 0; i < listIndexVertexSource.size(); i++)
{
  min_distance_temp = 100000.0;

  source_vertex_index = listIndexVertexSource[i];

  geodesic::SurfacePoint short_sources(&meshSP.vertices()[source_vertex_index]);

  std::vector<geodesic::SurfacePoint> all_sources(1,short_sources);
  dijkstra_algorithm->propagate(all_sources);

  for (unsigned j = 0; j < listIndexVertexTarget.size(); j++)
  //for(unsigned i=0; i<mesh->vertices().size(); ++i)
  {
    geodesic::SurfacePoint p(&meshSP.vertices()[listIndexVertexTarget[j]]);

    unsigned best_source = dijkstra_algorithm->best_source(p,distance_temp);   //for a given surface point, find closets source and distance to this source

    //on normalise la distance avec la proba
    distance_temp = (float)distance_temp * (
        1 - fabs(texSource[0].item(listIndexVertexSource[i]) - texTarget[0].item(listIndexVertexTarget[j]))
        );

//    if (distance_temp == 0)
//      {
//      distance_temp = fabs(texSource[0].item(listIndexVertexSource[i]) - texTarget[0].item(listIndexVertexTarget[j]));
//      j_min = j;
//      j = listIndexVertexTarget.size();
//      min_distance_temp = distance_temp;
//      }

    if (distance_temp < min_distance_temp )
    {
      min_distance_temp = distance_temp;
      j_min = j;
    }
  }

  std::cout << source_vertex_index << " " << listIndexVertexTarget[j_min] << " " << min_distance_temp << endl;
  min_distance_source_target[i] = min_distance_temp;
  max_distance_s_t = std::max(max_distance_s_t,min_distance_temp);
  moy_distance_s_t += min_distance_temp;

}

cout << "distance max = " << max_distance_s_t << endl;

moy_distance_s_t = (float)(moy_distance_s_t)/listIndexVertexSource.size();
cout << "distance mean = " << moy_distance_s_t << endl;

for (unsigned i = 0; i < listIndexVertexSource.size(); i++)
  ecart_type_distance += pow(min_distance_source_target[i]-moy_distance_s_t,2);

ecart_type_distance =  sqrt((float)(ecart_type_distance)/listIndexVertexSource.size());
cout << "distance standard deviation (ecart type) = " << ecart_type_distance << endl;

cout << "\ntarget -> source : " << endl;

std::vector<double> min_distance_target_source(listIndexVertexTarget.size(), 0.0);
double max_distance_t_s = 0;
double moy_distance_t_s = 0;
ecart_type_distance = 0;

for (unsigned i = 0; i < listIndexVertexTarget.size(); i++)
{
  min_distance_temp = 100000.0;

  source_vertex_index = listIndexVertexTarget[i];

  geodesic::SurfacePoint short_sources(&meshSP.vertices()[source_vertex_index]);

  std::vector<geodesic::SurfacePoint> all_sources(1,short_sources);
  dijkstra_algorithm->propagate(all_sources);

  for (unsigned j = 0; j < listIndexVertexSource.size(); j++)
  {
    geodesic::SurfacePoint p(&meshSP.vertices()[listIndexVertexSource[j]]);

    unsigned best_source = dijkstra_algorithm->best_source(p,distance_temp);   //for a given surface point, find closets source and distance to this source

    distance_temp = (float)distance_temp * (
    1- fabs(texSource[0].item(listIndexVertexSource[j]) - texTarget[0].item(listIndexVertexTarget[i]))
    );

//    if (distance_temp == 0)
//    {
//    distance_temp = fabs(texSource[0].item(listIndexVertexSource[j]) - texTarget[0].item(listIndexVertexTarget[i]));
//    j_min = j;
//    j = listIndexVertexSource.size();
//    min_distance_temp = distance_temp;
//    }

    if (distance_temp < min_distance_temp )
    {
      min_distance_temp = distance_temp;
      j_min = j;
    }
  }

  std::cout << source_vertex_index << " " << listIndexVertexSource[j_min] << " " << min_distance_temp << endl;
  min_distance_target_source[i] = min_distance_temp;
  max_distance_t_s = std::max(max_distance_t_s,min_distance_temp);
  moy_distance_t_s += min_distance_temp;

}
cout << "distance max = " << max_distance_t_s << endl;

moy_distance_t_s = (float)(moy_distance_t_s)/listIndexVertexTarget.size();
cout << "distance mean = " << moy_distance_t_s << endl;

for (unsigned i = 0; i < listIndexVertexTarget.size(); i++)
  ecart_type_distance += pow(min_distance_target_source[i]-moy_distance_t_s,2);

ecart_type_distance =  sqrt((float)(ecart_type_distance)/listIndexVertexTarget.size());
cout << "distance standard deviation (ecart type) = " << ecart_type_distance << endl;

cout << "distance hausdorff source -> target " << max(max_distance_s_t,max_distance_t_s) << endl;


//std::vector<geodesic::SurfacePoint> sources;
//sources.push_back(geodesic::SurfacePoint(&meshSP.vertices()[source]));
//
//std::vector<geodesic::SurfacePoint> targets;
//targets.push_back(geodesic::SurfacePoint(&meshSP.vertices()[target]));
//
//printf("indice source = %d target = %d \n",source, target);
//
//// clear path
//std::vector<geodesic::SurfacePoint> SPath;
//SPath.clear();
//
//
//// dijkstra method
//geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;
//dijkstra_algorithm = new geodesic::GeodesicAlgorithmDijkstra(&meshSP);
//
//std::vector<int> listIndexVertexPathSP;
//listIndexVertexPathSP.clear();
//
//geodesic::SurfacePoint short_sources(&meshSP.vertices()[source]);
//geodesic::SurfacePoint short_targets(&meshSP.vertices()[target]);
//
//dijkstra_algorithm->geodesic(short_sources,short_targets, SPath, listIndexVertexPathSP);
//
////std::vector<int>::iterator ite;
//reverse(listIndexVertexPathSP.begin(),listIndexVertexPathSP.end());
//listIndexVertexPathSP.push_back((int)target);
//
//cout << "shortest path (index vertex) = ";
//for (unsigned i = 0; i < listIndexVertexPathSP.size(); i++)
//cout << listIndexVertexPathSP[i] << " ";
//
//cout << endl;


//geodesic::GeodesicAlgorithmExact *exact_algorithm;
//exact_algorithm = new geodesic::GeodesicAlgorithmExact(&meshSP);
//exact_algorithm->propagate(sources); //cover the whole mesh
//exact_algorithm->print_statistics();
//exact_algorithm->trace_back(targets[0], SPath);
//
//int i;
//
//for (i = 0; i < SPath.size(); ++i)
//{
//  cout << "(" << SPath[i].x() << ',' << SPath[i].y() << ',' << SPath[i].z() << ")\n";
//}

cout << " done\n";

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


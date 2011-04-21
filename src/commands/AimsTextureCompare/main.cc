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

using namespace aims;
using namespace carto;
using namespace std;

int main( int argc, const char** argv )
{
try{
string meshFileIn;
std::string adrTexSource;
std::string adrTexTarget;
int strain = 3;
int value;

AimsApplication    app( argc, argv, "Compute the Hausdorff distance between two curves" );

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
//cout << "reading triangulation   : " << flush;
AimsSurfaceTriangle surface;
Reader<AimsSurfaceTriangle> triR( meshFileIn );
triR >> surface;
//cout << "done" << endl;

//cout << "reading input texture   : " << flush;

TimeTexture<float> texSource;
Reader < TimeTexture<float> > rits(adrTexSource);
rits.read( texSource );

TimeTexture<float> texTarget;
Reader < TimeTexture<float> > ritt(adrTexTarget);
ritt.read( texTarget );

//cout << "done" << endl;

GeodesicPath gp(surface,0,0);

std::vector<unsigned> listIndexVertexSource;

//cout << "source points" << endl;
for( uint i = 0; i < texSource[0].nItem(); i++)
{
  if (texSource[0].item(i)==value)
    {
    texSource[0].item(i) = (float)(texSource[0].item(i))/value;
    //cout << i << " " << texSource[0].item(i) << endl;
    listIndexVertexSource.push_back(i);
    }
}

std::vector<unsigned> listIndexVertexTarget;
//cout << "target points" << endl;

for( uint i = 0; i < texTarget[0].nItem(); i++)
{
  if (texTarget[0].item(i) == value)
  {
      listIndexVertexTarget.push_back(i);
      texTarget[0].item(i) = (float)(texTarget[0].item(i))/value;
      //cout << i << " " << texTarget[0].item(i) << endl;
  }
}

unsigned source_vertex_index;
unsigned target_vertex_index;

double distance_temp;
std::vector<double> min_distance_source_target(listIndexVertexSource.size(), 0.0);
double max_distance_s_t = 0;
double moy_distance_s_t = 0;
double ecart_type_distance = 0;

//cout << "\nsource -> target : " << endl;

for (unsigned i = 0; i < listIndexVertexSource.size(); i++)
{
  source_vertex_index = listIndexVertexSource[i];
  gp.shortestPath_1_N_ind(source_vertex_index,listIndexVertexTarget,&target_vertex_index,&distance_temp);

  min_distance_source_target[i] = distance_temp;
  max_distance_s_t = std::max(max_distance_s_t,distance_temp);
  moy_distance_s_t += distance_temp;
}

//cout << "\ntarget -> source : " << endl;
double max_distance_t_s = 0;
double moy_distance_t_s = 0;
ecart_type_distance = 0;
std::vector<double> min_distance_target_source(listIndexVertexTarget.size(), 0.0);

for (unsigned i = 0; i < listIndexVertexTarget.size(); i++)
{
  source_vertex_index = listIndexVertexTarget[i];
  gp.shortestPath_1_N_ind(source_vertex_index,listIndexVertexSource,&target_vertex_index,&distance_temp);

  min_distance_target_source[i] = distance_temp;
  max_distance_t_s = std::max(max_distance_t_s,distance_temp);
  moy_distance_t_s += distance_temp;
}

moy_distance_s_t = (float)(moy_distance_s_t)/listIndexVertexSource.size();

for (unsigned i = 0; i < listIndexVertexSource.size(); i++)
  ecart_type_distance += pow(min_distance_source_target[i]-moy_distance_s_t,2);

ecart_type_distance =  sqrt((float)(ecart_type_distance)/listIndexVertexSource.size());

moy_distance_t_s = (float)(moy_distance_t_s)/listIndexVertexTarget.size();
for (unsigned i = 0; i < listIndexVertexTarget.size(); i++)
  ecart_type_distance += pow(min_distance_target_source[i]-moy_distance_t_s,2);

ecart_type_distance =  sqrt((float)(ecart_type_distance)/listIndexVertexTarget.size());

//cout << "distance hausdorff " << endl;
cout << max(max_distance_s_t,max_distance_t_s) << endl;

//cout << " done\n";

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


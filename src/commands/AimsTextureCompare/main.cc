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

void compareCurves (map<int,vector<int> > &mapCurv1,map<int,vector<int> > &mapCurv2, AimsSurfaceTriangle &surface)
{
  vector<int>::iterator it1;
  vector<int>::iterator it2;

  map<int, vector<int> >::iterator mit1(mapCurv1.begin()),mend1(mapCurv1.end());
  map<int, vector<int> >::iterator mit2(mapCurv2.begin()),mend2(mapCurv2.end());

  GeodesicPath gp(surface,0,0);

  unsigned source_vertex_index;
  unsigned target_vertex_index;
  //
  //double distance_temp;
  //std::vector<double> min_distance_source_target(listIndexVertexSource.size(), 0.0);
  //double max_distance_s_t = 0;

  double somme_min_distance = 0.0;
  double nb_points = 0.0;

  //double ecart_type_distance = 0;
  //
  //cout << "\nsource -> target : " << endl;
  //
  //for (unsigned i = 0; i < listIndexVertexSource.size(); i++)
  //{
  //  source_vertex_index = listIndexVertexSource[i];
  //  gp.shortestPath_1_N_ind(source_vertex_index,listIndexVertexTarget,&target_vertex_index,&distance_temp);

  double length;
  double min_distance = 100000;
  int min_target;
  double max_source_target = 0.0;

  for (mit1 = mapCurv1.begin(); mit1 != mend1; mit1++)
  {

    vector <int> &c1 = mit1->second;
    //cout << "c1 = " << c1.size() << endl;

    nb_points+=c1.size();

    for (it1=c1.begin(); it1!=c1.end(); it1++)
    {
      min_distance = 100000;

      source_vertex_index = (*it1);

      for (mit2 = mapCurv2.begin(); mit2 != mend2; ++mit2)
      {
        vector <int> &c2 = mit2->second;

        vector<unsigned> listtemp(c2.begin(), c2.end());
//        for (it2=c2.begin(); it2!=c2.end(); it2++)
//          cout << "\n\n ******** " << *it2 << " ";

        gp.shortestPath_1_N_ind(source_vertex_index, listtemp, &target_vertex_index, &length);

        if (length < min_distance)
        {
          min_distance = length;
          min_target = target_vertex_index;
        }

       }

      //cout << source_vertex_index << "\t" << min_target << "\t" <<  min_distance << endl;
      somme_min_distance += min_distance;

      if (min_distance > max_source_target)
        max_source_target = min_distance;
    }


      //it1 = (mit1->second).begin();
    //it2 = (mit2->second).begin();


//    listIndexTemp.clear();
//    cout << "\nnb points cloud = " << (mit->second).size() << endl;
//
//    //pour chaque point du nuage, on cherche une extremit??
//    for (it=(mit->second).begin(); it!=(mit->second).end(); it++)
//    {
//      //cout << *it << " --> ";
//
//      set<uint> nei = _neigh[*it];
//      set<uint>::iterator neiIt = nei.begin();
//      nbNeigh = 0;
//      //on parcourt tous les voisins du sommet et on cherche les extremit??s
//      for (; neiIt != nei.end(); neiIt++)
//      {
//        itf=(mit->second).find(*neiIt);
//        if (itf != (mit->second).end())
//          nbNeigh++;
//      }
//
//      if (nbNeigh < 2)
//      {
//      listIndexTemp.push_back(*it);
//      listIndexTempSelect.insert(*it);
//      cout << *it << " ";
//      break; // on s'arr??te
//      }
//    }
//
//    int nbpoints = 0;
//    while (nbpoints < (mit->second).size()-1)
//    {
//      set<uint> nei = _neigh[listIndexTemp[nbpoints]];
//      set<uint>::iterator neiIt = nei.begin();
//      for (; neiIt != nei.end(); ++neiIt)
//      {
//        itf=(mit->second).find(*neiIt);
//        itf2=listIndexTempSelect.find(*neiIt);
//
//        if ( itf != (mit->second).end() && itf2 == listIndexTempSelect.end() )
//        {
//          listIndexTemp.push_back(*neiIt);
//          listIndexTempSelect.insert(*neiIt);
//          nbpoints++;
//          cout << *neiIt << " ";
//          //break;
//        }
//      }
//    }
//
//    cout << "\nnb points curv = " << listIndexTemp.size() << endl;
//    mapCurv.insert (pair<int, vector<int> >(nbCurv++, listIndexTemp));
  }

  cout << "Hav " << (float)somme_min_distance/nb_points  << " Hwor " << max_source_target << " ";
}

void cloud2curv (map<int,set<int> > &mapCloud,map<int,vector<int> > &mapCurv, AimsSurfaceTriangle &surface)
{
  int nbNeigh,value,nbCurv = 0;
  set<int>::iterator it;
  set<int>::iterator itf,itf2;
  vector<int> listIndexTemp;
  set<int> listIndexTempSelect;
  map<int, set<int> >::iterator mit(mapCloud.begin()),mend(mapCloud.end());
  std::vector<std::set<uint> > _neigh;

  _neigh = SurfaceManip::surfaceNeighbours( surface );

  //On parcourt les tous nuages de points
  for (; mit != mend; mit++)
  {
    listIndexTempSelect.clear();
    listIndexTemp.clear();

    //pour chaque point du nuage, on cherche une extremit??
    for (it=(mit->second).begin(); it!=(mit->second).end(); it++)
    {
      set<uint> nei = _neigh[*it];
      set<uint>::iterator neiIt = nei.begin();
      nbNeigh = 0;
      //on parcourt tous les voisins du sommet et on cherche les extremit??s
      for (; neiIt != nei.end(); neiIt++)
      {
        itf=(mit->second).find(*neiIt);
        if (itf != (mit->second).end())
          nbNeigh++;
      }

      if (nbNeigh < 2)
      {
      listIndexTemp.push_back(*it);
      listIndexTempSelect.insert(*it);
      //cout << *it << " " << endl;
      break; // on s'arr??te lorsque trouve une extremit??
      }
    }

    int nb_voisin,nbpoints = 0;
    while (nbpoints < (mit->second).size()-1)
    {
      set<uint> nei = _neigh[listIndexTemp[nbpoints]];
      set<uint>::iterator neiIt;
      nb_voisin = 0;

      for (neiIt = nei.begin(); neiIt != nei.end(); neiIt++)
      {
        //cout << *neiIt << "-";
        itf=(mit->second).find(*neiIt);
        itf2=listIndexTempSelect.find(*neiIt);

        if ( (itf != (mit->second).end()) && (itf2 == listIndexTempSelect.end()) )
        {
          listIndexTemp.push_back(*neiIt);
          listIndexTempSelect.insert(*neiIt);
          nbpoints++;
          nb_voisin++;
        }
      }

      if (nb_voisin > 1)
      {
      cout << "erreur : la courbe a plusieurs voisins au sommet " << listIndexTemp[nbpoints] << endl;
      return;
      }
    }

    mapCurv.insert (pair<int, vector<int> >(nbCurv++, listIndexTemp));
  }
}


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

vector<unsigned> listIndexVertexSource,listIndexVertexTarget;

TimeTexture<float> texSource;
Reader < TimeTexture<float> > rits(adrTexSource);
rits.read( texSource );

TimeTexture<float> texTarget;
Reader < TimeTexture<float> > ritt(adrTexTarget);
ritt.read( texTarget );

TimeTexture<short> texSourceShort(1,texSource[0].nItem());
for( uint i = 0; i < texSource[0].nItem(); i++)
{
  if (texSource[0].item(i)>= (float)value - 0.01 && texSource[0].item(i)<= (float)value + 0.01)
    {
    texSourceShort[0].item(i) = 1000;
    listIndexVertexSource.push_back(i);
    }
  else
    texSourceShort[0].item(i) = 0;
}

TimeTexture<short> texTargetShort(1,texTarget[0].nItem());
for( uint i = 0; i < texTarget[0].nItem(); i++)
{
  if (texTarget[0].item(i)>= (float)value - 0.01 && texTarget[0].item(i)<= (float)value + 0.01)
    {
    texTargetShort[0].item(i) = 1;
    listIndexVertexTarget.push_back(i);
    }
  else
    texTargetShort[0].item(i) = 0;
}

//string adr="";
//vector<float> p;
//SulcalLinesGeodesic slg(meshFileIn,adr,adr,adr,adr,adr,adr,adr,adr,1, 0, 3, p, adr, 0.0, adr, 0,1);
//
//map<int,set<int> > mapPointSource;
//map<int,vector<int> > mapCurvSource;
//slg.texConnectedComponent(texSourceShort, mapPointSource,1000);
////cout << "component number of source = " << mapPointSource.size() << endl;
//cloud2curv (mapPointSource,mapCurvSource,surface);
//
//map<int,vector<int> > mapCurvTarget;
//map<int,set<int> > mapPointTarget;
//slg.texConnectedComponent(texTargetShort, mapPointTarget,1000);
////cout << "component number of target = " << mapPointTarget.size() << endl;
//cloud2curv (mapPointTarget,mapCurvTarget,surface);
//
////cout << "source --> target\n";
//cout << value << " ";
//compareCurves (mapCurvSource,mapCurvTarget,surface);
//
////cout << "target --> source\n";
//compareCurves (mapCurvTarget,mapCurvSource,surface);
//
//cout << endl;

cout << value ;

GeodesicPath gp(surface,0,0);

unsigned source_vertex_index;
unsigned target_vertex_index;

double distance_temp;
std::vector<double> min_distance_source_target(listIndexVertexSource.size(), 0.0);
double max_distance_s_t = 0;
double moy_distance_s_t = 0;

for (unsigned i = 0; i < listIndexVertexSource.size(); i++)
{
  source_vertex_index = listIndexVertexSource[i];
  gp.shortestPath_1_N_ind(source_vertex_index,listIndexVertexTarget,&target_vertex_index,&distance_temp);
  min_distance_source_target[i] = distance_temp;
  max_distance_s_t = std::max(max_distance_s_t,distance_temp);
  moy_distance_s_t += distance_temp;
}

moy_distance_s_t = (float)(moy_distance_s_t)/listIndexVertexSource.size();

cout << " Hav(S,T) " << moy_distance_s_t << " Hwor(S,T) " << max_distance_s_t ;

std::vector<double> min_distance_target_source(listIndexVertexTarget.size(), 0.0);
double max_distance_t_s = 0;
double moy_distance_t_s = 0;

for (unsigned i = 0; i < listIndexVertexTarget.size(); i++)
{
  source_vertex_index = listIndexVertexTarget[i];
  gp.shortestPath_1_N_ind(source_vertex_index,listIndexVertexSource,&target_vertex_index,&distance_temp);
  min_distance_target_source[i] = distance_temp;
  max_distance_t_s = std::max(max_distance_t_s,distance_temp);
  moy_distance_t_s += distance_temp;
}

moy_distance_t_s = (float)(moy_distance_t_s)/listIndexVertexTarget.size();

cout << " Hav(T,S) " << moy_distance_t_s << " Hwor(T,S) " << max_distance_t_s << endl;


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

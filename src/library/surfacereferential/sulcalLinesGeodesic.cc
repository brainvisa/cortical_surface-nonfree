/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *  	CEA/DSV/SHFJ
 *  	4 place du General Leclerc
 *  	91401 Orsay cedex
 *  	France
 */

#include <cortical_surface/surfacereferential/sulcalLinesGeodesic.h>

#include <aims/mesh/surfaceOperation.h>
#include <aims/scalespace/meshDiffuse.h>
#include <aims/distancemap/meshmorphomat.h>

using namespace std;

SulcalLinesGeodesic::SulcalLinesGeodesic( std::string & adrMesh,std::string & adrGeodesicDepth, std::string & adrRootsLon,
    std::string & adrRootsLat, std::string & adrLines, std::string & adrBassins, std::string & adrLonGeodesicOut, std::string & adrLatGeodesicOut,int strain ) :
    _adrMesh(adrMesh),_adrGeodesicDepth(adrGeodesicDepth),_adrRootsLon(adrRootsLon),_adrRootsLat(adrRootsLat),_adrLines(adrLines),_adrBassins(adrBassins),
    _adrLonGeodesicOut(adrLonGeodesicOut),_adrLatGeodesicOut(adrLatGeodesicOut), _strain(strain)
{
  std::cout << "Read Mesh : ";
  Reader < AimsSurfaceTriangle > r(adrMesh);
  r.read( _mesh );
  std::cout << "done" << std::endl;

//  std::cout << "Read Geodesic Depth Texture : ";
//  Reader < TimeTexture<float> > rtDG(adrGeodesicDepth);
//  rtDG.read( _geoDepth );

  std::cout << "Read Roots Texture : ";
  Reader < TimeTexture<short> > rtLon(adrRootsLon);
  rtLon.read( _rootsLon );

  Reader < TimeTexture<short> > rtLat(adrRootsLat);
  rtLat.read( _rootsLat );

  std::cout << "done" << std::endl;

  std::cout << "compute neighbours : ";
  _neigh = SurfaceManip::surfaceNeighbours( _mesh );
  std::cout << "done " << std::endl;

  cout << "compute and save texture curvature : ";
  _texCurv = TimeTexture<float>(1, _mesh.vertex().size());
  _texCurv = AimsMeshCurvature(_mesh[0]);

  size_t found;
  found = _adrMesh.find_last_of(".");
  string adrCurv = _adrMesh.substr(0,found) + "_CurvMean.tex";

  Writer<TimeTexture<float> > texW(adrCurv);
  texW << _texCurv;

  cout << "done" << endl;

  computeGraphDijkstra(_mesh, 1 , strain);

}

SulcalLinesGeodesic::~SulcalLinesGeodesic()
{
}

void SulcalLinesGeodesic::run()
{
  //test path
  //computeShortestPathSulci(1000,3000);

  bassinsDetect();
}


void SulcalLinesGeodesic::floodFillIter(int indexVertex, float newTextureValue,
    float oldTextureValue)
{
  //cout << indexVertex << " " << newTextureValue <<  " : ";

  queue<int> stack;
  stack.push(indexVertex);

  int indexCurr;

  _listIndexVertexFill.clear();

  while (!stack.empty())
  {
    indexCurr = stack.front();

    _listIndexVertexFill.insert(indexCurr);
    _texBassins[0].item(indexCurr) = newTextureValue;

    stack.pop();

    set<uint> voisins = _neigh[indexCurr];
    set<uint>::iterator voisIt = voisins.begin();

    //on parcourt tous les voisins du sommet
    for (; voisIt != voisins.end(); voisIt++)
    {
      indexCurr = *voisIt;

      set<int>::const_iterator itef;
      itef = _listIndexVertexFill.find(indexCurr);

      if (itef != _listIndexVertexFill.end())
        continue;

      if ( _texBassins[0].item(indexCurr) == oldTextureValue)
      {
        _listIndexVertexFill.insert(indexCurr);
        stack.push(indexCurr);
        //cout << indexCurr << endl;
        _texBassins[0].item(indexCurr) = newTextureValue;
      }

    }
  }

  _mapBassins.insert (pair<int,set<int> >(newTextureValue, _listIndexVertexFill));
}

void SulcalLinesGeodesic::bassinsDetect()
{
  //writing path in the output texture
  _texBassins = TimeTexture<short>(1, _mesh.vertex().size());

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    if (_texCurv[0].item(i) < 0.0)
      _texBassins[0].item(i) = -1;
    else
      _texBassins[0].item(i) = 0;
  }

  TimeTexture<short> texBassinsErode(1, _mesh.vertex().size() );
  TimeTexture<short> texBassinsDilation(1, _mesh.vertex().size() );

//  texBassinsErode[0]=MeshErosion<short>( _mesh[0], texH0Erode[0], short(0), -1, 1, false);
//  texBassinsDilation[0]=MeshDilation<short>( _mesh[0], texBassinsErode[0], short(0), -1, 4, true);
//  texBassinsErode[0]=MeshErosion<short>( _mesh[0], texBassinsDilation[0], short(0), -1, 3, false);

//  TimeTexture<short> texH0Squel(1, _mesh.vertex().size() );
//
//  vector< list<unsigned> > neighbourso( _mesh.vertex().size());
//  neighbourso = AimsMeshOrderNode(_mesh[0]);
//  texH0Squel[0]=MeshSkeletization<short> ( _mesh[0], texH0Erode[0], short(1), short(0), neighbourso );
//
//  Writer<TimeTexture<short> > texW(_adrLines);
//  texW << texH0Squel;

  cout << "compute label of bassins : ";
  int j = 1;
  _mapBassins.clear();

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    if (_texBassins[0].item(i) == -1)
    {
    floodFillIter(i,j++,-1);
    }
    // pour intégrer les bassins de courbure H>0 (gyri)
//      if (_texBassins[0].item(i) == 0)
//      {
//      floodFillIter(i,j++,0);
//      }

    //on enlève les lat et lon projetés qui sont en dehors des bassins

    if (_texBassins[0].item(i) == 0)
    {
      _rootsLon[0].item(i) = 0;
      _rootsLat[0].item(i) = 0;
    }
  }

  cout << _mapBassins.size() << " bassins found" << endl;

  Writer<TimeTexture<short> > texWB(_adrBassins);
  texWB << _texBassins;

  cout << "sort constraints lat/lon by bassins : ";

  int value;
  int nb_voisins;

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    value = _rootsLon[0].item(i);
    if (value != 0)
      _listIndexLon.insert(i);

    value = _rootsLat[0].item(i);
    if (value != 0)
      _listIndexLat.insert(i);
  }

  map<int, set<int> >::const_iterator mit(_mapBassins.begin()),mend(_mapBassins.end());
  set<int>::iterator listVertexBassin;

  int lat,lon;
  set<int>::iterator itlon;
  set<int>::iterator itlat;
  set<int>::iterator ittemp;
  set<int>::iterator it;

  set<int> _listIndexTemp;

  int nbBassinsConstaintLat = 0;
  int nbBassinsConstaintLon = 0;

  //pour chaque bassin i
  for (; mit != mend; ++mit)
  {
    _listIndexTemp.clear();

    //cout << (int)mit->first << ": ";
    listVertexBassin = (mit->second).begin();

    //on parcourt la liste des contraintes lat
    if (!_listIndexLat.empty())
    for (itlat=_listIndexLat.begin(); itlat!=_listIndexLat.end(); itlat++)
    {
      it=(mit->second).find(*itlat);

      value = _rootsLat[0].item(*itlat);
      //si une lat est dans le bassin alors je l'enlève de la liste
      //et je l'ajoute à la liste temporaire des points du bassin i
      if (it != (mit->second).end())
      {
        if (value != 0)
        {
          set<uint> voisins = _neigh[*itlat];
          set<uint>::iterator voisIt = voisins.begin();
          nb_voisins = 0;


          //on parcourt tous les voisins du sommet
          for (; voisIt != voisins.end(); voisIt++)
          {
            if (_rootsLat[0].item(*voisIt) == value)
              nb_voisins++;
          }

          if (nb_voisins < 2)
            _listIndexTemp.insert(*it);
        }

        //_listIndexTemp.insert(*it);

        //        cout << *it << "," << _rootsLat.item(*it) << " ";
        // suprrimer l'élément copier pour accélérer le process: ça marche pas bizarre...?
        //_listIndexLat.erase(itlat);
      }
    }

    //pair<int,set<int> > temp;

    // on associe la liste temporaire à la map du bassin
    if (!_listIndexTemp.empty())
    {
      //temp = make_pair (_rootsLat.item(*it),_listIndexTemp);
      //cout << _rootsLat.item(*it) << " ";
      _mapConstraintLat.insert (pair<int, set<int> >(nbBassinsConstaintLat++, _listIndexTemp));
    }

    _listIndexTemp.clear();

    for (itlon=_listIndexLon.begin(); itlon!=_listIndexLon.end(); itlon++)
    {
      it=(mit->second).find(*itlon);

      value = _rootsLon[0].item(*itlon);
      //si une lat est dans le bassin alors je l'enlève de la liste
      //et je l'ajoute à la liste temporaire des points du bassin i
      if (it != (mit->second).end())
      {
        if (value != 0)
        {
          set<uint> voisins = _neigh[*itlon];
          set<uint>::iterator voisIt = voisins.begin();
          nb_voisins = 0;

          //on parcourt tous les voisins du sommet
          for (; voisIt != voisins.end(); voisIt++)
          {
            if (_rootsLon[0].item(*voisIt) == value)
              nb_voisins++;
          }

          if (nb_voisins < 2)
            _listIndexTemp.insert(*it);
        }
        // ça marche pas ?
        //_listIndexLat.erase(_listIndexLat.find(*it));
        //cout << *it << "#" << _rootsLon.item(*it) << " ";
      }
    }
    // on associe la liste temporaire à la map du bassin
    if (!_listIndexTemp.empty())
    {
      //temp = make_pair (_rootsLon.item(*it),_listIndexTemp);
      _mapConstraintLon.insert (pair<int, set<int> >(nbBassinsConstaintLon++, _listIndexTemp));
    }
  }

  cout << "done" << endl;

  //cout << _listIndexLat.size() << " points lat /" << _listIndexLon.size() << " points lon " << endl;
  cout << "nb bassin Lat= " << nbBassinsConstaintLat<< endl;

  map<int, set<int> >::const_iterator mclatit(_mapConstraintLat.begin()),mclatend(_mapConstraintLat.end());
  map<int, set<int> >::const_iterator mclonit(_mapConstraintLon.begin()),mclonend(_mapConstraintLon.end());

  set<int>::iterator itp1;

  vector<int> _vectorIndexTemp;

  //textures contenant les contraintes lat et lon
  TimeTexture<float> texOutLat(1, _mesh.vertex().size() );
  TimeTexture<float> texOutLon(1, _mesh.vertex().size() );

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    texOutLat[0].item(i) = 0.0;
    texOutLon[0].item(i) = 0.0;
  }

  //pour chaque bassin contenant des contraintes de latitude

  for (; mclatit != mclatend; ++mclatit)
  {
    _vectorIndexTemp.clear();

    cout << "bassin " << (int)mclatit->first << " : ";
    itp1 = (mclatit->second).begin();

    //on parcourt la liste des contraintes lat
    for (; itp1!=(mclatit->second).end(); itp1++)
    {
      //cout << *itp1 << "#" << _rootsLat.item(*itp1) << " ";
      _vectorIndexTemp.push_back (*itp1);
    }

    //cout << _vectorIndexTemp.size() << " - " ;
    int v1,v2;

    vector<int>::iterator itv1;

    int i;
    vector<int> _vectorIndexTempConstraint;
    if (_vectorIndexTemp.size() > 1)
    {
      while (!_vectorIndexTemp.empty())
      {
        itv1 = _vectorIndexTemp.begin();
        v1 = _rootsLat.item(*itv1);

        _vectorIndexTempConstraint.clear();

        for (i=0; i<_vectorIndexTemp.size(); i++)
        {
          //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
          v2 = _rootsLat.item(_vectorIndexTemp[i]);
          if (v1 == v2)
          {
            _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
            //cout << _vectorIndexTemp[i] << " ";
            _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
            i--;
          }
        }

        int source,target;

        if (_vectorIndexTempConstraint.size()>=2)
        {
          vector<int> listIndexVertexPathSP;

          listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);

          cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
          //listIndexVertexPathSP = computeShortestPathSulci(source,target);

          for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
            {
            //cout << listIndexVertexPathSP[t] << " ";
            texOutLat[0].item(listIndexVertexPathSP[t]) = (float) v1;
            }
        }
      }
      //cout << endl;
    }
  }

  Writer<TimeTexture<float> > texWLat(_adrLatGeodesicOut);
  texWLat << texOutLat;

  //pour chaque bassin contenant des contraintes de longitude

  cout << "nb bassin Lon= " << nbBassinsConstaintLon<< endl;

  for (; mclonit != mclonend; ++mclonit)
  {
    _vectorIndexTemp.clear();

    cout << "bassin " << (int)mclonit->first << " : ";
    itp1 = (mclonit->second).begin();

    //on parcourt la liste des contraintes lon
    for (; itp1!=(mclonit->second).end(); itp1++)
    {
      //cout << *itp1 << "#" << _rootsLon.item(*itp1) << " ";
      _vectorIndexTemp.push_back (*itp1);
    }

    //cout << _vectorIndexTemp.size() << " - " ;
    int v1,v2;

    vector<int>::iterator itv1;

    int i;
    vector<int> _vectorIndexTempConstraint;

    if (_vectorIndexTemp.size() > 1)
    {
      while (!_vectorIndexTemp.empty())
      {
        itv1 = _vectorIndexTemp.begin();
        v1 = _rootsLon.item(*itv1);

        _vectorIndexTempConstraint.clear();

        for (i=0; i<_vectorIndexTemp.size(); i++)
        {
          //cout << _vectorIndexTemp[i] << "#" << _rootsLon.item(_vectorIndexTemp[i]) << " ";
          v2 = _rootsLon.item(_vectorIndexTemp[i]);
          if (v1 == v2)
          {
            _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
            //cout << _vectorIndexTemp[i] << " ";
            _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
            i--;
          }
        }

        int source,target;

        vector<int> listIndexVertexPathSP;

        if (_vectorIndexTempConstraint.size()>=2)
        {
          listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);

          cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
          //listIndexVertexPathSP = computeShortestPathSulci(source,target);

          for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
            {
            //cout << listIndexVertexPathSP[t] << " ";
            texOutLon[0].item(listIndexVertexPathSP[t]) = (float) v1;
            }
        }
      }

      //cout << endl;
    }
  }

  Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
  texWLon << texOutLon;

  //on fusionne les textures lat et lon
  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
     if (texOutLon[0].item(i) != 0)
       texOutLat[0].item(i) = texOutLon[0].item(i);
  }

  Writer<TimeTexture<float> > texWLines(_adrLines);
  texWLines << texOutLat;
}

void SulcalLinesGeodesic::computeGraphDijkstra (AimsSurfaceTriangle surface, int constraintType,int strain)
{
  // compute and copy curvature

//  float *f = (float*) malloc (_texCurv[0].nItem() * sizeof(float));
//  for( uint i = 0; i < _texCurv[0].nItem(); i++)
//  {
//  f[i] = (float)(_texCurv[0].item(i));
//  }

  vector<float> &curv = _texCurv[0].data();

  // copy vertex and faces vector
  vector<Point3df> & vert = surface.vertex();
  vector<AimsVector<uint, 3> > & tri = surface.polygon();
  _pointsSP.resize(3*vert.size());
  _facesSP.resize(3*tri.size());

  for (uint j = 0; j < (int) vert.size(); j++)
  {
    _pointsSP[3*j] = vert[j][0];
    _pointsSP[3*j+1] = vert[j][1];
    _pointsSP[3*j+2] = vert[j][2];
  }
  for (uint j = 0; j < (int) tri.size(); j++)
  {
    _facesSP[3*j] = tri[j][0];
    _facesSP[3*j+1] = tri[j][1];
    _facesSP[3*j+2] = tri[j][2];
  }

  // compute adjacence graph

  cout << "compute adjacence graph : ";

  _meshSPc.initialize_mesh_data(_pointsSP,_facesSP, &curv[0],constraintType,strain);
  //_meshSP.initialize_mesh_data(_pointsSP,_facesSP, NULL ,0,0);
  dijkstra_algorithm = new geodesic::GeodesicAlgorithmDijkstra(&_meshSPc);

  cout << "done" << endl;


}

double SulcalLinesGeodesic::computeShortestPathSulci(unsigned source, unsigned target, vector<geodesic::SurfacePoint> & SPath, vector<int> &listIndexVertexPathSP)
{
  //vector<int> listIndexVertexPathSP;

  // compute shortest path
  //cout << "compute shortest path : ";

  vector<geodesic::SurfacePoint> sources;
  sources.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[source]));

  vector<geodesic::SurfacePoint> targets;
  targets.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[target]));

//  printf("indice source = %d target = %d \n",source, target);

  // clear path
  //vector<geodesic::SurfacePoint> SPath;
  SPath.clear();

  // dijkstra method
  //geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;

  listIndexVertexPathSP.clear();

  geodesic::SurfacePoint short_sources(&_meshSPc.vertices()[source]);
  geodesic::SurfacePoint short_targets(&_meshSPc.vertices()[target]);

  dijkstra_algorithm->geodesic(short_sources,short_targets, SPath, listIndexVertexPathSP);

  //std::vector<int>::iterator ite;
  reverse(listIndexVertexPathSP.begin(),listIndexVertexPathSP.end());
  listIndexVertexPathSP.push_back((int)target);

//  cout << "shortest path (index vertex) = ";
//  for (unsigned i = 0; i < listIndexVertexPathSP.size(); i++)
//    cout << listIndexVertexPathSP[i] << " " ;
//  cout << endl;

  return dijkstra_algorithm->length(SPath);
  //return listIndexVertexPathSP;

}

vector<int> SulcalLinesGeodesic::maxGeodesicDistance(vector<int> points, int constraint, int* s, int *d)
{
  int i,j;

  vector<int> listIndexVertexPathSP;
  vector<int> listIndexVertexPathSP_result;
  //cout << source << " " << dest << endl; ;
  // clear path
  vector<geodesic::SurfacePoint> SPath;
  SPath.clear();

  int index_max_i,index_max_j;
  double dist,dist_max;

  index_max_i = -1;
  index_max_i = -1;

  std::cout << endl;
  //dist_max_j = 0;

  int nb_combinaison_max = (points.size() * (points.size()-1))/2;
  int nb_combinaison = 0;

  dist_max = 0;

  for (j=0; j<points.size(); j++)
  {
    for(i=j+1; i<points.size(); i++)
    {
      nb_combinaison++;
      dist = computeShortestPathSulci( points[j], points[i],SPath,listIndexVertexPathSP);

      //cout << "(" << points[j] << "-" << points[i] << " " << dist << ") " << endl;
      //cout << " (" << i << "-" << j << ") " ;
      std::cout << "\r" << nb_combinaison << "/" << nb_combinaison_max << std::flush;

      if (dist > dist_max)
      {
        dist_max = dist;
        index_max_i = i;
        index_max_j = j;
        listIndexVertexPathSP_result = listIndexVertexPathSP;
      }
    }
  }

  *s = points[index_max_i];
  *d = points[index_max_j];

  return listIndexVertexPathSP_result;

  //  int index_max_i,index_max_j;
//  double dist,dist_max_i,dist_max_j;
//
//  index_max_i = -1;
//  index_max_i = -1;
//
//  std::cout << endl;
//  dist_max_j = 0;
//
//  for (j=0; j<points.size(); j++)
//  {
//    std::vector<geodesic::SurfacePoint> sources;
//    sources.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[points[j]]));
//    dijkstra_algorithm->propagate(sources);
//
//    //std::vector<double> max_distance(sources.size(), 0.0);    //distance to the furthest vertex that is covered by the given source
//    cout << endl << points[j] << endl;
//    dist_max_i = 0;
//
//    for(i=0; i<points.size(); i++)
//    {
//      unsigned best_source;
//      //geodesic::SurfacePoint p(&_meshSPc.vertices()[i]);
//      if (i!=index_max_j)
//      {
//        //geodesic::SurfacePoint p(&_meshSPc.vertices()[points[i]]);
//
//        //remplacer best_source par geodesic + lenght()
//        //best_source = dijkstra_algorithm->best_source(targets[i],dist);
//
//        cout << dist << " ";
//        //<< " (" << i << "-" << j << ")" ;
//        if (dist > dist_max_i)
//        {
//          dist_max_i = dist;
//          index_max_i = i;
//        }
//      }
//      //max_distance[best_source] = std::max(max_distance[best_source], distance);
//    }
//
//    if (dist_max_i > dist_max_j)
//    {
//      dist_max_j = dist_max_i;
//      index_max_j = j;
//    }
//
//    //std::cout << "distance " << dist_max_i << " " << index_max_i << endl;
//
//  }

//  std::cout << "max distance " << dist_max_j << " index s/d " <<  index_max_i << "/" << index_max_j << std::endl;
//
//  *s = points[index_max_i];
//  *d = points[index_max_j];


//    for(int i=0; i<max_distance.size(); ++i)
//    {
//      std::cout << "distance " << i
//          << " is " << max_distance[i]
//          << std::endl;
//    }
//  for(unsigned i=0; i<mesh->vertices().size(); ++i)
//  {
//    geodesic::SurfacePoint p(&(*mesh).vertices()[i]);
//
//    unsigned best_source = dijkstra_algorithm->best_source(p,distance_temp);   //for a given surface point, find closets source and distance to this source
//
//    distance[i] = distance_temp;
//    max_distance = std::max(max_distance, distance_temp);
//    //std::cout << i << " -" << distance[i] << "-";   //print geodesic distance for every vertex
//  }



}

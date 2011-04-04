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

SulcalLinesGeodesic::SulcalLinesGeodesic( std::string & adrMesh,std::string & adrCurv,std::string & adrGeodesicDepth,std::string & adrBassinsDepthNorm,std::string & adrRootsLon,
    std::string & adrRootsLat, std::string & adrLines, std::string & adrBassins, std::string & adrLonGeodesicOut, std::string & adrLatGeodesicOut,int strain ) :
    _adrMesh(adrMesh),_adrCurv(adrCurv),_adrGeodesicDepth(adrGeodesicDepth),_adrBassinsDepthNorm(adrBassinsDepthNorm),_adrRootsLon(adrRootsLon),_adrRootsLat(adrRootsLat),_adrLines(adrLines),_adrBassins(adrBassins),
    _adrLonGeodesicOut(adrLonGeodesicOut),_adrLatGeodesicOut(adrLatGeodesicOut), _strain(strain)
{
  std::cout << "Read Mesh : ";
  Reader < AimsSurfaceTriangle > r(adrMesh);
  r.read( _mesh );
  std::cout << "done" << std::endl;

  if (adrGeodesicDepth!="")
  {
  std::cout << "Read Geodesic Depth Texture : ";
  Reader < TimeTexture<float> > rtDG(adrGeodesicDepth);
  rtDG.read( _geoDepth );
  }

  std::cout << "done" << std::endl;

  std::cout << "Read Roots Texture : ";
  Reader < TimeTexture<short> > rtLon(adrRootsLon);
  rtLon.read( _rootsLon );

  Reader < TimeTexture<short> > rtLat(adrRootsLat);
  rtLat.read( _rootsLat );

  std::cout << "done" << std::endl;

  std::cout << "compute neighbours : ";
  _neigh = SurfaceManip::surfaceNeighbours( _mesh );
  std::cout << "done " << std::endl;

  if (adrCurv=="")
  {
  cout << "compute and save texture curvature : ";
  _texCurv = TimeTexture<float>(1, _mesh.vertex().size());
  _texCurv = AimsMeshCurvature(_mesh[0]);
  size_t found;
  found = _adrMesh.find_last_of(".");
  _adrCurv = _adrMesh.substr(0,found) + "_CurvMean.tex";
  Writer<TimeTexture<float> > texW(adrCurv);
  texW << _texCurv;
  cout << "done" << endl;
  }
  else
  {
  Reader < TimeTexture<float> > rtCurv(adrCurv);
  rtCurv.read( _texCurv );
  }


  if (_adrBassinsDepthNorm!="")
  {
  cout << "read bassins depth Norm\n";
  cout << _adrBassinsDepthNorm << endl;

  Reader < TimeTexture<float> > rtDepthNorm(_adrBassinsDepthNorm);
  rtDepthNorm.read( _texBassinsDepthNorm );
  }

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

  bassinsDetect3();
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

  cout << endl;

  ////////////////////////////////////////////////////
  map<int, set<int> >::const_iterator mit(_mapBassins.begin()),mend(_mapBassins.end());
   set<int>::iterator listVertexBassin;
  mit = _mapBassins.begin();

  cout << _mapBassins.size() << " bassins found" << endl;

  if (_adrBassins!="")
  {
    Writer<TimeTexture<short> > texWB(_adrBassins);
    texWB <<_texBassins;
  }

  cout << "sort constraints lat/lon by bassins : ";

  int value,nb_voisins;

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    value = _rootsLon[0].item(i);
    if (value != 0)
      _listIndexLon.insert(i);

    value = _rootsLat[0].item(i);
    if (value != 0)
      _listIndexLat.insert(i);
  }


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

  if (_adrLatGeodesicOut!="")
  {
    Writer<TimeTexture<float> > texWLat(_adrLatGeodesicOut);
    texWLat << texOutLat;
  }

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

  if (_adrLonGeodesicOut!="")
  {
    Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
    texWLon << texOutLon;
  }

  //on fusionne les textures lat et lon
  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
     if (texOutLon[0].item(i) != 0)
       texOutLat[0].item(i) = texOutLon[0].item(i);
  }

  Writer<TimeTexture<float> > texWLines(_adrLines);
  texWLines << texOutLat;
}

void SulcalLinesGeodesic::bassinsDetect2()
{
  //writing path in the output texture
  _texBassins = TimeTexture<short>(1, _mesh.vertex().size());

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    if (_texCurv[0].item(i) < 0.0)
      _texBassins[0].item(i) = 1;
    else
      _texBassins[0].item(i) = 0;
  }

  //idée 1 : ouverture des bassins (pour séparer les bassins connectés par de petits plis)
//  TimeTexture<short> texBassinsErode(1, _mesh.vertex().size() );
//  texBassinsErode[0]=MeshErosion<short>( _mesh[0], _texBassins[0], short(0), -1, 3 , false);
//  TimeTexture<short> texBassinsDil(1, _mesh.vertex().size() );
//  texBassinsDil[0]=MeshDilation<short>( _mesh[0], texBassinsErode[0], short(0), -1, 2, false);

  //pb avec une simple ouverture on a tendance à couper en deux le sillon central par exemple
  //mieux vaut garder de grands bassins (connectés par des plis de passage), l'étiquetage des sillons contribuera a séparer les lignes ...

  //idée2 : fermeture des bassins
  TimeTexture<short> texBassinsDil(1, _mesh.vertex().size() );
   texBassinsDil[0]=MeshDilation<short>( _mesh[0], _texBassins[0], short(0), -1, 1, true);
  TimeTexture<short> texBassinsErode(1, _mesh.vertex().size() );
  texBassinsErode[0]=MeshErosion<short>( _mesh[0], texBassinsDil[0], short(0), -1, 1 , true);



  //rajouter peut être en paramètre le facteur de dilatation ?

  cout << "compute label of bassins : ";
  int j = 2;
  _mapBassins.clear();

  for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
    if (texBassinsErode[0].item(i) > 0)
      _texBassins[0].item(i) = -1;
    else
      _texBassins[0].item(i) = 0;
    }

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    if (_texBassins[0].item(i) == -1)
    {
    floodFillIter(i,j++,-1);
    }
  }

  cout << "done\n";

  string toto0 = "texBassinsCloseEtiquette.tex";
  Writer<TimeTexture<short> > texW0(toto0);
  texW0 << _texBassins;

  cout << endl;

  //pour chaque bassin i
  _texBassinsDepthNorm = TimeTexture<float>(1, _mesh.vertex().size());
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//    _texBassinsDepthNorm[0].item(i) = 1;

  map<int, set<int> >::const_iterator mit(_mapBassins.begin()),mend(_mapBassins.end());
  set<int>::iterator listVertexBassin;
  int n_b = 0;
  float max_depth,val;

  //on normalise la profondeur dans chaque bassin

  if (_adrBassinsDepthNorm!="")
  {
  Reader < TimeTexture<float> > rtDepthNorm(_adrBassinsDepthNorm);
  rtDepthNorm.read( _texBassinsDepthNorm );
  cout << "read bassins depth Norm\n";
  }
  else
  {
    for (; mit != mend; ++mit)
    {
      max_depth = 0;

      //cout << "bassin " << n_b++ << " " << endl;;
      listVertexBassin = (mit->second).begin();
      for (; listVertexBassin!=(mit->second).end(); listVertexBassin++)
      {
      //cout << *listVertexBassin << " " ;
      val = _geoDepth[0].item(*listVertexBassin);
      if ( val > max_depth)
        max_depth = val;
      }
      //cout << max_depth << endl;

      listVertexBassin = (mit->second).begin();
      for (; listVertexBassin!=(mit->second).end(); listVertexBassin++)
      {
      if (max_depth!=0)
        _texBassinsDepthNorm[0].item(*listVertexBassin) =  (float)(_geoDepth[0].item(*listVertexBassin))/max_depth;
      else
        _texBassinsDepthNorm[0].item(*listVertexBassin) = 0 ;
      }
    }

    string toto1 = "texBassinsDepthNorm.tex";
    Writer<TimeTexture<float> > texW1(toto1);
    texW1 << _texBassinsDepthNorm;
  }

  //on binarise les bassins normalisés
//  TimeTexture<short> texBassinsDepthNormBin(1, _mesh.vertex().size() );
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//    if (_texBassinsDepthNorm[0].item(i) < 0.1)
//      texBassinsDepthNormBin[0].item(i) = 0;
//    else
//    texBassinsDepthNormBin[0].item(i) = 1;
//  }

//  string toto = "/home/arnaud/Bureau/texBassinsDepthNormBin.tex";
//  Writer<TimeTexture<short> > texW2(toto);
//  texW2 << texBassinsDepthNormBin;
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//    if (texBassinsDepthNormBin[0].item(i) == 1)
//      _texBassins[0].item(i) = -1;
//    else
//      _texBassins[0].item(i) = 0;
//  }
//
//  cout << "re-compute label of bassins norm : ";
//  j = 1;
//  _mapBassins.clear();
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//    if (_texBassins[0].item(i) == -1)
//    {
//    floodFillIter(i,j++,-1);
//    }
//  }
//
//  cout << endl;
//  string toto5 = "/home/arnaud/Bureau/texBassinsEtiquette.tex";
//  Writer<TimeTexture<short> > texW5(toto5);
//  texW5 << _texBassins;

  //on dilate les roots

   TimeTexture<short> texProjectionLatDil(1, _mesh.vertex().size() );
   texProjectionLatDil[0]=MeshDilation<short>( _mesh[0], _rootsLat[0], short(0), -1, 6, false);
   TimeTexture<short> texProjectionLonDil(1, _mesh.vertex().size() );
   texProjectionLonDil[0]=MeshDilation<short>( _mesh[0], _rootsLon[0], short(0), -1, 6, false);

   string lats = "texLatRootsDilate.tex";
   Writer<TimeTexture<short> > texWLats(lats);
   texWLats << texProjectionLatDil;

   string lons = "texLonRootsDilate.tex";
   Writer<TimeTexture<short> > texWLons(lons);
   texWLons << texProjectionLonDil;

   // on fait les intersections avec les bassins
   TimeTexture<short> texProjectionLatDilInter(1, _mesh.vertex().size() );
   TimeTexture<short> texProjectionLonDilInter(1, _mesh.vertex().size() );

   mit = _mapBassins.begin();
   mend = _mapBassins.end();

   for (; mit != mend; ++mit)
   {
     listVertexBassin = (mit->second).begin();
     for (; listVertexBassin!=(mit->second).end(); listVertexBassin++)
     {
     //cout << *listVertexBassin << " " ;
     //geoDepth[0].item(*listVertexBassin);
     }

//     listVertexBassin = (mit->second).begin();
//     for (; listVertexBassin!=(mit->second).end(); listVertexBassin++)
//     {
//     if (max_depth!=0)
//       _texBassinsDepthNorm[0].item(*listVertexBassin) =  (float)(_geoDepth[0].item(*listVertexBassin))/max_depth;
//     else
//       _texBassinsDepthNorm[0].item(*listVertexBassin) = 0 ;
//     }
   }

  if (_adrLines!="")
  {
    TimeTexture<short> texH0Squel(1, _mesh.vertex().size() );
    int nb_voisins;

    vector< list<unsigned> > neighbourso( _mesh.vertex().size());
    neighbourso = AimsMeshOrderNode(_mesh[0]);
    texH0Squel[0]=MeshSkeletization<short> ( _mesh[0], texBassinsDil[0], short(1), short(0), neighbourso );

    int value;

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      if (texH0Squel[0].item(i)==1)
        {
        value = 1;
        texH0Squel[0].item(i) = 180;
        }
      else
      value = _texBassins[0].item(i);

      //on marque les points de contour
      if (value > 0)
      {
        set<uint> voisins = _neigh[i];
        set<uint>::iterator voisIt = voisins.begin();
        nb_voisins = 0;

        //on parcourt tous les voisins du sommet
        if (value == 1)
        {
          for (; voisIt != voisins.end(); voisIt++)
          {
            if (texH0Squel[0].item(*voisIt) == 1 || texH0Squel[0].item(*voisIt) == 180 || texH0Squel[0].item(*voisIt) == 360 || texH0Squel[0].item(*voisIt) == 250)
              nb_voisins++;
          }

          if (nb_voisins > 2)
            texH0Squel[0].item(i) = 250;

          if (nb_voisins < 2)
            texH0Squel[0].item(i) = 360;

        }
        else
        {
          //on parcourt tous les voisins du sommet
          voisIt = voisins.begin();
          for (; voisIt != voisins.end(); voisIt++)
          {
            if (_texBassins[0].item(*voisIt) != value &&  _texBassins[0].item(*voisIt)!=20)
            {
              texH0Squel[0].item(i) = 20;
              continue;
            }
          }
        }
      }
    }




    TimeTexture<short> texExtremiteLat(1, _mesh.vertex().size() );
    TimeTexture<short> texExtremiteLon(1, _mesh.vertex().size() );


    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      if (texH0Squel[0].item(i)==20)
      {
        //si il y a une intersection avec la texture lat dilatée alors on attribue le label de la région
        if (texProjectionLatDil[0].item(i) > 0)
        {
        texExtremiteLat[0].item(i) = texProjectionLatDil[0].item(i);
        _listIndexLat.insert(i);
        }
        if (texProjectionLonDil[0].item(i) > 0)
        {
        texExtremiteLon[0].item(i) = texProjectionLonDil[0].item(i);
        _listIndexLon.insert(i);
        }
      }


    }


//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      if (texH0Squel[0].item(i)==250 || texH0Squel[0].item(i)==360 )
//      {
//        //si il y a une intersection avec la texture lat dilatée alors on attribue le label de la région
//        if (texProjectionLatDil[0].item(i) > 0)
//        {
//        texExtremiteLat[0].item(i) = texProjectionLatDil[0].item(i);
//        _listIndexLat.insert(i);
//        }
//        if (texProjectionLonDil[0].item(i) > 0)
//        {
//        texExtremiteLon[0].item(i) = texProjectionLonDil[0].item(i);
//        _listIndexLon.insert(i);
//        }
//      }
//    }

    //  cout << endl;
      string toto2 = "texBassinsContourSquel.tex";
      Writer<TimeTexture<short> > texW2(toto2);
      texW2 << texH0Squel;

      string toto3 = "texBassinsExtermiteLat.tex";
      Writer<TimeTexture<short> > texW3(toto3);
      texW3 << texExtremiteLat;

      string toto4 = "texBassinsExtermiteLon.tex";
      Writer<TimeTexture<short> > texW4(toto4);
      texW4 << texExtremiteLon;
//    Writer<TimeTexture<short> > texW(_adrLines);
//    texW << texH0Squel;

      mit = _mapBassins.begin();
      mend = _mapBassins.end();

      cout << "sort constraints lat/lon by bassins : ";

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

            value = texExtremiteLat[0].item(*itlat);
            //si une lat est dans le bassin alors je l'enlève de la liste
            //et je l'ajoute à la liste temporaire des points du bassin i
            if (it != (mit->second).end())
            {
              if (value != 0)
              {
                  _listIndexTemp.insert(*it);
              }

            }
          }

          // on associe la liste temporaire à la map du bassin
          if (!_listIndexTemp.empty())
          {
            //temp = make_pair (_rootsLat.item(*it),_listIndexTemp);
            //cout << _rootsLat.item(*it) << " ";
            _mapConstraintLat.insert (pair<int, set<int> >(nbBassinsConstaintLat++, _listIndexTemp));
          }

          _listIndexTemp.clear();

          if (!_listIndexLon.empty())
          for (itlon=_listIndexLon.begin(); itlon!=_listIndexLon.end(); itlon++)
          {
            it=(mit->second).find(*itlon);

            value = texExtremiteLon[0].item(*itlon);
            //si une lat est dans le bassin alors je l'enlève de la liste
            //et je l'ajoute à la liste temporaire des points du bassin i
            if (it != (mit->second).end())
            {
              if (value != 0)
              _listIndexTemp.insert(*it);
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


        myHistoLat.open ("../histoLat.txt");

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

        _texProbaPath = TimeTexture<short>(1, _mesh.vertex().size());

        //pour chaque bassin contenant des contraintes de latitude

        for (; mclatit != mclatend; ++mclatit)
        {
          _vectorIndexTemp.clear();



          cout << "bassin " << (int)mclatit->first << " : ";
          itp1 = (mclatit->second).begin();

          //ARN DEBUG
          //if (mclatit->first == 4)
          {
            myHistoLat << "bassin " << (int)mclatit->first << "\n";
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
                //v1 = _rootsLat.item(*itv1);
                v1 = texExtremiteLat.item(*itv1);

                _vectorIndexTempConstraint.clear();

                for (i=0; i<_vectorIndexTemp.size(); i++)
                {
                  //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
                  //v2 = _rootsLat.item(_vectorIndexTemp[i]);
                  v2 = texExtremiteLat.item(_vectorIndexTemp[i]);
                  if (v1 == v2 && (v1==45 || v1 == 5))
                  {
                    _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
                    //cout << _vectorIndexTemp[i] << " ";
                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
                    i--;
                  }
                  else
                    if (v1 == v2)
                    {
                    //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
                    //cout << _vectorIndexTemp[i] << " ";
                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
                    i--;
                    }
                }

                int source,target;

                if (_vectorIndexTempConstraint.size()>=2)
                {
                  vector<int> listIndexVertexPathSP;

                  //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
                  listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);

                  cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;

                  myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
                  saveHistoTemp (source,target);

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
          //ARN DEBUG
        }

        if (_adrLatGeodesicOut!="")
        {
          Writer<TimeTexture<float> > texWLat(_adrLatGeodesicOut);
          texWLat << texOutLat;
        }


        //pour chaque bassin contenant des contraintes de longitude

        cout << "nb bassin Lon= " << nbBassinsConstaintLon<< endl;

        for (; mclonit != mclonend; ++mclonit)
        {
          _vectorIndexTemp.clear();

          cout << "bassin " << (int)mclonit->first << " : ";
          itp1 = (mclonit->second).begin();

          //ARN DEBUG
          if (mclonit->first == 5000)
          {

            myHistoLat << "bassin " << (int)mclonit->first << "\n";
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
                v1 = texExtremiteLon.item(*itv1);

                _vectorIndexTempConstraint.clear();

                for (i=0; i<_vectorIndexTemp.size(); i++)
                {
                  //cout << _vectorIndexTemp[i] << "#" << _rootsLon.item(_vectorIndexTemp[i]) << " ";
                  v2 = texExtremiteLon.item(_vectorIndexTemp[i]);
                  if (v1 == v2 && v1==25)
                  {
                    _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
                    //cout << _vectorIndexTemp[i] << " ";
                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
                    i--;
                  }
                  else
                    if (v1 == v2)
                    {
                    //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
                    //cout << _vectorIndexTemp[i] << " ";
                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
                    i--;
                  }
                }

                int source,target;

                vector<int> listIndexVertexPathSP;

                if (_vectorIndexTempConstraint.size()>=2)
                {
                  //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
                  listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);

                  //listIndexVertexPathSP = computeShortestPathSulci(source,target);

                  //myHistoLat << "value" << v1 << "(" << source << "," << target << ")\n";

                  cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;

                  //saveHistoTemp (source,target);

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
          //ARN DEBUG
      }

    myHistoLat.close();



    string probpath = "texProbaPathNew.tex";
    Writer<TimeTexture<short> > texWProbaPath(probpath);
    //Reader<TimeTexture<short> > texWProbaPath(probpath);
    //texWProbaPath.read( _texProbaPath );
    texWProbaPath << _texProbaPath;

    TimeTexture<float> _texProbaPathNorm(1, _mesh.vertex().size() );

    map<int, set<int> >::const_iterator mit2(_mapBassins.begin()),mend2(_mapBassins.end());
    set<int>::iterator listVertexBassin2;
    int n_b2 = 0;
    float max_depth2,val2;

      //on normalise la profondeur dans chaque bassin

    for (; mit2 != mend2; ++mit2)
    {
      max_depth2 = 0;

      cout << "bassin " << n_b++ << " " << endl;;
      listVertexBassin2 = (mit2->second).begin();
      for (; listVertexBassin2!=(mit2->second).end(); listVertexBassin2++)
      {
      //cout << *listVertexBassin << " " ;
      val2 = _texProbaPath[0].item(*listVertexBassin2);
      if ( val2 > max_depth2)
        max_depth2 = val2;
      }

      cout << max_depth2 << endl;

      listVertexBassin2 = (mit2->second).begin();
      for (; listVertexBassin2!=(mit2->second).end(); listVertexBassin2++)
      {
      if (max_depth2!=0)
        _texProbaPathNorm[0].item(*listVertexBassin2) =  (float)(_texProbaPath[0].item(*listVertexBassin2))/max_depth2;
      else
        _texProbaPathNorm[0].item(*listVertexBassin2) = 0 ;
      }
    }

    string probpathnorm = "texProbaPathNewNorm.tex";
    Writer<TimeTexture<float> > texWProbaPathNorm(probpathnorm);
    texWProbaPathNorm << _texProbaPathNorm;



    if (_adrLonGeodesicOut!="")
    {
      Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
      texWLon << texOutLon;
    }

    //on fusionne les textures lat et lon
    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
       if (texOutLon[0].item(i) != 0)
         texOutLat[0].item(i) = texOutLon[0].item(i);
    }

    Writer<TimeTexture<float> > texWLines(_adrLines);
    texWLines << texOutLat;


  }


}

void SulcalLinesGeodesic::bassinsDetect3()
{

  myHistoLat.open ("./histoProba.txt");

  //writing path in the output texture
  _texBassins = TimeTexture<short>(1, _mesh.vertex().size());

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    if (_texCurv[0].item(i) < 0.0)
      _texBassins[0].item(i) = 1;
    else
      _texBassins[0].item(i) = 0;
  }


  //idée 1 : ouverture des bassins (pour séparer les bassins connectés par de petits plis)
//  TimeTexture<short> texBassinsErode(1, _mesh.vertex().size() );
//  texBassinsErode[0]=MeshErosion<short>( _mesh[0], _texBassins[0], short(0), -1, 3 , false);
//  TimeTexture<short> texBassinsDil(1, _mesh.vertex().size() );
//  texBassinsDil[0]=MeshDilation<short>( _mesh[0], texBassinsErode[0], short(0), -1, 2, false);

  //pb avec une simple ouverture on a tendance à couper en deux le sillon central par exemple
  //mieux vaut garder de grands bassins (connectés par des plis de passage), l'étiquetage des sillons contribuera a séparer les lignes ...

  //idée2 : fermeture des bassins
  TimeTexture<short> texBassinsDil(1, _mesh.vertex().size() );
   texBassinsDil[0]=MeshDilation<short>( _mesh[0], _texBassins[0], short(0), -1, 1, true);
  TimeTexture<short> texBassinsErode(1, _mesh.vertex().size() );
  texBassinsErode[0]=MeshErosion<short>( _mesh[0], texBassinsDil[0], short(0), -1, 1 , true);


  //rajouter peut être en paramètre le facteur de dilatation ?

  cout << "compute label of bassins : ";
  int j = 2;
  _mapBassins.clear();

  for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
    if (texBassinsErode[0].item(i) > 0)
      _texBassins[0].item(i) = -1;
    else
      _texBassins[0].item(i) = 0;
    }

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    if (_texBassins[0].item(i) == -1)
    {
    floodFillIter(i,j++,-1);
    }
  }

  TimeTexture<short> texBassinsSave(1, _mesh.vertex().size() );
  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
  texBassinsSave[0].item(i) = _texBassins[0].item(i);
  }

  cout << "done\n";

  string toto0 = "texBassinsCloseEtiquette.tex";
  Writer<TimeTexture<short> > texW0(toto0);
  texW0 << _texBassins;

  cout << endl;

  //pour chaque bassin i
  _texBassinsDepthNorm = TimeTexture<float>(1, _mesh.vertex().size());
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//    _texBassinsDepthNorm[0].item(i) = 1;

  map<int, set<int> >::const_iterator mit(_mapBassins.begin()),mend(_mapBassins.end());
  set<int>::iterator listVertexBassin;
  int n_b = 0;
  float max_depth,val;

  //on normalise la profondeur dans chaque bassin

  if (_adrBassinsDepthNorm=="")
  {
    for (; mit != mend; ++mit)
    {
      max_depth = 0;

      //cout << "bassin " << n_b++ << " " << endl;;
      listVertexBassin = (mit->second).begin();
      for (; listVertexBassin!=(mit->second).end(); listVertexBassin++)
      {
      //cout << *listVertexBassin << " " ;
      val = _geoDepth[0].item(*listVertexBassin);
      if ( val > max_depth)
        max_depth = val;
      }
      //cout << max_depth << endl;

      listVertexBassin = (mit->second).begin();
      for (; listVertexBassin!=(mit->second).end(); listVertexBassin++)
      {
      if (max_depth!=0)
        _texBassinsDepthNorm[0].item(*listVertexBassin) =  (float)(_geoDepth[0].item(*listVertexBassin))/max_depth;
      else
        _texBassinsDepthNorm[0].item(*listVertexBassin) = 0 ;
      }
    }

    string toto1 = "texBassinsDepthNorm.tex";
    Writer<TimeTexture<float> > texW1(toto1);
    texW1 << _texBassinsDepthNorm;
  }


  //on dilate les roots

   TimeTexture<short> texProjectionLatDil(1, _mesh.vertex().size() );
   texProjectionLatDil[0]=MeshDilation<short>( _mesh[0], _rootsLat[0], short(0), -1, 6, false);
   TimeTexture<short> texProjectionLonDil(1, _mesh.vertex().size() );
   texProjectionLonDil[0]=MeshDilation<short>( _mesh[0], _rootsLon[0], short(0), -1, 6, false);

   string lats = "texLatRootsDilate.tex";
   Writer<TimeTexture<short> > texWLats(lats);
   texWLats << texProjectionLatDil;

   string lons = "texLonRootsDilate.tex";
   Writer<TimeTexture<short> > texWLons(lons);
   texWLons << texProjectionLonDil;


   // on fait les intersections avec les bassins


    TimeTexture<short> texExtremiteLat(1, _mesh.vertex().size() );
    TimeTexture<short> texExtremiteLon(1, _mesh.vertex().size() );

    vector< list<unsigned> > neighbourso( _mesh.vertex().size());
    neighbourso = AimsMeshOrderNode(_mesh[0]);

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      if (_texBassins[0].item(i)!=0)
      {
        //si il y a une intersection avec la texture lat dilatée alors on attribue le label de la région
        if (texProjectionLatDil[0].item(i) > 0)
        {
        texExtremiteLat[0].item(i) = texProjectionLatDil[0].item(i);
        //_listIndexLat.insert(i);
        }
        if (texProjectionLonDil[0].item(i) > 0)
        {
        texExtremiteLon[0].item(i) = texProjectionLonDil[0].item(i);
        //_listIndexLon.insert(i);
        }
      }
    }



    cout << "re compute label of bassins with lat inter proj: ";
    j = 200;
    _mapBassins.clear();

    for (uint i = 0; i < _mesh.vertex().size(); i++)
      {
        if (texExtremiteLat[0].item(i) > 0)
          _texBassins[0].item(i) = texExtremiteLat[0].item(i);
        else
          _texBassins[0].item(i) = 0;
      }

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      if (_texBassins[0].item(i) > 0 && _texBassins[0].item(i) < 200)
      {
      floodFillIter(i,j++,_texBassins[0].item(i));
      }
    }

    _texBassinsLat = TimeTexture<short>(1, _mesh.vertex().size());
    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
    if (_texBassins[0].item(i)!=0)
      _texBassinsLat[0].item(i) = _texBassins[0].item(i) - 200;
    _texBassins[0].item(i) = texBassinsSave[0].item(i);
    }

    cout << "done\n";

    int nb_voisins;
    int value;

    TimeTexture<short> texContourLatDilInter(1, _mesh.vertex().size() );

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      value = texExtremiteLat[0].item(i);

      //on marque les points de contour
      if (value > 0)
      {
        set<uint> voisins = _neigh[i];
        set<uint>::iterator voisIt = voisins.begin();
        nb_voisins = 0;

        //on parcourt tous les voisins du sommet
        for (; voisIt != voisins.end(); voisIt++)
        {
          if (texExtremiteLat[0].item(*voisIt) != value)
            {
            nb_voisins++;
            texContourLatDilInter.item(i) = value;
            continue;
            }
        }
      }
    }

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      if (texContourLatDilInter[0].item(i)!=0)
      {
       _listIndexLat.insert(i);
      }
    }

    string latExt = "texLatExtInter.tex";
    Writer<TimeTexture<short> > texWLatExt(latExt);
    texWLatExt << texExtremiteLat;

    string latExt2 = "texLatContourExtInter.tex";
    Writer<TimeTexture<short> > texWLatExt2(latExt2);
    texWLatExt2 << texContourLatDilInter;



    string toto8 = "texBassinsInterLatCachiaEtiquette.tex";
    Writer<TimeTexture<short> > texW8(toto8);
    texW8 << _texBassinsLat;


    mit = _mapBassins.begin();
    mend = _mapBassins.end();

    cout << "sort constraints lat by bassins : ";

    cout << _mapBassins.size() << "\n";

    int lat;
    set<int>::iterator itlat;
    set<int>::iterator ittemp;
    set<int>::iterator it;

    set<int> _listIndexTemp;

    int nbBassinsConstaintLat = 0;

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

        value = texContourLatDilInter[0].item(*itlat);
        //si une lat est dans le bassin alors je l'enlève de la liste
        //et je l'ajoute à la liste temporaire des points du bassin i
        if (it != (mit->second).end())
        {
          if (value != 0)
          {
              _listIndexTemp.insert(*it);
          }

        }
      }

      // on associe la liste temporaire à la map du bassin
      if (!_listIndexTemp.empty())
      {
        //temp = make_pair (_rootsLat.item(*it),_listIndexTemp);
        //cout << _rootsLat.item(*it) << " ";
        _mapConstraintLat.insert (pair<int, set<int> >(nbBassinsConstaintLat++, _listIndexTemp));
      }
    }


    cout << "nb bassin Lat= " << nbBassinsConstaintLat<< endl;

    map<int, set<int> >::const_iterator mclatit(_mapConstraintLat.begin()),mclatend(_mapConstraintLat.end());

    set<int>::iterator itp1;
    vector<int> _vectorIndexTemp;

    //textures contenant les contraintes lat
    TimeTexture<float> texOutLat(1, _mesh.vertex().size() );

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
     texOutLat[0].item(i) = 0.0;
    }

    _texProbaPath = TimeTexture<short>(1, _mesh.vertex().size());

    //pour chaque bassin contenant des contraintes de latitude

    for (; mclatit != mclatend; ++mclatit)
    {
     _vectorIndexTemp.clear();

     cout << "bassin " << (int)mclatit->first << " : ";
     itp1 = (mclatit->second).begin();

     //ARN DEBUG
     //if (mclatit->first == 1)
     {
       //myHistoLat << "bassin " << (int)mclatit->first << "\n";
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
           //v1 = _rootsLat.item(*itv1);
           v1 = texContourLatDilInter.item(*itv1);

           _vectorIndexTempConstraint.clear();

           for (i=0; i<_vectorIndexTemp.size(); i++)
           {
             //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
             //v2 = _rootsLat.item(_vectorIndexTemp[i]);
             v2 = texContourLatDilInter.item(_vectorIndexTemp[i]);
             if (v1 == v2 && (v1==45))
             {
               _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
               //cout << _vectorIndexTemp[i] << " ";
               _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
               i--;
             }
              else
                if (v1 == v2)
                {
                //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
                //cout << _vectorIndexTemp[i] << " ";
                _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
                i--;
                }
           }

           int source,target;

           cout << _vectorIndexTempConstraint.size() << endl;
           if (_vectorIndexTempConstraint.size()>=30)
           {
             vector<int> listIndexVertexPathSP;


             //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
             listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);

             cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;

            // myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
             //saveHistoTemp (source,target);

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
     //ARN DEBUG
    }

    if (_adrLatGeodesicOut!="")
    {
     Writer<TimeTexture<float> > texWLat(_adrLatGeodesicOut);
     texWLat << texOutLat;
    }

    string probpathLat = "texProbaLat.tex";
    Writer<TimeTexture<short> > texWProbaLat(probpathLat);
    texWProbaLat << _texProbaPath;

    //on normalise la texture des probas des lat
    TimeTexture<float> _texProbaPathNormLat(1, _mesh.vertex().size() );

    map<int, set<int> >::const_iterator mit2(_mapBassins.begin()),mend2(_mapBassins.end());
    set<int>::iterator listVertexBassin2;
    int n_b2 = 0;
    float max_depth2,val2;

      //on normalise la profondeur dans chaque bassin lat

    for (; mit2 != mend2; ++mit2)
    {
      max_depth2 = 0.0;
      myHistoLat << "bassin lat" << n_b2++ << " " << endl;
      cout << "bassin " << n_b2-1 << " " << (mit2->second).size() << endl;;
      listVertexBassin2 = (mit2->second).begin();
      for (; listVertexBassin2!=(mit2->second).end(); listVertexBassin2++)
      {
      //cout << *listVertexBassin << " " ;
      val2 = _texProbaPath[0].item(*listVertexBassin2);
      if ( val2 > max_depth2)
        max_depth2 = val2;
      }

      if ( (mit2->second).size() > 150 )
      {
        cout << max_depth2 << endl;

        listVertexBassin2 = (mit2->second).begin();
        for (; listVertexBassin2!=(mit2->second).end(); listVertexBassin2++)
        {
        if (max_depth2!=0)
          _texProbaPathNormLat[0].item(*listVertexBassin2) =  (float)(_texProbaPath[0].item(*listVertexBassin2))/max_depth2;
        else
          _texProbaPathNormLat[0].item(*listVertexBassin2) = 0 ;

        myHistoLat << _texProbaPathNormLat[0].item(*listVertexBassin2) << "\t" << *listVertexBassin2 << "\n";
        }
      }
    }

    string probpathLatN = "texProbaLat_embc.tex";
    Writer<TimeTexture<float> > texWProbaLatN(probpathLatN);
    texWProbaLatN << _texProbaPathNormLat;

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      _texProbaPath[0].item(i) = 0;
    }

    cout << "re compute label of bassins with lon inter proj: \n";
    j = 200;
    _mapBassins.clear();
    _texBassinsLon = TimeTexture<short>(1, _mesh.vertex().size());

    for (uint i = 0; i < _mesh.vertex().size(); i++)
       _texBassinsLon[0].item(i) = _texBassins[0].item(i) - 200;

    for (uint i = 0; i < _mesh.vertex().size(); i++)
      {
        if (texExtremiteLon[0].item(i) > 0)
          _texBassins[0].item(i) = texExtremiteLon[0].item(i);
        else
          _texBassins[0].item(i) = 0;
      }

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      if (_texBassins[0].item(i) > 0 && _texBassins[0].item(i) < 200)
      {
      floodFillIter(i,j++,_texBassins[0].item(i));
      }
    }

    _texBassinsLon = TimeTexture<short>(1, _mesh.vertex().size());
    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
    if (_texBassins[0].item(i)!=0)
      _texBassinsLon[0].item(i) = _texBassins[0].item(i) - 200;
    _texBassins[0].item(i) = texBassinsSave[0].item(i);
    }

    cout << "done\n";

    string toto9 = "texBassinsInterLonCachiaEtiquette.tex";
    Writer<TimeTexture<short> > texW9(toto9);
    texW9 << _texBassinsLon;

    cout << endl;

    TimeTexture<short> texContourLonDilInter(1, _mesh.vertex().size() );

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      value = texExtremiteLon[0].item(i);

      //on marque les points de contour
      if (value > 0)
      {
        set<uint> voisins = _neigh[i];
        set<uint>::iterator voisIt = voisins.begin();
        nb_voisins = 0;

        //on parcourt tous les voisins du sommet
        for (; voisIt != voisins.end(); voisIt++)
        {
          if (texExtremiteLon[0].item(*voisIt) != value)
            {
            nb_voisins++;
            texContourLonDilInter.item(i) = value;
            continue;
            }
        }
      }
    }

     for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
      if (texContourLonDilInter[0].item(i)!=0)
      {
       _listIndexLon.insert(i);
      }
    }

    string lonExt = "texLonExtInter.tex";
    Writer<TimeTexture<short> > texWLonExt(lonExt);
    texWLonExt << texExtremiteLon;

    string lonExt2 = "texLonContourExtInter.tex";
    Writer<TimeTexture<short> > texWLonExt2(lonExt2);
    texWLonExt2 << texContourLonDilInter;

    cout << "done\n";



    mit = _mapBassins.begin();
    mend = _mapBassins.end();
//
    cout << "sort constraints lon by bassins : ";
//
    int lon;
    set<int>::iterator itlon;

    int nbBassinsConstaintLon = 0;
//
    //pour chaque bassin i
    for (; mit != mend; ++mit)
    {
      _listIndexTemp.clear();

      //cout << (int)mit->first << ": ";
      listVertexBassin = (mit->second).begin();

      //on parcourt la liste des contraintes lat
      if (!_listIndexLon.empty())
      for (itlon=_listIndexLon.begin(); itlon!=_listIndexLon.end(); itlon++)
      {
        it=(mit->second).find(*itlon);

        value = texContourLonDilInter[0].item(*itlon);
        //si une lat est dans le bassin alors je l'enlève de la liste
        //et je l'ajoute à la liste temporaire des points du bassin i
        if (it != (mit->second).end())
        {
          if (value != 0)
          {
              _listIndexTemp.insert(*it);
          }

        }
      }
    // on associe la liste temporaire à la map du bassin
      if (!_listIndexTemp.empty())
      {
        //temp = make_pair (_rootsLon.item(*it),_listIndexTemp);
        _mapConstraintLon.insert (pair<int, set<int> >(nbBassinsConstaintLon++, _listIndexTemp));
      }

    }


    cout << "nb bassin Lon= " << nbBassinsConstaintLon<< endl;

    map<int, set<int> >::const_iterator mclonit(_mapConstraintLon.begin()),mclonend(_mapConstraintLon.end());

    //textures contenant les contraintes lon
    TimeTexture<float> texOutLon(1, _mesh.vertex().size() );

    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
     texOutLon[0].item(i) = 0.0;
    }

    //pour chaque bassin contenant des contraintes de longitude

    for (; mclonit != mclonend; ++mclonit)
    {
     _vectorIndexTemp.clear();

     cout << "bassin " << (int)mclonit->first << " : ";
     itp1 = (mclonit->second).begin();

     //ARN DEBUG
    // if (mclonit->first == 100000)
     {
       //myHistoLat << "bassin " << (int)mclatit->first << "\n";
       //on parcourt la liste des contraintes lat
       for (; itp1!=(mclonit->second).end(); itp1++)
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
           //v1 = _rootsLat.item(*itv1);
           v1 = texContourLonDilInter.item(*itv1);

           _vectorIndexTempConstraint.clear();

           for (i=0; i<_vectorIndexTemp.size(); i++)
           {
             //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
             //v2 = _rootsLat.item(_vectorIndexTemp[i]);
             v2 = texContourLonDilInter.item(_vectorIndexTemp[i]);
             if (v1 == v2 && (v1==25))
             {
               _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
               //cout << _vectorIndexTemp[i] << " ";
               _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
               i--;
             }
            else
              if (v1 == v2)
              {
              //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
              //cout << _vectorIndexTemp[i] << " ";
              _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
              i--;
              }
           }

           int source,target;
           cout << _vectorIndexTempConstraint.size() << endl;
           if (_vectorIndexTempConstraint.size()>=30)
           {
             vector<int> listIndexVertexPathSP;

             //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
             listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);

             cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;

             //myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
             //saveHistoTemp (source,target);

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
     //ARN DEBUG
    }

    if (_adrLonGeodesicOut!="")
    {
     Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
     texWLon << texOutLon;
    }

    string probpathLon = "texProbaLon.tex";
    Writer<TimeTexture<short> > texWProbaLon(probpathLon);
    texWProbaLon << _texProbaPath;


    //on normalise la texture des probas des lon
    TimeTexture<float> _texProbaPathNormLon(1, _mesh.vertex().size() );

    mit2 =_mapBassins.begin();
    mend2 = _mapBassins.end();

    n_b2 = 0;

   //on normalise la profondeur dans chaque bassin lon

    for (; mit2 != mend2; ++mit2)
    {
      max_depth2 = 0.0;
      myHistoLat << "bassin lon" << n_b2++ << " " << endl;
      cout << "bassin " << n_b2-1 << " " << (mit2->second).size() << endl;;
      listVertexBassin2 = (mit2->second).begin();
      for (; listVertexBassin2!=(mit2->second).end(); listVertexBassin2++)
      {
      //cout << *listVertexBassin << " " ;
      val2 = _texProbaPath[0].item(*listVertexBassin2);
      if ( val2 > max_depth2)
        max_depth2 = val2;
      }

      cout << max_depth2 << endl;

      if ((mit2->second).size() > 150)
      {
        listVertexBassin2 = (mit2->second).begin();
        for (; listVertexBassin2!=(mit2->second).end(); listVertexBassin2++)
        {
        if (max_depth2!=0)
          _texProbaPathNormLon[0].item(*listVertexBassin2) =  (float)(_texProbaPath[0].item(*listVertexBassin2))/max_depth2;
        else
          _texProbaPathNormLon[0].item(*listVertexBassin2) = 0 ;

        myHistoLat << _texProbaPathNormLon[0].item(*listVertexBassin2) << "\t" << *listVertexBassin2 << "\n";
        }
      }
    }

    string probpathnormLon = "texProbaLon_embc.tex";
    Writer<TimeTexture<float> > texWProbaPathNormLon(probpathnormLon);
    texWProbaPathNormLon << _texProbaPathNormLon;

    if (_adrLonGeodesicOut!="")
    {
      Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
      texWLon << texOutLon;
    }

    //on fusionne les textures lat et lon
    for (uint i = 0; i < _mesh.vertex().size(); i++)
    {
       if (texOutLon[0].item(i) != 0)
         texOutLat[0].item(i) = texOutLon[0].item(i);
    }

    Writer<TimeTexture<float> > texWLines(_adrLines);
    texWLines << texOutLat;

    myHistoLat.close();

}
void SulcalLinesGeodesic::computeGraphDijkstra (AimsSurfaceTriangle surface, int constraintType,int strain)
{
  // compute and copy curvature

  cout << _texBassinsDepthNorm[0].nItem() << endl;

//  if (_adrBassinsDepthNorm=" "}
//    vector<float> &curv = _texCurv[0].data();
//  else
//
  vector<float> &curv = _texBassinsDepthNorm[0].data();

  //val = (float)fabs(p2 - p1)/(SPath[SPath.size()-2 - i].distance(&SPath[SP
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

double SulcalLinesGeodesic::computeDepthShortestPathSulci(unsigned source, unsigned target, vector<geodesic::SurfacePoint> & SPath, vector<int> &listIndexVertexPathSP)
{
  vector<geodesic::SurfacePoint> sources;
  sources.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[source]));

  vector<geodesic::SurfacePoint> targets;
  targets.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[target]));

  //printf("indice source = %d target = %d \n",source, target);

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
//  cout << endl;

  double length = 0;
  float p1,p2,val,val2;
//    for (unsigned i = 0; i < listIndexVertexPathSP.size() - 1; i++)
//      cout << listIndexVertexPathSP[i+1] << " " ;

  if(!SPath.empty())
  {
    for(unsigned i=0; i<SPath.size()-1; ++i)
    {
      p1 = _texBassinsDepthNorm[0].item(listIndexVertexPathSP[i]);
      p2 = _texBassinsDepthNorm[0].item(listIndexVertexPathSP[i+1]);

//      val = (float)fabs(p2 - p1)/(SPath[SPath.size()-2 - i].distance(&SPath[SPath.size()-i-1]));
//      val2 = (float)(SPath[SPath.size()-2 - i].distance(&SPath[SPath.size()-i-1]))* (1 - 5*sqrt(val));
//
//      length += val2;

      length += SPath[i].distance(&SPath[i+1]);

    }
  }

  for(unsigned i=0; i<SPath.size(); ++i)
   {
    _texProbaPath[0].item(listIndexVertexPathSP[i])++;
   }

  return length;

}

double SulcalLinesGeodesic::saveHistoTemp(unsigned source, unsigned target)
{
  vector<int> listIndexVertexPathSP;
  vector<geodesic::SurfacePoint> SPath;

  SPath.clear();
  vector<geodesic::SurfacePoint> sources;
  sources.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[source]));

  vector<geodesic::SurfacePoint> targets;
  targets.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[target]));

  //printf("indice source = %d target = %d \n",source, target);

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
//  cout << endl;

  double length = 0;
  float p1,p2,val;

  myHistoLat << SPath.size() <<"\n";
  myHistoLat << "\nhisto\n";

  if(!SPath.empty())
  {
    for(unsigned i=0; i<SPath.size(); ++i)
    {
      p1 = _texBassinsDepthNorm[0].item(listIndexVertexPathSP[i]);
      myHistoLat << 1-p1 << "\t";
    }
  }

  myHistoLat << "\ngradient\n";
  if(!SPath.empty())
  {
    for(unsigned i=0; i<SPath.size()-1; ++i)
    {
      p1 = _texBassinsDepthNorm[0].item(listIndexVertexPathSP[i]);
      p2 = _texBassinsDepthNorm[0].item(listIndexVertexPathSP[i+1]);

      //val = (float)fabs(p2 - p1)/(SPath[SPath.size()-2 - i].distance(&SPath[SPath.size()-i-1]));
      val = (float)fabs(p2 - p1);

      myHistoLat << val << "\t";
    }
  }

  myHistoLat << "\npath\n";
  if(!SPath.empty())
  {
    for(unsigned i=0; i<SPath.size()-1; ++i)
    {
      myHistoLat << SPath[i].distance(&SPath[i+1]) << "\t";
    }
  }
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

  dist_max = -100000;

  for (j=0; j<points.size(); j++)
  {
    for(i=j+1; i<points.size(); i++)
    {
      nb_combinaison++;
      dist = computeShortestPathSulci( points[j], points[i],SPath,listIndexVertexPathSP);

      //cout << "(" << points[j] << "-" << points[i] << " " << dist << ") " << endl;
      //cout << " (" << i << "-" << j << ") " ;
      std::cout << "\r" << nb_combinaison << "/" << nb_combinaison_max << std::flush;

      if (dist >= dist_max)
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

vector<int> SulcalLinesGeodesic::maxGeodesicDistanceDepthStable(vector<int> points, int constraint, int* s, int *d)
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

  dist_max = -100000;

  for (j=0; j<points.size(); j++)
  {

    for(i=j+1; i<points.size(); i++)
    {
      //myHistoLat <<  "\t" << points[j] << "\t" << points[i] <<"\t";
      nb_combinaison++;
      //dist = computeShortestPathSulci( points[j], points[i],SPath,listIndexVertexPathSP);

      dist = computeDepthShortestPathSulci( points[j], points[i],SPath,listIndexVertexPathSP);
      //cout << "(" << points[j] << "-" << points[i] << " " << dist << ") ";
      //cout << " (" << i << "-" << j << ") " ;

      //myHistoLat << dist << "\t" << points[j] << "\t" << points[i] << "\t";

//      myHistoLat << "\n" << dist << " (" << points[j] << "," << points[i] << ")\n";
      //saveHistoTemp (points[j],points[i]);

      std::cout << "\r" << i << " " << j << " "<< nb_combinaison << "/" << nb_combinaison_max << std::flush;

      if (dist >= dist_max)
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

  //myHistoLat << "\ndist = " << dist_max << "\n";
  return listIndexVertexPathSP_result;
}

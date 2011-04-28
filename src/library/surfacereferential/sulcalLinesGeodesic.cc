/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *    CEA/DSV/SHFJ
 *    4 place du General Leclerc
 *    91401 Orsay cedex
 *    France
 */

#include <cortical_surface/surfacereferential/sulcalLinesGeodesic.h>

#include <aims/mesh/surfaceOperation.h>
#include <aims/scalespace/meshDiffuse.h>
#include <aims/distancemap/meshmorphomat.h>

using namespace std;

SulcalLinesGeodesic::SulcalLinesGeodesic(  string & adrMesh,string & adrCurv, string & adrGeodesicDepth,
    string & adrRootsLon, string & adrRootsLat, int extremeties_method, int segmentation_method, int strain ) :
    _adrMesh(adrMesh),_adrCurv(adrCurv),_adrGeodesicDepth(adrGeodesicDepth),
    _adrRootsLon(adrRootsLon),_adrRootsLat(adrRootsLat),_segmentation_method(segmentation_method),_extremeties_method(extremeties_method), _strain(strain)
{
  std::cout << "Read mesh : ";
  Reader < AimsSurfaceTriangle > r(adrMesh);
  r.read( _mesh );
  std::cout << "done" << std::endl;

  std::cout << "Read roots texture : ";
  Reader < TimeTexture<short> > rtLon(adrRootsLon);
  rtLon.read( _rootsLon );
  Reader < TimeTexture<short> > rtLat(adrRootsLat);
  rtLat.read( _rootsLat );
  std::cout << "done" << std::endl;

  if (adrCurv=="")
  {
    cout << "compute and save curvature texture : ";
    _texCurv = TimeTexture<float>(1, _mesh.vertex().size());
    _texCurv = AimsMeshCurvature(_mesh[0]);

    writeFloatTexture("curv_new.tex",_texCurv);
  }
  else
  {
    std::cout << "Read curvature texture : ";
    Reader < TimeTexture<float> > rtCurv(adrCurv);
    rtCurv.read( _texCurv );
    std::cout << "done" << std::endl;
  }

  if (adrGeodesicDepth!="")
  {
    std::cout << "Read Geodesic Depth Texture : ";
    Reader < TimeTexture<float> > rtDG(adrGeodesicDepth);
    rtDG.read( _geoDepth );
    std::cout << "done" << std::endl;
  }

  std::cout << "compute neighbours : ";
  _neigh = SurfaceManip::surfaceNeighbours( _mesh );
  std::cout << "done " << std::endl;

}

SulcalLinesGeodesic::~SulcalLinesGeodesic()
{
}

void SulcalLinesGeodesic::run()
{
  std::cout << "START : ";

  std::cout << "Normalize Geodesic Depth Texture with sulcal basins " << endl;

  TimeTexture<short> texBasins(1, _mesh.vertex().size());
  map<int,set<int> > mapBasins;

  segmentationSulcalBasins (_texCurv, texBasins, mapBasins,0.0,1);

  _geoDepthNorm = TimeTexture<float>(1, _mesh.vertex().size());
  normalizeDepthMap (_geoDepth, _geoDepthNorm, mapBasins);

  std::cout << "Save Normalize Geodesic Depth texture : ";
  writeFloatTexture("depth_norm.tex",_geoDepthNorm);

  switch (_extremeties_method)
  {
  case 1 :
    cout << "extraction of extremities method : projection crop by basins" << endl;
    sulcalLinesExtract_projection();
    break;
  case 2 :
    cout << "extraction of extremities method : map of probability" << endl ;
    sulcalLinesExtract_probability();
    break;
  }

}

void SulcalLinesGeodesic::writeShortTexture (string name,TimeTexture<short> &out)
{
  size_t found;
  found = _adrRootsLon.find_last_of(".");
  string adrBasins = _adrRootsLon.substr(0,found-9) + name;
  Writer<TimeTexture<short> > texW(adrBasins);
  texW << out;
  cout << "write " << adrBasins << " done" << endl;
}

void SulcalLinesGeodesic::writeFloatTexture (string name,TimeTexture<float> &out)
{
  size_t found;
  found = _adrRootsLon.find_last_of(".");
  string adrBasins = _adrRootsLon.substr(0,found-9) + name;
  Writer<TimeTexture<float> > texW(adrBasins);
  texW << out;
  cout << "write " << adrBasins << " done" << endl;
}

void SulcalLinesGeodesic::floodFillIter(int indexVertex, float newTextureValue,float oldTextureValue,
    TimeTexture<short> &texBasinsTemp, map<int,set<int> > &mapBasins)
{
  queue<int> stack;
  stack.push(indexVertex);

  int indexCurr;

  set<int> listIndexVertexFill;
  listIndexVertexFill.clear();

  while (!stack.empty())
  {
    indexCurr = stack.front();

    listIndexVertexFill.insert(indexCurr);
    texBasinsTemp[0].item(indexCurr) = newTextureValue;

    stack.pop();

    set<uint> voisins = _neigh[indexCurr];
    set<uint>::iterator voisIt = voisins.begin();

    //on parcourt tous les voisins du sommet
    for (; voisIt != voisins.end(); voisIt++)
    {
      indexCurr = *voisIt;

      set<int>::const_iterator itef;
      itef = listIndexVertexFill.find(indexCurr);

      if (itef != listIndexVertexFill.end())
        continue;

      if ( texBasinsTemp[0].item(indexCurr) == oldTextureValue)
      {
        listIndexVertexFill.insert(indexCurr);
        stack.push(indexCurr);
        //cout << indexCurr << endl;
        texBasinsTemp[0].item(indexCurr) = newTextureValue;
      }

    }
  }

  mapBasins.insert (pair<int,set<int> >(newTextureValue, listIndexVertexFill));

}

void SulcalLinesGeodesic::texBinarizeF2S(TimeTexture<float> &texIn, TimeTexture<short> &texOut, float threshold,int inf,int sup)
{
  for (uint i = 0; i < texIn[0].nItem(); i++)
  {
    if (texIn[0].item(i) < threshold)
      texOut[0].item(i) = inf;
    else
      texOut[0].item(i) = sup;
  }
}

void SulcalLinesGeodesic::texBinarizeS2S(TimeTexture<short> &texIn, TimeTexture<short> &texOut, int threshold,int inf,int sup)
{
  for (uint i = 0; i < texIn[0].nItem(); i++)
  {
    if (texIn[0].item(i) < threshold)
      texOut[0].item(i) = inf;
    else
      texOut[0].item(i) = sup;
  }
}


void SulcalLinesGeodesic::texConnectedComponent(TimeTexture<short> &texBasins, map<int,set<int> > &mapBasins)
{
  int j = 1;
  mapBasins.clear();

  for (uint i = 0; i < texBasins[0].nItem(); i++)
  {
    if (texBasins[0].item(i) == -1)
      floodFillIter(i,j++,-1,texBasins,mapBasins);
  }
}

void SulcalLinesGeodesic::listRootsProjections(TimeTexture<short> &texBasins,set<int> &listIndexLat,set<int> &listIndexLon)
{
  int value;

  for (uint i = 0; i < texBasins[0].nItem(); i++)
  {
    value = texBasins[0].item(i);
    if (value != 0)
    {
      if (_rootsLat[0].item(i)!=0)
        listIndexLat.insert(i);

      if (_rootsLon[0].item(i)!=0)
        listIndexLon.insert(i);
    }
    else
    {
      //on efface les projections qui sont en dehors des bassins
      _rootsLat[0].item(i) = 0;
      _rootsLon[0].item(i) = 0;
    }
  }

}

void SulcalLinesGeodesic::computeListLabelProjectionsBasins (TimeTexture<short> &roots, map<int,set<int> > &mapBasins,set<int> &listIndex, map<int,set<int> > &mapConstraint)
{
  int nbNeigh,value;
  set<int>::iterator itl;
  set<int>::iterator it;
  set<int> listIndexTemp;
  int nbConstraint = 0;

  set<int> listLabelBasins;
  set<int>::iterator itLabel;

  map<int, set<int> >::const_iterator mit(mapBasins.begin()),mend(mapBasins.end());

  //On parcourt les tous bassins
  for (; mit != mend; ++mit)
  {
    listIndexTemp.clear();
    //cout << "\nbasin " << mit->first << endl;

    //recherche des étiquettes présentes dans le bassin
    listLabelBasins.clear();
    for (itl=listIndex.begin(); itl!=listIndex.end(); itl++)
    {
      it=(mit->second).find(*itl);
      value = roots[0].item(*itl);
      if (it != (mit->second).end())
        listLabelBasins.insert(value);
    }
    //pour chaque etiquette du bassin
    for (itLabel=listLabelBasins.begin(); itLabel!=listLabelBasins.end(); itLabel++)
    {
      //cout << *itLabel << " --> ";
      for (itl=listIndex.begin(); itl!=listIndex.end(); itl++)
      {
        it=(mit->second).find(*itl);
        value = roots[0].item(*itl);

        //si une lat est dans le bassin alors on l'ajoute à la liste temporaire des points du bassin de valeur itLabel
        if (it != (mit->second).end() && (value==*itLabel) )
        {
          set<uint> nei = _neigh[*itl];
          set<uint>::iterator neiIt = nei.begin();
          nbNeigh = 0;
          //on parcourt tous les voisins du sommet
          for (; neiIt != nei.end(); neiIt++)
          {
            if (roots[0].item(*neiIt) == value)
              nbNeigh++;
          }

          if (nbNeigh < 2)
            {
            listIndexTemp.insert(*it);
            //cout << *it << " ";
            }
        }
      }
      // on insert la liste temporaire à la map lat du bassin
      if (!listIndexTemp.empty())
        mapConstraint.insert (pair<int, set<int> >(nbConstraint++, listIndexTemp));
    }
  }

}

void SulcalLinesGeodesic::computeLongestPathBasins (TimeTexture<short> &roots, TimeTexture<short> &out, map<int,set<int> > &mapConstraint)
{
  map<int, set<int> >::const_iterator mit(mapConstraint.begin()),mend(mapConstraint.end());

  set<int>::iterator it;
  vector<int> pathTemp,indexTemp;

  for (uint i = 0; i < _mesh.vertex().size(); i++)
    out[0].item(i) = 0.0;

  GeodesicPath sp(_mesh,1,_strain);

  int source,target;
  int constraintValue;

  for (; mit != mend; ++mit)
  {
    indexTemp.clear();

    cout << "\nbasin " << (int)mit->first << " : \n";
    it = (mit->second).begin();

    //on copie la liste des index de sommets set --> vector
    for (; it!=(mit->second).end(); it++)
      indexTemp.push_back (*it);

    //on récupère la valeur de la contrainte
    if (!indexTemp.empty())
    {
      constraintValue = roots.item( *(indexTemp.begin()) );

      //on calcule le plus long chemin
      pathTemp = sp.longestPath_ind(indexTemp, &source, &target);

      for (int i = 0; i < pathTemp.size(); i++)
        out[0].item(pathTemp[i]) = constraintValue;
    }

  }

//  vector<int> pathIndex;
//  pathIndex = sp.shortestPath_1_1_ind(source,target);
//  for (int i = 0; i < pathIndex.size(); i++)
//  cout << pathIndex[i] << " ";
//  vector<Point3df> path3D;
//  path3D = sp.shortestPath_1_1_xyz(source,target);
//  for (int i = 0; i < path3D.size(); i++)
//  cout << path3D[i][0] << " " << path3D[i][1] << " "<< path3D[i][2] << "\n";
//
//  double len = sp.shortestPath_1_1_len (source,target);
//
//  vector<unsigned> targetList;
//  targetList.push_back(40);
//  targetList.push_back(400);
//  targetList.push_back(340);
//  targetList.push_back(403);
//  sp.shortestPath_1_N_ind(source,targetList,&target,&len);
//  cout << "best target = " << target << " length = " << len << endl;
//
}



void SulcalLinesGeodesic::normalizeDepthMap (TimeTexture<float> &depth, TimeTexture<float> &depthNorm, map<int,set<int> > &mapBasins)
{
   map<int, set<int> >::const_iterator mit(mapBasins.begin()),mend(mapBasins.end());
   set<int>::iterator listVertexbasin;
   int n_b = 0;
   float max_depth,val;

   //on normalise la profondeur dans chaque bassin entre 0 et 1
   for (; mit != mend; ++mit)
   {
     max_depth = 0;

     listVertexbasin = (mit->second).begin();
     for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
     {
     val = depth[0].item(*listVertexbasin);
     if ( val > max_depth)
       max_depth = val;
     }

     listVertexbasin = (mit->second).begin();
     for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
     {
     if (max_depth!=0)
       depthNorm[0].item(*listVertexbasin) =  (float)(depth[0].item(*listVertexbasin))/max_depth;
     else
       depthNorm[0].item(*listVertexbasin) = 0 ;
     }
   }

 }

void SulcalLinesGeodesic::segmentationSulcalBasins (TimeTexture<float> &texIn,TimeTexture<short> &texBasins,map<int,set<int> > &mapBasins,float threshold, int ESsize)
{
//binarisation texture de courbure < 0.0
texBinarizeF2S(texIn, texBasins, threshold , -1, 0);

//fermeture des bassins
TimeTexture<short> texBasinsDil(1, _mesh.vertex().size() );
texBasinsDil[0]=MeshDilation<short>( _mesh[0], texBasins[0], short(0), -1, ESsize, true);
TimeTexture<short> texBasinsErode(1, _mesh.vertex().size() );
texBasinsErode[0]=MeshErosion<short>( _mesh[0], texBasinsDil[0], short(0), -1, ESsize , true);

texBinarizeS2S(texBasinsErode, texBasins, 0 ,-1 ,0);

// étiquetage des composantes connexes
texConnectedComponent(texBasins, mapBasins);
cout << mapBasins.size() << " Basins found" << endl;
}


void SulcalLinesGeodesic::sulcalLinesExtract_projection()
{
  TimeTexture<short> texBasins(1, _mesh.vertex().size());
  map<int,set<int> > mapBasins;

  segmentationSulcalBasins (_texCurv, texBasins, mapBasins,0.0,1);

  std::cout << "Save connected components texture : ";
  writeShortTexture("basins.tex",texBasins);

  //liste les projections roots et conservent seulement celles inclues dans les bassins
  set<int> listIndexLon,listIndexLat;
  listRootsProjections(texBasins,listIndexLat,listIndexLon);

  //textures contenant les contraintes lat et lon
  TimeTexture<short> texSulcalinesLat(1, _mesh.vertex().size() );
  TimeTexture<short> texSulcalinesLon(1, _mesh.vertex().size() );

  map<int, set<int> > mapConstraintLat,mapConstraintLon;

  cout << "Sort Lat/Lon by basin : "<< endl;

  //groupe et étiquette dans des map les projections lat et lon de chaque bassin
  //conserve seulement les projections qui sont potentiellement des extremités des lignes (points qui ont au plus un voisin)

  computeListLabelProjectionsBasins (_rootsLat,mapBasins,listIndexLat,mapConstraintLat);
  cout << endl << mapConstraintLat.size() << " Basins Latitude extracted";
  computeLongestPathBasins (_rootsLat, texSulcalinesLat, mapConstraintLat);

  computeListLabelProjectionsBasins (_rootsLon,mapBasins,listIndexLon,mapConstraintLon);
  cout << endl << mapConstraintLon.size() << " Basins Longitude extracted";
  computeLongestPathBasins (_rootsLon, texSulcalinesLon, mapConstraintLon);

  std::cout << "\nSave sulcal lines texture : " << endl;

  writeShortTexture("lat_sulcal_lines.tex",texSulcalinesLat);
  writeShortTexture("lon_sulcal_lines.tex",texSulcalinesLon);
}

void SulcalLinesGeodesic::sulcalLinesExtract_probability()
{
  TimeTexture<short> texBasins(1, _mesh.vertex().size());
  map<int,set<int> > mapBasins;

  segmentationSulcalBasins (_texCurv, texBasins, mapBasins,0.0,1);

  std::cout << "Save connected components texture : ";
  size_t found;
  found = _adrRootsLon.find_last_of(".");
  string adrBasins = _adrRootsLon.substr(0,found-9) + "basins.tex";
  Writer<TimeTexture<short> > texWB(adrBasins);
  texWB << texBasins;
  cout << adrBasins << " done" << endl;

  //
  //
  //  //on dilate les roots
  //
  //   TimeTexture<short> texProjectionLatDil(1, _mesh.vertex().size() );
  //   texProjectionLatDil[0]=MeshDilation<short>( _mesh[0], _rootsLat[0], short(0), -1, 6, false);
  //   TimeTexture<short> texProjectionLonDil(1, _mesh.vertex().size() );
  //   texProjectionLonDil[0]=MeshDilation<short>( _mesh[0], _rootsLon[0], short(0), -1, 6, false);
  //
  //   string lats = "texLatRootsDilate.tex";
  //   Writer<TimeTexture<short> > texWLats(lats);
  //   texWLats << texProjectionLatDil;
  //
  //   string lons = "texLonRootsDilate.tex";
  //   Writer<TimeTexture<short> > texWLons(lons);
  //   texWLons << texProjectionLonDil;
  //
  //
  //   // on fait les intersections avec les basins
  //
  //
  //    TimeTexture<short> texExtremiteLat(1, _mesh.vertex().size() );
  //    TimeTexture<short> texExtremiteLon(1, _mesh.vertex().size() );
  //
  //    vector< list<unsigned> > neighbourso( _mesh.vertex().size());
  //    neighbourso = AimsMeshOrderNode(_mesh[0]);
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      if (_texbasins[0].item(i)!=0)
  //      {
  //        //si il y a une intersection avec la texture lat dilatée alors on attribue le label de la région
  //        if (texProjectionLatDil[0].item(i) > 0)
  //        {
  //        texExtremiteLat[0].item(i) = texProjectionLatDil[0].item(i);
  //        //_listIndexLat.insert(i);
  //        }
  //        if (texProjectionLonDil[0].item(i) > 0)
  //        {
  //        texExtremiteLon[0].item(i) = texProjectionLonDil[0].item(i);
  //        //_listIndexLon.insert(i);
  //        }
  //      }
  //    }
  //
  //
  //
  //    cout << "re compute label of basins with lat inter proj: ";
  //    j = 200;
  //    _mapbasins.clear();
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //      {
  //        if (texExtremiteLat[0].item(i) > 0)
  //          _texbasins[0].item(i) = texExtremiteLat[0].item(i);
  //        else
  //          _texbasins[0].item(i) = 0;
  //      }
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      if (_texbasins[0].item(i) > 0 && _texbasins[0].item(i) < 200)
  //      {
  //      floodFillIter(i,j++,_texbasins[0].item(i));
  //      }
  //    }
  //
  //    _texbasinsLat = TimeTexture<short>(1, _mesh.vertex().size());
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //    if (_texbasins[0].item(i)!=0)
  //      _texbasinsLat[0].item(i) = _texbasins[0].item(i) - 200;
  //    _texbasins[0].item(i) = texbasinsSave[0].item(i);
  //    }
  //
  //    cout << "done\n";
  //
  //    int nb_voisins;
  //    int value;
  //
  //    TimeTexture<short> texContourLatDilInter(1, _mesh.vertex().size() );
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      value = texExtremiteLat[0].item(i);
  //
  //      //on marque les points de contour
  //      if (value > 0)
  //      {
  //        set<uint> voisins = _neigh[i];
  //        set<uint>::iterator voisIt = voisins.begin();
  //        nb_voisins = 0;
  //
  //        //on parcourt tous les voisins du sommet
  //        for (; voisIt != voisins.end(); voisIt++)
  //        {
  //          if (texExtremiteLat[0].item(*voisIt) != value)
  //            {
  //            nb_voisins++;
  //            texContourLatDilInter.item(i) = value;
  //            continue;
  //            }
  //        }
  //      }
  //    }
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      if (texContourLatDilInter[0].item(i)!=0)
  //      {
  //       _listIndexLat.insert(i);
  //      }
  //    }
  //
  //    string latExt = "texLatExtInter.tex";
  //    Writer<TimeTexture<short> > texWLatExt(latExt);
  //    texWLatExt << texExtremiteLat;
  //
  //    string latExt2 = "texLatContourExtInter.tex";
  //    Writer<TimeTexture<short> > texWLatExt2(latExt2);
  //    texWLatExt2 << texContourLatDilInter;
  //
  //
  //
  //    string toto8 = "texbasinsInterLatCachiaEtiquette.tex";
  //    Writer<TimeTexture<short> > texW8(toto8);
  //    texW8 << _texbasinsLat;
  //
  //
  //    mit = _mapbasins.begin();
  //    mend = _mapbasins.end();
  //
  //    cout << "sort constraints lat by basins : ";
  //
  //    cout << _mapbasins.size() << "\n";
  //
  //    int lat;
  //    set<int>::iterator itlat;
  //    set<int>::iterator ittemp;
  //    set<int>::iterator it;
  //
  //    set<int> _listIndexTemp;
  //
  //    int nbbasinsConstaintLat = 0;
  //
  //    //pour chaque basin i
  //    for (; mit != mend; ++mit)
  //    {
  //      _listIndexTemp.clear();
  //
  //      //cout << (int)mit->first << ": ";
  //      listVertexbasin = (mit->second).begin();
  //
  //      //on parcourt la liste des contraintes lat
  //      if (!_listIndexLat.empty())
  //      for (itlat=_listIndexLat.begin(); itlat!=_listIndexLat.end(); itlat++)
  //      {
  //        it=(mit->second).find(*itlat);
  //
  //        value = texContourLatDilInter[0].item(*itlat);
  //        //si une lat est dans le basin alors je l'enlève de la liste
  //        //et je l'ajoute à la liste temporaire des points du basin i
  //        if (it != (mit->second).end())
  //        {
  //          if (value != 0)
  //          {
  //              _listIndexTemp.insert(*it);
  //          }
  //
  //        }
  //      }
  //
  //      // on associe la liste temporaire à la map du basin
  //      if (!_listIndexTemp.empty())
  //      {
  //        //temp = make_pair (_rootsLat.item(*it),_listIndexTemp);
  //        //cout << _rootsLat.item(*it) << " ";
  //        _mapConstraintLat.insert (pair<int, set<int> >(nbbasinsConstaintLat++, _listIndexTemp));
  //      }
  //    }
  //
  //
  //    cout << "nb basin Lat= " << nbbasinsConstaintLat<< endl;
  //
  //    map<int, set<int> >::const_iterator mclatit(_mapConstraintLat.begin()),mclatend(_mapConstraintLat.end());
  //
  //    set<int>::iterator itp1;
  //    vector<int> _vectorIndexTemp;
  //
  //    //textures contenant les contraintes lat
  //    TimeTexture<float> texOutLat(1, _mesh.vertex().size() );
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //     texOutLat[0].item(i) = 0.0;
  //    }
  //
  //    _texProbaPath = TimeTexture<short>(1, _mesh.vertex().size());
  //
  //    //pour chaque basin contenant des contraintes de latitude
  //
  //    for (; mclatit != mclatend; ++mclatit)
  //    {
  //     _vectorIndexTemp.clear();
  //
  //     cout << "basin " << (int)mclatit->first << " : ";
  //     itp1 = (mclatit->second).begin();
  //
  //     //ARN DEBUG
  //     //if (mclatit->first == 10000000)
  //     {
  //       //myHistoLat << "basin " << (int)mclatit->first << "\n";
  //       //on parcourt la liste des contraintes lat
  //       for (; itp1!=(mclatit->second).end(); itp1++)
  //       {
  //         //cout << *itp1 << "#" << _rootsLat.item(*itp1) << " ";
  //         _vectorIndexTemp.push_back (*itp1);
  //       }
  //
  //       //cout << _vectorIndexTemp.size() << " - " ;
  //       int v1,v2;
  //
  //       vector<int>::iterator itv1;
  //
  //       int i;
  //       vector<int> _vectorIndexTempConstraint;
  //       if (_vectorIndexTemp.size() > 1)
  //       {
  //         while (!_vectorIndexTemp.empty())
  //         {
  //           itv1 = _vectorIndexTemp.begin();
  //           //v1 = _rootsLat.item(*itv1);
  //           v1 = texContourLatDilInter.item(*itv1);
  //
  //           _vectorIndexTempConstraint.clear();
  //
  //           for (i=0; i<_vectorIndexTemp.size(); i++)
  //           {
  //             //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
  //             //v2 = _rootsLat.item(_vectorIndexTemp[i]);
  //             v2 = texContourLatDilInter.item(_vectorIndexTemp[i]);
  //             if (v1 == v2 && (v1==45))
  //             {
  //               _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
  //               //cout << _vectorIndexTemp[i] << " ";
  //               _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
  //               i--;
  //             }
  //              else
  //                if (v1 == v2)
  //                {
  //                //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
  //                //cout << _vectorIndexTemp[i] << " ";
  //                _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
  //                i--;
  //                }
  //           }
  //
  //           int source,target;
  //
  //           cout << _vectorIndexTempConstraint.size() << endl;
  //           if (_vectorIndexTempConstraint.size()>=30)
  //           {
  //             vector<int> listIndexVertexPathSP;
  //
  //
  //             //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
  //             listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);
  //
  //             cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
  //
  //            // myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
  //             //saveHistoTemp (source,target);
  //
  //             //listIndexVertexPathSP = computeShortestPathSulci(source,target);
  //
  //             for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
  //               {
  //               //cout << listIndexVertexPathSP[t] << " ";
  //               texOutLat[0].item(listIndexVertexPathSP[t]) = (float) v1;
  //               }
  //           }
  //         }
  //         //cout << endl;
  //       }
  //     }
  //     //ARN DEBUG
  //    }
  //
  //    if (_adrLatGeodesicOut!="")
  //    {
  //     Writer<TimeTexture<float> > texWLat(_adrLatGeodesicOut);
  //     texWLat << texOutLat;
  //    }
  //
  //    string probpathLat = "texProbaLat.tex";
  //    Writer<TimeTexture<short> > texWProbaLat(probpathLat);
  //    texWProbaLat << _texProbaPath;
  //
  //    //on normalise la texture des probas des lat
  //    TimeTexture<float> _texProbaPathNormLat(1, _mesh.vertex().size() );
  //
  //    map<int, set<int> >::const_iterator mit2(_mapbasins.begin()),mend2(_mapbasins.end());
  //    set<int>::iterator listVertexbasin2;
  //    int n_b2 = 0;
  //    float max_depth2,val2;
  //
  //      //on normalise la profondeur dans chaque basin lat
  //
  //    for (; mit2 != mend2; ++mit2)
  //    {
  //      max_depth2 = 0.0;
  //      myHistoLat << "basin lat" << n_b2++ << " " << endl;
  //      cout << "basin " << n_b2-1 << " " << (mit2->second).size() << endl;;
  //      listVertexbasin2 = (mit2->second).begin();
  //      for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //      {
  //      //cout << *listVertexbasin << " " ;
  //      val2 = _texProbaPath[0].item(*listVertexbasin2);
  //      if ( val2 > max_depth2)
  //        max_depth2 = val2;
  //      }
  //
  //      if ( (mit2->second).size() > 150 )
  //      {
  //        cout << max_depth2 << endl;
  //
  //        listVertexbasin2 = (mit2->second).begin();
  //        for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //        {
  //        if (max_depth2!=0)
  //          _texProbaPathNormLat[0].item(*listVertexbasin2) =  (float)(_texProbaPath[0].item(*listVertexbasin2))/max_depth2;
  //        else
  //          _texProbaPathNormLat[0].item(*listVertexbasin2) = 0 ;
  //
  //        myHistoLat << _texProbaPathNormLat[0].item(*listVertexbasin2) << "\t" << *listVertexbasin2 << "\n";
  //        }
  //      }
  //    }
  //
  //    string probpathLatN = "texProbaLat_embc.tex";
  //    Writer<TimeTexture<float> > texWProbaLatN(probpathLatN);
  //    texWProbaLatN << _texProbaPathNormLat;
  //
  //
  //
  //
  //
  //    TimeTexture<float> _texProbaPathNormLatThresh(1, _mesh.vertex().size() );
  //
  //        for (uint i = 0; i < _mesh.vertex().size(); i++)
  //        {
  //           if (_texProbaPathNormLat[0].item(i) > 0.4)
  //             _texProbaPathNormLatThresh[0].item(i) = 1;
  //           else
  //             _texProbaPathNormLatThresh[0].item(i) = 0;
  //        }
  //
  //
  //
  //
  //
  //
  //
  //        TimeTexture<float> _texProbaPathNormLat7(1, _mesh.vertex().size() );
  //
  //        //on cherche le plus long chemin
  //
  //        GeodesicPath sp(_mesh,_texbasinsDepthNorm,1,-3);
  //        std::vector<int> listIndexVertexTarget;
  //
  //
  //        mit2 =_mapbasins.begin();
  //        mend2 = _mapbasins.end();
  //
  //        for (; mit2 != mend2; ++mit2)
  //        {
  //          listIndexVertexTarget.clear();
  //
  //          //listVertexbasin2 = (mit2->second).begin();
  //          //for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //
  //          if ((mit2->second).size() > 150)
  //          {
  //            listVertexbasin2 = (mit2->second).begin();
  //            for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //            {
  //              //listVertexbasin2++;
  //              if (_texProbaPathNormLatThresh[0].item(*listVertexbasin2)==1)
  //                listIndexVertexTarget.push_back(*listVertexbasin2);
  //            }
  //
  //            vector<int> pathIndex;
  //            int s,d;
  //
  //            if (listIndexVertexTarget.size() > 10)
  //              {
  //              pathIndex = sp.longestPath_ind(listIndexVertexTarget, &s, &d);
  //
  //              for (int i = 0; i < pathIndex.size(); i++)
  //                _texProbaPathNormLat7[0].item(pathIndex[i]) = 45;
  //              }
  //          }
  //        }
  //
  //        string probpathnormLatBest = "texProbaLat_embc_best_depth.tex";
  //        Writer<TimeTexture<float> > texWProbaPathNormLatBest(probpathnormLatBest);
  //        texWProbaPathNormLatBest << _texProbaPathNormLat7;
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      _texProbaPath[0].item(i) = 0;
  //    }
  //
  //    cout << "re compute label of basins with lon inter proj: \n";
  //    j = 200;
  //    _mapbasins.clear();
  //    _texbasinsLon = TimeTexture<short>(1, _mesh.vertex().size());
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //       _texbasinsLon[0].item(i) = _texbasins[0].item(i) - 200;
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //      {
  //        if (texExtremiteLon[0].item(i) > 0)
  //          _texbasins[0].item(i) = texExtremiteLon[0].item(i);
  //        else
  //          _texbasins[0].item(i) = 0;
  //      }
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      if (_texbasins[0].item(i) > 0 && _texbasins[0].item(i) < 200)
  //      {
  //      floodFillIter(i,j++,_texbasins[0].item(i));
  //      }
  //    }
  //
  //    _texbasinsLon = TimeTexture<short>(1, _mesh.vertex().size());
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //    if (_texbasins[0].item(i)!=0)
  //      _texbasinsLon[0].item(i) = _texbasins[0].item(i) - 200;
  //    _texbasins[0].item(i) = texbasinsSave[0].item(i);
  //    }
  //
  //    cout << "done\n";
  //
  //    string toto9 = "texbasinsInterLonCachiaEtiquette.tex";
  //    Writer<TimeTexture<short> > texW9(toto9);
  //    texW9 << _texbasinsLon;
  //
  //    cout << endl;
  //
  //    TimeTexture<short> texContourLonDilInter(1, _mesh.vertex().size() );
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      value = texExtremiteLon[0].item(i);
  //
  //      //on marque les points de contour
  //      if (value > 0)
  //      {
  //        set<uint> voisins = _neigh[i];
  //        set<uint>::iterator voisIt = voisins.begin();
  //        nb_voisins = 0;
  //
  //        //on parcourt tous les voisins du sommet
  //        for (; voisIt != voisins.end(); voisIt++)
  //        {
  //          if (texExtremiteLon[0].item(*voisIt) != value)
  //            {
  //            nb_voisins++;
  //            texContourLonDilInter.item(i) = value;
  //            continue;
  //            }
  //        }
  //      }
  //    }
  //
  //     for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //      if (texContourLonDilInter[0].item(i)!=0)
  //      {
  //       _listIndexLon.insert(i);
  //      }
  //    }
  //
  //    string lonExt = "texLonExtInter.tex";
  //    Writer<TimeTexture<short> > texWLonExt(lonExt);
  //    texWLonExt << texExtremiteLon;
  //
  //    string lonExt2 = "texLonContourExtInter.tex";
  //    Writer<TimeTexture<short> > texWLonExt2(lonExt2);
  //    texWLonExt2 << texContourLonDilInter;
  //
  //    cout << "done\n";
  //
  //
  //
  //    mit = _mapbasins.begin();
  //    mend = _mapbasins.end();
  ////
  //    cout << "sort constraints lon by basins : ";
  ////
  //    int lon;
  //    set<int>::iterator itlon;
  //
  //    int nbbasinsConstaintLon = 0;
  ////
  //    //pour chaque basin i
  //    for (; mit != mend; ++mit)
  //    {
  //      _listIndexTemp.clear();
  //
  //      //cout << (int)mit->first << ": ";
  //      listVertexbasin = (mit->second).begin();
  //
  //      //on parcourt la liste des contraintes lat
  //      if (!_listIndexLon.empty())
  //      for (itlon=_listIndexLon.begin(); itlon!=_listIndexLon.end(); itlon++)
  //      {
  //        it=(mit->second).find(*itlon);
  //
  //        value = texContourLonDilInter[0].item(*itlon);
  //        //si une lat est dans le basin alors je l'enlève de la liste
  //        //et je l'ajoute à la liste temporaire des points du basin i
  //        if (it != (mit->second).end())
  //        {
  //          if (value != 0)
  //          {
  //              _listIndexTemp.insert(*it);
  //          }
  //
  //        }
  //      }
  //    // on associe la liste temporaire à la map du basin
  //      if (!_listIndexTemp.empty())
  //      {
  //        //temp = make_pair (_rootsLon.item(*it),_listIndexTemp);
  //        _mapConstraintLon.insert (pair<int, set<int> >(nbbasinsConstaintLon++, _listIndexTemp));
  //      }
  //
  //    }
  //
  //
  //    cout << "nb basin Lon= " << nbbasinsConstaintLon<< endl;
  //
  //    map<int, set<int> >::const_iterator mclonit(_mapConstraintLon.begin()),mclonend(_mapConstraintLon.end());
  //
  //    //textures contenant les contraintes lon
  //    TimeTexture<float> texOutLon(1, _mesh.vertex().size() );
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //     texOutLon[0].item(i) = 0.0;
  //    }
  //
  //    //pour chaque basin contenant des contraintes de longitude
  //
  //    for (; mclonit != mclonend; ++mclonit)
  //    {
  //     _vectorIndexTemp.clear();
  //
  //     cout << "basin " << (int)mclonit->first << " : ";
  //     itp1 = (mclonit->second).begin();
  //
  //     //ARN DEBUG
  //    // if (mclonit->first == 100000)
  //     {
  //       //myHistoLat << "basin " << (int)mclatit->first << "\n";
  //       //on parcourt la liste des contraintes lat
  //       for (; itp1!=(mclonit->second).end(); itp1++)
  //       {
  //         //cout << *itp1 << "#" << _rootsLat.item(*itp1) << " ";
  //         _vectorIndexTemp.push_back (*itp1);
  //       }
  //
  //       //cout << _vectorIndexTemp.size() << " - " ;
  //       int v1,v2;
  //
  //       vector<int>::iterator itv1;
  //
  //       int i;
  //       vector<int> _vectorIndexTempConstraint;
  //       if (_vectorIndexTemp.size() > 1)
  //       {
  //         while (!_vectorIndexTemp.empty())
  //         {
  //           itv1 = _vectorIndexTemp.begin();
  //           //v1 = _rootsLat.item(*itv1);
  //           v1 = texContourLonDilInter.item(*itv1);
  //
  //           _vectorIndexTempConstraint.clear();
  //
  //           for (i=0; i<_vectorIndexTemp.size(); i++)
  //           {
  //             //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
  //             //v2 = _rootsLat.item(_vectorIndexTemp[i]);
  //             v2 = texContourLonDilInter.item(_vectorIndexTemp[i]);
  //             if (v1 == v2 && (v1==25))
  //             {
  //               _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
  //               //cout << _vectorIndexTemp[i] << " ";
  //               _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
  //               i--;
  //             }
  //            else
  //              if (v1 == v2)
  //              {
  //              //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
  //              //cout << _vectorIndexTemp[i] << " ";
  //              _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
  //              i--;
  //              }
  //           }
  //
  //           int source,target;
  //           cout << _vectorIndexTempConstraint.size() << endl;
  //           if (_vectorIndexTempConstraint.size()>=30)
  //           {
  //             vector<int> listIndexVertexPathSP;
  //
  //             //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
  //             listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);
  //
  //             cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
  //
  //             //myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
  //             //saveHistoTemp (source,target);
  //
  //             //listIndexVertexPathSP = computeShortestPathSulci(source,target);
  //
  //             for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
  //               {
  //               //cout << listIndexVertexPathSP[t] << " ";
  //               texOutLon[0].item(listIndexVertexPathSP[t]) = (float) v1;
  //               }
  //           }
  //         }
  //         //cout << endl;
  //       }
  //     }
  //     //ARN DEBUG
  //    }
  //
  //    if (_adrLonGeodesicOut!="")
  //    {
  //     Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
  //     texWLon << texOutLon;
  //    }
  //
  //    string probpathLon = "texProbaLon.tex";
  //    Writer<TimeTexture<short> > texWProbaLon(probpathLon);
  //    texWProbaLon << _texProbaPath;
  //
  //
  //    //on normalise la texture des probas des lon
  //    TimeTexture<float> _texProbaPathNormLon(1, _mesh.vertex().size() );
  //
  //    mit2 =_mapbasins.begin();
  //    mend2 = _mapbasins.end();
  //
  //    n_b2 = 0;
  //
  //    for (; mit2 != mend2; ++mit2)
  //    {
  //      max_depth2 = 0.0;
  //      //myHistoLat << "basin lon" << n_b2++ << " " << endl;
  //      //cout << "basin " << n_b2-1 << " " << (mit2->second).size() << endl;;
  //      listVertexbasin2 = (mit2->second).begin();
  //      for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //      {
  //      //cout << *listVertexbasin << " " ;
  //      val2 = _texProbaPath[0].item(*listVertexbasin2);
  //      if ( val2 > max_depth2)
  //        max_depth2 = val2;
  //      }
  //
  //      //cout << max_depth2 << endl;
  //
  //      if ((mit2->second).size() > 150)
  //      {
  //        listVertexbasin2 = (mit2->second).begin();
  //        for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //        {
  //        if (max_depth2!=0)
  //          _texProbaPathNormLon[0].item(*listVertexbasin2) =  (float)(_texProbaPath[0].item(*listVertexbasin2))/max_depth2;
  //        else
  //          _texProbaPathNormLon[0].item(*listVertexbasin2) = 0 ;
  //
  //        myHistoLat << _texProbaPathNormLon[0].item(*listVertexbasin2) << "\t" << *listVertexbasin2 << "\n";
  //        }
  //      }
  //    }
  //
  //    string probpathnormLon = "texProbaLon_embc.tex";
  //    Writer<TimeTexture<float> > texWProbaPathNormLon(probpathnormLon);
  //    texWProbaPathNormLon << _texProbaPathNormLon;
  //
  //    TimeTexture<float> _texProbaPathNormLonThresh(1, _mesh.vertex().size() );
  //
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //       if (_texProbaPathNormLon[0].item(i) > 0.3)
  //         _texProbaPathNormLonThresh[0].item(i) = 1;
  //       else
  //         _texProbaPathNormLonThresh[0].item(i) = 0;
  //    }
  //
  //
  //    TimeTexture<float> _texProbaPathNormLon7(1, _mesh.vertex().size() );
  //
  //    //on cherche le plus long chemin
  //
  //    mit2 =_mapbasins.begin();
  //    mend2 = _mapbasins.end();
  //
  //    for (; mit2 != mend2; ++mit2)
  //    {
  //      listIndexVertexTarget.clear();
  //
  //      //listVertexbasin2 = (mit2->second).begin();
  //      //for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //
  //      if ((mit2->second).size() > 150)
  //      {
  //        listVertexbasin2 = (mit2->second).begin();
  //        for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
  //        {
  //          //listVertexbasin2++;
  //          if (_texProbaPathNormLonThresh[0].item(*listVertexbasin2)==1)
  //            listIndexVertexTarget.push_back(*listVertexbasin2);
  //        }
  //
  //        vector<int> pathIndex;
  //        int s,d;
  //
  //        if (listIndexVertexTarget.size() > 20)
  //          {
  //          pathIndex = sp.longestPath_ind(listIndexVertexTarget, &s, &d);
  //
  //          for (int i = 0; i < pathIndex.size(); i++)
  //            _texProbaPathNormLon7[0].item(pathIndex[i]) = 25;
  //          }
  //      }
  //    }
  //
  //    string probpathnormLonBest = "texProbaLon_embc_best_depth.tex";
  //    Writer<TimeTexture<float> > texWProbaPathNormLonBest(probpathnormLonBest);
  //    texWProbaPathNormLonBest << _texProbaPathNormLon7;
  //
  //
  //
  //
  //
  //
  //    if (_adrLonGeodesicOut!="")
  //    {
  //      Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
  //      texWLon << texOutLon;
  //    }
  //
  //    //on fusionne les textures lat et lon
  //    for (uint i = 0; i < _mesh.vertex().size(); i++)
  //    {
  //       if (texOutLon[0].item(i) != 0)
  //         texOutLat[0].item(i) = texOutLon[0].item(i);
  //    }
  //
  //    Writer<TimeTexture<float> > texWLines(_adrLines);
  //    texWLines << texOutLat;
  //
  //    myHistoLat.close();
}

//
//void SulcalLinesGeodesic::basinsDetect2()
//{
//  //writing path in the output texture
//  _texbasins = TimeTexture<short>(1, _mesh.vertex().size());
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//    if (_texCurv[0].item(i) < 0.0)
//      _texbasins[0].item(i) = 1;
//    else
//      _texbasins[0].item(i) = 0;
//  }
//
//  //idée 1 : ouverture des basins (pour séparer les basins connectés par de petits plis)
////  TimeTexture<short> texbasinsErode(1, _mesh.vertex().size() );
////  texbasinsErode[0]=MeshErosion<short>( _mesh[0], _texbasins[0], short(0), -1, 3 , false);
////  TimeTexture<short> texbasinsDil(1, _mesh.vertex().size() );
////  texbasinsDil[0]=MeshDilation<short>( _mesh[0], texbasinsErode[0], short(0), -1, 2, false);
//
//  //pb avec une simple ouverture on a tendance à couper en deux le sillon central par exemple
//  //mieux vaut garder de grands basins (connectés par des plis de passage), l'étiquetage des sillons contribuera a séparer les lignes ...
//
//  //idée2 : fermeture des basins
//  TimeTexture<short> texbasinsDil(1, _mesh.vertex().size() );
//   texbasinsDil[0]=MeshDilation<short>( _mesh[0], _texbasins[0], short(0), -1, 1, true);
//  TimeTexture<short> texbasinsErode(1, _mesh.vertex().size() );
//  texbasinsErode[0]=MeshErosion<short>( _mesh[0], texbasinsDil[0], short(0), -1, 1 , true);
//
//
//
//  //rajouter peut être en paramètre le facteur de dilatation ?
//
//  cout << "compute label of basins : ";
//  int j = 2;
//  _mapbasins.clear();
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//    if (texbasinsErode[0].item(i) > 0)
//      _texbasins[0].item(i) = -1;
//    else
//      _texbasins[0].item(i) = 0;
//    }
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//    if (_texbasins[0].item(i) == -1)
//    {
//    floodFillIter(i,j++,-1);
//    }
//  }
//
//  cout << "done\n";
//
//  string toto0 = "texbasinsCloseEtiquette.tex";
//  Writer<TimeTexture<short> > texW0(toto0);
//  texW0 << _texbasins;
//
//  cout << endl;
//
//  //pour chaque basin i
//  _texbasinsDepthNorm = TimeTexture<float>(1, _mesh.vertex().size());
////  for (uint i = 0; i < _mesh.vertex().size(); i++)
////    _texbasinsDepthNorm[0].item(i) = 1;
//
//  map<int, set<int> >::const_iterator mit(_mapbasins.begin()),mend(_mapbasins.end());
//  set<int>::iterator listVertexbasin;
//  int n_b = 0;
//  float max_depth,val;
//
//  //on normalise la profondeur dans chaque basin
//
//  if (_adrbasinsDepthNorm!="")
//  {
//  Reader < TimeTexture<float> > rtDepthNorm(_adrbasinsDepthNorm);
//  rtDepthNorm.read( _texbasinsDepthNorm );
//  cout << "read basins depth Norm\n";
//  }
//  else
//  {
//    for (; mit != mend; ++mit)
//    {
//      max_depth = 0;
//
//      //cout << "basin " << n_b++ << " " << endl;;
//      listVertexbasin = (mit->second).begin();
//      for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
//      {
//      //cout << *listVertexbasin << " " ;
//      val = _geoDepth[0].item(*listVertexbasin);
//      if ( val > max_depth)
//        max_depth = val;
//      }
//      //cout << max_depth << endl;
//
//      listVertexbasin = (mit->second).begin();
//      for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
//      {
//      if (max_depth!=0)
//        _texbasinsDepthNorm[0].item(*listVertexbasin) =  (float)(_geoDepth[0].item(*listVertexbasin))/max_depth;
//      else
//        _texbasinsDepthNorm[0].item(*listVertexbasin) = 0 ;
//      }
//    }
//
//    string toto1 = "texbasinsDepthNorm.tex";
//    Writer<TimeTexture<float> > texW1(toto1);
//    texW1 << _texbasinsDepthNorm;
//  }
//
//  //on binarise les basins normalisés
////  TimeTexture<short> texbasinsDepthNormBin(1, _mesh.vertex().size() );
////  for (uint i = 0; i < _mesh.vertex().size(); i++)
////  {
////    if (_texbasinsDepthNorm[0].item(i) < 0.1)
////      texbasinsDepthNormBin[0].item(i) = 0;
////    else
////    texbasinsDepthNormBin[0].item(i) = 1;
////  }
//
////  string toto = "/home/arnaud/Bureau/texbasinsDepthNormBin.tex";
////  Writer<TimeTexture<short> > texW2(toto);
////  texW2 << texbasinsDepthNormBin;
////
////  for (uint i = 0; i < _mesh.vertex().size(); i++)
////  {
////    if (texbasinsDepthNormBin[0].item(i) == 1)
////      _texbasins[0].item(i) = -1;
////    else
////      _texbasins[0].item(i) = 0;
////  }
////
////  cout << "re-compute label of basins norm : ";
////  j = 1;
////  _mapbasins.clear();
////
////  for (uint i = 0; i < _mesh.vertex().size(); i++)
////  {
////    if (_texbasins[0].item(i) == -1)
////    {
////    floodFillIter(i,j++,-1);
////    }
////  }
////
////  cout << endl;
////  string toto5 = "/home/arnaud/Bureau/texbasinsEtiquette.tex";
////  Writer<TimeTexture<short> > texW5(toto5);
////  texW5 << _texbasins;
//
//  //on dilate les roots
//
////   TimeTexture<short> texProjectionLatDil(1, _mesh.vertex().size() );
////   texProjectionLatDil[0]=MeshDilation<short>( _mesh[0], _rootsLat[0], short(0), -1, 6, false);
////   TimeTexture<short> texProjectionLonDil(1, _mesh.vertex().size() );
////   texProjectionLonDil[0]=MeshDilation<short>( _mesh[0], _rootsLon[0], short(0), -1, 6, false);
////
////   string lats = "texLatRootsDilate.tex";
////   Writer<TimeTexture<short> > texWLats(lats);
////   texWLats << texProjectionLatDil;
////
////   string lons = "texLonRootsDilate.tex";
////   Writer<TimeTexture<short> > texWLons(lons);
////   texWLons << texProjectionLonDil;
////
////   // on fait les intersections avec les basins
////   TimeTexture<short> texProjectionLatDilInter(1, _mesh.vertex().size() );
////   TimeTexture<short> texProjectionLonDilInter(1, _mesh.vertex().size() );
////
////   mit = _mapbasins.begin();
////   mend = _mapbasins.end();
////
////   for (; mit != mend; ++mit)
////   {
////     listVertexbasin = (mit->second).begin();
////     for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
////     {
////     //cout << *listVertexbasin << " " ;
////     //geoDepth[0].item(*listVertexbasin);
////     }
////
//////     listVertexbasin = (mit->second).begin();
//////     for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
//////     {
//////     if (max_depth!=0)
//////       _texbasinsDepthNorm[0].item(*listVertexbasin) =  (float)(_geoDepth[0].item(*listVertexbasin))/max_depth;
//////     else
//////       _texbasinsDepthNorm[0].item(*listVertexbasin) = 0 ;
//////     }
////   }
////
////  if (_adrLines!="")
////  {
////    TimeTexture<short> texH0Squel(1, _mesh.vertex().size() );
////    int nb_voisins;
////
////    vector< list<unsigned> > neighbourso( _mesh.vertex().size());
////    neighbourso = AimsMeshOrderNode(_mesh[0]);
////    texH0Squel[0]=MeshSkeletization<short> ( _mesh[0], texbasinsDil[0], short(1), short(0), neighbourso );
////
////    int value;
////
////    for (uint i = 0; i < _mesh.vertex().size(); i++)
////    {
////      if (texH0Squel[0].item(i)==1)
////        {
////        value = 1;
////        texH0Squel[0].item(i) = 180;
////        }
////      else
////      value = _texbasins[0].item(i);
////
////      //on marque les points de contour
////      if (value > 0)
////      {
////        set<uint> voisins = _neigh[i];
////        set<uint>::iterator voisIt = voisins.begin();
////        nb_voisins = 0;
////
////        //on parcourt tous les voisins du sommet
////        if (value == 1)
////        {
////          for (; voisIt != voisins.end(); voisIt++)
////          {
////            if (texH0Squel[0].item(*voisIt) == 1 || texH0Squel[0].item(*voisIt) == 180 || texH0Squel[0].item(*voisIt) == 360 || texH0Squel[0].item(*voisIt) == 250)
////              nb_voisins++;
////          }
////
////          if (nb_voisins > 2)
////            texH0Squel[0].item(i) = 250;
////
////          if (nb_voisins < 2)
////            texH0Squel[0].item(i) = 360;
////
////        }
////        else
////        {
////          //on parcourt tous les voisins du sommet
////          voisIt = voisins.begin();
////          for (; voisIt != voisins.end(); voisIt++)
////          {
////            if (_texbasins[0].item(*voisIt) != value &&  _texbasins[0].item(*voisIt)!=20)
////            {
////              texH0Squel[0].item(i) = 20;
////              continue;
////            }
////          }
////        }
////      }
////    }
////
////
////
////
////    TimeTexture<short> texExtremiteLat(1, _mesh.vertex().size() );
////    TimeTexture<short> texExtremiteLon(1, _mesh.vertex().size() );
////
////
////    for (uint i = 0; i < _mesh.vertex().size(); i++)
////    {
////      if (texH0Squel[0].item(i)==20)
////      {
////        //si il y a une intersection avec la texture lat dilatée alors on attribue le label de la région
////        if (texProjectionLatDil[0].item(i) > 0)
////        {
////        texExtremiteLat[0].item(i) = texProjectionLatDil[0].item(i);
////        _listIndexLat.insert(i);
////        }
////        if (texProjectionLonDil[0].item(i) > 0)
////        {
////        texExtremiteLon[0].item(i) = texProjectionLonDil[0].item(i);
////        _listIndexLon.insert(i);
////        }
////      }
////
////
////    }
////
////
//////    for (uint i = 0; i < _mesh.vertex().size(); i++)
//////    {
//////      if (texH0Squel[0].item(i)==250 || texH0Squel[0].item(i)==360 )
//////      {
//////        //si il y a une intersection avec la texture lat dilatée alors on attribue le label de la région
//////        if (texProjectionLatDil[0].item(i) > 0)
//////        {
//////        texExtremiteLat[0].item(i) = texProjectionLatDil[0].item(i);
//////        _listIndexLat.insert(i);
//////        }
//////        if (texProjectionLonDil[0].item(i) > 0)
//////        {
//////        texExtremiteLon[0].item(i) = texProjectionLonDil[0].item(i);
//////        _listIndexLon.insert(i);
//////        }
//////      }
//////    }
////
////    //  cout << endl;
////      string toto2 = "texbasinsContourSquel.tex";
////      Writer<TimeTexture<short> > texW2(toto2);
////      texW2 << texH0Squel;
////
////      string toto3 = "texbasinsExtermiteLat.tex";
////      Writer<TimeTexture<short> > texW3(toto3);
////      texW3 << texExtremiteLat;
////
////      string toto4 = "texbasinsExtermiteLon.tex";
////      Writer<TimeTexture<short> > texW4(toto4);
////      texW4 << texExtremiteLon;
//////    Writer<TimeTexture<short> > texW(_adrLines);
//////    texW << texH0Squel;
////
////      mit = _mapbasins.begin();
////      mend = _mapbasins.end();
////
////      cout << "sort constraints lat/lon by basins : ";
////
////        int lat,lon;
////        set<int>::iterator itlon;
////        set<int>::iterator itlat;
////        set<int>::iterator ittemp;
////        set<int>::iterator it;
////
////        set<int> _listIndexTemp;
////
////        int nbbasinsConstaintLat = 0;
////        int nbbasinsConstaintLon = 0;
////
////        //pour chaque basin i
////        for (; mit != mend; ++mit)
////        {
////          _listIndexTemp.clear();
////
////          //cout << (int)mit->first << ": ";
////          listVertexbasin = (mit->second).begin();
////
////          //on parcourt la liste des contraintes lat
////          if (!_listIndexLat.empty())
////          for (itlat=_listIndexLat.begin(); itlat!=_listIndexLat.end(); itlat++)
////          {
////            it=(mit->second).find(*itlat);
////
////            value = texExtremiteLat[0].item(*itlat);
////            //si une lat est dans le basin alors je l'enlève de la liste
////            //et je l'ajoute à la liste temporaire des points du basin i
////            if (it != (mit->second).end())
////            {
////              if (value != 0)
////              {
////                  _listIndexTemp.insert(*it);
////              }
////
////            }
////          }
////
////          // on associe la liste temporaire à la map du basin
////          if (!_listIndexTemp.empty())
////          {
////            //temp = make_pair (_rootsLat.item(*it),_listIndexTemp);
////            //cout << _rootsLat.item(*it) << " ";
////            _mapConstraintLat.insert (pair<int, set<int> >(nbbasinsConstaintLat++, _listIndexTemp));
////          }
////
////          _listIndexTemp.clear();
////
////          if (!_listIndexLon.empty())
////          for (itlon=_listIndexLon.begin(); itlon!=_listIndexLon.end(); itlon++)
////          {
////            it=(mit->second).find(*itlon);
////
////            value = texExtremiteLon[0].item(*itlon);
////            //si une lat est dans le basin alors je l'enlève de la liste
////            //et je l'ajoute à la liste temporaire des points du basin i
////            if (it != (mit->second).end())
////            {
////              if (value != 0)
////              _listIndexTemp.insert(*it);
////            }
////          }
////          // on associe la liste temporaire à la map du basin
////          if (!_listIndexTemp.empty())
////          {
////            //temp = make_pair (_rootsLon.item(*it),_listIndexTemp);
////            _mapConstraintLon.insert (pair<int, set<int> >(nbbasinsConstaintLon++, _listIndexTemp));
////          }
////          }
////
////        cout << "done" << endl;
////
////
////        myHistoLat.open ("../histoLat.txt");
////
////        //cout << _listIndexLat.size() << " points lat /" << _listIndexLon.size() << " points lon " << endl;
////        cout << "nb basin Lat= " << nbbasinsConstaintLat<< endl;
////
////        map<int, set<int> >::const_iterator mclatit(_mapConstraintLat.begin()),mclatend(_mapConstraintLat.end());
////        map<int, set<int> >::const_iterator mclonit(_mapConstraintLon.begin()),mclonend(_mapConstraintLon.end());
////
////        set<int>::iterator itp1;
////
////        vector<int> _vectorIndexTemp;
////
////        //textures contenant les contraintes lat et lon
////        TimeTexture<float> texOutLat(1, _mesh.vertex().size() );
////        TimeTexture<float> texOutLon(1, _mesh.vertex().size() );
////
////        for (uint i = 0; i < _mesh.vertex().size(); i++)
////        {
////          texOutLat[0].item(i) = 0.0;
////          texOutLon[0].item(i) = 0.0;
////        }
////
////        _texProbaPath = TimeTexture<short>(1, _mesh.vertex().size());
////
////        //pour chaque basin contenant des contraintes de latitude
////
////        for (; mclatit != mclatend; ++mclatit)
////        {
////          _vectorIndexTemp.clear();
////
////
////
////          cout << "basin " << (int)mclatit->first << " : ";
////          itp1 = (mclatit->second).begin();
////
////          //ARN DEBUG
////          //if (mclatit->first == 4)
////          {
////            myHistoLat << "basin " << (int)mclatit->first << "\n";
////            //on parcourt la liste des contraintes lat
////            for (; itp1!=(mclatit->second).end(); itp1++)
////            {
////              //cout << *itp1 << "#" << _rootsLat.item(*itp1) << " ";
////              _vectorIndexTemp.push_back (*itp1);
////            }
////
////            //cout << _vectorIndexTemp.size() << " - " ;
////            int v1,v2;
////
////            vector<int>::iterator itv1;
////
////            int i;
////            vector<int> _vectorIndexTempConstraint;
////            if (_vectorIndexTemp.size() > 1)
////            {
////              while (!_vectorIndexTemp.empty())
////              {
////                itv1 = _vectorIndexTemp.begin();
////                //v1 = _rootsLat.item(*itv1);
////                v1 = texExtremiteLat.item(*itv1);
////
////                _vectorIndexTempConstraint.clear();
////
////                for (i=0; i<_vectorIndexTemp.size(); i++)
////                {
////                  //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
////                  //v2 = _rootsLat.item(_vectorIndexTemp[i]);
////                  v2 = texExtremiteLat.item(_vectorIndexTemp[i]);
////                  if (v1 == v2 && (v1==45 || v1 == 5))
////                  {
////                    _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
////                    //cout << _vectorIndexTemp[i] << " ";
////                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
////                    i--;
////                  }
////                  else
////                    if (v1 == v2)
////                    {
////                    //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
////                    //cout << _vectorIndexTemp[i] << " ";
////                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
////                    i--;
////                    }
////                }
////
////                int source,target;
////
////                if (_vectorIndexTempConstraint.size()>=2)
////                {
////                  vector<int> listIndexVertexPathSP;
////
////                  //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
////                  listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);
////
////                  cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
////
////                  myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
////                  saveHistoTemp (source,target);
////
////                  //listIndexVertexPathSP = computeShortestPathSulci(source,target);
////
////                  for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
////                    {
////                    //cout << listIndexVertexPathSP[t] << " ";
////                    texOutLat[0].item(listIndexVertexPathSP[t]) = (float) v1;
////                    }
////                }
////              }
////              //cout << endl;
////            }
////          }
////          //ARN DEBUG
////        }
////
////        if (_adrLatGeodesicOut!="")
////        {
////          Writer<TimeTexture<float> > texWLat(_adrLatGeodesicOut);
////          texWLat << texOutLat;
////        }
////
////
////        //pour chaque basin contenant des contraintes de longitude
////
////        cout << "nb basin Lon= " << nbbasinsConstaintLon<< endl;
////
////        for (; mclonit != mclonend; ++mclonit)
////        {
////          _vectorIndexTemp.clear();
////
////          cout << "basin " << (int)mclonit->first << " : ";
////          itp1 = (mclonit->second).begin();
////
////          //ARN DEBUG
////          if (mclonit->first == 5000)
////          {
////
////            myHistoLat << "basin " << (int)mclonit->first << "\n";
////            //on parcourt la liste des contraintes lon
////            for (; itp1!=(mclonit->second).end(); itp1++)
////            {
////              //cout << *itp1 << "#" << _rootsLon.item(*itp1) << " ";
////              _vectorIndexTemp.push_back (*itp1);
////            }
////
////            //cout << _vectorIndexTemp.size() << " - " ;
////            int v1,v2;
////
////            vector<int>::iterator itv1;
////
////            int i;
////            vector<int> _vectorIndexTempConstraint;
////
////            if (_vectorIndexTemp.size() > 1)
////            {
////              while (!_vectorIndexTemp.empty())
////              {
////                itv1 = _vectorIndexTemp.begin();
////                v1 = texExtremiteLon.item(*itv1);
////
////                _vectorIndexTempConstraint.clear();
////
////                for (i=0; i<_vectorIndexTemp.size(); i++)
////                {
////                  //cout << _vectorIndexTemp[i] << "#" << _rootsLon.item(_vectorIndexTemp[i]) << " ";
////                  v2 = texExtremiteLon.item(_vectorIndexTemp[i]);
////                  if (v1 == v2 && v1==25)
////                  {
////                    _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
////                    //cout << _vectorIndexTemp[i] << " ";
////                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
////                    i--;
////                  }
////                  else
////                    if (v1 == v2)
////                    {
////                    //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
////                    //cout << _vectorIndexTemp[i] << " ";
////                    _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
////                    i--;
////                  }
////                }
////
////                int source,target;
////
////                vector<int> listIndexVertexPathSP;
////
////                if (_vectorIndexTempConstraint.size()>=2)
////                {
////                  //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
////                  listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);
////
////                  //listIndexVertexPathSP = computeShortestPathSulci(source,target);
////
////                  //myHistoLat << "value" << v1 << "(" << source << "," << target << ")\n";
////
////                  cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
////
////                  //saveHistoTemp (source,target);
////
////                  for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
////                    {
////                    //cout << listIndexVertexPathSP[t] << " ";
////                    texOutLon[0].item(listIndexVertexPathSP[t]) = (float) v1;
////                    }
////                }
////              }
////
////              //cout << endl;
////            }
////          }
////          //ARN DEBUG
////      }
////
////    myHistoLat.close();
////
////
////
////    string probpath = "texProbaPathNew.tex";
////    Writer<TimeTexture<short> > texWProbaPath(probpath);
////    //Reader<TimeTexture<short> > texWProbaPath(probpath);
////    //texWProbaPath.read( _texProbaPath );
////    texWProbaPath << _texProbaPath;
////
////    TimeTexture<float> _texProbaPathNorm(1, _mesh.vertex().size() );
////
////    map<int, set<int> >::const_iterator mit2(_mapbasins.begin()),mend2(_mapbasins.end());
////    set<int>::iterator listVertexbasin2;
////    int n_b2 = 0;
////    float max_depth2,val2;
////
////      //on normalise la profondeur dans chaque basin
////
////    for (; mit2 != mend2; ++mit2)
////    {
////      max_depth2 = 0;
////
////      cout << "basin " << n_b++ << " " << endl;;
////      listVertexbasin2 = (mit2->second).begin();
////      for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
////      {
////      //cout << *listVertexbasin << " " ;
////      val2 = _texProbaPath[0].item(*listVertexbasin2);
////      if ( val2 > max_depth2)
////        max_depth2 = val2;
////      }
////
////      cout << max_depth2 << endl;
////
////      listVertexbasin2 = (mit2->second).begin();
////      for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
////      {
////      if (max_depth2!=0)
////        _texProbaPathNorm[0].item(*listVertexbasin2) =  (float)(_texProbaPath[0].item(*listVertexbasin2))/max_depth2;
////      else
////        _texProbaPathNorm[0].item(*listVertexbasin2) = 0 ;
////      }
////    }
////
////    string probpathnorm = "texProbaPathNewNorm.tex";
////    Writer<TimeTexture<float> > texWProbaPathNorm(probpathnorm);
////    texWProbaPathNorm << _texProbaPathNorm;
////
////
////
////    if (_adrLonGeodesicOut!="")
////    {
////      Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
////      texWLon << texOutLon;
////    }
////
////    //on fusionne les textures lat et lon
////    for (uint i = 0; i < _mesh.vertex().size(); i++)
////    {
////       if (texOutLon[0].item(i) != 0)
////         texOutLat[0].item(i) = texOutLon[0].item(i);
////    }
////
////    Writer<TimeTexture<float> > texWLines(_adrLines);
////    texWLines << texOutLat;
////
////
////  }
//
//
//}
//
//
//void SulcalLinesGeodesic::basinsDetect3()
//{
//
//  myHistoLat.open ("./histoProba.txt");
//
//  //writing path in the output texture
//  _texbasins = TimeTexture<short>(1, _mesh.vertex().size());
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//    if (_texCurv[0].item(i) < 0.0)
//      _texbasins[0].item(i) = 1;
//    else
//      _texbasins[0].item(i) = 0;
//  }
//
//
//  //idée 1 : ouverture des basins (pour séparer les basins connectés par de petits plis)
////  TimeTexture<short> texbasinsErode(1, _mesh.vertex().size() );
////  texbasinsErode[0]=MeshErosion<short>( _mesh[0], _texbasins[0], short(0), -1, 3 , false);
////  TimeTexture<short> texbasinsDil(1, _mesh.vertex().size() );
////  texbasinsDil[0]=MeshDilation<short>( _mesh[0], texbasinsErode[0], short(0), -1, 2, false);
//
//  //pb avec une simple ouverture on a tendance à couper en deux le sillon central par exemple
//  //mieux vaut garder de grands basins (connectés par des plis de passage), l'étiquetage des sillons contribuera a séparer les lignes ...
//
//  //idée2 : fermeture des basins
//  TimeTexture<short> texbasinsDil(1, _mesh.vertex().size() );
//   texbasinsDil[0]=MeshDilation<short>( _mesh[0], _texbasins[0], short(0), -1, 1, true);
//  TimeTexture<short> texbasinsErode(1, _mesh.vertex().size() );
//  texbasinsErode[0]=MeshErosion<short>( _mesh[0], texbasinsDil[0], short(0), -1, 1 , true);
//
//
//  //rajouter peut être en paramètre le facteur de dilatation ?
//
//  cout << "compute label of basins : ";
//  int j = 2;
//  _mapbasins.clear();
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//    if (texbasinsErode[0].item(i) > 0)
//      _texbasins[0].item(i) = -1;
//    else
//      _texbasins[0].item(i) = 0;
//    }
//
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//    if (_texbasins[0].item(i) == -1)
//    {
//    floodFillIter(i,j++,-1);
//    }
//  }
//
//  TimeTexture<short> texbasinsSave(1, _mesh.vertex().size() );
//  for (uint i = 0; i < _mesh.vertex().size(); i++)
//  {
//  texbasinsSave[0].item(i) = _texbasins[0].item(i);
//  }
//
//  cout << "done\n";
//
//  string toto0 = "texbasinsCloseEtiquette.tex";
//  Writer<TimeTexture<short> > texW0(toto0);
//  texW0 << _texbasins;
//
//  cout << endl;
//
//  //pour chaque basin i
//  //_texbasinsDepthNorm = TimeTexture<float>(1, _mesh.vertex().size());
////  for (uint i = 0; i < _mesh.vertex().size(); i++)
////    _texbasinsDepthNorm[0].item(i) = 1;
//
//  map<int, set<int> >::const_iterator mit(_mapbasins.begin()),mend(_mapbasins.end());
//  set<int>::iterator listVertexbasin;
//  int n_b = 0;
//  float max_depth,val;
//
//  //on normalise la profondeur dans chaque basin
//
//  if (_adrbasinsDepthNorm=="")
//  {
//    for (; mit != mend; ++mit)
//    {
//      max_depth = 0;
//
//      //cout << "basin " << n_b++ << " " << endl;;
//      listVertexbasin = (mit->second).begin();
//      for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
//      {
//      //cout << *listVertexbasin << " " ;
//      val = _geoDepth[0].item(*listVertexbasin);
//      if ( val > max_depth)
//        max_depth = val;
//      }
//      //cout << max_depth << endl;
//
//      listVertexbasin = (mit->second).begin();
//      for (; listVertexbasin!=(mit->second).end(); listVertexbasin++)
//      {
//      if (max_depth!=0)
//        _texbasinsDepthNorm[0].item(*listVertexbasin) =  (float)(_geoDepth[0].item(*listVertexbasin))/max_depth;
//      else
//        _texbasinsDepthNorm[0].item(*listVertexbasin) = 0 ;
//      }
//    }
//
//    string toto1 = "texbasinsDepthNorm_e10.tex";
//    Writer<TimeTexture<float> > texW1(toto1);
//    texW1 << _texbasinsDepthNorm;
//  }
//
//
//  //on dilate les roots
//
//   TimeTexture<short> texProjectionLatDil(1, _mesh.vertex().size() );
//   texProjectionLatDil[0]=MeshDilation<short>( _mesh[0], _rootsLat[0], short(0), -1, 6, false);
//   TimeTexture<short> texProjectionLonDil(1, _mesh.vertex().size() );
//   texProjectionLonDil[0]=MeshDilation<short>( _mesh[0], _rootsLon[0], short(0), -1, 6, false);
//
//   string lats = "texLatRootsDilate.tex";
//   Writer<TimeTexture<short> > texWLats(lats);
//   texWLats << texProjectionLatDil;
//
//   string lons = "texLonRootsDilate.tex";
//   Writer<TimeTexture<short> > texWLons(lons);
//   texWLons << texProjectionLonDil;
//
//
//   // on fait les intersections avec les basins
//
//
//    TimeTexture<short> texExtremiteLat(1, _mesh.vertex().size() );
//    TimeTexture<short> texExtremiteLon(1, _mesh.vertex().size() );
//
//    vector< list<unsigned> > neighbourso( _mesh.vertex().size());
//    neighbourso = AimsMeshOrderNode(_mesh[0]);
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      if (_texbasins[0].item(i)!=0)
//      {
//        //si il y a une intersection avec la texture lat dilatée alors on attribue le label de la région
//        if (texProjectionLatDil[0].item(i) > 0)
//        {
//        texExtremiteLat[0].item(i) = texProjectionLatDil[0].item(i);
//        //_listIndexLat.insert(i);
//        }
//        if (texProjectionLonDil[0].item(i) > 0)
//        {
//        texExtremiteLon[0].item(i) = texProjectionLonDil[0].item(i);
//        //_listIndexLon.insert(i);
//        }
//      }
//    }
//
//
//
//    cout << "re compute label of basins with lat inter proj: ";
//    j = 200;
//    _mapbasins.clear();
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//      {
//        if (texExtremiteLat[0].item(i) > 0)
//          _texbasins[0].item(i) = texExtremiteLat[0].item(i);
//        else
//          _texbasins[0].item(i) = 0;
//      }
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      if (_texbasins[0].item(i) > 0 && _texbasins[0].item(i) < 200)
//      {
//      floodFillIter(i,j++,_texbasins[0].item(i));
//      }
//    }
//
//    _texbasinsLat = TimeTexture<short>(1, _mesh.vertex().size());
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//    if (_texbasins[0].item(i)!=0)
//      _texbasinsLat[0].item(i) = _texbasins[0].item(i) - 200;
//    _texbasins[0].item(i) = texbasinsSave[0].item(i);
//    }
//
//    cout << "done\n";
//
//    int nb_voisins;
//    int value;
//
//    TimeTexture<short> texContourLatDilInter(1, _mesh.vertex().size() );
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      value = texExtremiteLat[0].item(i);
//
//      //on marque les points de contour
//      if (value > 0)
//      {
//        set<uint> voisins = _neigh[i];
//        set<uint>::iterator voisIt = voisins.begin();
//        nb_voisins = 0;
//
//        //on parcourt tous les voisins du sommet
//        for (; voisIt != voisins.end(); voisIt++)
//        {
//          if (texExtremiteLat[0].item(*voisIt) != value)
//            {
//            nb_voisins++;
//            texContourLatDilInter.item(i) = value;
//            continue;
//            }
//        }
//      }
//    }
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      if (texContourLatDilInter[0].item(i)!=0)
//      {
//       _listIndexLat.insert(i);
//      }
//    }
//
//    string latExt = "texLatExtInter.tex";
//    Writer<TimeTexture<short> > texWLatExt(latExt);
//    texWLatExt << texExtremiteLat;
//
//    string latExt2 = "texLatContourExtInter.tex";
//    Writer<TimeTexture<short> > texWLatExt2(latExt2);
//    texWLatExt2 << texContourLatDilInter;
//
//
//
//    string toto8 = "texbasinsInterLatCachiaEtiquette.tex";
//    Writer<TimeTexture<short> > texW8(toto8);
//    texW8 << _texbasinsLat;
//
//
//    mit = _mapbasins.begin();
//    mend = _mapbasins.end();
//
//    cout << "sort constraints lat by basins : ";
//
//    cout << _mapbasins.size() << "\n";
//
//    int lat;
//    set<int>::iterator itlat;
//    set<int>::iterator ittemp;
//    set<int>::iterator it;
//
//    set<int> _listIndexTemp;
//
//    int nbbasinsConstaintLat = 0;
//
//    //pour chaque basin i
//    for (; mit != mend; ++mit)
//    {
//      _listIndexTemp.clear();
//
//      //cout << (int)mit->first << ": ";
//      listVertexbasin = (mit->second).begin();
//
//      //on parcourt la liste des contraintes lat
//      if (!_listIndexLat.empty())
//      for (itlat=_listIndexLat.begin(); itlat!=_listIndexLat.end(); itlat++)
//      {
//        it=(mit->second).find(*itlat);
//
//        value = texContourLatDilInter[0].item(*itlat);
//        //si une lat est dans le basin alors je l'enlève de la liste
//        //et je l'ajoute à la liste temporaire des points du basin i
//        if (it != (mit->second).end())
//        {
//          if (value != 0)
//          {
//              _listIndexTemp.insert(*it);
//          }
//
//        }
//      }
//
//      // on associe la liste temporaire à la map du basin
//      if (!_listIndexTemp.empty())
//      {
//        //temp = make_pair (_rootsLat.item(*it),_listIndexTemp);
//        //cout << _rootsLat.item(*it) << " ";
//        _mapConstraintLat.insert (pair<int, set<int> >(nbbasinsConstaintLat++, _listIndexTemp));
//      }
//    }
//
//
//    cout << "nb basin Lat= " << nbbasinsConstaintLat<< endl;
//
//    map<int, set<int> >::const_iterator mclatit(_mapConstraintLat.begin()),mclatend(_mapConstraintLat.end());
//
//    set<int>::iterator itp1;
//    vector<int> _vectorIndexTemp;
//
//    //textures contenant les contraintes lat
//    TimeTexture<float> texOutLat(1, _mesh.vertex().size() );
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//     texOutLat[0].item(i) = 0.0;
//    }
//
//    _texProbaPath = TimeTexture<short>(1, _mesh.vertex().size());
//
//    //pour chaque basin contenant des contraintes de latitude
//
//    for (; mclatit != mclatend; ++mclatit)
//    {
//     _vectorIndexTemp.clear();
//
//     cout << "basin " << (int)mclatit->first << " : ";
//     itp1 = (mclatit->second).begin();
//
//     //ARN DEBUG
//     //if (mclatit->first == 10000000)
//     {
//       //myHistoLat << "basin " << (int)mclatit->first << "\n";
//       //on parcourt la liste des contraintes lat
//       for (; itp1!=(mclatit->second).end(); itp1++)
//       {
//         //cout << *itp1 << "#" << _rootsLat.item(*itp1) << " ";
//         _vectorIndexTemp.push_back (*itp1);
//       }
//
//       //cout << _vectorIndexTemp.size() << " - " ;
//       int v1,v2;
//
//       vector<int>::iterator itv1;
//
//       int i;
//       vector<int> _vectorIndexTempConstraint;
//       if (_vectorIndexTemp.size() > 1)
//       {
//         while (!_vectorIndexTemp.empty())
//         {
//           itv1 = _vectorIndexTemp.begin();
//           //v1 = _rootsLat.item(*itv1);
//           v1 = texContourLatDilInter.item(*itv1);
//
//           _vectorIndexTempConstraint.clear();
//
//           for (i=0; i<_vectorIndexTemp.size(); i++)
//           {
//             //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
//             //v2 = _rootsLat.item(_vectorIndexTemp[i]);
//             v2 = texContourLatDilInter.item(_vectorIndexTemp[i]);
//             if (v1 == v2 && (v1==45))
//             {
//               _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
//               //cout << _vectorIndexTemp[i] << " ";
//               _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
//               i--;
//             }
//              else
//                if (v1 == v2)
//                {
//                //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
//                //cout << _vectorIndexTemp[i] << " ";
//                _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
//                i--;
//                }
//           }
//
//           int source,target;
//
//           cout << _vectorIndexTempConstraint.size() << endl;
//           if (_vectorIndexTempConstraint.size()>=30)
//           {
//             vector<int> listIndexVertexPathSP;
//
//
//             //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
//             listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);
//
//             cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
//
//            // myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
//             //saveHistoTemp (source,target);
//
//             //listIndexVertexPathSP = computeShortestPathSulci(source,target);
//
//             for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
//               {
//               //cout << listIndexVertexPathSP[t] << " ";
//               texOutLat[0].item(listIndexVertexPathSP[t]) = (float) v1;
//               }
//           }
//         }
//         //cout << endl;
//       }
//     }
//     //ARN DEBUG
//    }
//
//    if (_adrLatGeodesicOut!="")
//    {
//     Writer<TimeTexture<float> > texWLat(_adrLatGeodesicOut);
//     texWLat << texOutLat;
//    }
//
//    string probpathLat = "texProbaLat.tex";
//    Writer<TimeTexture<short> > texWProbaLat(probpathLat);
//    texWProbaLat << _texProbaPath;
//
//    //on normalise la texture des probas des lat
//    TimeTexture<float> _texProbaPathNormLat(1, _mesh.vertex().size() );
//
//    map<int, set<int> >::const_iterator mit2(_mapbasins.begin()),mend2(_mapbasins.end());
//    set<int>::iterator listVertexbasin2;
//    int n_b2 = 0;
//    float max_depth2,val2;
//
//      //on normalise la profondeur dans chaque basin lat
//
//    for (; mit2 != mend2; ++mit2)
//    {
//      max_depth2 = 0.0;
//      myHistoLat << "basin lat" << n_b2++ << " " << endl;
//      cout << "basin " << n_b2-1 << " " << (mit2->second).size() << endl;;
//      listVertexbasin2 = (mit2->second).begin();
//      for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//      {
//      //cout << *listVertexbasin << " " ;
//      val2 = _texProbaPath[0].item(*listVertexbasin2);
//      if ( val2 > max_depth2)
//        max_depth2 = val2;
//      }
//
//      if ( (mit2->second).size() > 150 )
//      {
//        cout << max_depth2 << endl;
//
//        listVertexbasin2 = (mit2->second).begin();
//        for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//        {
//        if (max_depth2!=0)
//          _texProbaPathNormLat[0].item(*listVertexbasin2) =  (float)(_texProbaPath[0].item(*listVertexbasin2))/max_depth2;
//        else
//          _texProbaPathNormLat[0].item(*listVertexbasin2) = 0 ;
//
//        myHistoLat << _texProbaPathNormLat[0].item(*listVertexbasin2) << "\t" << *listVertexbasin2 << "\n";
//        }
//      }
//    }
//
//    string probpathLatN = "texProbaLat_embc.tex";
//    Writer<TimeTexture<float> > texWProbaLatN(probpathLatN);
//    texWProbaLatN << _texProbaPathNormLat;
//
//
//
//
//
//    TimeTexture<float> _texProbaPathNormLatThresh(1, _mesh.vertex().size() );
//
//        for (uint i = 0; i < _mesh.vertex().size(); i++)
//        {
//           if (_texProbaPathNormLat[0].item(i) > 0.4)
//             _texProbaPathNormLatThresh[0].item(i) = 1;
//           else
//             _texProbaPathNormLatThresh[0].item(i) = 0;
//        }
//
//
//
//
//
//
//
//        TimeTexture<float> _texProbaPathNormLat7(1, _mesh.vertex().size() );
//
//        //on cherche le plus long chemin
//
//        GeodesicPath sp(_mesh,_texbasinsDepthNorm,1,-3);
//        std::vector<int> listIndexVertexTarget;
//
//
//        mit2 =_mapbasins.begin();
//        mend2 = _mapbasins.end();
//
//        for (; mit2 != mend2; ++mit2)
//        {
//          listIndexVertexTarget.clear();
//
//          //listVertexbasin2 = (mit2->second).begin();
//          //for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//
//          if ((mit2->second).size() > 150)
//          {
//            listVertexbasin2 = (mit2->second).begin();
//            for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//            {
//              //listVertexbasin2++;
//              if (_texProbaPathNormLatThresh[0].item(*listVertexbasin2)==1)
//                listIndexVertexTarget.push_back(*listVertexbasin2);
//            }
//
//            vector<int> pathIndex;
//            int s,d;
//
//            if (listIndexVertexTarget.size() > 10)
//              {
//              pathIndex = sp.longestPath_ind(listIndexVertexTarget, &s, &d);
//
//              for (int i = 0; i < pathIndex.size(); i++)
//                _texProbaPathNormLat7[0].item(pathIndex[i]) = 45;
//              }
//          }
//        }
//
//        string probpathnormLatBest = "texProbaLat_embc_best_depth.tex";
//        Writer<TimeTexture<float> > texWProbaPathNormLatBest(probpathnormLatBest);
//        texWProbaPathNormLatBest << _texProbaPathNormLat7;
//
//
//
//
//
//
//
//
//
//
//
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      _texProbaPath[0].item(i) = 0;
//    }
//
//    cout << "re compute label of basins with lon inter proj: \n";
//    j = 200;
//    _mapbasins.clear();
//    _texbasinsLon = TimeTexture<short>(1, _mesh.vertex().size());
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//       _texbasinsLon[0].item(i) = _texbasins[0].item(i) - 200;
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//      {
//        if (texExtremiteLon[0].item(i) > 0)
//          _texbasins[0].item(i) = texExtremiteLon[0].item(i);
//        else
//          _texbasins[0].item(i) = 0;
//      }
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      if (_texbasins[0].item(i) > 0 && _texbasins[0].item(i) < 200)
//      {
//      floodFillIter(i,j++,_texbasins[0].item(i));
//      }
//    }
//
//    _texbasinsLon = TimeTexture<short>(1, _mesh.vertex().size());
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//    if (_texbasins[0].item(i)!=0)
//      _texbasinsLon[0].item(i) = _texbasins[0].item(i) - 200;
//    _texbasins[0].item(i) = texbasinsSave[0].item(i);
//    }
//
//    cout << "done\n";
//
//    string toto9 = "texbasinsInterLonCachiaEtiquette.tex";
//    Writer<TimeTexture<short> > texW9(toto9);
//    texW9 << _texbasinsLon;
//
//    cout << endl;
//
//    TimeTexture<short> texContourLonDilInter(1, _mesh.vertex().size() );
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      value = texExtremiteLon[0].item(i);
//
//      //on marque les points de contour
//      if (value > 0)
//      {
//        set<uint> voisins = _neigh[i];
//        set<uint>::iterator voisIt = voisins.begin();
//        nb_voisins = 0;
//
//        //on parcourt tous les voisins du sommet
//        for (; voisIt != voisins.end(); voisIt++)
//        {
//          if (texExtremiteLon[0].item(*voisIt) != value)
//            {
//            nb_voisins++;
//            texContourLonDilInter.item(i) = value;
//            continue;
//            }
//        }
//      }
//    }
//
//     for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//      if (texContourLonDilInter[0].item(i)!=0)
//      {
//       _listIndexLon.insert(i);
//      }
//    }
//
//    string lonExt = "texLonExtInter.tex";
//    Writer<TimeTexture<short> > texWLonExt(lonExt);
//    texWLonExt << texExtremiteLon;
//
//    string lonExt2 = "texLonContourExtInter.tex";
//    Writer<TimeTexture<short> > texWLonExt2(lonExt2);
//    texWLonExt2 << texContourLonDilInter;
//
//    cout << "done\n";
//
//
//
//    mit = _mapbasins.begin();
//    mend = _mapbasins.end();
////
//    cout << "sort constraints lon by basins : ";
////
//    int lon;
//    set<int>::iterator itlon;
//
//    int nbbasinsConstaintLon = 0;
////
//    //pour chaque basin i
//    for (; mit != mend; ++mit)
//    {
//      _listIndexTemp.clear();
//
//      //cout << (int)mit->first << ": ";
//      listVertexbasin = (mit->second).begin();
//
//      //on parcourt la liste des contraintes lat
//      if (!_listIndexLon.empty())
//      for (itlon=_listIndexLon.begin(); itlon!=_listIndexLon.end(); itlon++)
//      {
//        it=(mit->second).find(*itlon);
//
//        value = texContourLonDilInter[0].item(*itlon);
//        //si une lat est dans le basin alors je l'enlève de la liste
//        //et je l'ajoute à la liste temporaire des points du basin i
//        if (it != (mit->second).end())
//        {
//          if (value != 0)
//          {
//              _listIndexTemp.insert(*it);
//          }
//
//        }
//      }
//    // on associe la liste temporaire à la map du basin
//      if (!_listIndexTemp.empty())
//      {
//        //temp = make_pair (_rootsLon.item(*it),_listIndexTemp);
//        _mapConstraintLon.insert (pair<int, set<int> >(nbbasinsConstaintLon++, _listIndexTemp));
//      }
//
//    }
//
//
//    cout << "nb basin Lon= " << nbbasinsConstaintLon<< endl;
//
//    map<int, set<int> >::const_iterator mclonit(_mapConstraintLon.begin()),mclonend(_mapConstraintLon.end());
//
//    //textures contenant les contraintes lon
//    TimeTexture<float> texOutLon(1, _mesh.vertex().size() );
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//     texOutLon[0].item(i) = 0.0;
//    }
//
//    //pour chaque basin contenant des contraintes de longitude
//
//    for (; mclonit != mclonend; ++mclonit)
//    {
//     _vectorIndexTemp.clear();
//
//     cout << "basin " << (int)mclonit->first << " : ";
//     itp1 = (mclonit->second).begin();
//
//     //ARN DEBUG
//    // if (mclonit->first == 100000)
//     {
//       //myHistoLat << "basin " << (int)mclatit->first << "\n";
//       //on parcourt la liste des contraintes lat
//       for (; itp1!=(mclonit->second).end(); itp1++)
//       {
//         //cout << *itp1 << "#" << _rootsLat.item(*itp1) << " ";
//         _vectorIndexTemp.push_back (*itp1);
//       }
//
//       //cout << _vectorIndexTemp.size() << " - " ;
//       int v1,v2;
//
//       vector<int>::iterator itv1;
//
//       int i;
//       vector<int> _vectorIndexTempConstraint;
//       if (_vectorIndexTemp.size() > 1)
//       {
//         while (!_vectorIndexTemp.empty())
//         {
//           itv1 = _vectorIndexTemp.begin();
//           //v1 = _rootsLat.item(*itv1);
//           v1 = texContourLonDilInter.item(*itv1);
//
//           _vectorIndexTempConstraint.clear();
//
//           for (i=0; i<_vectorIndexTemp.size(); i++)
//           {
//             //cout << _vectorIndexTemp[i] << "#" << _rootsLat.item(_vectorIndexTemp[i]) << " ";
//             //v2 = _rootsLat.item(_vectorIndexTemp[i]);
//             v2 = texContourLonDilInter.item(_vectorIndexTemp[i]);
//             if (v1 == v2 && (v1==25))
//             {
//               _vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
//               //cout << _vectorIndexTemp[i] << " ";
//               _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
//               i--;
//             }
//            else
//              if (v1 == v2)
//              {
//              //_vectorIndexTempConstraint.push_back (_vectorIndexTemp[i]);
//              //cout << _vectorIndexTemp[i] << " ";
//              _vectorIndexTemp.erase (_vectorIndexTemp.begin()+i);
//              i--;
//              }
//           }
//
//           int source,target;
//           cout << _vectorIndexTempConstraint.size() << endl;
//           if (_vectorIndexTempConstraint.size()>=30)
//           {
//             vector<int> listIndexVertexPathSP;
//
//             //listIndexVertexPathSP = maxGeodesicDistance (_vectorIndexTempConstraint,v1,&source,&target);
//             listIndexVertexPathSP = maxGeodesicDistanceDepthStable (_vectorIndexTempConstraint,v1,&source,&target);
//
//             cout << " longest shortestpath (" << source << "," << target << ") --> value = " << v1 << endl;
//
//             //myHistoLat << "\n" << source << "," << target << ") --> value = " << v1 << endl;
//             //saveHistoTemp (source,target);
//
//             //listIndexVertexPathSP = computeShortestPathSulci(source,target);
//
//             for (unsigned t = 0; t < listIndexVertexPathSP.size(); t++)
//               {
//               //cout << listIndexVertexPathSP[t] << " ";
//               texOutLon[0].item(listIndexVertexPathSP[t]) = (float) v1;
//               }
//           }
//         }
//         //cout << endl;
//       }
//     }
//     //ARN DEBUG
//    }
//
//    if (_adrLonGeodesicOut!="")
//    {
//     Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
//     texWLon << texOutLon;
//    }
//
//    string probpathLon = "texProbaLon.tex";
//    Writer<TimeTexture<short> > texWProbaLon(probpathLon);
//    texWProbaLon << _texProbaPath;
//
//
//    //on normalise la texture des probas des lon
//    TimeTexture<float> _texProbaPathNormLon(1, _mesh.vertex().size() );
//
//    mit2 =_mapbasins.begin();
//    mend2 = _mapbasins.end();
//
//    n_b2 = 0;
//
//    for (; mit2 != mend2; ++mit2)
//    {
//      max_depth2 = 0.0;
//      //myHistoLat << "basin lon" << n_b2++ << " " << endl;
//      //cout << "basin " << n_b2-1 << " " << (mit2->second).size() << endl;;
//      listVertexbasin2 = (mit2->second).begin();
//      for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//      {
//      //cout << *listVertexbasin << " " ;
//      val2 = _texProbaPath[0].item(*listVertexbasin2);
//      if ( val2 > max_depth2)
//        max_depth2 = val2;
//      }
//
//      //cout << max_depth2 << endl;
//
//      if ((mit2->second).size() > 150)
//      {
//        listVertexbasin2 = (mit2->second).begin();
//        for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//        {
//        if (max_depth2!=0)
//          _texProbaPathNormLon[0].item(*listVertexbasin2) =  (float)(_texProbaPath[0].item(*listVertexbasin2))/max_depth2;
//        else
//          _texProbaPathNormLon[0].item(*listVertexbasin2) = 0 ;
//
//        myHistoLat << _texProbaPathNormLon[0].item(*listVertexbasin2) << "\t" << *listVertexbasin2 << "\n";
//        }
//      }
//    }
//
//    string probpathnormLon = "texProbaLon_embc.tex";
//    Writer<TimeTexture<float> > texWProbaPathNormLon(probpathnormLon);
//    texWProbaPathNormLon << _texProbaPathNormLon;
//
//    TimeTexture<float> _texProbaPathNormLonThresh(1, _mesh.vertex().size() );
//
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//       if (_texProbaPathNormLon[0].item(i) > 0.3)
//         _texProbaPathNormLonThresh[0].item(i) = 1;
//       else
//         _texProbaPathNormLonThresh[0].item(i) = 0;
//    }
//
//
//    TimeTexture<float> _texProbaPathNormLon7(1, _mesh.vertex().size() );
//
//    //on cherche le plus long chemin
//
//    mit2 =_mapbasins.begin();
//    mend2 = _mapbasins.end();
//
//    for (; mit2 != mend2; ++mit2)
//    {
//      listIndexVertexTarget.clear();
//
//      //listVertexbasin2 = (mit2->second).begin();
//      //for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//
//      if ((mit2->second).size() > 150)
//      {
//        listVertexbasin2 = (mit2->second).begin();
//        for (; listVertexbasin2!=(mit2->second).end(); listVertexbasin2++)
//        {
//          //listVertexbasin2++;
//          if (_texProbaPathNormLonThresh[0].item(*listVertexbasin2)==1)
//            listIndexVertexTarget.push_back(*listVertexbasin2);
//        }
//
//        vector<int> pathIndex;
//        int s,d;
//
//        if (listIndexVertexTarget.size() > 20)
//          {
//          pathIndex = sp.longestPath_ind(listIndexVertexTarget, &s, &d);
//
//          for (int i = 0; i < pathIndex.size(); i++)
//            _texProbaPathNormLon7[0].item(pathIndex[i]) = 25;
//          }
//      }
//    }
//
//    string probpathnormLonBest = "texProbaLon_embc_best_depth.tex";
//    Writer<TimeTexture<float> > texWProbaPathNormLonBest(probpathnormLonBest);
//    texWProbaPathNormLonBest << _texProbaPathNormLon7;
//
//
//
//
//
//
//    if (_adrLonGeodesicOut!="")
//    {
//      Writer<TimeTexture<float> > texWLon(_adrLonGeodesicOut);
//      texWLon << texOutLon;
//    }
//
//    //on fusionne les textures lat et lon
//    for (uint i = 0; i < _mesh.vertex().size(); i++)
//    {
//       if (texOutLon[0].item(i) != 0)
//         texOutLat[0].item(i) = texOutLon[0].item(i);
//    }
//
//    Writer<TimeTexture<float> > texWLines(_adrLines);
//    texWLines << texOutLat;
//
//    myHistoLat.close();
//
//}

//void SulcalLinesGeodesic::computeGraphDijkstra (AimsSurfaceTriangle surface, int constraintType,int strain)
//{
//  // compute and copy curvature
//
//  cout << _texbasinsDepthNorm[0].nItem() << endl;
//
////  if (_adrbasinsDepthNorm=" "}
////    vector<float> &curv = _texCurv[0].data();
////  else
////
//  vector<float> &curv = _texbasinsDepthNorm[0].data();
//  //vector<float> &curv = _texCurv[0].data();
//
//  //val = (float)fabs(p2 - p1)/(SPath[SPath.size()-2 - i].distance(&SPath[SP
//  // copy vertex and faces vector
//  vector<Point3df> & vert = surface.vertex();
//  vector<AimsVector<uint, 3> > & tri = surface.polygon();
//  _pointsSP.resize(3*vert.size());
//  _facesSP.resize(3*tri.size());
//
//  for (uint j = 0; j < (int) vert.size(); j++)
//  {
//    _pointsSP[3*j] = vert[j][0];
//    _pointsSP[3*j+1] = vert[j][1];
//    _pointsSP[3*j+2] = vert[j][2];
//  }
//  for (uint j = 0; j < (int) tri.size(); j++)
//  {
//    _facesSP[3*j] = tri[j][0];
//    _facesSP[3*j+1] = tri[j][1];
//    _facesSP[3*j+2] = tri[j][2];
//  }
//
//  // compute adjacence graph
//
//  cout << "compute adjacence graph : ";
//
//  _meshSPc.initialize_mesh_data(_pointsSP,_facesSP, &curv[0],constraintType,strain);
//  //_meshSP.initialize_mesh_data(_pointsSP,_facesSP, NULL ,0,0);
//  dijkstra_algorithm = new geodesic::GeodesicAlgorithmDijkstra(&_meshSPc);
//
//  cout << "done" << endl;
//
//
//}
//
//double SulcalLinesGeodesic::computeDepthShortestPathSulci(unsigned source, unsigned target, vector<geodesic::SurfacePoint> & SPath, vector<int> &listIndexVertexPathSP)
//{
//  vector<geodesic::SurfacePoint> sources;
//  sources.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[source]));
//
//  vector<geodesic::SurfacePoint> targets;
//  targets.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[target]));
//
//  //printf("indice source = %d target = %d \n",source, target);
//
//  // clear path
//  //vector<geodesic::SurfacePoint> SPath;
//  SPath.clear();
//
//  // dijkstra method
//  //geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;
//
//  listIndexVertexPathSP.clear();
//
//  geodesic::SurfacePoint short_sources(&_meshSPc.vertices()[source]);
//  geodesic::SurfacePoint short_targets(&_meshSPc.vertices()[target]);
//
//  dijkstra_algorithm->geodesic(short_sources,short_targets, SPath, listIndexVertexPathSP);
//
//  //std::vector<int>::iterator ite;
//  reverse(listIndexVertexPathSP.begin(),listIndexVertexPathSP.end());
//  listIndexVertexPathSP.push_back((int)target);
//
////  cout << "shortest path (index vertex) = ";
////  cout << endl;
//
//  double length = 0;
//  float p1,p2,val,val2;
////    for (unsigned i = 0; i < listIndexVertexPathSP.size() - 1; i++)
////      cout << listIndexVertexPathSP[i+1] << " " ;
//
//  if(!SPath.empty())
//  {
//    for(unsigned i=0; i<SPath.size()-1; ++i)
//    {
//      p1 = _texbasinsDepthNorm[0].item(listIndexVertexPathSP[i]);
//      p2 = _texbasinsDepthNorm[0].item(listIndexVertexPathSP[i+1]);
//
////      val = (float)fabs(p2 - p1)/(SPath[SPath.size()-2 - i].distance(&SPath[SPath.size()-i-1]));
////      val2 = (float)(SPath[SPath.size()-2 - i].distance(&SPath[SPath.size()-i-1]))* (1 - 5*sqrt(val));
////
////      length += val2;
//
//      length += SPath[i].distance(&SPath[i+1]);
//
//    }
//  }
//
//  for(unsigned i=0; i<SPath.size(); ++i)
//   {
//    _texProbaPath[0].item(listIndexVertexPathSP[i])++;
//   }
//
//  return length;
//
//}


//
//double SulcalLinesGeodesic::saveHistoTemp(unsigned source, unsigned target)
//{
//  vector<int> listIndexVertexPathSP;
//  vector<geodesic::SurfacePoint> SPath;
//
//  SPath.clear();
//  vector<geodesic::SurfacePoint> sources;
//  sources.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[source]));
//
//  vector<geodesic::SurfacePoint> targets;
//  targets.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[target]));
//
//  //printf("indice source = %d target = %d \n",source, target);
//
//  // clear path
//  //vector<geodesic::SurfacePoint> SPath;
//  SPath.clear();
//
//  // dijkstra method
//  //geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;
//
//  listIndexVertexPathSP.clear();
//
//  geodesic::SurfacePoint short_sources(&_meshSPc.vertices()[source]);
//  geodesic::SurfacePoint short_targets(&_meshSPc.vertices()[target]);
//
//  dijkstra_algorithm->geodesic(short_sources,short_targets, SPath, listIndexVertexPathSP);
//
//  //std::vector<int>::iterator ite;
//  reverse(listIndexVertexPathSP.begin(),listIndexVertexPathSP.end());
//  listIndexVertexPathSP.push_back((int)target);
//
////  cout << "shortest path (index vertex) = ";
////  cout << endl;
//
//  double length = 0;
//  float p1,p2,val;
//
//  myHistoLat << SPath.size() <<"\n";
//  myHistoLat << "\nhisto\n";
//
//  if(!SPath.empty())
//  {
//    for(unsigned i=0; i<SPath.size(); ++i)
//    {
//      p1 = _texbasinsDepthNorm[0].item(listIndexVertexPathSP[i]);
//      myHistoLat << 1-p1 << "\t";
//    }
//  }
//
//  myHistoLat << "\ngradient\n";
//  if(!SPath.empty())
//  {
//    for(unsigned i=0; i<SPath.size()-1; ++i)
//    {
//      p1 = _texbasinsDepthNorm[0].item(listIndexVertexPathSP[i]);
//      p2 = _texbasinsDepthNorm[0].item(listIndexVertexPathSP[i+1]);
//
//      //val = (float)fabs(p2 - p1)/(SPath[SPath.size()-2 - i].distance(&SPath[SPath.size()-i-1]));
//      val = (float)fabs(p2 - p1);
//
//      myHistoLat << val << "\t";
//    }
//  }
//
//  myHistoLat << "\npath\n";
//  if(!SPath.empty())
//  {
//    for(unsigned i=0; i<SPath.size()-1; ++i)
//    {
//      myHistoLat << SPath[i].distance(&SPath[i+1]) << "\t";
//    }
//  }
//}

//double SulcalLinesGeodesic::computeShortestPathSulci(unsigned source, unsigned target, vector<geodesic::SurfacePoint> & SPath, vector<int> &listIndexVertexPathSP)
//{
//  //vector<int> listIndexVertexPathSP;
//
//  // compute shortest path
//  //cout << "compute shortest path : ";
//
//  vector<geodesic::SurfacePoint> sources;
//  sources.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[source]));
//
//  vector<geodesic::SurfacePoint> targets;
//  targets.push_back(geodesic::SurfacePoint(&_meshSPc.vertices()[target]));
//
////  printf("indice source = %d target = %d \n",source, target);
//
//  // clear path
//  //vector<geodesic::SurfacePoint> SPath;
//  SPath.clear();
//
//  // dijkstra method
//  //geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;
//
//  listIndexVertexPathSP.clear();
//
//  geodesic::SurfacePoint short_sources(&_meshSPc.vertices()[source]);
//  geodesic::SurfacePoint short_targets(&_meshSPc.vertices()[target]);
//
//  dijkstra_algorithm->geodesic(short_sources,short_targets, SPath, listIndexVertexPathSP);
//
//  //std::vector<int>::iterator ite;
//  reverse(listIndexVertexPathSP.begin(),listIndexVertexPathSP.end());
//  listIndexVertexPathSP.push_back((int)target);
//
////  cout << "shortest path (index vertex) = ";
////  for (unsigned i = 0; i < listIndexVertexPathSP.size(); i++)
////    cout << listIndexVertexPathSP[i] << " " ;
////  cout << endl;
//
//  return dijkstra_algorithm->length(SPath);
//  //return listIndexVertexPathSP;
//
//}
//
//vector<int> SulcalLinesGeodesic::maxGeodesicDistance(vector<int> points, int constraint, int* s, int *d)
//{
//  int i,j;
//
//  vector<int> listIndexVertexPathSP;
//  vector<int> listIndexVertexPathSP_result;
//  //cout << source << " " << dest << endl; ;
//  // clear path
//  vector<geodesic::SurfacePoint> SPath;
//  SPath.clear();
//
//  int index_max_i,index_max_j;
//  double dist,dist_max;
//
//  index_max_i = -1;
//  index_max_i = -1;
//
//  std::cout << endl;
//  //dist_max_j = 0;
//
//  int nb_combinaison_max = (points.size() * (points.size()-1))/2;
//  int nb_combinaison = 0;
//
//  dist_max = -100000;
//
//  for (j=0; j<points.size(); j++)
//  {
//    for(i=j+1; i<points.size(); i++)
//    {
//      nb_combinaison++;
//      dist = computeShortestPathSulci( points[j], points[i],SPath,listIndexVertexPathSP);
//
//      //cout << "(" << points[j] << "-" << points[i] << " " << dist << ") " << endl;
//      //cout << " (" << i << "-" << j << ") " ;
//      std::cout << "\r" << nb_combinaison << "/" << nb_combinaison_max << std::flush;
//
//      if (dist >= dist_max)
//      {
//        dist_max = dist;
//        index_max_i = i;
//        index_max_j = j;
//        listIndexVertexPathSP_result = listIndexVertexPathSP;
//      }
//    }
//  }
//
//  *s = points[index_max_i];
//  *d = points[index_max_j];
//
//  return listIndexVertexPathSP_result;
//
//}
//
//vector<int> SulcalLinesGeodesic::maxGeodesicDistanceDepthStable(vector<int> points, int constraint, int* s, int *d)
//{
//  int i,j;
//
//  vector<int> listIndexVertexPathSP;
//  vector<int> listIndexVertexPathSP_result;
//  //cout << source << " " << dest << endl; ;
//  // clear path
//  vector<geodesic::SurfacePoint> SPath;
//  SPath.clear();
//
//  int index_max_i,index_max_j;
//  double dist,dist_max;
//
//  index_max_i = -1;
//  index_max_i = -1;
//
//  std::cout << endl;
//  //dist_max_j = 0;
//
//  int nb_combinaison_max = (points.size() * (points.size()-1))/2;
//  int nb_combinaison = 0;
//
//  dist_max = -100000;
//
//  for (j=0; j<points.size(); j++)
//  {
//
//    for(i=j+1; i<points.size(); i++)
//    {
//      //myHistoLat <<  "\t" << points[j] << "\t" << points[i] <<"\t";
//      nb_combinaison++;
//      //dist = computeShortestPathSulci( points[j], points[i],SPath,listIndexVertexPathSP);
//
//      dist = computeDepthShortestPathSulci( points[j], points[i],SPath,listIndexVertexPathSP);
//      //cout << "(" << points[j] << "-" << points[i] << " " << dist << ") ";
//      //cout << " (" << i << "-" << j << ") " ;
//
//      //myHistoLat << dist << "\t" << points[j] << "\t" << points[i] << "\t";
//
////      myHistoLat << "\n" << dist << " (" << points[j] << "," << points[i] << ")\n";
//      //saveHistoTemp (points[j],points[i]);
//
//      std::cout << "\r" << i << " " << j << " "<< nb_combinaison << "/" << nb_combinaison_max << std::flush;
//
//      if (dist >= dist_max)
//      {
//        dist_max = dist;
//        index_max_i = i;
//        index_max_j = j;
//        listIndexVertexPathSP_result = listIndexVertexPathSP;
//      }
//    }
//  }
//
//  *s = points[index_max_i];
//  *d = points[index_max_j];
//
//  //myHistoLat << "\ndist = " << dist_max << "\n";
//  return listIndexVertexPathSP_result;
//}

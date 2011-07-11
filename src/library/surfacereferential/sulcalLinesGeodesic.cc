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
    string & adrRootsLon, string & adrRootsLat, int extremeties_method, int constraint_type, int strain, vector<float> proba, string saveFolder, float curv_thresh) :
    _adrMesh(adrMesh),_adrCurv(adrCurv),_adrGeodesicDepth(adrGeodesicDepth),
    _adrRootsLon(adrRootsLon),_adrRootsLat(adrRootsLat),_constraint_type(constraint_type),_extremeties_method(extremeties_method), _strain(strain),
    _proba(proba), _adrSaveFolder(saveFolder), _curv_thresh(curv_thresh)
{
  std::cout << "Read mesh : ";
  Reader < AimsSurfaceTriangle > r(adrMesh);
  r.read( _mesh );
  std::cout << "done" << std::endl;

  if (_adrSaveFolder!="")
    _save = true;
  else
    _save = false;

  if (adrRootsLon!="")
  {
  std::cout << "Read roots texture : ";
  Reader < TimeTexture<short> > rtLon(adrRootsLon);
  rtLon.read( _rootsLon );
  }

  if (adrRootsLat!="")
  {
  Reader < TimeTexture<short> > rtLat(adrRootsLat);
  rtLat.read( _rootsLat );
  std::cout << "done" << std::endl;
  }

  if (adrCurv=="")
  {
    cout << "compute and save curvature texture : ";
    _texCurv = TimeTexture<float>(1, _mesh.vertex().size());
    _texCurv = AimsMeshCurvature(_mesh[0]);

    if (_save)
    {
      writeFloatTexture("curv_new.tex",_texCurv);
      cout << "done" << endl;
    }
  }
  else
  {
    std::cout << "Read curvature texture : ";
    Reader < TimeTexture<float> > rtCurv(adrCurv);
    rtCurv.read( _texCurv );
    cout << "done" << endl;
  }

  if (adrGeodesicDepth!="")
  {
    std::cout << "Read Geodesic Depth Texture : ";
    Reader < TimeTexture<float> > rtDG(adrGeodesicDepth);
    rtDG.read( _geoDepth );
    std::cout << "done" << std::endl;
  }
  else
    constraint_type = 1;

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

  std::cout << "Sulcal basins segmentation " << endl;

  TimeTexture<short> texBasins(1, _mesh.vertex().size());

  map<int,set<int> > mapBasins;

  segmentationSulcalBasins (_texCurv, texBasins, mapBasins,_curv_thresh,1,2);

  if (_save)
  {
    std::cout << "Save connected components texture : ";
    writeShortTexture("basins.tex",texBasins);
  }

  if (_adrGeodesicDepth!="")
  {
    cout << "Normalize Geodesic Depth Texture with sulcal basins " << endl;
    _geoDepthNorm = TimeTexture<float>(1, _mesh.vertex().size());
    normalizeDepthMap (_geoDepth, _geoDepthNorm, mapBasins);

    if (_save)
    {
      cout << "Save Normalize Geodesic Depth texture : ";
      writeFloatTexture("depth_norm.tex",_geoDepthNorm);
    }
  }

  switch (_extremeties_method)
  {
  case 1 :
    cout << "extraction of extremities method : projection crop by basins" << endl;
    sulcalLinesExtract_projection(mapBasins,texBasins);
    break;
  case 2 :
    cout << "extraction of extremities method : map of probability" << endl ;
    sulcalLinesExtract_probability(mapBasins,texBasins);
    break;
  }

}
void SulcalLinesGeodesic::probaMap()
{
  std::cout << "START : ";

  std::cout << "reading sulcal basins in roots lat/lon texture " << endl;

  map<int,set<int> > mapBasins;

  TimeTexture<float> texBasinsF;

  if (_adrCurv!="")
  {
  std::cout << "Read basins texture : ";
  Reader < TimeTexture<float> > rtLon(_adrCurv);
  rtLon.read( texBasinsF );
  }

  TimeTexture<short> texBasins(1, _mesh.vertex().size());

  texBinarizeF2S(texBasinsF, texBasins, 0 ,0 ,1);
  texConnectedComponent(texBasins, mapBasins,1000);

  if (_adrGeodesicDepth!="")
  {
    cout << "Normalize Geodesic Depth Texture with sulcal basins " << endl;
    _geoDepthNorm = TimeTexture<float>(1, _mesh.vertex().size());
    normalizeDepthMap (_geoDepth, _geoDepthNorm, mapBasins);

    if (_save)
    {
      cout << "Save Normalize Geodesic Depth texture : ";
      writeFloatTexture("depth_norm_new.tex",_geoDepthNorm);
    }
  }

  cout << mapBasins.size() << " Basins Lat found" << endl;

  TimeTexture<short> texContourBasins(1, _mesh.vertex().size() );
  map<int,set<int> > mapContourBasins;

  contourBasins(mapBasins,texBasins,mapContourBasins, texContourBasins);

  if (_save)
  {
    writeShortTexture("contour_basins_new.tex",texContourBasins);
  }

  cout << "compute probability map" << endl;
  TimeTexture<float> texProba(1, _mesh.vertex().size() );
  TimeTexture<float> texProbaNorm(1, _mesh.vertex().size() );
  computeProbabiltyMap(mapContourBasins,texContourBasins,texProba);
  normalizeProbabiltyMap(mapBasins,mapContourBasins,texContourBasins,texProba,texProbaNorm);
  if (_save)
  {
    writeFloatTexture("proba_new.tex",texProba);
    writeFloatTexture("proba_norm_new.tex",texProbaNorm);
  }

  cout << "done " << endl;

}

void SulcalLinesGeodesic::writeShortTexture (string name,TimeTexture<short> &out)
{
//  size_t found;
//  found = _adrRootsLon.find_last_of(".");
//  string adrBasins = _adrRootsLon.substr(0,found-9) + name;
//  Writer<TimeTexture<short> > texW(adrBasins);
//  texW << out;
//  cout << "write " << adrBasins << " done" << endl;
//

  string adrBasins = _adrSaveFolder + name;
  Writer<TimeTexture<short> > texW(adrBasins);
  texW << out;
  cout << "write " << adrBasins << " done" << endl;
}

void SulcalLinesGeodesic::writeFloatTexture (string name,TimeTexture<float> &out)
{
//  size_t found;
//  found = _adrRootsLon.find_last_of(".");
//  string adrBasins = _adrRootsLon.substr(0,found-9) + name;
//  Writer<TimeTexture<float> > texW(adrBasins);
//  texW << out;
//  cout << "write " << adrBasins << " done" << endl;
//
  string adrBasins = _adrSaveFolder + name;
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

      if ( texBasinsTemp[0].item(indexCurr) <= (oldTextureValue + 0.001) &&
          texBasinsTemp[0].item(indexCurr) >= (oldTextureValue - 0.001) )
      //if ( texBasinsTemp[0].item(indexCurr) == oldTextureValue)
      {
        listIndexVertexFill.insert(indexCurr);
        stack.push(indexCurr);
        //cout << indexCurr << " " <<  newTextureValue << endl;
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
    if (texIn[0].item(i) <= threshold)
      texOut[0].item(i) = inf;
    else
      texOut[0].item(i) = sup;
  }
}

void SulcalLinesGeodesic::texBinarizeS2S(TimeTexture<short> &texIn, TimeTexture<short> &texOut, int threshold,int inf,int sup)
{
  for (uint i = 0; i < texIn[0].nItem(); i++)
  {
    if (texIn[0].item(i) <= threshold)
      texOut[0].item(i) = inf;
    else
      texOut[0].item(i) = sup;
  }
}

TimeTexture<short> SulcalLinesGeodesic::texConnectedComponent(TimeTexture<short> &texBasins, map<int,set<int> > &mapBasins, int offset)
{
  int j = 1;
  int val;
  mapBasins.clear();

  TimeTexture<short> texTemp(1,texBasins[0].nItem());

  for (uint i = 0; i < texBasins[0].nItem(); i++)
    texTemp[0].item(i) = texBasins[0].item(i) + offset;

  for (uint i = 0; i < texBasins[0].nItem(); i++)
  {
    if (texTemp[0].item(i) > offset  )
      floodFillIter(i,j++,texTemp[0].item(i),texTemp,mapBasins);
  }

  return texTemp;
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
  vector<unsigned> pathTemp,indexTemp;

  for (uint i = 0; i < _mesh.vertex().size(); i++)
    out[0].item(i) = 0.0;

  TimeTexture<float> texConstraint(1, _mesh.vertex().size() );
  int method;

  // contrainte sur la courbure
  if (_constraint_type == 1)
  {
    texConstraint = _texCurv;
    method = 1; // valeur pour les sulci
  }

  // contrainte sur la profondeur normalisée (par les bassins décrits sur les infos de courbures)
  if (_constraint_type == 2)
  {
    texConstraint = _geoDepthNorm;
    method = 2; // valeur pour les gyri (car la depth map est inversé : les maxima sont dans les vallées)
  }

  GeodesicPath sp(_mesh,texConstraint,method,_strain);

  int source,target;
  int constraintValue;

  for (; mit != mend; mit++)
  {
    indexTemp.clear();

    cout << "\nbasin " << (int)mit->first << " : \n";
    it = (mit->second).begin();

    //on copie la liste des index de sommets set --> vector (non nul)
    for (; it!=(mit->second).end(); it++)
    {
      constraintValue = roots.item( *it );
      if (constraintValue != 0)
        indexTemp.push_back (*it);
    }

    //si il y a plus d'un point dans le bassin
    if (indexTemp.size()>1)
    {
      //on récupère la valeur de la contrainte
      constraintValue = roots.item( *(indexTemp.begin()) );

      cout << "value " << constraintValue << endl;
      //on calcule le plus long chemin
      pathTemp.clear();
      double len;

      pathTemp = sp.longestPath_N_N_ind(indexTemp, &source, &target,&len, 1);

      for (int i = 0; i < pathTemp.size(); i++)
        out[0].item(pathTemp[i]) = constraintValue;
    }
    else
      cout << "not enough points ! " << endl;
  }
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

void SulcalLinesGeodesic::segmentationSulcalBasins (TimeTexture<float> &texIn,TimeTexture<short> &texBasins,map<int,set<int> > &mapBasins,float threshold, int close, int open)
{
//binarisation texture de courbure < 0.0
texBinarizeF2S(texIn, texBasins, threshold , 1, 0);

TimeTexture<short> texBasinsDil(1, _mesh.vertex().size() );
TimeTexture<short> texBasinsErode(1, _mesh.vertex().size() );

//fermeture des bassins
texBasinsDil[0]=MeshDilation<short>( _mesh[0], texBasins[0], short(0), -1, close, true);
texBasinsErode[0]=MeshErosion<short>( _mesh[0], texBasinsDil[0], short(0), -1, close , true);

//ouverture des bassins
texBasinsDil[0]=MeshErosion<short>( _mesh[0], texBasinsErode[0], short(0), -1, open , true);
texBasinsErode[0]=MeshDilation<short>( _mesh[0], texBasinsDil[0], short(0), -1, open, true);

texBinarizeS2S(texBasinsErode, texBasins, 0 ,0 ,1);

// étiquetage des composantes connexes
texConnectedComponent(texBasins, mapBasins,1000);
cout << mapBasins.size() << " Basins found" << endl;
}


void SulcalLinesGeodesic::sulcalLinesExtract_projection(map<int,set<int> > &mapBasins, TimeTexture<short> &texBasins)
{
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

void SulcalLinesGeodesic::dilationRoots(TimeTexture<short> &texLatDil,TimeTexture<short> &texLonDil,int size)
{
  texLatDil[0]=MeshDilation<short>( _mesh[0], _rootsLat[0], short(0), -1, size, false);
  texLonDil[0]=MeshDilation<short>( _mesh[0], _rootsLon[0], short(0), -1, size, false);
}

void SulcalLinesGeodesic::interRootsDilBasins(TimeTexture<short> &texBasins,TimeTexture<short> &texDil,TimeTexture<short> &texInter)
{
  for (uint i = 0; i < _mesh.vertex().size(); i++)
    if (texBasins[0].item(i)!=0 && texDil[0].item(i) > 0)
      texInter[0].item(i) = texDil[0].item(i);
}

void SulcalLinesGeodesic::cleanBasins(map<int,set<int> > &mapBasins,TimeTexture<short> &texBasins,int nbPoint)
{
  map<int,set<int> >::iterator it;
  set<int>::iterator its;

  for ( it=mapBasins.begin() ; it != mapBasins.end(); it++ )
  {
  its = ((*it).second).begin();

  if ( ((*it).second).size() < nbPoint)
    {
    //on efface tous les vertex de la texture texBasins groupés dans des bassins qui contiennent moins de nbPoint
    for ( ; its != ((*it).second).end(); its++ )
      texBasins[0].item(*its) = 0;

    // on l'enlève ensuite de la map des bassins
    mapBasins.erase ((*it).first);
    }
  }
}

void SulcalLinesGeodesic::contourBasins(map<int,set<int> > &mapBasins,TimeTexture<short> &texBasins,map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins)
{
  map<int,set<int> >::iterator it;
  set<int>::iterator its;
  set<int> listIndexTemp;
  int value;
  int nb = 0;

  mapContourBasins.clear();

  //pour tous les bassins
  for ( it=mapBasins.begin() ; it != mapBasins.end(); it++ )
  {
    listIndexTemp.clear();

    //pour tous les points its du bassin it
    for (its = ((*it).second).begin() ; its != ((*it).second).end(); its++ )
    {
      //on récupère l'étiquette du bassin
      value = texBasins[0].item(*its);
      set<uint> voisins = _neigh[*its];
      set<uint>::iterator voisIt = voisins.begin();

      //on parcourt tous les voisins du sommet, si un voisin a une valeur differente alors le point est un contour
      for (; voisIt != voisins.end(); voisIt++)
      {
        if (texBasins[0].item(*voisIt) != value)
        {
          texContourBasins[0].item(*its) = value;
          listIndexTemp.insert(*its);
          continue;
        }
      }
    }
    mapContourBasins.insert (pair<int,set<int> >(nb++, listIndexTemp));
  }
}

void SulcalLinesGeodesic::computeProbabiltyMap(map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins,TimeTexture<float> &texProba)
{
  map<int,set<int> >::iterator itc;
  set<int>::iterator its;
  int value,nb_c,nb_v;
  vector<unsigned> listIndexTarget;
  vector<vector<unsigned> >indices;
  vector<vector<unsigned> >::iterator it_vv;
  vector<unsigned>::iterator it_v;

  for (uint i = 0; i < texProba[0].nItem(); i++)
    texProba[0].item(i) = 1.0;

  TimeTexture<float> texConstraint(1, _mesh.vertex().size() );
  int method;

  // contrainte sur la courbure
  if (_constraint_type == 1)
  {
    texConstraint = _texCurv;
    method = 1; // valeur pour les sulci
  }

  // contrainte sur la profondeur normalisée (par les bassins décrits sur les infos de courbures)
  if (_constraint_type == 2)
  {
    texConstraint = _geoDepthNorm;
    method = 2; // valeur pour les gyri (car la depth map est inversé : les maxima sont dans les vallées)
  }

  GeodesicPath sp(_mesh,texConstraint,method,_strain);

  //pour tous les bassins
  for ( itc=mapContourBasins.begin(); itc != mapContourBasins.end(); itc++)
  {
    listIndexTarget.clear();

    //on récupère l'étiquette du bassin (avec le premier point)
    its = ((*itc).second).begin();

    value = texContourBasins[0].item(*its);
    cout << "value = " << value << endl;

    //pour tous les points its du bassin it
    for ( ; its != ((*itc).second).end(); its++ )
      listIndexTarget.push_back(*its);

    //nombre de points de contours
    nb_c = 0;
    //nombre de sommets parcourus
    nb_v = 0;

    for (its = ((*itc).second).begin() ; its != ((*itc).second).end(); its++ )
    {
      indices.clear();

      sp.shortestPath_1_N_All_ind(*its, listIndexTarget, indices);

      nb_v += indices.size();
      cout << "\r\033[K" <<  "extremities : " << ++nb_c << "/" << listIndexTarget.size()<< " vextex : " << nb_v << flush;

      //on incrémente les valeurs de la texture proba pour chaque point ayant été parcourue
      for ( it_vv=indices.begin() ; it_vv != indices.end(); it_vv++ )
        for ( it_v=(*it_vv).begin() ; it_v != (*it_vv).end(); it_v++ )
        {
          texProba[0].item(*it_v)++;
        }
    }

    cout << endl;
  }
}

void SulcalLinesGeodesic::normalizeProbabiltyMap(map<int,set<int> > &mapBasins, map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins,TimeTexture<float> &texProba,TimeTexture<float> &texProbaNorm)
{
  map<int,set<int> >::iterator itc;
  map<int,set<int> >::iterator itb;
  set<int>::iterator its;

  int value;
  float accu,accu_max;

  for (uint i = 0; i < texProba[0].nItem(); i++)
    texProbaNorm[0].item(i) = 0.0;

  //pour tous les bassins
  for ( itc=mapContourBasins.begin(), itb=mapBasins.begin(); itc != mapContourBasins.end(), itb != mapBasins.end(); itc++ ,itb++)
  {
    accu_max = 0;

    //on récupère l'étiquette du bassin (avec le premier point)
    its = ((*itb).second).begin();

    value = texContourBasins[0].item(*its);

    //pour tous les points its du bassin it
    for ( ; its != ((*itb).second).end(); its++ )
      accu_max = max (accu_max,texProba[0].item(*its));

    //on normalise entre 0 et 1
    for (its = ((*itb).second).begin() ; its != ((*itb).second).end(); its++ )
       texProbaNorm[0].item(*its) = (float)(texProba[0].item(*its)/accu_max);

    //on met tous les points de contour du bassins à 0
    //      for (its = ((*itc).second).begin() ; its != ((*itc).second).end(); its++ )
    //        texProbaNorm[0].item(*its) = 0;


  }
}

void SulcalLinesGeodesic::textureBin2Label(TimeTexture<short> &texLabel, TimeTexture<short> &texIn, TimeTexture<short> &texOut)
{
  for (uint i = 0; i < texIn[0].nItem(); i++)
    {
    if (texIn[0].item(i) == 1)
      texOut[0].item(i) = texLabel[0].item(i);
    else
      texOut[0].item(i) = 0;
    }
}


void SulcalLinesGeodesic::sulcalLinesExtract_probability(map<int,set<int> > &mapBasins, TimeTexture<short> &texBasins)
{
  TimeTexture<short> texTemp(1, _mesh.vertex().size() );

  //on dilate les roots
  TimeTexture<short> texProjectionLatDil(1, _mesh.vertex().size() );
  TimeTexture<short> texProjectionLonDil(1, _mesh.vertex().size() );
  dilationRoots(texProjectionLatDil,texProjectionLonDil,5);

  if (_save)
  {
    writeShortTexture("lat_roots_dil.tex",texProjectionLatDil);
    writeShortTexture("lon_roots_dil.tex",texProjectionLonDil);
  }

  // on calcule les intersections avec les basins
  TimeTexture<short> texInterRootsBasinsLat(1, _mesh.vertex().size() );
  TimeTexture<short> texInterRootsBasinsLon(1, _mesh.vertex().size() );

  map<int,set<int> > mapBasinsLat;
  map<int,set<int> > mapBasinsLon;

  interRootsDilBasins(texBasins,texProjectionLatDil,texInterRootsBasinsLat);
  interRootsDilBasins(texBasins,texProjectionLonDil,texInterRootsBasinsLon);

  texConnectedComponent(texInterRootsBasinsLat, mapBasinsLat,1000);
  texConnectedComponent(texInterRootsBasinsLon, mapBasinsLon,1000);

  cout << mapBasinsLat.size() << " Basins Lat found" << endl;
  cout << mapBasinsLon.size() << " Basins Lon found" << endl;
  cleanBasins(mapBasinsLat,texInterRootsBasinsLat,15);
  cleanBasins(mapBasinsLon,texInterRootsBasinsLon,15);
  cout << "after cleaning" << endl;
  cout << mapBasinsLat.size() << " Basins Lat found" << endl;
  cout << mapBasinsLon.size() << " Basins Lon found" << endl;

  if (_save)
  {
    writeShortTexture("lat_inter_roots_dil.tex",texInterRootsBasinsLat);
    writeShortTexture("lon_inter_roots_dil.tex",texInterRootsBasinsLon);
  }

  TimeTexture<short> texContourBasinsLat(1, _mesh.vertex().size() );
  TimeTexture<short> texContourBasinsLon(1, _mesh.vertex().size() );
  map<int,set<int> > mapContourBasinsLat;
  map<int,set<int> > mapContourBasinsLon;
  contourBasins(mapBasinsLat,texInterRootsBasinsLat,mapContourBasinsLat, texContourBasinsLat);
  contourBasins(mapBasinsLon,texInterRootsBasinsLon,mapContourBasinsLon, texContourBasinsLon);
  if (_save)
  {
    writeShortTexture("lat_contour_basins.tex",texContourBasinsLat);
    writeShortTexture("lon_contour_basins.tex",texContourBasinsLon);
  }

  cout << "compute probability map" << endl;
  cout << "latitude :" << endl;
  TimeTexture<float> texProbaLat(1, _mesh.vertex().size() );
  TimeTexture<float> texProbaNormLat(1, _mesh.vertex().size() );
  computeProbabiltyMap(mapContourBasinsLat,texContourBasinsLat,texProbaLat);
  normalizeProbabiltyMap(mapBasinsLat,mapContourBasinsLat,texContourBasinsLat,texProbaLat,texProbaNormLat);
  if (_save)
  {
    writeFloatTexture("lat_proba.tex",texProbaLat);
    writeFloatTexture("lat_proba_norm.tex",texProbaNormLat);
  }
  cout << "done " << endl;
  cout << "longitude :" << endl;
  TimeTexture<float> texProbaLon(1, _mesh.vertex().size() );
  TimeTexture<float> texProbaNormLon(1, _mesh.vertex().size() );
  computeProbabiltyMap(mapContourBasinsLon,texContourBasinsLon,texProbaLon);
  normalizeProbabiltyMap(mapBasinsLon,mapContourBasinsLon,texContourBasinsLon,texProbaLon,texProbaNormLon);
  if (_save)
  {
    writeFloatTexture("lon_proba.tex",texProbaLon);
    writeFloatTexture("lon_proba_norm.tex",texProbaNormLon);
  }
  cout << "done " << endl;

  vector<float>::iterator itp;
  for ( itp=_proba.begin(); itp != _proba.end(); itp++)
  {
    cout << "proba value = " << *itp << endl;

    cout << "threshold probability map" << endl;
    TimeTexture<short> texProbaThreshLat(1, _mesh.vertex().size() );
    texBinarizeF2S(texProbaNormLat, texProbaThreshLat, *itp ,0 ,1);
    TimeTexture<short> texProbaThreshLon(1, _mesh.vertex().size() );
    texBinarizeF2S(texProbaNormLon, texProbaThreshLon, *itp ,0 ,1);
    cout << "done " << endl;

    std::ostringstream buff;
    buff<<(*itp);
    size_t found;
    string convert = buff.str();
    found=convert.find(".");
    convert.replace(found,1,",");
    string texname;

    if (_save)
    {
      texname = "lat_proba_thresh_" + convert + ".tex";
      writeShortTexture(texname.c_str(),texProbaThreshLat);

      texname = "lon_proba_thresh_" + convert + ".tex";
      writeShortTexture(texname.c_str(),texProbaThreshLon);
    }

    cout << "compute the longest path in threshold probability map" << endl;
    cout << "latitude :" << endl;
    TimeTexture<short> texProbaSulcalinesLat(1, _mesh.vertex().size() );
    textureBin2Label(texInterRootsBasinsLat,texProbaThreshLat,texTemp);
    if (_save)
      {
      texname = "lat_proba_thresh_" + convert + ".tex";
      writeShortTexture(texname.c_str(),texTemp);
      }

    computeLongestPathBasins (texTemp, texProbaSulcalinesLat, mapBasinsLat);

    cout << "done " << endl;
    cout << "longitude :" << endl;
    TimeTexture<short> texProbaSulcalinesLon(1, _mesh.vertex().size() );
    textureBin2Label(texInterRootsBasinsLon,texProbaThreshLon,texTemp);
    if (_save)
      {
      texname = "lon_proba_thresh_" + convert + ".tex";
      writeShortTexture(texname.c_str(),texTemp);
      }

    computeLongestPathBasins (texTemp, texProbaSulcalinesLon, mapBasinsLon);
    cout << "done " << endl;

    texname = "lat_proba_lines_" + convert + ".tex";
    writeShortTexture(texname.c_str(),texProbaSulcalinesLat);
    texname = "lon_proba_lines_" + convert + ".tex";
    writeShortTexture(texname.c_str(),texProbaSulcalinesLon);
  }
}

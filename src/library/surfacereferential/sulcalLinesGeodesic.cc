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
#include <fstream>

using namespace std;

SulcalLinesGeodesic::SulcalLinesGeodesic(string & adrMesh, string & adrCurv,
    string & adrGeodesicDepth, string & adrRootsLon, string & adrRootsLat,string & adrRootsBottom,string & adrModelParam,
    int extremeties_method, int constraint_type, int strain,
    vector<float> proba, string saveFolder, float curv_thresh) :
  _adrMesh(adrMesh), _adrCurv(adrCurv), _adrGeodesicDepth(adrGeodesicDepth),
      _adrRootsLon(adrRootsLon), _adrRootsLat(adrRootsLat),_adrRootsBottom(adrRootsBottom),_adrModelParam(adrModelParam),
      _constraint_type(constraint_type), _extremeties_method(extremeties_method), _strain(
          strain), _proba(proba), _adrSaveFolder(saveFolder), _curv_thresh(
          curv_thresh)
{
//  std::cout << "Read mesh : ";
  Reader<AimsSurfaceTriangle> r(adrMesh);
  r.read(_mesh);
//  std::cout << "done" << std::endl;

  if (_adrSaveFolder != "")
    _save = true;
  else
    _save = false;

  if (adrRootsLon != "")
  {
    std::cout << "Read roots texture : ";
    Reader<TimeTexture<short> > rtLon(adrRootsLon);
    rtLon.read(_rootsLon);
  }

  if (adrRootsLat != "")
  {
    Reader<TimeTexture<short> > rtLat(adrRootsLat);
    rtLat.read(_rootsLat);
    std::cout << "done" << std::endl;
  }

  if (adrCurv == "")
  {
//    cout << "compute and save curvature texture : ";
    _texCurv = TimeTexture<float> (1, _mesh.vertex().size());
    _texCurv = AimsMeshCurvature(_mesh[0]);

    if (_save)
    {
      writeFloatTexture("curv_new.tex", _texCurv);
//      cout << "done" << endl;
    }
  }
  else
  {
    std::cout << "Read curvature texture : ";
    Reader<TimeTexture<float> > rtCurv(adrCurv);
    rtCurv.read(_texCurv);
    cout << "done" << endl;
  }

  if (adrGeodesicDepth != "")
  {
    std::cout << "Read Geodesic Depth Texture : ";
    Reader<TimeTexture<float> > rtDG(adrGeodesicDepth);
    rtDG.read(_geoDepth);
    std::cout << "done" << std::endl;
  }
  else
    constraint_type = 1;

//  std::cout << "compute neighbours : ";
  _neigh = SurfaceManip::surfaceNeighbours(_mesh);
//  std::cout << "done " << std::endl;
}

SulcalLinesGeodesic::~SulcalLinesGeodesic()
{

}

void SulcalLinesGeodesic::run()
{
  std::cout << "START : ";

  std::cout << "Sulcal basins segmentation " << endl;

  TimeTexture<short> texBasins(1, _mesh.vertex().size());

  map<int, set<int> > mapBasins;

  segmentationSulcalBasins(_texCurv, texBasins, mapBasins, _curv_thresh, 1, 2);

  if (_save)
  {
    std::cout << "Save connected components texture : ";
    writeShortTexture("basins.tex", texBasins);
  }

  if (_adrGeodesicDepth != "")
  {
    cout << "Normalize Geodesic Depth Texture with sulcal basins " << endl;
    _geoDepthNorm = TimeTexture<float> (1, _mesh.vertex().size());
    normalizeDepthMap(_geoDepth, _geoDepthNorm, mapBasins);

    if (_save)
    {
      cout << "Save Normalize Geodesic Depth texture : ";
      writeFloatTexture("depth_norm.tex", _geoDepthNorm);
    }
  }

  switch (_extremeties_method)
  {
  case 1:
    cout << "extraction of extremities method : projection crop by basins"
        << endl;
    sulcalLinesExtract_projection(mapBasins, texBasins);
    break;
  case 2:
    cout << "extraction of extremities method (EMBC11): map of probability" << endl;
    sulcalLinesExtract_probability(mapBasins, texBasins);
    break;
  case 3:

    cout << "extraction of extremities method (NEUROIMAGE): map of density" << endl;
    sulcalLinesExtract_density(mapBasins, texBasins);
    break;


  }

}

void SulcalLinesGeodesic::probaMap()
{
  std::cout << "START : ";

  std::cout << "reading sulcal basins in the curvature texture file " << endl;

  map<int, set<int> > mapBasins;

  TimeTexture<float> texBasinsF;

  if (_adrCurv != "")
  {
    std::cout << "Read basins texture : ";
    Reader<TimeTexture<float> > rtB(_adrCurv);
    rtB.read(texBasinsF);
  }

  TimeTexture<short> texBasins(1, _mesh.vertex().size());
  for (uint i = 0; i < texBasins[0].nItem(); i++)
     texBasins[0].item(i) = (int)(texBasinsF[0].item(i)+ 0.5);

  //texBinarizeF2S(texBasinsF, texBasins, 0, 0, 1);
  texConnectedComponent(texBasins, mapBasins, 1000);

  if (_adrGeodesicDepth != "")
  {
    cout << "Normalize Geodesic Depth Texture with sulcal basins " << endl;
    _geoDepthNorm = TimeTexture<float> (1, _mesh.vertex().size());
    normalizeDepthMap(_geoDepth, _geoDepthNorm, mapBasins);

    if (_save)
    {
      cout << "Save Normalize Geodesic Depth texture : ";
      writeFloatTexture("depth_norm_new.tex", _geoDepthNorm);
    }
  }

  cout << mapBasins.size() << " Basins found" << endl;

  TimeTexture<short> texContourBasins(1, _mesh.vertex().size());
  map<int, set<int> > mapContourBasins;

  contourBasins(mapBasins, texBasins, mapContourBasins, texContourBasins);

  cout << mapContourBasins.size() << " Contour found" << endl;

  if (_save)
  {
    writeShortTexture("contour_basins_new.tex", texContourBasins);
  }

  cout << "compute probability map" << endl;
  TimeTexture<float> texProba(1, _mesh.vertex().size());
  TimeTexture<float> texProbaNorm(1, _mesh.vertex().size());
  computeProbabiltyMap(mapContourBasins, texContourBasins, texProba);
  normalizeProbabiltyMap(mapBasins, mapContourBasins, texContourBasins,
      texProba, texProbaNorm);
  if (_save)
  {
    TimeTexture<short> texAutoThreshold(1, _mesh.vertex().size());
    automaticThresholdDensityMap(mapBasins,texBasins,texProbaNorm,texAutoThreshold,50);
    if (_save)
      writeShortTexture("threshold_auto_ud.tex", texAutoThreshold);

    TimeTexture<short> texSulcalines(1, _mesh.vertex().size() );
    computeLongestPathBasins (texAutoThreshold, texSulcalines, mapBasins);
    if (_save)
      writeShortTexture("sulcalines_ud.tex", texSulcalines);

    writeFloatTexture("proba_ud.tex", texProba);
    writeFloatTexture("proba_norm_ud.tex", texProbaNorm);
  }

  cout << "done " << endl;
}

void SulcalLinesGeodesic::writeShortTexture(string name,
    TimeTexture<short> &out)
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

void SulcalLinesGeodesic::writeFloatTexture(string name,
    TimeTexture<float> &out)
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

void SulcalLinesGeodesic::floodFillIter(int indexVertex, float newTextureValue,
    float oldTextureValue, TimeTexture<short> &texBasinsTemp,
    map<int, set<int> > &mapBasins)
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

      if (texBasinsTemp[0].item(indexCurr) <= (oldTextureValue + 0.001)
          && texBasinsTemp[0].item(indexCurr) >= (oldTextureValue - 0.001))
      //if ( texBasinsTemp[0].item(indexCurr) == oldTextureValue)
      {
        listIndexVertexFill.insert(indexCurr);
        stack.push(indexCurr);
        //cout << indexCurr << " " <<  newTextureValue << endl;
        texBasinsTemp[0].item(indexCurr) = newTextureValue;
      }

    }
  }

  mapBasins.insert(pair<int, set<int> > (newTextureValue, listIndexVertexFill));

}

void SulcalLinesGeodesic::texBinarizeF2S(TimeTexture<float> &texIn,
    TimeTexture<short> &texOut, float threshold, int inf, int sup)
{
  for (uint i = 0; i < texIn[0].nItem(); i++)
  {
    if (texIn[0].item(i) <= threshold)
      texOut[0].item(i) = inf;
    else
      texOut[0].item(i) = sup;
  }
}

void SulcalLinesGeodesic::texBinarizeS2S(TimeTexture<short> &texIn,
    TimeTexture<short> &texOut, int threshold, int inf, int sup)
{
  for (uint i = 0; i < texIn[0].nItem(); i++)
  {
    if (texIn[0].item(i) <= threshold)
      texOut[0].item(i) = inf;
    else
      texOut[0].item(i) = sup;
  }
}

TimeTexture<short> SulcalLinesGeodesic::texConnectedComponent(
    TimeTexture<short> &texBasins, map<int, set<int> > &mapBasins, int offset)
{
  int j = 1;
  int val;
  mapBasins.clear();

  TimeTexture<short> texTemp(1, texBasins[0].nItem());

  for (uint i = 0; i < texBasins[0].nItem(); i++)
    texTemp[0].item(i) = texBasins[0].item(i) + offset;

  for (uint i = 0; i < texBasins[0].nItem(); i++)
  {
    if (texTemp[0].item(i) > offset)
      floodFillIter(i, j++, texTemp[0].item(i), texTemp, mapBasins);
  }

  return texTemp;
}

void SulcalLinesGeodesic::constraintListOpen(set<int> & constraintListValues)
{
  //QString filt = "(*.txt)" ;
  //QString capt = "Open constraints file" ;

  //QString filename;

  //_adrModelParam
  //filename = QFileDialog::getOpenFileName( QString::null, filt, 0, 0, capt );

  string line;
  ifstream myfile(_adrModelParam.c_str());

  cout << _adrModelParam << endl;

  int value;

  if (myfile.is_open())
  {
    cout << "open " << _adrModelParam.c_str() << endl;

    while (myfile.good())
    {
      getline(myfile, line);
      if (line.length() != 0)
        {
        //constraintListNames.insert(line);
        //istringstream iss(line);
        //iss >> value;
        int position = line.find_last_of('\t');
        std::istringstream strin(line.substr(position + 1));
        strin >> value;
        //cout << value << endl;
        constraintListValues.insert(value);
        }

    }
    myfile.close();
  }

}

void SulcalLinesGeodesic::computeRootsBottomMap(TimeTexture<short> &texBasins,TimeTexture<short> &texRootsBottom, float dist_max)
{
  cout << "reading bottom volume  : " << _adrRootsBottom << flush;
  AimsData<short> bottom;
  Reader<AimsData<short> > bottomR(_adrRootsBottom);
  bottomR >> bottom;
  cout << "done" << endl;

  set<int> constraintListValues;

  constraintListOpen(constraintListValues);

  int x,y,z,sx,sy,sz;
  float dx, dy, dz;

  sx=bottom.dimX(); sy=bottom.dimY(); sz=bottom.dimZ();
  dx=bottom.sizeX(); dy=bottom.sizeY(); dz=bottom.sizeZ();

  std::cout << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << endl;
  std::cout << "sx=" << sx << ", sy=" << sy << ", sz=" << sz << endl;

  const vector<Point3df>      & vert = _mesh.vertex() ;
  vector<pair<Point3df,int> >  cloud;
  set<int>::iterator it;

  for (int z=0; z<sz; z++)
   for (int y=0; y<sy; y++)
     for (int x=0; x<sx; x++)
     {
       if (bottom(x,y,z)> 0)
         cloud.push_back(pair<Point3df,int>(Point3df(dx*x,dy*y,dz*z), bottom(x,y,z) ));
     }

  vector<pair<Point3df,int> >::const_iterator itb(cloud.begin()), ite(cloud.end());

  cout << "nb bottom point = " << cloud.size() << endl;

  const vector<Point3df> & norm = _mesh.normal();

  int imin;
  float min,dist,xx,yy,zz,no;
  float vvx,vvy,vvz,angle;

  float testx,testy,testz;

  AimsSurfaceTriangle *tmpMeshOut,meshOut;
  tmpMeshOut = new AimsSurfaceTriangle;

//  for (itb = cloud.begin(); itb != cloud.end(); ++itb)
//  {
//    tmpMeshOut = SurfaceGenerator::sphere((itb->first), 0.05 ,10 );
//    SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
//  }

  for (uint i = 0; i < _mesh.vertex().size(); i++)
  {
    texRootsBottom[0].item(i) = 0;
    min = 10000000.;


    if (texBasins[0].item(i) > 0)
    {
      for (itb = cloud.begin(); itb != cloud.end(); ++itb)
      {
        xx = ((itb->first)[0] - vert[i][0]);
        yy = ((itb->first)[1] - vert[i][1]);
        zz = ((itb->first)[2] - vert[i][2]);

//        Point3df tt((itb->first)[0],(itb->first)[1],(itb->first)[2]);

//        tmpMeshOut = SurfaceGenerator::sphere(tt, 0.05 ,10 );
//        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );

        dist = (xx*xx + yy*yy+ zz*zz);

        double scalar;

        //on calcule le produit scalaire p entre le vecteur vv [vertex,voxel] et la normale au vertex
        //si p est négatif alors on attribue le label -1 sinon on attribue le label du voxel
        scalar = xx*norm[i][0] + yy*norm[i][1] + zz*norm[i][2];
        no = sqrt(xx*xx + yy*yy + zz*zz);
        angle=acos((float)scalar/no)*(180./M_PI);

        //if (scalar > 0 )
        if (angle < 45)
        {
        if (dist < min)
          {
          min = dist;
          imin = itb->second;
          vvx = xx;
          vvy = yy;
          vvz = zz;
          }
        }
//        else
//        {
//          if (dist < min)
//          {
//            no = sqrt(xx*xx + yy*yy + zz*zz);
//            angle=acos((float)p/no)*(180./M_PI);
//            if (angle < 120)
//            {
//              min = dist;
//              imin = itb->second;
//              vvx = xx;
//              vvy = yy;
//              vvz = zz;
//            }
//          }
//
//
//        }

        //cout << itb->second << " " << (itb->first)[0]<< " " << (itb->first)[1]<< " " << (itb->first)[2] << endl;
      }

      it=constraintListValues.find(imin);
      if (it!=constraintListValues.end())
      {
        if (min < dist_max)
          texRootsBottom[0].item(i) = imin;
        else
          {
          texRootsBottom[0].item(i) = -1;

          //cout << i << " " << min << endl;
          tmpMeshOut = SurfaceGenerator::sphere(vert[i], 0.1 ,10 );
          SurfaceManip::meshMerge( meshOut, *tmpMeshOut );

          tmpMeshOut = SurfaceGenerator::cylinder(vert[i], vert[i]+ norm[i],0.05, 0.05, 12, false, true);
          SurfaceManip::meshMerge( meshOut, *tmpMeshOut );

          tmpMeshOut = SurfaceGenerator::sphere(vert[i]+ norm[i], 0.1 ,10 );
          SurfaceManip::meshMerge( meshOut, *tmpMeshOut );

          Point3df test(vvx,vvy,vvz);
          tmpMeshOut = SurfaceGenerator::sphere(vert[i]+ test, 0.1 ,10 );

          SurfaceManip::meshMerge( meshOut, *tmpMeshOut );

          tmpMeshOut = SurfaceGenerator::cylinder(vert[i], vert[i]+ test,0.05, 0.05, 12, false, true);
          SurfaceManip::meshMerge( meshOut, *tmpMeshOut );

          }
      }


      //on calcule le produit scalaire p entre le vecteur vv [vertex,voxel] et la normale au vertex
      //si p est négatif alors on attribue le label -1 sinon on attribue le label du voxel
//      double p;
//      //Point3df vv(vvx,vvy,vvz);
//
//      p = vvx*norm[i][0] + vvy*norm[i][1] + vvz*norm[i][2];




//      if (i == 87440)
//        cout << " scal = " << p << " " << norm[i][0]  << " " << norm[i][1] << " " << norm[i][2] << endl;
//
//      if (i == 50544)
//        cout << " scal = " << p << " " << norm[i][0]  << " " << norm[i][1] << " " << norm[i][2] << endl;

//      if (p < 0 )
//        {
//        no = sqrt(vvx*vvx + vvy*vvy + vvz*vvz);
//
//        angle=acos((float)p/no)*(180./M_PI);
//
////        cout << "angle = " << angle <<  endl;
////
////        cout << " scal = " << p << " " << norm[i][0]  << " " << norm[i][1] << " " << norm[i][2] << endl;
////
//        testx = vvx;
//        testy = vvy;
//        testz = vvz;
//
//        tmpMeshOut = SurfaceGenerator::sphere(vert[i], 0.1 ,10 );
//        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
//
//        tmpMeshOut = SurfaceGenerator::cylinder(vert[i], vert[i]+ norm[i],0.05, 0.05, 12, false, true);
//        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
//
//        tmpMeshOut = SurfaceGenerator::sphere(vert[i]+ norm[i], 0.1 ,10 );
//        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
//
//
//        Point3df test(testx,testy,testz);
//
//        tmpMeshOut = SurfaceGenerator::sphere(vert[i]+ test, 0.1 ,10 );
//
//        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
//
//        tmpMeshOut = SurfaceGenerator::cylinder(vert[i], vert[i]+ test,0.05, 0.05, 12, false, true);
//        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
//        }
//
//      it=constraintListValues.find(imin);
//      if (it!=constraintListValues.end())
//      {
//        if (p >= 0)
//          texRootsBottom[0].item(i) = imin;
//        else
//          texRootsBottom[0].item(i) = 500;
//      }
//
    }
    cout << "%\r\033[K" << (int)(100*(float)i/_mesh.vertex().size());

  }






  Writer<AimsSurfaceTriangle> wm("test.mesh");
  wm.write(meshOut);


}

void SulcalLinesGeodesic::listRootsProjections(TimeTexture<short> &texBasins,
    set<int> &listIndexLat, set<int> &listIndexLon)
{
  int value;

  for (uint i = 0; i < texBasins[0].nItem(); i++)
  {
    value = texBasins[0].item(i);
    if (value != 0)
    {
      if (_rootsLat[0].item(i) != 0)
        listIndexLat.insert(i);

      if (_rootsLon[0].item(i) != 0)
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

void SulcalLinesGeodesic::computeListLabelProjectionsBasins(
    TimeTexture<short> &roots, map<int, set<int> > &mapBasins,
    set<int> &listIndex, map<int, set<int> > &mapConstraint)
{
  int nbNeigh, value;
  set<int>::iterator itl;
  set<int>::iterator it;
  set<int> listIndexTemp;
  int nbConstraint = 0;

  set<int> listLabelBasins;
  set<int>::iterator itLabel;

  map<int, set<int> >::const_iterator mit(mapBasins.begin()), mend(
      mapBasins.end());

  //On parcourt les tous bassins
  for (; mit != mend; ++mit)
  {
    listIndexTemp.clear();
    //cout << "\nbasin " << mit->first << endl;

    //recherche des étiquettes présentes dans le bassin
    listLabelBasins.clear();
    for (itl = listIndex.begin(); itl != listIndex.end(); itl++)
    {
      it = (mit->second).find(*itl);
      value = roots[0].item(*itl);
      if (it != (mit->second).end())
        listLabelBasins.insert(value);
    }
    //pour chaque etiquette du bassin
    for (itLabel = listLabelBasins.begin(); itLabel != listLabelBasins.end(); itLabel++)
    {
      //cout << *itLabel << " --> ";
      for (itl = listIndex.begin(); itl != listIndex.end(); itl++)
      {
        it = (mit->second).find(*itl);
        value = roots[0].item(*itl);

        //si une lat est dans le bassin alors on l'ajoute à la liste temporaire des points du bassin de valeur itLabel
        if (it != (mit->second).end() && (value == *itLabel))
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
        mapConstraint.insert(pair<int, set<int> > (nbConstraint++,
            listIndexTemp));
    }
  }

}

void SulcalLinesGeodesic::computeLongestPathBasins(TimeTexture<short> &roots,
    TimeTexture<short> &out, map<int, set<int> > &mapConstraint)
{
  map<int, set<int> >::const_iterator mit(mapConstraint.begin()), mend(
      mapConstraint.end());

  set<int>::iterator it;
  vector<unsigned> pathTemp, indexTemp;

  for (uint i = 0; i < _mesh.vertex().size(); i++)
    out[0].item(i) = 0.0;

  TimeTexture<float> texConstraint(1, _mesh.vertex().size());
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

  GeodesicPath sp(_mesh, texConstraint, method, _strain);

  int source, target;
  int constraintValue;

  for (; mit != mend; mit++)
  {
    indexTemp.clear();

    cout << "\nbasin " << (int) mit->first << " : \n";
    it = (mit->second).begin();

    //on copie la liste des index de sommets set --> vector (non nul)
    for (; it != (mit->second).end(); it++)
    {
      constraintValue = roots.item(*it);
      if (constraintValue != 0)
        indexTemp.push_back(*it);
    }

    //si il y a plus d'un point dans le bassin
    if (indexTemp.size() > 1)
    {
      //on récupère la valeur de la contrainte
      constraintValue = roots.item(*(indexTemp.begin()));

      cout << "value " << constraintValue << endl;
      //on calcule le plus long chemin
      pathTemp.clear();
      double len;

      pathTemp = sp.longestPath_N_N_ind(indexTemp, &source, &target, &len, 1);

      for (int i = 0; i < pathTemp.size(); i++)
        out[0].item(pathTemp[i]) = constraintValue;
    }
    else
      cout << "not enough points ! " << endl;
  }
}

void SulcalLinesGeodesic::normalizeDepthMap(TimeTexture<float> &depth,
    TimeTexture<float> &depthNorm, map<int, set<int> > &mapBasins)
{
  map<int, set<int> >::const_iterator mit(mapBasins.begin()), mend(
      mapBasins.end());
  set<int>::iterator listVertexbasin;
  int n_b = 0;
  float max_depth, val;

  //on normalise la profondeur dans chaque bassin entre 0 et 1
  for (; mit != mend; ++mit)
  {
    max_depth = 0;

    listVertexbasin = (mit->second).begin();
    for (; listVertexbasin != (mit->second).end(); listVertexbasin++)
    {
      val = depth[0].item(*listVertexbasin);
      if (val > max_depth)
        max_depth = val;
    }

    listVertexbasin = (mit->second).begin();
    for (; listVertexbasin != (mit->second).end(); listVertexbasin++)
    {
      if (max_depth != 0)
        depthNorm[0].item(*listVertexbasin) = (float) (depth[0].item(
            *listVertexbasin)) / max_depth;
      else
        depthNorm[0].item(*listVertexbasin) = 0;
    }
  }
}

void SulcalLinesGeodesic::segmentationSulcalBasins(TimeTexture<float> &texIn,
    TimeTexture<short> &texBasins, map<int, set<int> > &mapBasins,
    float threshold, int close, int open)
{
  //binarisation texture de courbure < 0.0
  texBinarizeF2S(texIn, texBasins, threshold, 1, 0);

  TimeTexture<short> texBasinsDil(1, _mesh.vertex().size());
  TimeTexture<short> texBasinsErode(1, _mesh.vertex().size());

  //fermeture des bassins
  texBasinsDil[0] = MeshDilation<short> (_mesh[0], texBasins[0], short(0), -1,
      close, true);
  texBasinsErode[0] = MeshErosion<short> (_mesh[0], texBasinsDil[0], short(0),
      -1, close, true);

  //ouverture des bassins
  texBasinsDil[0] = MeshErosion<short> (_mesh[0], texBasinsErode[0], short(0),
      -1, open, true);
  texBasinsErode[0] = MeshDilation<short> (_mesh[0], texBasinsDil[0], short(0),
      -1, open, true);

  texBinarizeS2S(texBasinsErode, texBasins, 0, 0, 1);

  // étiquetage des composantes connexes
  texConnectedComponent(texBasins, mapBasins, 1000);
  cout << mapBasins.size() << " Basins found" << endl;
}

void SulcalLinesGeodesic::sulcalLinesExtract_projection(
    map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins)
{
  //liste les projections roots et conservent seulement celles inclues dans les bassins
  set<int> listIndexLon, listIndexLat;
  listRootsProjections(texBasins, listIndexLat, listIndexLon);

  //textures contenant les contraintes lat et lon
  TimeTexture<short> texSulcalinesLat(1, _mesh.vertex().size());
  TimeTexture<short> texSulcalinesLon(1, _mesh.vertex().size());

  map<int, set<int> > mapConstraintLat, mapConstraintLon;

  cout << "Sort Lat/Lon by basin : " << endl;

  //groupe et étiquette dans des map les projections lat et lon de chaque bassin
  //conserve seulement les projections qui sont potentiellement des extremités des lignes (points qui ont au plus un voisin)

  computeListLabelProjectionsBasins(_rootsLat, mapBasins, listIndexLat,
      mapConstraintLat);
  cout << endl << mapConstraintLat.size() << " Basins Latitude extracted";
  computeLongestPathBasins(_rootsLat, texSulcalinesLat, mapConstraintLat);

  computeListLabelProjectionsBasins(_rootsLon, mapBasins, listIndexLon,
      mapConstraintLon);
  cout << endl << mapConstraintLon.size() << " Basins Longitude extracted";
  computeLongestPathBasins(_rootsLon, texSulcalinesLon, mapConstraintLon);

  std::cout << "\nSave sulcal lines texture : " << endl;

  writeShortTexture("lat_sulcal_lines.tex", texSulcalinesLat);
  writeShortTexture("lon_sulcal_lines.tex", texSulcalinesLon);
}

void SulcalLinesGeodesic::dilationRoots(TimeTexture<short> &texLatDil,
    TimeTexture<short> &texLonDil, int size)
{
  texLatDil[0] = MeshDilation<short> (_mesh[0], _rootsLat[0], short(0), -1,
      size, false);
  texLonDil[0] = MeshDilation<short> (_mesh[0], _rootsLon[0], short(0), -1,
      size, false);
}

void SulcalLinesGeodesic::interRootsDilBasins(TimeTexture<short> &texBasins,
    TimeTexture<short> &texDil, TimeTexture<short> &texInter)
{
  for (uint i = 0; i < _mesh.vertex().size(); i++)
    if (texBasins[0].item(i) != 0 && texDil[0].item(i) > 0)
      texInter[0].item(i) = texDil[0].item(i);
}

void SulcalLinesGeodesic::cleanBasins(map<int, set<int> > &mapBasins,
    TimeTexture<short> &texBasins, int nbPoint)
{
  map<int, set<int> >::iterator it;
  set<int>::iterator its;

  map<int, set<int> > mapBasinsClean;

  for (it = mapBasins.begin(); it != mapBasins.end(); ++it)
  {

    its = ((*it).second).begin();

    //cout << "bassin " << (*it).first << " " << (*it).second.size() << "\n";

    //calcul de l'aire du bassin A = 1/2||a^b||

//    n.x = (a.y * b.z) - (a.z * b.y)
//    n.y = (a.z * b.x) - (a.x * b.z)
//    n.z = (a.x * b.y) - (a.y * b.x)

    if (((*it).second).size() < nbPoint)
    {
      //on efface tous les vertex de la texture texBasins groupés dans des bassins qui contiennent moins de nbPoint
      for (; its != ((*it).second).end(); its++)
        texBasins[0].item(*its) = 0;

      //cout << " del\n";
      // on l'enlève ensuite de la map des bassins
      //mapBasins.erase (it);
    }
    else
      mapBasinsClean.insert(*it);
  }

  //cout << "cleaning\n" << endl;
  for (it = mapBasinsClean.begin(); it != mapBasinsClean.end(); ++it)
    cout << "bassin " << (*it).first << " " << (*it).second.size() << "\n";

  mapBasins = mapBasinsClean;
}

void SulcalLinesGeodesic::contourBasins(map<int, set<int> > &mapBasins,
    TimeTexture<short> &texBasins, map<int, set<int> > &mapContourBasins,
    TimeTexture<short> &texContourBasins)
{
  map<int, set<int> >::iterator it;
  set<int>::iterator its;
  set<int> listIndexTemp;
  int value;
  int nb = 0;

  mapContourBasins.clear();

  //pour tous les bassins
  for (it = mapBasins.begin(); it != mapBasins.end(); it++)
  {
    listIndexTemp.clear();

    //pour tous les points its du bassin it
    for (its = ((*it).second).begin(); its != ((*it).second).end(); its++)
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
    mapContourBasins.insert(pair<int, set<int> > (nb++, listIndexTemp));
  }
}

void SulcalLinesGeodesic::computeProbabiltyMap(
    map<int, set<int> > &mapContourBasins,
    TimeTexture<short> &texContourBasins, TimeTexture<float> &texProba)
{
  map<int, set<int> >::iterator itc;
  set<int>::iterator its;
  int value, nb_c, nb_v;
  vector<unsigned> listIndexTarget;
  vector<vector<unsigned> > indices;
  vector<vector<unsigned> >::iterator it_vv;
  vector<unsigned>::iterator it_v;

  for (uint i = 0; i < texProba[0].nItem(); i++)
    texProba[0].item(i) = 1.0;

  TimeTexture<float> texConstraint(1, _mesh.vertex().size());
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

  GeodesicPath sp(_mesh, texConstraint, method, _strain);

  //pour tous les bassins
  for (itc = mapContourBasins.begin(); itc != mapContourBasins.end(); itc++)
  {
    listIndexTarget.clear();

    //on récupère l'étiquette du bassin (avec le premier point)
    its = ((*itc).second).begin();

    value = texContourBasins[0].item(*its);
    cout << "value = " << value << endl;

    //pour tous les points its du bassin it
    for (; its != ((*itc).second).end(); its++)
      listIndexTarget.push_back(*its);

    //nombre de points de contours
    nb_c = 0;
    //nombre de sommets parcourus
    nb_v = 0;

    for (its = ((*itc).second).begin(); its != ((*itc).second).end(); its++)
    {
      indices.clear();

      sp.shortestPath_1_N_All_ind(*its, listIndexTarget, indices);

      nb_v += indices.size();
      cout << "\r\033[K" << "extremities : " << ++nb_c << "/"
          << listIndexTarget.size() << " vertex : " << nb_v << flush;

      //on incrémente les valeurs de la texture proba pour chaque point ayant été parcourue
      for (it_vv = indices.begin(); it_vv != indices.end(); it_vv++)
        for (it_v = (*it_vv).begin(); it_v != (*it_vv).end(); it_v++)
        {
          texProba[0].item(*it_v)++;
        }
    }

    cout << endl;
  }
}

void SulcalLinesGeodesic::normalizeProbabiltyMap(
    map<int, set<int> > &mapBasins, map<int, set<int> > &mapContourBasins,
    TimeTexture<short> &texContourBasins, TimeTexture<float> &texProba,
    TimeTexture<float> &texProbaNorm)
{
  map<int, set<int> >::iterator itc;
  map<int, set<int> >::iterator itb;
  set<int>::iterator its;

  int value;
  float accu, accu_max;

  for (uint i = 0; i < texProba[0].nItem(); i++)
    texProbaNorm[0].item(i) = 0.0;

  //pour tous les bassins
  for (itc = mapContourBasins.begin(), itb = mapBasins.begin(); itc
      != mapContourBasins.end(), itb != mapBasins.end(); itc++, itb++)
  {
    accu_max = 0;

    //on récupère l'étiquette du bassin (avec le premier point)
    its = ((*itb).second).begin();

    value = texContourBasins[0].item(*its);

    //pour tous les points its du bassin it
    for (; its != ((*itb).second).end(); its++)
      accu_max = max(accu_max, texProba[0].item(*its));

    //on normalise entre 0 et 1
    for (its = ((*itb).second).begin(); its != ((*itb).second).end(); its++)
      texProbaNorm[0].item(*its) = (float) (texProba[0].item(*its) / accu_max);

    //on met tous les points de contour du bassins à 0
    //      for (its = ((*itc).second).begin() ; its != ((*itc).second).end(); its++ )
    //        texProbaNorm[0].item(*its) = 0;


  }
}

void SulcalLinesGeodesic::automaticThresholdDensityMap(map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins,
    TimeTexture<float> &texProbaNorm, TimeTexture<short> &texAutoThreshold,int nb_bin)
{
  map<int, set<int> >::iterator itb;
  set<int>::iterator its;

  int value;
  map<uint, float> histo, smooth, sig;

  int i;

  ofstream myfile;

  if (_save)
  {
    myfile.open("infos.txt");
  }

  //pour tous les bassins
  for (itb = mapBasins.begin(); itb != mapBasins.end(); itb++)
  {
    for (i = 0; i <= nb_bin; i++)
    {
      histo[i] = 0.0;
      smooth[i] = 0.0;
      sig[i] = 0.0;
    }

    //on récupère l'étiquette du bassin (avec le premier point)
    its = ((*itb).second).begin();

    value = texBasins[0].item(*its);

    cout << "\nbassin " << (*itb).first << " --> " << value << " "
        << (*itb).second.size() << "\n";

    if (_save)
      myfile << "\nbassin " << (*itb).first << " --> " << value << " " << (*itb).second.size() << "\n";
    //pour tous les points its du bassin it
    for (; its != ((*itb).second).end(); its++)
    {
      //cout << " " << texProbaNorm[0].item(*its);
      //myfile << texProbaNorm[0].item(*its) << "\n";
      histo[(texProbaNorm[0].item(*its)) * nb_bin]++;
    }

    map<uint, float>::iterator ith;


    //cout << "smoothing histo" << "\n";

    for (i = 0; i <= nb_bin; i++)
    {
      smooth[i] = histo[i];
      //smooth.insert(std::pair<uint, float>(i, curvM[i]));
      //sig.insert(std::pair<uint, float>(i, curvM[i]));
    }

    float lapl;
    for (uint t = 0; t < 10; t++)
    {
      for (i = 0; i <= nb_bin; i++)
        sig[i] = smooth[i];

      for (i = 0; i <= nb_bin; i++)
      {
        if (i == 0)
          lapl = sig[1] - sig[0];
        else if (i == nb_bin)
          lapl = sig[nb_bin - 1] - sig[nb_bin];
        else
          lapl = sig[i - 1] - 2 * sig[i] + sig[i + 1];

        smooth[i] = sig[i] + 0.20 * lapl * 0.5;
      }
    }

    if (_save)
      for (ith = smooth.begin(); ith != smooth.end(); ith++)
        myfile << ith->first << "\t" << ith->second << endl;

    if (_save)
      myfile << "proba" << "\t" << "nb points" << "\t" << "nb extremites" << "\n";

    // histo pour déterminer le seuil automatique
    // on compte le nombre de vertex "point extrémité" (un seul voisin)
    vector<float>::iterator itp;
    vector<int> extremities;

    for (int j = 0 ; j <= 100 ; j+=(100/nb_bin))
    {
      if (_save)
        myfile << (float)j/100. << "\t";

      TimeTexture<short> texProbaThresh(1, _mesh.vertex().size());
      texBinarizeF2S(texProbaNorm, texProbaThresh, (float)j/100. , 0, 1);

      int nb_points;
      nb_points = 0;
      int nb_ext;
      nb_ext = 0;

      set<int>::iterator its;

      //pour tous les points its du bassin it
      for (its = ((*itb).second).begin(); its != ((*itb).second).end(); its++)
      {
        //on récupère l'étiquette du bassin
        int value_thresh = texProbaThresh[0].item(*its);

        if (value_thresh == 1)
        {
          nb_points++;

          set<uint> voisins = _neigh[*its];
          set<uint>::iterator voisIt = voisins.begin();
          int nb_vois;
          nb_vois = 0;
          //on parcourt tous les voisins du sommet, si un voisin a une valeur differente alors le point est un contour
          for (; voisIt != voisins.end(); voisIt++)
          {
            if (texProbaThresh[0].item(*voisIt) == value_thresh)
              nb_vois++;
          }

          if (nb_vois <= 1)
            nb_ext++;
        }
      }

      extremities.push_back(nb_ext);

      if (_save)
      myfile << nb_points << "\t" << nb_ext << endl;
    }

    vector<float> dx;
    vector<uint> gz, lz, lo;
    //recherche des minima locaux
    for (i = 0; i < nb_bin; i++)
    {
      dx.push_back(smooth[i + 1] - smooth[i]);
      if (dx[i] > 0)
      {
        gz.push_back(1);
        lz.push_back(0);
      }
      else
      {
        gz.push_back(0);
        lz.push_back(1);
      }

    }

    for (i = 0; i < nb_bin - 1; i++)
    {
      if (gz[i + 1] && lz[i])
        lo.push_back(i + 1);
    }

    if (_save)
    myfile << "minima locaux\n" << endl;

    vector<uint>::iterator itlo;

    if (_save)
      {
      for (itlo = lo.begin(); itlo < lo.end(); itlo++)
        myfile << *itlo << " ";
      myfile << endl;
      }

    // extraction automatique du seuil
    int val_ext;
    int min_ext;
    int seuil = 10000;
    min_ext = 10000;


//    for (i = 0; i < nb_bin; i++)
//    //for (itext = extremities.begin(); itext != extremities.end(); itext++)
//    {
//      val_ext = extremities[i];
//      //cout << val_ext << " ";
//
//      if (val_ext < min_ext && val_ext > 1)
//      {
//        min_ext = val_ext;
//        seuil = i;
//      }
//    }

    for (itlo = lo.begin(); itlo < lo.end(); itlo++)
    {
      val_ext = extremities[*itlo];
      if (val_ext < min_ext && val_ext > 1)
      {
        min_ext = val_ext;
        seuil = *itlo;
      }
    }

    if (_save)
    myfile << "seuil select = " << seuil << " nb_ext = " << min_ext << endl;

    cout << "seuil select = " << seuil << " nb_ext = " << min_ext << endl;

    while (extremities[seuil] <= min_ext)
      seuil--;

    seuil++;

    if (_save)
    myfile << "descente seuil auto = " << (float)2*seuil/100. << "\n";

    cout << "descente seuil auto = " << (float)2*seuil/100. << "\n";

    TimeTexture<short> texProbaThresh2(1, _mesh.vertex().size());
    texBinarizeF2S(texProbaNorm, texProbaThresh2, (float)2*seuil/100. , 0, 1);

    //pour tous les points its du bassin it
    for (its = ((*itb).second).begin(); its != ((*itb).second).end(); its++)
    //for (uint i = 0; i < texProbaThresh2[0].nItem(); i++)
      if (texProbaThresh2[0].item(*its) == 1)
        texAutoThreshold[0].item(*its) = value;

  }

  if (_save)
    myfile.close();
}

void SulcalLinesGeodesic::textureBin2Label(TimeTexture<short> &texLabel,
    TimeTexture<short> &texIn, TimeTexture<short> &texOut)
{
  for (uint i = 0; i < texIn[0].nItem(); i++)
  {
    if (texIn[0].item(i) == 1)
      texOut[0].item(i) = texLabel[0].item(i);
    else
      texOut[0].item(i) = 0;
  }
}

void SulcalLinesGeodesic::sulcalLinesExtract_probability(
    map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins)
{
  TimeTexture<short> texTemp(1, _mesh.vertex().size());

  // dilatation sur les maillages white_fine de 7, sur les maillages white 5
  int dilation_size = 7;

  //on dilate les roots
  TimeTexture<short> texProjectionLatDil(1, _mesh.vertex().size());
  TimeTexture<short> texProjectionLonDil(1, _mesh.vertex().size());
  dilationRoots(texProjectionLatDil, texProjectionLonDil, dilation_size);

  if (_save)
  {
    writeShortTexture("lat_roots_dil.tex", texProjectionLatDil);
    writeShortTexture("lon_roots_dil.tex", texProjectionLonDil);
  }

  // on calcule les intersections avec les basins
  TimeTexture<short> texInterRootsBasinsLat(1, _mesh.vertex().size());
  TimeTexture<short> texInterRootsBasinsLon(1, _mesh.vertex().size());

  map<int, set<int> > mapBasinsLat;
  map<int, set<int> > mapBasinsLon;

  interRootsDilBasins(texBasins, texProjectionLatDil, texInterRootsBasinsLat);
  interRootsDilBasins(texBasins, texProjectionLonDil, texInterRootsBasinsLon);

  texConnectedComponent(texInterRootsBasinsLat, mapBasinsLat, 1000);
  texConnectedComponent(texInterRootsBasinsLon, mapBasinsLon, 1000);

  cout << mapBasinsLat.size() << " Basins Lat found" << endl;
  cout << mapBasinsLon.size() << " Basins Lon found" << endl;

  if (_save)
  {
    writeShortTexture("lat_inter_roots_dil_before.tex", texInterRootsBasinsLat);
    writeShortTexture("lon_inter_roots_dil_before.tex", texInterRootsBasinsLon);
  }

  //on ne conserve que les bassins qui au plus de "clean_size" vertex
  int clean_size = 250;

  cleanBasins(mapBasinsLat, texInterRootsBasinsLat, clean_size);
  cleanBasins(mapBasinsLon, texInterRootsBasinsLon, clean_size);
  cout << "after cleaning" << endl;
  cout << mapBasinsLat.size() << " Basins Lat found" << endl;
  cout << mapBasinsLon.size() << " Basins Lon found" << endl;

  if (_save)
  {
    writeShortTexture("lat_inter_roots_dil.tex", texInterRootsBasinsLat);
    writeShortTexture("lon_inter_roots_dil.tex", texInterRootsBasinsLon);
  }

  TimeTexture<short> texContourBasinsLat(1, _mesh.vertex().size());
  TimeTexture<short> texContourBasinsLon(1, _mesh.vertex().size());
  map<int, set<int> > mapContourBasinsLat;
  map<int, set<int> > mapContourBasinsLon;
  contourBasins(mapBasinsLat, texInterRootsBasinsLat, mapContourBasinsLat,
      texContourBasinsLat);
  contourBasins(mapBasinsLon, texInterRootsBasinsLon, mapContourBasinsLon,
      texContourBasinsLon);
  if (_save)
  {
    writeShortTexture("lat_contour_basins.tex", texContourBasinsLat);
    writeShortTexture("lon_contour_basins.tex", texContourBasinsLon);
  }

  cout << "compute probability map" << endl;
  cout << "latitude :" << endl;
  TimeTexture<float> texProbaLat(1, _mesh.vertex().size());
  TimeTexture<float> texProbaNormLat(1, _mesh.vertex().size());
  computeProbabiltyMap(mapContourBasinsLat, texContourBasinsLat, texProbaLat);
  normalizeProbabiltyMap(mapBasinsLat, mapContourBasinsLat,
      texContourBasinsLat, texProbaLat, texProbaNormLat);

  if (_save)
  {
    writeFloatTexture("lat_proba.tex", texProbaLat);
    writeFloatTexture("lat_proba_norm.tex", texProbaNormLat);
  }
  cout << "done " << endl;
  cout << "longitude :" << endl;
  TimeTexture<float> texProbaLon(1, _mesh.vertex().size());
  TimeTexture<float> texProbaNormLon(1, _mesh.vertex().size());
  computeProbabiltyMap(mapContourBasinsLon, texContourBasinsLon, texProbaLon);
  normalizeProbabiltyMap(mapBasinsLon, mapContourBasinsLon,
      texContourBasinsLon, texProbaLon, texProbaNormLon);

  if (_save)
  {
    writeFloatTexture("lon_proba.tex", texProbaLon);
    writeFloatTexture("lon_proba_norm.tex", texProbaNormLon);
  }
  cout << "done " << endl;

  //
  //  vector<float>::iterator itp;
  //  for ( itp=_proba.begin(); itp != _proba.end(); itp++)
  //  {
  //    cout << "proba value = " << *itp << endl;
  //
  //    cout << "threshold probability map" << endl;
  //    TimeTexture<short> texProbaThreshLat(1, _mesh.vertex().size() );
  //    texBinarizeF2S(texProbaNormLat, texProbaThreshLat, *itp ,0 ,1);
  //    TimeTexture<short> texProbaThreshLon(1, _mesh.vertex().size() );
  //    texBinarizeF2S(texProbaNormLon, texProbaThreshLon, *itp ,0 ,1);
  //    cout << "done " << endl;
  //
  //    std::ostringstream buff;
  //    buff<<(*itp);
  //    size_t found;
  //    string convert = buff.str();
  //    found=convert.find(".");
  //    convert.replace(found,1,",");
  //    string texname;
  //
  //    if (_save)
  //    {
  //      texname = "lat_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texProbaThreshLat);
  //
  //      texname = "lon_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texProbaThreshLon);
  //    }
  //
  //    cout << "compute the longest path in threshold probability map" << endl;
  //    cout << "latitude :" << endl;
  //    TimeTexture<short> texProbaSulcalinesLat(1, _mesh.vertex().size() );
  //    textureBin2Label(texInterRootsBasinsLat,texProbaThreshLat,texTemp);
  //    if (_save)
  //      {
  //      texname = "lat_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texTemp);
  //      }
  //
  //    computeLongestPathBasins (texTemp, texProbaSulcalinesLat, mapBasinsLat);
  //
  //    cout << "done " << endl;
  //    cout << "longitude :" << endl;
  //    TimeTexture<short> texProbaSulcalinesLon(1, _mesh.vertex().size() );
  //    textureBin2Label(texInterRootsBasinsLon,texProbaThreshLon,texTemp);
  //    if (_save)
  //      {
  //      texname = "lon_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texTemp);
  //      }
  //
  //    computeLongestPathBasins (texTemp, texProbaSulcalinesLon, mapBasinsLon);
  //    cout << "done " << endl;
  //
  //    texname = "lat_proba_lines_" + convert + ".tex";
  //    writeShortTexture(texname.c_str(),texProbaSulcalinesLat);
  //    texname = "lon_proba_lines_" + convert + ".tex";
  //    writeShortTexture(texname.c_str(),texProbaSulcalinesLon);
  //  }
}


void SulcalLinesGeodesic::sulcalLinesExtract_density(
    map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins)
{

  TimeTexture<short> texRootsBottom(1, _mesh.vertex().size());

  // carré de la distance max (entre un vertex et le plus proche voxel) pour l'attribution du label du voxel au vertex
  float dist_max = 100;
  computeRootsBottomMap(texBasins,texRootsBottom,dist_max);

  if (_save)
    {
      cout << "Save Roots Bottom texture : ";
      writeShortTexture("roots_bottom.tex", texRootsBottom);
    }

  map<int, set<int> > mapBasinsRoots;

  texConnectedComponent(texRootsBottom, mapBasinsRoots, 1000);

  cout << mapBasinsRoots.size() << " Basins Roots found" << endl;

  //on ne conserve que les bassins qui au plus de "clean_size" vertex

  //utiliser une mesure en mm2 et pas un nombre de sommets

  int clean_size = 50;

  TimeTexture<short> texRootsBottomClean(1, _mesh.vertex().size());

  for (uint i = 0; i < texRootsBottom[0].nItem(); i++)
    texRootsBottomClean[0].item(i) = texRootsBottom[0].item(i);

  cleanBasins(mapBasinsRoots, texRootsBottomClean, clean_size);
  cout << "after cleaning" << endl;
  cout << mapBasinsRoots.size() << " Basins Roots found" << endl;

  if (_save)
  {
    writeShortTexture("roots_bottom_clean.tex", texRootsBottomClean);
  }

  TimeTexture<short> texContourBasins(1, _mesh.vertex().size());

  map<int, set<int> > mapContourBasins;
  contourBasins(mapBasinsRoots, texRootsBottomClean, mapContourBasins,texContourBasins);
  if (_save)
  {
    writeShortTexture("contour_basins_roots.tex", texContourBasins);
  }

  cout << "compute probability map" << endl;
  TimeTexture<float> texProba(1, _mesh.vertex().size());
  TimeTexture<float> texProbaNorm(1, _mesh.vertex().size());
  computeProbabiltyMap(mapContourBasins, texContourBasins, texProba);
  normalizeProbabiltyMap(mapBasinsRoots, mapContourBasins,
      texContourBasins, texProba, texProbaNorm);

  if (_save)
  {
    writeFloatTexture("proba_roots.tex", texProba);
    writeFloatTexture("proba_roots_norm.tex", texProbaNorm);

    TimeTexture<short> texAutoThreshold(1, _mesh.vertex().size());
    //saveHistoProbabiltyMap(mapBasinsRoots, texRootsBottomClean,"histoRoots.txt");
    automaticThresholdDensityMap(mapBasinsRoots,texRootsBottomClean,texProbaNorm,texAutoThreshold,50);
    if (_save)
      writeShortTexture("threshold_auto.tex", texAutoThreshold);

    TimeTexture<short> texSulcalines(1, _mesh.vertex().size() );
    computeLongestPathBasins (texAutoThreshold, texSulcalines, mapBasinsRoots);
    if (_save)
      writeShortTexture("sulcalines.tex", texSulcalines);


  }
  cout << "done " << endl;


  //
  //  vector<float>::iterator itp;
  //  for ( itp=_proba.begin(); itp != _proba.end(); itp++)
  //  {
  //    cout << "proba value = " << *itp << endl;
  //
  //    cout << "threshold probability map" << endl;
  //    TimeTexture<short> texProbaThreshLat(1, _mesh.vertex().size() );
  //    texBinarizeF2S(texProbaNormLat, texProbaThreshLat, *itp ,0 ,1);
  //    TimeTexture<short> texProbaThreshLon(1, _mesh.vertex().size() );
  //    texBinarizeF2S(texProbaNormLon, texProbaThreshLon, *itp ,0 ,1);
  //    cout << "done " << endl;
  //
  //    std::ostringstream buff;
  //    buff<<(*itp);
  //    size_t found;
  //    string convert = buff.str();
  //    found=convert.find(".");
  //    convert.replace(found,1,",");
  //    string texname;
  //
  //    if (_save)
  //    {
  //      texname = "lat_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texProbaThreshLat);
  //
  //      texname = "lon_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texProbaThreshLon);
  //    }
  //
  //    cout << "compute the longest path in threshold probability map" << endl;
  //    cout << "latitude :" << endl;
  //    TimeTexture<short> texProbaSulcalinesLat(1, _mesh.vertex().size() );
  //    textureBin2Label(texInterRootsBasinsLat,texProbaThreshLat,texTemp);
  //    if (_save)
  //      {
  //      texname = "lat_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texTemp);
  //      }
  //
  //    computeLongestPathBasins (texTemp, texProbaSulcalinesLat, mapBasinsLat);
  //
  //    cout << "done " << endl;
  //    cout << "longitude :" << endl;
  //    TimeTexture<short> texProbaSulcalinesLon(1, _mesh.vertex().size() );
  //    textureBin2Label(texInterRootsBasinsLon,texProbaThreshLon,texTemp);
  //    if (_save)
  //      {
  //      texname = "lon_proba_thresh_" + convert + ".tex";
  //      writeShortTexture(texname.c_str(),texTemp);
  //      }
  //
  //    computeLongestPathBasins (texTemp, texProbaSulcalinesLon, mapBasinsLon);
  //    cout << "done " << endl;
  //
  //    texname = "lat_proba_lines_" + convert + ".tex";
  //    writeShortTexture(texname.c_str(),texProbaSulcalinesLat);
  //    texname = "lon_proba_lines_" + convert + ".tex";
  //    writeShortTexture(texname.c_str(),texProbaSulcalinesLon);
  //  }
}

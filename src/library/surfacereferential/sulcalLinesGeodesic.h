/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *    CEA/DSV/SHFJ
 *    4 place du General Leclerc
 *    91401 Orsay cedex
 *    France
 */



#ifndef AIMS_SULCUS_CORTICAL_GEO_H
#define AIMS_SULCUS_CORTICAL_GEO_H

#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfacegen.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshvoronoi.h>
#include <aims/scalespace/meshDiffuse.h>
#include <aims/distancemap/meshmorphomat.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/geodesicpath/geodesicPath.h>

#include <iostream>
#include <fstream>
#include <queue>
#include <map>

using namespace aims;
using namespace aims::meshdistance;
using namespace std;

class SulcalLinesGeodesic
{
  public:

    //Parametres pour l'execution
    string _adrMesh;

    string _adrRootsLat;
    string _adrRootsLon;

    string _adrCurv;
    string _adrGeodesicDepth;
    //string _adrBassinsDepthNorm;

    int _strain;
    int _extremeties_method;
    int _segmentation_method;

    AimsSurfaceTriangle _mesh;

//    string _adrBassins;
//       string _adrLines;
//
//       string _adrLatGeodesicOut;
//       string _adrLonGeodesicOut;
//
//
    std::vector<std::set<uint> > _neigh;
//
//    TimeTexture<short> _texProbaPath;

    TimeTexture<short> _rootsLon;
    TimeTexture<short> _rootsLat;
    TimeTexture<float> _texCurv;
    TimeTexture<float> _geoDepth;
    TimeTexture<float> _geoDepthNorm;



//    TimeTexture<float> _texCurvSeuil;
//
//    TimeTexture<short> _texBassins;
//    TimeTexture<float> _texBassinsDepthNorm;
//
//    TimeTexture<short> _texBassinsLat;
//    TimeTexture<short> _texBassinsLon;

//      geodesic::Mesh _meshSPc;
//      geodesic::Mesh _meshSP;
//
//      geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;
//
//      vector<double> _pointsSP;
//      vector<unsigned> _facesSP;

//    set<int> _listIndexVertexFill;
//    map<int,set<int> > _mapBassins;
//
//    set<int> _listIndexLon;
//    set<int> _listIndexLat;
//
//    map<int, set<int> > _mapConstraintLat;
//    map<int, set<int> > _mapConstraintLon;

    //Constructor
    SulcalLinesGeodesic( string & adrMesh,string & adrCurv, string & adrGeodesicDepth,
        string & adrRootsLon, string & adrRootsLat, int extremeties_method, int segmentation_method, int strain );

    ~SulcalLinesGeodesic();

    //public methods
    void run();

  private :
    void segmentationSulcalBasins (TimeTexture<float> &texIn,TimeTexture<short> &texBasins,map<int,set<int> > &mapBasins,float threshold, int close, int open);

    void writeShortTexture (string name,TimeTexture<short> &out);
    void writeFloatTexture (string name,TimeTexture<float> &out);

    void floodFillIter(int indexVertex, float newTextureValue,float oldTextureValue,TimeTexture<short> &texBasinsTemp, map<int,set<int> > &mapBasins);
    void texConnectedComponent(TimeTexture<short> &texBasins, map<int,set<int> > &mapBasins);
    void texBinarizeF2S(TimeTexture<float> &texIn, TimeTexture<short> &texOut, float threshold,int inf,int sup);
    void texBinarizeS2S(TimeTexture<short> &texIn, TimeTexture<short> &texOut, int threshold,int inf,int sup);

    void listRootsProjections(TimeTexture<short> &texBasins,set<int> &listIndexLat,set<int> &listIndexLon);
    void computeListLabelProjectionsBasins (TimeTexture<short> &roots, map<int,set<int> > &mapBasins,set<int> &listIndex, map<int,set<int> > &mapConstraint);
    void computeLongestPathBasins (TimeTexture<short> &roots, TimeTexture<short> &out, map<int,set<int> > &mapConstraint);

    void normalizeDepthMap (TimeTexture<float> &depth, TimeTexture<float> &depthNorm, map<int,set<int> > &mapBasins);

    void sulcalLinesExtract_projection(map<int,set<int> > &mapBasins, TimeTexture<short> &texBasins);

    void dilationRoots(TimeTexture<short> &texLatDil,TimeTexture<short> &texLonDil,int size);
    void interRootsDilBasins(TimeTexture<short> &texBasins,TimeTexture<short> &texDil,TimeTexture<short> &texInter);
    void cleanBasins(map<int,set<int> > &mapBasins,TimeTexture<short> &texBasins,int nbPoint);
    void contourBasins(map<int,set<int> > &mapBasins,TimeTexture<short> &texBasins,map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins);

    void sulcalLinesExtract_probability(map<int,set<int> > &mapBasins, TimeTexture<short> &texBasins);

};


#endif



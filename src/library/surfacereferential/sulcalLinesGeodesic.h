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

typedef std::vector<double> double_vec_t;

class Comp
{

  double_vec_t& _v;
 public:
   Comp(double_vec_t& v) : _v(v) {}
   bool operator()(double i, double j){
         return _v[i] < _v[j];
   }
};

class SulcalLinesGeodesic
{
  public:

    //Parametres pour l'execution
    string _adrMesh;
    string _adrRootsLat;
    string _adrRootsLon;
    string _adrCurv;
    string _adrGeodesicDepth;
    string _adrSaveFolder;
    string _adrRootsBottom;
    string _adrLabelBasins;
    string _adrLabelSulcalines;
    string _adrSulcalines;
    string _side;

    int _strain;
    int _extremeties_method;
    int _constraint_type;
    vector<float> _proba;
    bool _save;
    float _curv_thresh;
    float _clean_size;
    int _constraintValue;
    int _max_extremities;

    AimsSurfaceTriangle _mesh;
    std::vector<std::set<uint> > _neigh;

    TimeTexture<short> _rootsLon;
    TimeTexture<short> _rootsLat;
    TimeTexture<float> _texCurv;
    TimeTexture<float> _geoDepth;
    TimeTexture<float> _geoDepthNorm;
    TimeTexture<float> _geoCurvDepthNorm;

    //Constructor
    SulcalLinesGeodesic(string & adrMesh, string & adrCurv,
    string & adrGeodesicDepth, string & adrRootsLon, string & adrRootsLat,string & adrRootsBottom, string & adrLabelBasins, string & adrLabelSulcalines, string & adrSulcalines,
    int extremeties_method, int constraint_type, int strain,
    vector<float> proba, string saveFolder, float curv_thresh, string side, float clean_size, int constraintValue,int max_extremities);

    ~SulcalLinesGeodesic();

    //public methods
    void run();
    void probaMap();

    void writeShortTexture (string name,TimeTexture<short> &out);
    void writeFloatTexture (string name,TimeTexture<float> &out);
    void floodFillIter(int indexVertex, float newTextureValue,float oldTextureValue,TimeTexture<short> &texBasinsTemp, map<int,set<int> > &mapBasins);
    TimeTexture<short> texConnectedComponent(TimeTexture<short> &texBasins, map<int,set<int> > &mapBasins, int offset);

    void texBinarizeF2S(TimeTexture<float> &texIn, TimeTexture<short> &texOut, float threshold,int inf,int sup);
    void texBinarizeS2S(TimeTexture<short> &texIn, TimeTexture<short> &texOut, int threshold,int inf,int sup);
    void computeListLabelProjectionsBasins (TimeTexture<short> &roots, map<int,set<int> > &mapBasins,set<int> &listIndex, map<int,set<int> > &mapConstraint);

  private :
    void segmentationSulcalBasins (TimeTexture<float> &texIn,TimeTexture<short> &texBasins,map<int,set<int> > &mapBasins,float threshold, int close, int open);

    void computeRootsBottomMap(TimeTexture<short> &texBasins,TimeTexture<short> &texRootsBottom, float dist_max);

    void computeConstraintList(map<int,int> & listValues);

    void listRootsProjections(TimeTexture<short> &texBasins,set<int> &listIndexLat,set<int> &listIndexLon);
    void computeLongestPathBasins (TimeTexture<short> &roots, TimeTexture<short> &out, map<int,set<int> > &mapConstraint);

    void normalizeDepthMap (TimeTexture<float> &depth, TimeTexture<float> &depthNorm, map<int,set<int> > &mapBasins);

    void sulcalLinesExtract_projection(map<int,set<int> > &mapBasins, TimeTexture<short> &texBasins);

    void dilationRoots(TimeTexture<short> &texLatDil,TimeTexture<short> &texLonDil,int size);
    void interRootsDilBasins(TimeTexture<short> &texBasins,TimeTexture<short> &texDil,TimeTexture<short> &texInter);
    void cleanBasins(map<int,set<int> > &mapBasins,TimeTexture<short> &texBasins,map<int, vector<int> > &mapPolygonSetBasins,float min_area_size);
    void contourBasins(map<int,set<int> > &mapBasins,TimeTexture<short> &texBasins,map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins);
    void computeProbabiltyMap(map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins,TimeTexture<short> &texBasins, TimeTexture<float> &texProba);
    void computeMaximalProbabiltyMap(map<int, set<int> > &mapContourBasins, TimeTexture<short> &texContourBasins, TimeTexture<short> &texBasins, TimeTexture<float> &texProba);
    void normalizeProbabiltyMap(map<int,set<int> > &mapBasins, map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins,TimeTexture<float> &texProba,TimeTexture<float> &texProbaNorm);
    void normalizeMaximalProbabiltyMap(map<int,set<int> > &mapBasins, map<int,set<int> > &mapContourBasins,TimeTexture<short> &texContourBasins,TimeTexture<float> &texProba,TimeTexture<float> &texProbaNorm);
    void textureBin2Label(TimeTexture<short> &texLabel, TimeTexture<short> &texIn, TimeTexture<short> &texOut);

    //embc11
    void sulcalLinesExtract_probability(map<int,set<int> > &mapBasins, TimeTexture<short> &texBasins);
    //neuroimage11
    void sulcalLinesExtract_density(map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins);
    //miccai12
    void sulcalLinesExtract_maximal_density(map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins);

    void automaticThresholdDensityMap(map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins, TimeTexture<float> &texProbaNorm,TimeTexture<short> &texAutoThreshold,int nb_bin);
    void automaticThresholdMaximalDensityMap(map<int, set<int> > &mapBasins, TimeTexture<short> &texBasins, TimeTexture<float> &texProbaNorm,TimeTexture<short> &texAutoThreshold_begin,TimeTexture<short> &texAutoThreshold_middle,TimeTexture<short> &texAutoThreshold_end,int nb_bin);

    void vertexmap2polygonMap(map<int, set<int> > &mapVertexSetBasins, map<int, vector<int> > &mapPolygonSetBasins);

};


#endif

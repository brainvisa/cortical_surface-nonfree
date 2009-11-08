#ifndef SURF_BLOBS_H
#define SURF_BLOBS_H
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace carto;
using namespace std;

class SubjectData{
  public :
    string subject;
    AimsSurfaceTriangle mesh;
    TimeTexture<float> tex;
    TimeTexture<float> lat;
    TimeTexture<float> lon;    
};

class Blob{
  public :
    uint index;
//     uint group_index;
    uint parent;
//     string subject;
    float t;
    float scale;
    Point3df boundingbox_max;
    Point3df boundingbox_min;
    set<int> nodes_set;
    vector<int> nodes_list;
};

class SSBlob{
  public :
    uint index;
    uint graph_index;
    float t;
//     string label;
    string subject;
    float tmin;
    float tmax;
    set<int> representation;
    set<Blob *> blobs;
};

class SSBClique{
  public :
    SSBlob *ssb1;
    SSBlob *ssb2;
    float similarity;
    SSBClique(SSBlob *s1, SSBlob *s2, float sim){ssb1=s1; ssb2=s2; similarity=sim;}

};


#endif 

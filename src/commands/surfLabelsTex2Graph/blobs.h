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
    uint parent;
//     string subject;
    float t;
    float scale;
    Point3df boundingbox_max;
    Point3df boundingbox_min;
    set<uint> nodes_list;
};

class SSBlob{
  public :
    uint index;
    float t;
    string subject;
    float tmin;
    float tmax;
    set<uint> representation;
    set<Blob *> blobs;
};

// vector<Blob *> construireBlobs(Graph &sketch);
#endif 

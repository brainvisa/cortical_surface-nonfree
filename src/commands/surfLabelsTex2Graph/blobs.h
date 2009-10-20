#ifndef SURF_BLOBS_H
#define SURF_BLOBS_H
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace carto;
using namespace std;


class Blob{
  public :
    uint index;
    uint parent;
    string subject;
    float t;
    float scale;
    int rank;
    Point3df boundingbox_max;
    Point3df boundingbox_min;
    set<uint> nodes_list;
};

class SSBlob{
  public :
    uint index;
    set<Blob *> blobs;
};


#endif 

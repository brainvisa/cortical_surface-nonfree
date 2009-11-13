#ifndef SURF_BLOBS_H
#define SURF_BLOBS_H
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace carto;
using namespace std;



namespace surf{
  
  // Of course, there is no info about the corresponding mesh but let's say the user is supposed to know that...

  class Blob{
    public: 
    set<int> nodes;
  };
  
  class GreyLevelBlob: public Blob{
    public :
      uint index;
  //     uint group_index;
      uint parent;
  //     string subject;
      float t;
      float scale;      
      Point3df boundingbox_max;
      Point3df boundingbox_min;
  };
  
  class ScaleSpaceBlob: public Blob{
    public :
      uint index;
      uint graph_index;
      float t;
      int label;
      string subject;
      float tmin;
      float tmax;
      set<GreyLevelBlob *> blobs;
  };
  
  class SSBClique{
    public :
      ScaleSpaceBlob *ssb1;
      ScaleSpaceBlob *ssb2;
      float similarity;
      SSBClique(ScaleSpaceBlob *s1, ScaleSpaceBlob *s2, float sim){ssb1=s1; ssb2=s2; similarity=sim;}
  
  };
  
}


#endif 

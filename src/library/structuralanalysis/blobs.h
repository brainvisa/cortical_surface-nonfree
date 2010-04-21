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
            int index;
            int label;
            set<int> nodes;
            map<int, vector<float> > coordinates;

            AimsSurface<3, Void> getAimsMeshPatch ( AimsSurface<3, Void> &mesh, set<int> &nodes_list );
            AimsSurface<3, Void> getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh,
                Texture<float> &lat,
                Texture<float> &lon,
                float radius,
                set<int> &nodes_list );
            AimsSurface<3, Void> getAimsPatchOnAPlane  ( AimsSurface<3, Void> &mesh,
                Texture<float> &lat,
                Texture<float> &lon,
                float height,
                set<int> &nodes_list );
    };
    class ScaleSpaceBlob;
  

    class GreyLevelBlob: public Blob{
        public :
            float t;
            float scale;
            Point3df boundingbox_max;
            Point3df boundingbox_min;
            ScaleSpaceBlob *ssb_parent;

            AimsSurface<3, Void> getAimsMeshPatch ( AimsSurface<3, Void> &mesh, set<int> &nodes_list );
            AimsSurface<3, Void> getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh,
                                                        Texture<float> &lat,
                                                        Texture<float> &lon,
                                                        set<int> &nodes_list );
            AimsSurface<3, Void> getAimsPatchOnAPlane  ( AimsSurface<3, Void> &mesh,
                                                        Texture<float> &lat,
                                                        Texture<float> &lon,
                                                        set<int> &nodes_list );
            Point3df getBlobBarycenterOnASphere( );

    };

    class ScaleSpaceBlob: public Blob{
        public :
            float t;
            int label;
            string subject;
            float tmin;
            float tmax;
            set<GreyLevelBlob *> blobs;
            set<ScaleSpaceBlob *> topBlobs, bottomBlobs;
            AimsSurface<3, Void> getAimsMeshPatch ( AimsSurface<3, Void> &mesh, set<int> &nodes_list );
            AimsSurface<3, Void> getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh,
                Texture<float> &lat,
                Texture<float> &lon,
                set<int> &nodes_list );
                AimsSurface<3, Void> getAimsPatchOnAPlane  ( AimsSurface<3, Void> &mesh,
                    Texture<float> &lat,
                    Texture<float> &lon,
                    set<int> &nodes_list );

    };

    class SSBClique{
        public :
            ScaleSpaceBlob *ssb1;
            ScaleSpaceBlob *ssb2;
            float similarity;
            SSBClique(ScaleSpaceBlob *s1, ScaleSpaceBlob *s2, float sim){ssb1=s1; ssb2=s2; similarity=sim;}

    };

    class SSBBifurcation{
        public :
            set<ScaleSpaceBlob *> topBlobs;
            set<ScaleSpaceBlob *> bottomBlobs;
            string type;
            SSBBifurcation ( set<ScaleSpaceBlob *> &s1, set< ScaleSpaceBlob *> &s2, string _type){topBlobs = set<ScaleSpaceBlob *>(s1); bottomBlobs = set<ScaleSpaceBlob *>(s2); type = _type;}
          
    };
}

//##############################################################################

double getOverlapMeasure( Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap );

void filteringBlobs (  vector<surf::ScaleSpaceBlob *> & ssblobs,
                       vector<surf::GreyLevelBlob *> &filteredBlobs,
                       vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
                       Point3df bbmin2,
                       Point3df bbmax2 );

//##############################################################################

pair<Point2df, Point2df> getBoundingBox(set<int> &nodes_list, TimeTexture<float> &lat, TimeTexture<float> &lon);

pair<Point2df, Point2df> getBoundingBox(set<int> &nodes_list, map<int, float> &lat, map<int, float> &lon);

pair<Point2df, Point2df> getBoundingBox ( set<int> &nodes_list,
                                          map<int, vector<float> > &coordinates );

float compareBlobsScales(const surf::GreyLevelBlob *s1, const surf::GreyLevelBlob *s2);

struct ltBlobs
{
  bool operator()(const surf::GreyLevelBlob * s1, const surf::GreyLevelBlob * s2) const
  {
    return compareBlobsScales(s1, s2) < 0.0;
  }
};
#endif

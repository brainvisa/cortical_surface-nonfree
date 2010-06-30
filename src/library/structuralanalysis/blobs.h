#ifndef SURF_BLOBS_H
#define SURF_BLOBS_H
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace carto;
using namespace std;


enum typesRepresentation {
    RAW, SPHERE, FLAT, NODES_BARYCENTERS
};

namespace surf{

    // Of course, there is no info about the corresponding mesh but let's say the user is supposed to know that...

    class Blob{
        public:
            int index;
            set<int> nodes;
            map<int, vector<float> > coordinates;
            map<int, vector<float> > raw_coordinates;
            AimsSurface<3, Void> mesh;

            pair<Point2df, Point2df> get2DBoundingBox ( );

            void getAimsMesh (  AimsSurface<3, Void> &mesh,
                                float radius,
                                int representation_mode );
            void getAimsEllipsoid ( float abscissa, float depth, float height, float area );

            void getAimsMeshPatch ( AimsSurface<3, Void> &mesh,
                                    float radius = 1.0 );
            void getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh,
                float radius);
            void getAimsPatchOnAPlane  ( AimsSurface<3, Void> &mesh,
                float height );
            Blob(){}
            ~Blob(){}
    };
    class ScaleSpaceBlob;
  

    class GreyLevelBlob: public Blob {
        public :
            float t;
            float scale;
//             float area;
            // Spatial dispersion : to describe if the blob is isolated are 
//             float dispersion;
            // Coordinate defined along an anatomical (possibly defined by user) axis.
            // It is used to measure the alignment of structures along axes.
//             map<int, float> x_along_axis;

            Point3df boundingbox_max;
            Point3df boundingbox_min;
            ScaleSpaceBlob *ssb_parent;

            void getAimsMesh (  AimsSurface<3, Void> &mesh,
                                int representation_mode );
            void getAimsEllipsoid ( void );

            void getAimsMeshPatch ( AimsSurface<3, Void> &mesh );
            void getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh );
            void getAimsPatchOnAPlane  ( AimsSurface<3, Void> &mesh );
            Point3df getBlobBarycenterOnASphere( );
            Point3df getBlobBarycenter( );
            Point3df getBlobBarycenterFromMesh( );
            Point3df getBlobBarycenterOnAPlane( );

            pair<Point2df, Point2df> get2DBoundingBox ( );
            
            GreyLevelBlob(){}
            ~GreyLevelBlob(){}
            GreyLevelBlob( GreyLevelBlob *glb ) {
//                 area = glb->area;
                index = glb->index;
                nodes = set<int>(glb->nodes);
//                 x_along_axis = map<int, float>(glb->x_along_axis);
                coordinates = map<int, vector<float> >(glb->coordinates);
                raw_coordinates = map<int, vector<float> >(glb->raw_coordinates);
                ssb_parent = glb->ssb_parent;
                mesh = glb->mesh;
                t = glb->t;
                scale = glb->scale;
                for ( uint i = 0 ; i < 3 ; i ++ ) {
                    boundingbox_min[i] = glb->boundingbox_min[i];
                    boundingbox_max[i] = glb->boundingbox_max[i];
                }
            }
                

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

            ScaleSpaceBlob(){}
            ~ScaleSpaceBlob(){}
            ScaleSpaceBlob( ScaleSpaceBlob *ssb ) {
                index = ssb->index;
                label = ssb->label;
                subject = ssb->subject;
                nodes = ssb->nodes;
                coordinates = ssb->coordinates;
                raw_coordinates = ssb->raw_coordinates;
                t = ssb->t;
                tmin = ssb->tmin;
                tmax = ssb->tmax;
                blobs = ssb->blobs;
                topBlobs = ssb->topBlobs;
                bottomBlobs = ssb->bottomBlobs;
            }

            pair<Point2df, Point2df> get2DBoundingBox ( );


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

void computeBlobsDispersion( vector<surf::ScaleSpaceBlob *> & ssblobs );

double getOverlapMeasure( Point2df bbmin1, Point2df bbmax1, Point2df bbmin2, Point2df bbmax2, uint *no_overlap );

bool isInside2DBox( Point2df p1, Point2df bbmin, Point2df bbmax);

// void filteringBlobs (  vector<surf::ScaleSpaceBlob *> & ssblobs,
//                        vector<surf::GreyLevelBlob *> &filteredBlobs,
//                        vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
//                        Point2df bbmin2,
//                        Point2df bbmax2 );

void filteringBlobs (  vector<surf::ScaleSpaceBlob *> & ssblobs,
                        vector<surf::GreyLevelBlob *> &filteredBlobs,
                        vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
                        set< int > &nodes );

//##############################################################################


// pair<Point2df, Point2df> getBoundingBox(set<int> &nodes_list, map<int, float> &lat, map<int, float> &lon);

float compareBlobsScales(const surf::GreyLevelBlob *s1, const surf::GreyLevelBlob *s2);

struct ltBlobs
{
  bool operator()(const surf::GreyLevelBlob * s1, const surf::GreyLevelBlob * s2) const
  {
    return compareBlobsScales(s1, s2) < 0.0;
  }
};
#endif

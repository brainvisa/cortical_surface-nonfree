#ifndef SURF_BLOBS_H
#define SURF_BLOBS_H
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>



enum typesRepresentation {
    RAW, SPHERE, FLAT, NODES_BARYCENTERS
};

namespace surf{

    class Blob{
        public:
            int index;
            std::set<int> nodes;
            std::map<int, std::vector<float> > coordinates;
            std::map<int, std::vector<float> > raw_coordinates;
            AimsSurface<3, Void> mesh;

            std::pair<Point2df, Point2df> get2DBoundingBox ( );
            int getMaximumNode( Texture<float> &tex );

            void getAimsMesh (  AimsSurface<3, Void> &mesh );
            void getAimsSphereAtMaxNode (  Texture<float> &tex, float radius = 0.1 );
            void getAimsEllipsoid ( float abscissa, float height, float depth, float area );
            void moveMeshToSphericalAtlas ( float radius ) ;
            void moveMeshToPlaneAtlas ( float height ) ;
            void getNodesFromBlob( surf::Blob * blob);

            Point3df getBlobBarycenterOnASphere( );
            Point3df getBlobBarycenter( );
            Point3df getBlobBarycenterFromMesh( );
            Point3df getBlobBarycenterOnAPlane( );

            Blob(){}
            ~Blob(){}
    };
    class ScaleSpaceBlob;


    class GreyLevelBlob: public Blob {
        public :
            float t;
            float scale;
            int label;

            surf::ScaleSpaceBlob *ssb_parent;

//            void getAimsMesh (  AimsSurface<3, Void> &mesh );
            void getAimsEllipsoid ( void );
            void getAimsEllipsoidAtMaxNode (  Texture<float> &tex ) ;
            void moveMeshToSphericalAtlas ( void ) ;
            void moveMeshToPlaneAtlas ( void ) ;

            std::pair<Point2df, Point2df> get2DBoundingBox ( );

            GreyLevelBlob(){}
            ~GreyLevelBlob(){}
            GreyLevelBlob( GreyLevelBlob *glb ) {
                index = glb->index;
                nodes = std::set<int>(glb->nodes);
//                polygons = set::set<Point3di>(glb->polygons);
                coordinates = std::map<int, std::vector<float> >(glb->coordinates);
                raw_coordinates = std::map<int, std::vector<float> >(glb->raw_coordinates);
                ssb_parent = glb->ssb_parent;
                mesh = glb->mesh;
                t = glb->t;
                scale = glb->scale;
            }

    };

    class ScaleSpaceBlob: public Blob{
        public :
            float t;
            int label;
            std::string subject;
            float tmin;
            float tmax;
            std::set<GreyLevelBlob *> blobs;
            std::set<ScaleSpaceBlob *> topBlobs, bottomBlobs;

//            void getAimsMesh ( AimsSurface<3, Void> &mesh );


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
                mesh = ssb->mesh;
            }

            std::pair<Point2df, Point2df> get2DBoundingBox ( );


    };

    class SSBClique{
        public :
            ScaleSpaceBlob *ssb1;
            ScaleSpaceBlob *ssb2;
            float similarity;
            float distance;
            SSBClique(ScaleSpaceBlob *s1, ScaleSpaceBlob *s2, float dist ) {
                ssb1 = s1;
                ssb2 = s2;
                distance = dist;
                similarity = dist;
                //if (distance < threshold)
                    //similarity = 1.0;
                //else
                    //similarity = - 0.5 * (distance-3.0) +1.0;}
                    //similarity = - 1.0 / (threshold - alpha) * (distance - alpha) + 1.0;
                    //similarity = 
            }

    };

    class SSBBifurcation{
        public :
            std::set<ScaleSpaceBlob *> topBlobs;
            std::set<ScaleSpaceBlob *> bottomBlobs;
            std::string type;
            SSBBifurcation ( std::set<ScaleSpaceBlob *> &s1, std::set< ScaleSpaceBlob *> &s2, std::string _type){topBlobs = std::set<ScaleSpaceBlob *>(s1); bottomBlobs = std::set<ScaleSpaceBlob *>(s2); type = _type;}

    };
}

//##############################################################################

void computeBlobsDispersion( std::vector<surf::ScaleSpaceBlob *> & ssblobs );

double getOverlapMeasure( Point2df bbmin1, Point2df bbmax1, Point2df bbmin2, Point2df bbmax2, uint *no_overlap );

bool isInside2DBox( Point2df p1, Point2df bbmin, Point2df bbmax);

void filteringBlobs (  std::vector<surf::ScaleSpaceBlob *> & ssblobs,
        std::vector<surf::GreyLevelBlob *> &filteredBlobs,
        std::vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
        std::set< int > &nodes );

//##############################################################################



float compareBlobsScales(const surf::GreyLevelBlob *s1, const surf::GreyLevelBlob *s2);

struct ltBlobs
{
  bool operator()(const surf::GreyLevelBlob * s1, const surf::GreyLevelBlob * s2) const
  {
    return compareBlobsScales(s1, s2) < 0.0;
  }
};
#endif

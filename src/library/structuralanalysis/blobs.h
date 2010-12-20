#ifndef SURF_BLOBS_H
#define SURF_BLOBS_H
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <cortical_surface/structuralanalysis/subjectdata.h>




enum typesRepresentation {
    RAW, SPHERE, FLAT, NODES_BARYCENTERS
};

namespace surf {

    class Blob {
        public:
            int index;
            std::set<int> nodes;
            std::string subject;
            float t;

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
            float scale;
            //int label;
            surf::ScaleSpaceBlob *ssb_parent;

            void getAimsEllipsoid ( void );
            void getAimsEllipsoidAtMaxNode (  Texture<float> &tex ) ;
            void moveMeshToSphericalAtlas ( void ) ;
            void moveMeshToPlaneAtlas ( void ) ;

            std::pair<Point2df, Point2df> get2DBoundingBox ( );

            GreyLevelBlob(){}
            ~GreyLevelBlob(){}
            GreyLevelBlob( GreyLevelBlob *glb ) {
                index = glb->index;
                //label = glb->label;
                nodes = std::set<int>(glb->nodes);
                subject = glb->subject;
                nodes = std::set<int>(glb->nodes);
                coordinates = std::map<int, std::vector<float> >(glb->coordinates);
                raw_coordinates = std::map<int, std::vector<float> >(glb->raw_coordinates);
                ssb_parent = glb->ssb_parent;
                mesh = glb->mesh;
                t = glb->t;
                scale = glb->scale;
            }

    };

    class ScaleSpaceBlob: public Blob {
        public :
            float t;
            int label;
            std::set<float> scales;
            float tmin;
            float tmax;
            std::set<GreyLevelBlob *> blobs;
            std::set<ScaleSpaceBlob *> topBlobs, bottomBlobs;

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

    class Clique {
        public :
            ScaleSpaceBlob *ssb1;
            ScaleSpaceBlob *ssb2;
            float similarity;
            float distance;
            Clique ( ScaleSpaceBlob *s1, ScaleSpaceBlob *s2, float _distance, float _similarity ) {
                ssb1 = s1;
                ssb2 = s2;
                distance = _distance;
                similarity = _similarity;
            }

    };

    class Bifurcation {
        public :
            std::set<ScaleSpaceBlob *> topBlobs;
            std::set<ScaleSpaceBlob *> bottomBlobs;
            std::string type;
            Bifurcation ( std::set<ScaleSpaceBlob *> &s1, std::set< ScaleSpaceBlob *> &s2, std::string _type){topBlobs = std::set<ScaleSpaceBlob *>(s1); bottomBlobs = std::set<ScaleSpaceBlob *>(s2); type = _type;}

    };


}

//##############################################################################




#endif


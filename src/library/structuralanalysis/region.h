#ifndef REGION_H_
#define REGION_H_

#include <cortical_surface/structuralanalysis/subjectdata.h>
#include <cortical_surface/structuralanalysis/blobs.h>


namespace surf {
    class Region {
        public:
//            AimsSurface<3,Void> *mesh;
//            Texture<float> *lat, *tex, *lon;
            SubjectData *subject;
            std::vector<uint> nodes; //gyrusVertices
            std::vector<uint> corres;
            std::vector<uint> triangles;
            std::map<unsigned, std::set<std::pair<unsigned,float> > > weightLapl;

            AimsSurface<3,Void> regionMesh;
            Texture<float> regionTex, regionLat, regionLon;

            std::vector<surf::GreyLevelBlob *> blobs;
            std::vector<surf::ScaleSpaceBlob *> ssblobs;

            Region () {}
            Region ( SubjectData *subject,
                     std::vector<uint> &_nodes ) ;

            void getGlobalFromLocalNodes( surf::Blob *blob );

            Texture<float> getGlobalTexture( float background_value );
            Texture<float> getGlobalFromLocalTexture( Texture<float> &corticalTex, float background_value = 0.0 );
            Texture<float> getLocalFromGlobalTexture ( Texture<float> &corticalTex );
    };
}


#endif /*REGION_H_*/

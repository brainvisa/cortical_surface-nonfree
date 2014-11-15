#ifndef SURF_PRIMALSKETCH_OPERATIONS_H_
#define SURF_PRIMALSKETCH_OPERATIONS_H_


#include <cortical_surface/structuralanalysis/texturetoblobs.h>

#include <aims/primalsketch/scalespace.h>
#include <aims/primalsketch/primalSketch.h>

#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <aims/primalsketch/primalSketchUtil.h>

namespace TextureToBlobs
{

    void PrimalSketchRegionMode ( //std::vector<surf::GreyLevelBlob *> &blobs,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
            surf::Region &region,
            SubjectData &regionData,
            std::string scaleSpacePath,
            std::string blobsPath,
            bool recover = false,
            float scale_max = -1.0 );

//    void PrimalSketchRegionMode ( std::vector<surf::ScaleSpaceBlob *> &ssblobs,
//            surf::Region &region,
//            SubjectData &regionData,
//            std::string scaleSpacePath,
//            std::string blobsPath,
//            bool recover = false,
//            float scale_max = -1.0 ) {
//        //std::vector<surf::GreyLevelBlob *> blobs;
//        PrimalSketchRegionMode ( ssblobs, region, regionData, scaleSpacePath, blobsPath, recover, scale_max );
//    }

    void PrimalSketchGlobalMode ( //std::vector<surf::GreyLevelBlob *> &blobs,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
            SubjectData &subject,
            std::string scaleSpacePath,
            std::string blobsPath,
            bool recover = false,
            float scale_max = -1.0 );

//    void PrimalSketchGlobalMode ( std::vector<surf::ScaleSpaceBlob *> &ssblobs,
//            SubjectData &subject,
//            std::string scaleSpacePath,
//            std::string blobsPath,
//            bool recover = false,
//            float scale_max = -1.0 ) {
//        //std::vector<surf::GreyLevelBlob *> blobs;
//        PrimalSketchGlobalMode ( ssblobs, subject, scaleSpacePath, blobsPath, recover, scale_max );
//    }

    void PrimalSketch ( //std::vector<surf::GreyLevelBlob *> &blobs,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
            SubjectData &subject,
            std::string scaleSpacePath,
            std::string blobsPath,
            bool recover = false,
            float scale_max = -1.0 ) {
        PrimalSketchGlobalMode ( ssblobs, subject, scaleSpacePath, blobsPath, recover, scale_max );
    }

//    void PrimalSketch ( std::vector<surf::ScaleSpaceBlob *> &ssblobs,
//            SubjectData &subject,
//            std::string scaleSpacePath,
//            std::string blobsPath,
//            bool recover = false,
//            float scale_max = -1.0 ) {
//        //std::vector<surf::GreyLevelBlob *> blobs;
//        PrimalSketchGlobalMode ( ssblobs, subject, scaleSpacePath, blobsPath, recover, scale_max );
//    }

    void PrimalSketch ( SubjectData &subject,
                        //std::vector<surf::GreyLevelBlob *> &blobs,
                        std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                        aims::ScaleSpace<AimsSurface<3, Void>, Texture<float> > *ss,
                        TimeTexture<float> &blobs_texture,
                        float scale_max = -1.0);

    void GreyLevelBlobsFromTexture ( SubjectData &subject,
            //std::vector<surf::GreyLevelBlob *> &blobs,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
            std::string blobsPath ) ;

    void getBlobsFromPrimalSketch ( SubjectData & subject,
            aims::PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch,
            //std::vector<surf::GreyLevelBlob *> &blobs,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs ) ;

}

#endif /*SURF_PRIMALSKETCH_OPERATIONS_H_*/

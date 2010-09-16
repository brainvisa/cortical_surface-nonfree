#ifndef TEXTURETOBLOBS_H_
#define TEXTURETOBLOBS_H_

#include <aims/primalsketch/scalespace.h>
#include <aims/primalsketch/primalSketch.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>

#include <cortical_surface/structuralanalysis/blobs.h>
#include <cortical_surface/structuralanalysis/region.h>



enum filteringMode {
    GYRUS, AXIS, NO_FILTER, NO_REPRESENTATION
};

enum blobsMode {
    PRIMALSKETCH, LABELS
};


class TextureToBlobs {

    public:

        TextureToBlobs () { }

        static void PrimalSketchRegionMode ( vector<surf::GreyLevelBlob *> &blobs,
                vector<surf::ScaleSpaceBlob *> &ssblobs,
                surf::Region &region,
                SubjectData &regionData,
                string scaleSpacePath,
                string blobsPath,
                bool recover = false,
                float scale_max = -1.0 );

        static void PrimalSketchGlobalMode ( vector<surf::GreyLevelBlob *> &blobs,
                vector<surf::ScaleSpaceBlob *> &ssblobs,
                SubjectData &subject,
                string scaleSpacePath,
                string blobsPath,
                bool recover = false,
                float scale_max = -1.0 );

        static void PrimalSketch ( vector<surf::GreyLevelBlob *> &blobs,
                vector<surf::ScaleSpaceBlob *> &ssblobs,
                SubjectData &subject,
                string scaleSpacePath,
                string blobsPath,
                bool recover = false,
                float scale_max = -1.0 ) {
           PrimalSketchGlobalMode ( blobs, ssblobs, subject, scaleSpacePath, blobsPath, recover, scale_max );
        }

        static void PrimalSketch ( SubjectData &subject,
                            std::vector<surf::GreyLevelBlob *> &blobs,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            ScaleSpace<AimsSurface<3, Void>, Texture<float> > *ss,
                            TimeTexture<float> &blobs_texture,
                            float scale_max = -1.0);

        static void GreyLevelBlobsFromTexture ( SubjectData &subject,
                std::vector<surf::GreyLevelBlob *> &blobs,
                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                std::string blobsPath ) ;

        static void getBlobsFromPrimalSketch ( SubjectData & subject,
                aims::PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch,
                std::vector<surf::GreyLevelBlob *> &blobs,
                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                bool initNull = true ) ;

        static std::set<int> getFilteringNodes( SubjectData & subject );

        static void filterBlobs ( SubjectData & subject,
        		std::vector<surf::ScaleSpaceBlob *> & ssblobs,
                std::vector<surf::GreyLevelBlob *> & filteredBlobs,
                std::vector<surf::ScaleSpaceBlob *> & filteredSsblobs ) ;

//        static void RecoverBlobsFromIndivGraph( Graph *graph,
//                                    std::vector<surf::GreyLevelBlob *> &blobs,
//                                    std::vector<surf::ScaleSpaceBlob *> &ssblobs,
//                                    bool initNull = true);

        static void getGreyLevelBlobsFromIndividualGraph ( Graph *graph,
                                    SubjectData &subject,
                                    std::vector<surf::GreyLevelBlob *> &blobs,
                                    bool initNull );

        static void getScaleSpaceBlobsFromIndividualGraph ( Graph *graph,
//                                    std::string subject_id,
                                    std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                    std::map<int, std::set<int> > &listGLBindices,
                                    bool initNull );

        static void RecoverBlobsFromIndividualGraph( Graph *graph,
                SubjectData &subject,
                vector<surf::GreyLevelBlob *> &blobs,
                vector<surf::ScaleSpaceBlob *> &ssblobs,
                bool initNull = true);

        static void RecoverBlobsFromGLBOnly( Graph *graph,
                                    SubjectData &subject,
                                    vector<surf::GreyLevelBlob *> &blobs,
                                    vector<surf::ScaleSpaceBlob *> &ssblobs,
                                    bool initNull = true);
        
        static void AimsGraph( Graph *graph,
						SubjectData & subject,
						std::vector<surf::Blob *> & blobs ) ;
        
        static void AimsGraph( Graph *graph,
                                SubjectData & subject,
                                std::vector<surf::GreyLevelBlob *> & blobs,
                                std::vector<surf::ScaleSpaceBlob *> & ssblobs ) ;

        static void AimsGroupGraph( Graph *graph,
                                std::map<std::string, SubjectData *> data,
                                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                std::vector<surf::SSBClique> &cliques );

        static void ReadAimsGroupGraph ( Graph &graph,
                GroupData &data,
                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                std::vector<surf::SSBClique> &cliques,
                std::vector<Vertex *> &listVertex );



};

#endif /*TEXTURETOBLOBS_H_*/

#ifndef TEXTURETOBLOBS_H_
#define TEXTURETOBLOBS_H_

#include <cortical_surface/structuralanalysis/region.h>

enum representationMode {
    CORTICAL_PATCHES, SPHERES, NONE
};

enum graphMode {
    DEFAULT, NO_SCALESPACEBLOBS_MESHES_AND_NO_RELATIONS_MESHES, NO_SCALESPACEBLOBS_MESHES
};

namespace TextureToBlobs {

    std::set<int> getFilteringNodes( SubjectData & subject );

    void getGreyLevelBlobsFromGraph ( Graph *graph,
                                SubjectData &subject,
                                std::vector<surf::GreyLevelBlob *> &blobs,
                                bool initNull );

    void getScaleSpaceBlobsFromGraph ( Graph *graph,
                                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                std::map< std::string, std::map<int, std::set<int> > > &listGLBindices,
                                bool initNull );

    void BlobsFromLabelTexture ( std::vector<surf::Blob *> &blobs,
                                 SubjectData &subject );

    void BlobsFromLabelTextureGlobalMode ( std::vector<surf::Blob *> &blobs,
                                           SubjectData &subject );

    void BlobsFromLabelTextureRegionMode ( std::vector<surf::Blob *> &blobs,
                                           surf::Region &region,
                                           SubjectData &regionData );

    surf::GreyLevelBlob * findBlob ( const std::vector<surf::GreyLevelBlob *> &blobs,
                                    std::string subject_id,
                                    int index ) ;
    surf::ScaleSpaceBlob * findBlob ( const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                    std::string subject_id,
                                    int index ) ;
    int findBlobIndex ( const std::vector<surf::GreyLevelBlob *> &blobs,
                                    std::string subject_id,
                                    int index ) ;
    int findBlobIndex ( const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                    std::string subject_id,
                                    int index ) ;    
    
    void RecoverBlobsFromGraph( Graph *graph,
            SubjectData &subject,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
            bool initNull = true);

    void RecoverBlobsFromGLBOnly( Graph *graph,
                                SubjectData &subject,
                                std::vector<surf::GreyLevelBlob *> &blobs,
                                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                bool initNull = true);

    void AimsGraph( Graph *graph,
					SubjectData & subject,
					const std::vector<surf::Blob *> & blobs,
					int representation_mode = CORTICAL_PATCHES) ;


    void AimsGraph ( Graph *graph,
                     SubjectData & subject,
                     const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                     int graph_mode = DEFAULT,
                     int representation_mode = NONE );

    void AimsGroupGraph ( Graph *graph,
                            std::map<std::string, SubjectData *> data,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            std::vector<surf::Clique> &cliques,
                            bool buildAndStoreRelationsMeshes = false );

    void ReadAimsGroupGraph ( Graph &graph,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
            std::vector<surf::Clique> &cliques );

    void DestroyBlobs ( std::vector<surf::ScaleSpaceBlob *> &ssblobs ) ;
    void DestroyBlobs ( std::vector<surf::Blob *> &blobs ) ;

    std::vector<uint> getClustersListsFromGLB ( std::vector< surf::GreyLevelBlob *> &blobs,
                                                  GroupData &data,
                                                  float clustering_distance_threshold );

    void buildBlobsFromClustersLists ( std::vector< surf::GreyLevelBlob *> &blobs,
                                              GroupData & data,
                                              std::vector<uint> &clusters,
                                              std::vector<surf::ScaleSpaceBlob *> &clusteredSsblobs,
                                              float clustering_distance_threshold = -1.0,
                                              std::string outputTextFile = "/tmp/blobsCountTable.py",
                                              bool uniqueGLB = false,
                                              int representation_mode = SPHERES ) ;

    double getOverlapMeasure( Point2df bbmin1, Point2df bbmax1, Point2df bbmin2, Point2df bbmax2, uint *no_overlap );

    bool isInside2DBox( Point2df p1, Point2df bbmin, Point2df bbmax);

    void filteringBlobs (  std::vector<surf::ScaleSpaceBlob *> & ssblobs,
            std::vector<surf::GreyLevelBlob *> &filteredBlobs,
            std::vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
            std::set< int > &nodes );
}

#endif /*TEXTURETOBLOBS_H_*/

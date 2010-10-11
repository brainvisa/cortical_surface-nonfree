#ifndef TEXTURETOBLOBS_H_
#define TEXTURETOBLOBS_H_

#include <cortical_surface/structuralanalysis/region.h>

enum filteringMode {
    GYRUS, AXIS, NO_FILTER, NO_REPRESENTATION
};

enum blobsMode {
    PRIMALSKETCH, LABELS
};


namespace TextureToBlobs {

    std::set<int> getFilteringNodes( SubjectData & subject );

    void filterBlobs ( SubjectData & subject,
    		std::vector<surf::ScaleSpaceBlob *> & ssblobs,
            std::vector<surf::GreyLevelBlob *> & filteredBlobs,
            std::vector<surf::ScaleSpaceBlob *> & filteredSsblobs ) ;




    
    void getGreyLevelBlobsFromIndividualGraph ( Graph *graph,
                                SubjectData &subject,
                                std::vector<surf::GreyLevelBlob *> &blobs,
                                bool initNull );

    void getScaleSpaceBlobsFromIndividualGraph ( Graph *graph,
                                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                std::map<int, std::set<int> > &listGLBindices,
                                bool initNull );

    void RecoverBlobsFromIndividualGraph( Graph *graph,
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
					const std::vector<surf::Blob *> & blobs ) ;
    
    void AimsGraph ( Graph *graph,
                            SubjectData & subject,
                            //const std::vector<surf::GreyLevelBlob *> &blobs,
                            const std::vector<surf::ScaleSpaceBlob *> &ssblobs );

    void AimsGroupGraph( Graph *graph,
                            std::map<std::string, SubjectData *> data,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            std::vector<surf::SSBClique> &cliques );

    void ReadAimsGroupGraph ( Graph &graph,
            GroupData &data,
            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
            std::vector<surf::SSBClique> &cliques,
            std::vector<Vertex *> &listVertex );

    void DestroyBlobs ( std::vector<surf::ScaleSpaceBlob *> &ssblobs ) ;
    
    std::vector<uint> getClustersListsFromGLB ( std::vector< surf::GreyLevelBlob *> &blobs, 
                                                  GroupData &data,
                                                  float clustering_distance_threshold );
    
    void buildBlobsFromClustersLists ( std::vector< surf::GreyLevelBlob *> &blobs, 
                                              GroupData & data,
                                              std::vector<uint> &clusters,
                                              std::vector<surf::ScaleSpaceBlob *> &clusteredSsblobs,
                                              float clustering_distance_threshold = -1.0) ;

    void computeBlobsDispersion( std::vector<surf::ScaleSpaceBlob *> & ssblobs );

    double getOverlapMeasure( Point2df bbmin1, Point2df bbmax1, Point2df bbmin2, Point2df bbmax2, uint *no_overlap );

    bool isInside2DBox( Point2df p1, Point2df bbmin, Point2df bbmax);

    void filteringBlobs (  std::vector<surf::ScaleSpaceBlob *> & ssblobs,
            std::vector<surf::GreyLevelBlob *> &filteredBlobs,
            std::vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
            std::set< int > &nodes );
}

#endif /*TEXTURETOBLOBS_H_*/

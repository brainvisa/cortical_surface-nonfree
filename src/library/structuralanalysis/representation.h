#ifndef SURF_REPRESENTATION_H
#define SURF_REPRESENTATION_H

#include <aims/mesh/texture.h>
#include <cortical_surface/structuralanalysis/blobs.h>


float compareSurfBlobsScales ( const surf::GreyLevelBlob *s1, const surf::GreyLevelBlob *s2 );

struct ltSurfBlobs
{
    bool operator()(const surf::GreyLevelBlob * s1, const surf::GreyLevelBlob * s2) const
    {
        return compareSurfBlobsScales(s1, s2) < 0.0;
    }
};

//AimsSurfaceTriangle getFlatMap(vector<vector<int> > &nodes_lists, TimeTexture<float> &lat, TimeTexture<float> &lon, TimeTexture<float> &tex);

AimsSurfaceTriangle getLabelObjectsOnASphere( TimeTexture<short> &tex,
                                AimsSurface<3,Void> &mesh,
                                Texture<float> &lat,
                                Texture<float> &lon,
                                std::vector< std::set<int> > &nodes_lists);

//AimsSurfaceTriangle getLabelObjectsOnAMesh( TimeTexture<short> &tex,
//                                AimsSurface<3,Void> &mesh,
//                                vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getG2GRelationsMeshes ( std::vector< std::pair<surf::GreyLevelBlob *, surf::GreyLevelBlob *> > &blobsPairs,
                                            int representation_mode = SPHERE );

AimsSurfaceTriangle getBifurcationRelationsMeshes ( std::vector< std::pair< surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> > & bifurcPairs,
                                                    int representation_mode = SPHERE );

AimsSurfaceTriangle getB2BRelationsMeshes ( std::vector<surf::Clique> &cliques, int representation_mode = SPHERE  ) ;

#endif


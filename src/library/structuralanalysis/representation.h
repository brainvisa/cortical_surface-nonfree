#ifndef SURF_REPRESENTATION_H
#define SURF_REPRESENTATION_H
#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>
#include <cortical_surface/structuralanalysis/blobs.h>




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

AimsSurfaceTriangle getB2BRelationsMeshes ( std::vector<surf::SSBClique> &cliques, int representation_mode = SPHERE  ) ;

#endif


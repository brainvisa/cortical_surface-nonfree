#ifndef SURF_REPRESENTATION_H
#define SURF_REPRESENTATION_H
#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>
#include <cortical_surface/structuralanalysis/blobs.h>



using namespace aims;
using namespace carto;
using namespace std;


AimsSurfaceTriangle getBlobsMeshesEllipsoid ( vector<surf::GreyLevelBlob *> &blobs );

AimsSurfaceTriangle getBlobsMeshesSphericalAtlas ( vector<surf::GreyLevelBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh );

//AimsSurfaceTriangle getBlobsMeshesSpherical ( vector<surf::ScaleSpaceBlob *> &blobs,
//                                     AimsSurface<3, Void> &mesh );

AimsSurfaceTriangle getBlobsMeshes2DAtlas ( vector<surf::GreyLevelBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh );

AimsSurfaceTriangle getBlobsMeshesFromMesh( vector<surf::GreyLevelBlob *> &blobs,
                                    AimsSurface<3, Void> &mesh );

//AimsSurfaceTriangle getBlobsMeshesFromMesh( vector<surf::ScaleSpaceBlob *> &blobs,
//                                    AimsSurface<3, Void> &mesh );

//AimsSurfaceTriangle getFlatMap(vector<vector<int> > &nodes_lists, TimeTexture<float> &lat, TimeTexture<float> &lon, TimeTexture<float> &tex);

AimsSurfaceTriangle getLabelObjectsOnASphere( TimeTexture<short> &tex,
                                AimsSurface<3,Void> &mesh,
                                Texture<float> &lat,
                                Texture<float> &lon,
                                vector<set<int> > &nodes_lists);

//AimsSurfaceTriangle getLabelObjectsOnAMesh( TimeTexture<short> &tex,
//                                AimsSurface<3,Void> &mesh,
//                                vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getG2GRelationsMeshes ( vector< pair<surf::GreyLevelBlob *, surf::GreyLevelBlob *> > &blobsPairs,
                                            int representation_mode = SPHERE );

AimsSurfaceTriangle getBifurcationRelationsMeshes ( vector< pair< surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> > & bifurcPairs,
                                                    int representation_mode = SPHERE );

#endif


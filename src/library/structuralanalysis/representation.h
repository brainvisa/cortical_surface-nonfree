#ifndef SURF_REPRESENTATION_H
#define SURF_REPRESENTATION_H
#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>
#include <cortical_surface/structuralanalysis/blobs.h>



using namespace aims;
using namespace carto;
using namespace std;


AimsSurfaceTriangle getBlobsSphericalMeshes ( vector<surf::GreyLevelBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     Texture<float> &lat,
                                     Texture<float> &lon,
                                     vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getBlobsSphericalMeshes ( vector<surf::ScaleSpaceBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     Texture<float> &lat,
                                     Texture<float> &lon,
                                     vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getBlobs2DMeshes ( vector<surf::GreyLevelBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     Texture<float> &lat,
                                     Texture<float> &lon,
                                     vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getBlobsMeshes( vector<surf::GreyLevelBlob *> &blobs,
                                    AimsSurface<3, Void> &mesh,
                                    vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getBlobsMeshes( vector<surf::ScaleSpaceBlob *> &blobs,
                                    AimsSurface<3, Void> &mesh,
                                    vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getFlatMap(vector<vector<int> > &nodes_lists, TimeTexture<float> &lat, TimeTexture<float> &lon, TimeTexture<float> &tex);

AimsSurfaceTriangle getLabelObjectsOnASphere( TimeTexture<short> &tex,
                                AimsSurface<3,Void> &mesh,
                                Texture<float> &lat,
                                Texture<float> &lon,
                                vector<set<int> > &nodes_lists);

AimsSurfaceTriangle getLabelObjectsOnAMesh( TimeTexture<short> &tex,
                                AimsSurface<3,Void> &mesh,
                                vector<set<int> > &nodes_lists);                                

#endif


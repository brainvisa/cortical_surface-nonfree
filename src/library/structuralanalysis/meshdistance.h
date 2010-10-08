#ifndef AIMS_RELAXMESHDISTANCE_H
#define AIMS_RELAXMESHDISTANCE_H


#include <aims/mesh/surface.h>


//vector<map<uint,float> > CalculeCarteDistances(AimsSurfaceTriangle mesh, set<uint> nodes,float dist_thresh);
// vector<map<uint,float> > CalculeDistancesBlob(AimsSurfaceTriangle mesh, set<uint> nodes, vector<uint> &sites);
//map<float, vector<pair<float, uint> > > getAlternateMesh(AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &longit);
//map<uint,float> getDistMap( AimsSurfaceTriangle *mesh,  map<unsigned, set<unsigned> >    &neighbours,  int dep, float dist_thresh);

std::map<uint, float> LocalMeshDistanceMap ( AimsSurface<3,Void> *mesh,
                                             const std::vector< std::set<unsigned> >    &neighbours,
                                             uint max_node,
                                             float ind ) ;

#endif


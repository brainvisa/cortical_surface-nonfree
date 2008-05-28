#ifndef AIMS_POINT_DISTANCE_H
#define AIMS_POINT_DISTANCE_H

#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>

namespace aims
{

// goedesic distance between points  on a mesh

class MeshPointDistance
{
public:
     
     MeshPointDistance(AimsSurfaceTriangle mesh) : 
     _mesh(mesh) {computeNeighbours();}
     inline void computeNeighbours() {_neigh= SurfaceManip::surfaceNeighbours(_mesh);}
     float compute(uint p1, uint p2);  // distance between nodes i and j
     
private:
     AimsSurfaceTriangle _mesh;
     std::vector<std::set<uint> >  _neigh; 

};  //fin du namespace

}

#endif

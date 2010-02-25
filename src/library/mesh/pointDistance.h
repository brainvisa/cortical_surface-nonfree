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
};

class MeshPointNeighborhoodFromDistance
{
public:

	 MeshPointNeighborhoodFromDistance(AimsSurfaceTriangle mesh) :
     _mesh(mesh) {computeNeighbours();}
     inline void computeNeighbours() {_neigh= SurfaceManip::surfaceNeighbours(_mesh);}
     std::set<uint> compute(uint node, float distance);

     void includeNeighbors(uint ind);

private:
     AimsSurfaceTriangle _mesh;
     std::vector<std::set<uint> >  _neigh;
     std::set<uint> _liste;
     uint _start;
     float _distance;
};

//fin du namespace

}

#endif

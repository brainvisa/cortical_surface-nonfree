/*
 *  This software and supporting documentation were developed by
 *  	CNRS, labo LSIS, UMR6168, Marseille, France
 *      olivier.coulon@univmed.fr
 */



#ifndef AIMS_CONFORMALMAPPING_H
#define AIMS_CONFORMALMAPPING_H

#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>


// This class implement the conformal mapping of a mesh to 
// the unit sphere as defined in Gu et al., IEEE Tr on Med Im.
//, vol23, no8, pp 949-958, 2004
//
// the idea is to minimize the harmonic energy and add constraints
// such as zero-mass center and maybe landmarks


namespace aims
{
     class ConformalMapping
     {
          public:
               ConformalMapping(std::string adr_mesh, float dt, float dE) : 
                    _dt(dt), _dE(dE)
                    {
                         Reader < AimsSurfaceTriangle > r(adr_mesh);
                         r.read( _mesh );
                         _nv=_mesh.vertex().size();
                         _tuetteMap=std::vector<Point3df>(_nv);
                         _conformalMap=std::vector<Point3df>(_nv);
                         _neigh= SurfaceManip::surfaceNeighbours(_mesh);
                         ComputeEdges();
                         ComputeHarmonicCoefficients();
                    }
               ConformalMapping(AimsSurfaceTriangle mesh, float dt, float dE) : 
                    _mesh(mesh), _dt(dt), _dE(dE)
                    {
                         _nv=_mesh.vertex().size();
                         _tuetteMap=std::vector<Point3df>(_nv);
                         _conformalMap=std::vector<Point3df>(_nv);
                         _neigh= SurfaceManip::surfaceNeighbours(_mesh);
                         ComputeEdges();
                         ComputeHarmonicCoefficients();
                    }
                    
               void ComputeEdges();
               void ComputeHarmonicCoefficients();
               void TuetteMap();
               void ConformalMap(std::string adr_tuette);
               
               AimsSurfaceTriangle GetTuetteMapping();
               AimsSurfaceTriangle GetConformalMapping(std::string adr_tuette);
               float TuetteEnergy();
               float HarmonicEnergy();
               void WriteConformalCoordinates(std::string adr_coord);

          private:
               AimsSurfaceTriangle _mesh;
               float _dt;
               float _dE;
               
               uint _nv;
               std::vector<std::set<uint> >  _neigh;
               std::vector<AimsVector<uint,2> > _edges;
               std::vector<std::set<std::pair<uint, uint> > > _edgeMap;
               
               std::vector<float> _harmonicCoef;
               std::vector<Point3df> _tuetteMap;
               std::vector<Point3df> _conformalMap;
               std::vector<float> _conformalLon;
               std::vector<float> _conformalLat;
               
               void buildEdgeMap();
               uint whichEdge(uint v1, uint v2);
               Point3df getCenterOfMass(std::vector<Point3df> tmpMap);
               void ComputeConformalCoordinates();
     };
}

#endif

/*
 *  This software and supporting documentation were developed by
 *  	CNRS, labo LSIS, UMR6168, Marseille, France
 *      olivier.coulon@univmed.fr
 */



#include <aims/mesh/texture.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surfaceOperation.h>
#include <cortical_surface/surfacereferential/conformalMapping.h>


namespace aims
{

//--- a simple utility that computes edges of the mesh ---

     void ComputeEdges()
     {
          uint i,j;
          std::vector<std::set<uint> >  neigh= SurfaceManip::surfaceNeighbours(_mesh);
          for (i=0; i< _nv; i++)
          {
               std::set<uint> vois=neigh[i];
               std::set<uint>::iterator itVois=vois.begin();
               for ( ; itVois!=vois.end(); ++itVois)
               {
                    j=(*itVois);
                    if (j>i)
                    {
                         AimsVector<uint, 2> edge(i,j);
                         _edges.push_back(edge);
                    }
          }
     }

//------------------- Tuette mapping ---------------------------
// this is the first step to computing the conformal mapping
//--------------------------------------------------------------

     void ConformalMapping::TuetteMap()
     {
          uint i;
          float E;
          _mesh.updateNormals();
          std::vector(Point3df) norm=_mesh.normal();
          for (i=0; i<_nv; i++)
          {
               _tuetteMap[i]=norm[i];
          }
          deltaE=10000;
          while (deltaE >= _dE) do
          {
               //iteration sur les edges
          
          
          
          
          }

     }

//------------------------- Conformal Mapping ----------------------
// the actual conformal mapping with zero gravisty center constraint
// initialized with the Tuette Mapping
//------------------------------------------------------------------

     void ConformalMapping::ConformalMap()
     {
     }
     
     



//------------------------------------------------------------------

     AimsSurfaceTriangle GetConformalMapping()
     {
          AimsSurfaceTriangle result;
          result.vertex()=_conformalMap;
          result.polygon()=_mesh.polygon();
          return(result);
     }

     AimsSurfaceTriangle GetTuetteMapping()
     {
          AimsSurfaceTriangle result;
          result.vertex()=_tuetteMap;
          result.polygon()=_mesh.polygon();
          return(result);
     }



}
#ifndef AIMS_SHORTESTPATH_H
#define AIMS_SHORTESTPATH_H

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <limits>

// #include "cathier/aims_wrap.h"
// #include "cathier/triangle_mesh_geodesic_map.h"
// #include "cathier/math_functions.h"

namespace aims
{


// cette classe calcule le plus court chemin entre deux points d'un maillage,
// contraint à passer dans un ensemble de noeuds défini par une texture et une valeur
// L'algo utilise une carte de distance par fast marching sur la surface et un backtracking tout con
// a l'intérieur de l'ensemble de noeuds.

template<typename Val> class GraphPath
{

     public:
     
     GraphPath() {longueur =0;}
     float getLongueur(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr);
     TimeTexture<Val> process(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr);
     
     private:

     float longueur;
};

//------------------------------------------------------------------------------
// calcul du plus court chemin
//------------------------------------------------------------------------------

template<typename Val>
TimeTexture<Val> GraphPath<Val>::process(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr)
{


     uint ns=initmesh.vertex().size();


     std::cerr << "Shortest Path In" << std::endl;

     // ici on construit la carte de distance avec le fast marching à la cathier
//      typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
//      til::Mesh_N meshCathier;
//      til::Mesh1 mesh0;
//      til::convert(mesh0, initmesh);
//      meshCathier = addNeighborsToMesh(mesh0);
// 
//      shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(meshCathier), getFaceIndices(meshCathier));

// /*     til::Triangle_mesh_geodesic_map<til::Mesh_N::VertexCollection, CNeighborhoods, double> geomap(getVertices(meshCathier), *pneighc);*/
//      til::Graph_distance_map<til::Mesh_N::VertexCollection, CNeighborhoods, double> geomap(getVertices(meshCathier), *pneighc);
// 
//      std::cerr << "Shortest Path In" << std::endl;
// 
//      std::vector<std::size_t> startPoints;
//      std::vector<double> dist;
//      startPoints.push_back(dep);
//      dist.push_back(0.0);*/

     // la boucle ci dessous intitialise la carte de distance à 0 pour le départ
     // et l'infini pour les points hors de l'ensemble de recherche, pour avoir une distance
     // limitée à cet ensemble.
     
//      for (uint l=0; l<ns; l++)
//      {
//           if (tex[0].item(l) != value)
//           {
//                startPoints.push_back(l);
//                dist.push_back( std::numeric_limits<double>::infinity() );
//           }
//      }
// 
//      geomap.init(startPoints, dist);
//      geomap.process();
//      shared_ptr<til::sparse_vector<double> > sres = geomap.distanceMap(); 
//      std::map<uint,float> res;
//      til::sparse_vector<double>::const_iterator iRes2 = sres->begin();
//      uint i=0;
// 
//      for (; iRes2 != sres->end();++iRes2,i++)
//      {
//           res[i] = (float) *iRes2;
//      }


     // ici on utilise finalement le bon vieil algo de distanceMap d'aimsalgo

     

     TimeTexture<float> distanceMap;
     Texture<short> initMap(ns);
     for (uint l=0; l<ns; l++)
     {
          if (tex[0].item(l)==value)
               initMap.item(l)=(short) 0; // domaine de calcul de la distance
          else initMap.item(l)=(short) -1; // domaine interdit
     }
     initMap.item(dep)=(short) 100; // point de départ
     
     distanceMap[0]=meshdistance::MeshDistance( initmesh[0] , initMap , true);
     std::map<uint,float> res;
     uint i=0;
     for (uint l=0; l<ns; l++)
     {
          res[l]=(float) distanceMap[0].item(l);
     }


     std::vector<std::set<uint> > neigh=SurfaceManip::surfaceNeighbours(initmesh);
     TimeTexture<Val> result(1,ns);
     i=arr;
     uint j, inext, iprevious=i;
     double distmin;
     result[0].item(i)=(Val) value;
     distmin=10000.0;
     TimeTexture<Val> debugTex(tex);


     // ici on fait le backtracking

     do
     {
          std::set<uint> vois=neigh[i];
          std::set<uint>::iterator itVois=vois.begin();
          inext=i; // distmin=10000.0;
          for (; itVois != vois.end(); ++itVois)
          {
               j=(*itVois);
               if ((tex[0].item(j)==value) && (j != iprevious))
               {
                    if (res[j]<distmin)
                    {
                         inext=j; distmin=res[j];
                    }
               }
          }
          debugTex[0].item(inext)=value*2;
          if (inext==i)
          {
               std::cerr << "ShortestPath->GraphPath<Val>::process : problem. There is no path between start and end included in the provided set" << std::endl;
               Writer<TimeTexture<Val> > debugW("~/debugTex.tex");
               debugTex[0].item(dep)=value*3;
               debugTex[0].item(arr)=value*3;
               debugW.write(debugTex);
               exit(EXIT_FAILURE);
               
          }
          longueur+=(initmesh.vertex()[i] - initmesh.vertex()[inext]).norm();
          iprevious=i; i=inext;
          result[0].item(i)=(Val) value;
     }
     while (i!=dep);
     std::cerr << "Shortest Path Out" << std::endl;

     return(result);
     
}

//------------------------------------------------------------------------------
// renvoie seulement la longueur du plus court chemin
//------------------------------------------------------------------------------

template<typename Val>
float GraphPath<Val>::getLongueur(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr)
{
     process(tex, initmesh, value, dep, arr);
     return longueur;
}

}



#endif

/*
 *  This software and supporting documentation were developed by
 *  	CNRS, labo LSIS, UMR6168, Marseille, France
 *      olivier.coulon@univmed.fr
 */



#include <cstdlib>
#include <aims/mesh/texture.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surfaceOperation.h>
#include <cortical_surface/surfacereferential/conformalMapping.h>


namespace aims
{

//--- a simple utility that computes edges of the mesh ---

     void ConformalMapping::ComputeEdges()
     {
          uint i,j;
          std::cout << "Computing edges" << std::endl;
          for (i=0; i< _nv; i++)
          {
               std::set<uint> vois=_neigh[i];
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
          buildEdgeMap();
     }
     
     
//--- building a vertex, edge association map ------------------
//--- a bit redundant but easy ---------------------------------
     void ConformalMapping::buildEdgeMap()
     { 
          uint i, ne=_edges.size();
          std::cout<< "Building edge map" << std::endl;
          for (i=0; i<ne; i++)
          {
               std::set<std::pair<uint, uint> > a;
               _edgeMap.push_back(a);
          }
          for (i=0; i<ne; i++)
          {
               uint u=_edges[i][0];
               uint v=_edges[i][1];
               (_edgeMap[u]).insert(std::pair<uint, uint>(v, i));
               (_edgeMap[v]).insert(std::pair<uint, uint>(u, i));
          }
     }



//--- another utility that returns the index of an edge --------

     uint ConformalMapping::whichEdge(uint v1, uint v2)
     {
/*          uint i, ne=_edges.size();*/
          std::set<std::pair<uint, uint> > vois;
          std::set<std::pair<uint, uint> >::iterator voisIt;
          vois=_edgeMap[v1];
          for (voisIt=vois.begin(); voisIt!=vois.end(); ++voisIt)
          {
               if ((*voisIt).first==v2) return((*voisIt).second);
          }

/*          
          for (i=0; i<ne; i++)
               if (((v1==_edges[i][0]) || (v1==_edges[i][1])) && ((v2==_edges[i][0]) || (v2==_edges[i][1])))
                    return(i);*/
          std::cerr << "Pb with whichEdge(" << v1 << "," << v2 << ") : edge does not exist" << std::endl;
          exit(EXIT_FAILURE);
     }
     
//--- compute the center of mass to constrain the conformal mapping ---

     Point3df ConformalMapping::getCenterOfMass(std::vector<Point3df> tmpMap)
     {
          std::vector< AimsVector<uint,3> > poly=_mesh.polygon();
          std::vector< AimsVector<uint,3> >::iterator polyIt;
          Point3df v1, v2, v3;
          Point3df center, bary;
          double area=0.0;

          for (polyIt=poly.begin(); polyIt!=poly.end(); ++polyIt)
          {
               v1=tmpMap[(*polyIt)[0]];
               v2=tmpMap[(*polyIt)[1]];
               v3=tmpMap[(*polyIt)[2]];
               bary=(v1+v2+v3)/3.0;
               double aire;
               Point3df cross=vectProduct( v2-v1, v3-v1);
               aire=cross.dnorm()/2.0;
               center=center + (bary*aire);
               area += aire;
          }
          center = center / area;
          return(center);
          
          
//           float cx=cy=cz=0.0;
//           for (i=0; i<_nv; i++)
//           {
//                cx+=tmpMap[i][0];
//                cy+=tmpMap[i][1];
//                cz+=tmpMap[i][2];
//           }
//           cx/=(float)_nv; cy/=(float)_nv; cz/=(float)_nv;


     }

//------------------- Tuette mapping ---------------------------
// this is the first step to computing the conformal mapping
//--------------------------------------------------------------

     float ConformalMapping::TuetteEnergy()
     {
          std::vector<AimsVector<uint,2> >::iterator edgeIt;
          float E=0.0;
          for (edgeIt=_edges.begin(); edgeIt!=_edges.end(); ++edgeIt)
          {
               uint i,j;
               i=(*edgeIt)[0];
               j=(*edgeIt)[1];
               Point3df u=(_tuetteMap[i] - _tuetteMap[j]);
               E += u.norm2();
          }
          return(E);
     }

     //------------------------------------------------------------------

     void ConformalMapping::TuetteMap()
     {
          uint i;
          float E, E0, deltaE;
          std::vector<Point3df> tmpMap(_nv);
          
/*          AimsSurfaceTriangle tmpResult;*/
          std::cout << "Computing Tuette mapping with dt=" << _dt << std::endl;
          _mesh.updateNormals();
          std::vector<Point3df> norm=_mesh.normal();
          for (i=0; i<_nv; i++)
          {
               _tuetteMap[i]=norm[i];
               tmpMap[i]=norm[i];
          }
          deltaE=10000;
          E0=TuetteEnergy(); 
          uint iter=0;
//           uint time=0;
          
          while (deltaE >= _dE)
          {
               for (i=0; i<_nv; i++)
               {
               // question : comment calculer le piecewise laplacian en i?
                    float Dx=0.0, Dy=0.0, Dz=0.0; // les 3 composantes du laplacien
                    float ortho0, ortho1, ortho2, deriv0, deriv1, deriv2;
                    
                    // Calcul du Laplacien au point i 
                    std::set<uint> vois=_neigh[i];
                    std::set<uint>::iterator voisIt;
                    for (voisIt=vois.begin(); voisIt!=vois.end(); ++voisIt)
                    {
                         uint j=*voisIt;
                         Dx+=(_tuetteMap[i][0]-_tuetteMap[j][0]);
                         Dy+=(_tuetteMap[i][1]-_tuetteMap[j][1]);
                         Dz+=(_tuetteMap[i][2]-_tuetteMap[j][2]);
                    }
                    
               // here the normal to the target surface is the normal to the unit sphere
               // i.e. n(f(v))=f(v) with f(v) the Tuette mapping.
               // actually we normalize systematically n(f(v)) for stability.
                    float norm=sqrt((_tuetteMap[i][0]*_tuetteMap[i][0]) + (_tuetteMap[i][1]*_tuetteMap[i][1]) + (_tuetteMap[i][2]*_tuetteMap[i][2]));
                    float nx=_tuetteMap[i][0] / norm;
                    float ny=_tuetteMap[i][1] / norm;
                    float nz=_tuetteMap[i][2] / norm;
                    ortho0=(nx*Dx + ny*Dy + nz*Dz)*nx;
                    ortho1=(nx*Dx + ny*Dy + nz*Dz)*ny;
                    ortho2=(nx*Dx + ny*Dy + nz*Dz)*nz;
                    deriv0=Dx-ortho0;
                    deriv1=Dy-ortho1;
                    deriv2=Dz-ortho2;
                    tmpMap[i][0]=_tuetteMap[i][0] - _dt*deriv0;
                    tmpMap[i][1]=_tuetteMap[i][1] - _dt*deriv1;
                    tmpMap[i][2]=_tuetteMap[i][2] - _dt*deriv2;
               }
               for (i=0; i<_nv; i++)
               {
                    float norm=sqrt((tmpMap[i][0]*tmpMap[i][0]) + (tmpMap[i][1]*tmpMap[i][1]) + (tmpMap[i][2]*tmpMap[i][2]));
                    _tuetteMap[i][0]=tmpMap[i][0]/norm;
                    _tuetteMap[i][1]=tmpMap[i][1]/norm;
                    _tuetteMap[i][2]=tmpMap[i][2]/norm;
               }
               E=TuetteEnergy();
               deltaE=fabs(E-E0);
/*               std::cout << "Iteration " << iter << " : E=" << E << ", deltaE=" << deltaE << std::endl;*/
               if (((iter%100)==0) /*&& (iter<10000)*/)
               {
                    std::cout << "Iteration " << iter << " : E=" << E << ", deltaE=" << deltaE << std::endl;
//                     tmpResult[time].vertex()=_tuetteMap;
//                     tmpResult[time].polygon()=_mesh.polygon();
//                     time++;
               }
               iter++;
               E0=E;
          }
          Writer<AimsSurfaceTriangle> meshW("/home/olivier/tuette.mesh");
          AimsSurfaceTriangle tmpT;
          tmpT.vertex()=_tuetteMap;
          tmpT.polygon()=_mesh.polygon();
          tmpT.updateNormals();
          meshW.write(tmpT);
//           Writer<AimsSurfaceTriangle> meshW("evolutionConformal.mesh");
//           meshW.write(tmpResult);

     }

//------------------------- Conformal Mapping ----------------------
// the actual conformal mapping with zero gravity center constraint
// initialized with the Tuette Mapping
//------------------------------------------------------------------

     void ConformalMapping::ComputeHarmonicCoefficients()
     {
          
          std::vector< AimsVector<uint,3> > poly=_mesh.polygon();
          std::vector<Point3df> vert=_mesh.vertex();
          std::vector< AimsVector<uint,3> >::iterator polyIt;
          uint ne=_edges.size();
          uint i;
          uint v1,v2,v3;
          float a1, a2, a3;
          float x1, y1, z1;
          float x2, y2, z2;
          float x3, y3, z3;
          std::cout << "Computing harmonic coeeficients" << std::endl;

          for (i=0; i<ne; i++)
               _harmonicCoef.push_back(0);


          for (polyIt=poly.begin(); polyIt!=poly.end(); ++polyIt)
          {
               v1=(*polyIt)[0]; v2=(*polyIt)[1]; v3=(*polyIt)[2];
               x1=vert[v1][0]; y1=vert[v1][1]; z1=vert[v1][2];
               x2=vert[v2][0]; y2=vert[v2][1]; z2=vert[v2][2];
               x3=vert[v3][0]; y3=vert[v3][1]; z3=vert[v1][2];
               
               Point3df cross;
               
               cross=vectProduct( vert[v2]-vert[v1], vert[v3]-vert[v1]);
               a1=0.5*((x2-x1)*(x3-x1) + (y2-y1)*(y3-y1) + (z2-z1)*(z3-z1))/cross.norm();
               cross=vectProduct( vert[v1]-vert[v2], vert[v3]-vert[v2]);
               a2=0.5*((x1-x2)*(x3-x2) + (y1-y2)*(y3-y2) + (z1-z2)*(z3-z2))/cross.norm();
               cross=vectProduct( vert[v1]-vert[v3], vert[v2]-vert[v3]);
               a3=0.5*((x1-x3)*(x2-x3) + (y1-y3)*(y2-y3) + (z1-z3)*(z2-z3))/cross.norm();
               
               i=whichEdge(v2,v3); _harmonicCoef[i]=_harmonicCoef[i]+a1;
               i=whichEdge(v1,v3); _harmonicCoef[i]=_harmonicCoef[i]+a2;
               i=whichEdge(v1,v2); _harmonicCoef[i]=_harmonicCoef[i]+a3;
          }
     }
     
     //------------------------------------------------------------------

     float ConformalMapping::HarmonicEnergy()
     {
          uint i, ne=_edges.size();
          float E=0.0;
          for (i=0; i<ne; i++)
          {
               uint u,v;
               u=(_edges[i])[0];
               v=(_edges[i])[1];
               Point3df p=(_conformalMap[u] - _conformalMap[v]);
               E += (_harmonicCoef[i]*p.norm2());
          }
          return(E);
     }
     
     //------------------------------------------------------------------

     void ConformalMapping::ConformalMap(std::string adr_tuette="nofile")
     {
     //Here we assume that Tuette mapping has been computed or is given as a file
     // name (optional argument)
          uint i;
          float E, E0, deltaE;
          std::vector<Point3df> tmpMap(_nv);
          float cx, cy, cz;
          float x, y, z;
          AimsSurfaceTriangle tmpResult;
          

          if (adr_tuette=="nofile") TuetteMap();
          else
          {
               Reader < AimsSurfaceTriangle > r(adr_tuette);
               AimsSurfaceTriangle meshT;
               r.read( meshT );
               uint nvT=meshT.vertex().size();
               if (nvT==_nv)
               {
                    _tuetteMap=std::vector<Point3df>(_nv);
                    std::vector<Point3df> vert=meshT.vertex();
                    for (i=0; i<_nv; i++)
                    {
                         _tuetteMap[i]=vert[i];
                    }
               }
               else
               {
                    std::cerr << "Pb : Tuette file has wrong number of vertex" << std::endl;
               }
          }
          // TEMP !!!!!!!!!! 
          _dt=0.05;
          //je sais, c'est mal...
          
          std::cout << "Computing conformal mapping with _dt=" << _dt << std::endl;

          for (i=0; i<_nv; i++)
          {
               _conformalMap[i]=_tuetteMap[i];
               tmpMap[i]=_tuetteMap[i];
          }
          deltaE=10000;
          E0=HarmonicEnergy(); 
          uint iter=0;
          uint time=0;
          std::cout << "\t\t\tStarting with E0=" << E0 << std::endl;

          while (deltaE >= _dE)
          {
                for (i=0; i<_nv; i++)
               {
                    float Dx=0.0, Dy=0.0, Dz=0.0; // les 3 composantes du laplacien
                    float ortho0, ortho1, ortho2, deriv0, deriv1, deriv2;
                    
                    // Calcul du Laplacien au point i 
                    std::set<uint> vois=_neigh[i];
                    std::set<uint>::iterator voisIt;
                    for (voisIt=vois.begin(); voisIt!=vois.end(); ++voisIt)
                    {
                         uint j=*voisIt;
                         uint edge=whichEdge(i,j);
                         float h=_harmonicCoef[edge];
                         Dx+=h*(_conformalMap[i][0]-_conformalMap[j][0]);
                         Dy+=h*(_conformalMap[i][1]-_conformalMap[j][1]);
                         Dz+=h*(_conformalMap[i][2]-_conformalMap[j][2]);
                    }
                    
               // here the normal to the target surface is the normal to the unit sphere
               // i.e. n(f(v))=f(v) with f(v) the conrmal mapping.
               // actually we normalize systematically n(f(v)) for stability.
                    float norm=sqrt((_conformalMap[i][0]*_conformalMap[i][0]) + (_conformalMap[i][1]*_conformalMap[i][1]) + (_conformalMap[i][2]*_conformalMap[i][2]));
                    float nx=_conformalMap[i][0] / norm;
                    float ny=_conformalMap[i][1] / norm;
                    float nz=_conformalMap[i][2] / norm;
                    ortho0=(nx*Dx + ny*Dy + nz*Dz)*nx;
                    ortho1=(nx*Dx + ny*Dy + nz*Dz)*ny;
                    ortho2=(nx*Dx + ny*Dy + nz*Dz)*nz;
                    deriv0=Dx-ortho0;
                    deriv1=Dy-ortho1;
                    deriv2=Dz-ortho2;
                    tmpMap[i][0]=_conformalMap[i][0] - _dt*deriv0;
                    tmpMap[i][1]=_conformalMap[i][1] - _dt*deriv1;
                    tmpMap[i][2]=_conformalMap[i][2] - _dt*deriv2;
               }
               // zero mass-center condition  - VERIFIER CETTE BOUCLE
//                cx=cy=cz=0.0;
//                for (i=0; i<_nv; i++)
//                {
//                     cx+=tmpMap[i][0];
//                     cy+=tmpMap[i][1];
//                     cz+=tmpMap[i][2];
//                }
//                cx/=(float)_nv; cy/=(float)_nv; cz/=(float)_nv;

               Point3df center=getCenterOfMass(tmpMap);
               cx=center[0]; cy=center[1]; cz=center[2];

               for (i=0; i<_nv; i++)
               {
                    x=(tmpMap[i][0]-cx); y=(tmpMap[i][1]-cy); z=(tmpMap[i][2]-cz);
                    float norm=sqrt(x*x + y*y + z*z);
                    _conformalMap[i][0]=x/norm;
                    _conformalMap[i][1]=y/norm;
                    _conformalMap[i][2]=z/norm;
               }
               E=HarmonicEnergy();
               deltaE=fabs(E-E0);
//               std::cout << "Iteration " << iter << " : E=" << E << ", deltaE=" << deltaE << std::endl;
               if (((iter%10)==0) /*&& (iter<10000)*/)
               {
                    std::cout << "Iteration " << iter << " : E=" << E << ", deltaE=" << deltaE << std::endl;
                    std::cout << "\tcx=" << cx << ", cy=" << cy << ", cz=" << cz << std::endl;
                    tmpResult[time].vertex()=_conformalMap;
                    tmpResult[time].polygon()=_mesh.polygon();
                    time++;
               }
               iter++;
               E0=E;
          }
          tmpResult[time].vertex()=_conformalMap;
          tmpResult[time].polygon()=_mesh.polygon();
          Writer<AimsSurfaceTriangle> meshW("evolutionConformal.mesh");
          meshW.write(tmpResult);

     }

     //------------------------------------------------------------------

     void ConformalMapping::ComputeConformalCoordinates()
     {
          Point3df v;
          float alpha, sigma;
          float x, y, z;
          
          _conformalLon=std::vector<float>(_nv);
          _conformalLat=std::vector<float>(_nv);
     
          for (uint i=0; i<_nv; i++)
          {
               v=_conformalMap[i];
               x=v[0]; y=v[1]; z=v[2];
               sigma=acosf(z);
               _conformalLat[i]=sigma*180.0/M_PI;
               if ((z<1.0) && (z>-1.0))
               {
                    alpha=acosf(x/sinf(sigma));
                    if (y>=0) _conformalLon[i]=alpha*180.0/M_PI;
                    else _conformalLon[i]=360.0-(alpha*180.0/M_PI);
               }
               else
                _conformalLat[i]=0;
          }
     }

//-------------------------------------------------------------------
// function to access and compute the mapping meshes and coordinates
//-------------------------------------------------------------------

     AimsSurfaceTriangle ConformalMapping::GetConformalMapping(std::string adr_tuette="nofile")
     {
          AimsSurfaceTriangle result;
          
          ConformalMap(adr_tuette);
          result.vertex()=_conformalMap;
          result.polygon()=_mesh.polygon();
          ComputeConformalCoordinates();
          return(result);
     }

     AimsSurfaceTriangle ConformalMapping::GetTuetteMapping()
     {
          AimsSurfaceTriangle result;
          TuetteMap();
          result.vertex()=_tuetteMap;
          result.polygon()=_mesh.polygon();
          return(result);
     }
     
     void ConformalMapping::WriteConformalCoordinates(std::string adr_coord)
     {
          TimeTexture<float> lat(1, _nv), lon(1, _nv);
          uint i;
          for (i=0; i<_nv; i++)
          {
               lat[0].item(i)=_conformalLat[i];
               lon[0].item(i)=_conformalLon[i];
          } 
          std::string strlon=adr_coord+"_lon.tex";
          Writer<TimeTexture<float> > lonW(strlon);
          lonW.write(lon);
          std::string strlat=adr_coord+"_lat.tex";
          Writer<TimeTexture<float> > latW(strlat);
          latW.write(lat);
     }
}

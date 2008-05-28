#ifndef AIMS_CONNECT_PATH_H
#define AIMS_CONNECT_PATH_H


#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <cortical_surface/mesh/pointDistance.h>

namespace aims
{

// When a line on a mesh is disconnected in two connected components, this class
// connect it in a single component with value val 

template<typename T> class ConnectMeshPath
{


public:
     
     ConnectMeshPath(AimsSurfaceTriangle mesh, TimeTexture<T> tex, T label1, T label2) : 
     _mesh(mesh), _tex(tex), _l1(label1), _l2(label2) {}
     Texture<T> run(T val);
     
private:
     AimsSurfaceTriangle _mesh;
     TimeTexture<T> _tex;
     T _l1;
     T _l2;
};


template<typename T> Texture<T> ConnectMeshPath<T>::run(T val)
{
     std::vector<std::set<uint> >  neigh = SurfaceManip::surfaceNeighbours(_mesh);
     std::set<uint>::iterator neighIt;
     
     std::set<uint> set1, set2;
     std::set<uint> end1, end2;
     std::set<uint>::iterator setIt, setIt2;
     uint ext1, ext2;
     
     MeshPointDistance distance(_mesh);
     
     int ns=_tex[0].nItem(), i;
     
     Texture<T> texOut(ns);
     std::vector<Point3df> vert=_mesh.vertex();
     Point3df pt1, pt2;

     std::cout << "Connecting Paths" << std::endl;
     std::cout << "\tBuilding the two sets" << std::endl;
     // building the two sets to connect
     for (i=0; i<ns; i++)
     {
          if (fabs((double)(_tex[0].item(i) - _l1)) < 0.001)
               set1.insert(i);
          else if (fabs((double)(_tex[0].item(i) - _l2)) < 0.001)
               set2.insert(i);
     }
     // finding extremities of both sets
     std::cout << "\tFinding their extremities" << std::endl;

     for (setIt=set1.begin(); setIt!=set1.end(); ++setIt)
     {
          std::set<uint> voisin=neigh[*setIt];
          std::set<uint>::iterator vIt;
          int count=0; 
          for (vIt=voisin.begin(); vIt!=voisin.end(); ++vIt)
               if (set1.find(*vIt) != set1.end())
                    count++;
          if (count < 2)
               end1.insert(*setIt); 
     }
     for (setIt=set2.begin(); setIt!=set2.end(); ++setIt)
     {
          std::set<uint> voisin=neigh[*setIt];
          std::set<uint>::iterator vIt;
          int count=0; 
          for (vIt=voisin.begin(); vIt!=voisin.end(); ++vIt)
               if (set2.find(*vIt) != set2.end())
                    count++;
          if (count < 2)
               end2.insert(*setIt); 
     }
     
     // the two extremities to join are the closest to eachother
     // at the moment I use the 3D distance 
     
     float dist, distMin=10000.0;
     for (setIt=end1.begin(); setIt!=end1.end(); ++setIt)
          for (setIt2=end2.begin(); setIt2!=end2.end(); ++setIt2)
          {
               pt1=vert[*setIt]; pt2=vert[*setIt2];
               dist=sqrt( (pt1[0]-pt2[0])*(pt1[0]-pt2[0])
                        + (pt1[1]-pt2[1])*(pt1[1]-pt2[1])
                        + (pt1[2]-pt2[2])*(pt1[2]-pt2[2]) );
               if (dist<distMin)
               {
                    distMin=dist;
                    ext1=*setIt;
                    ext2=*setIt2;
               }
          }
     
     // now let's join ext1 and ext2; 
     // this time we have no other choice than using geodesic distance
     std::cout << "\tJoining them" << std::endl;

     for (i=0; i<ns; i++)
     {
          if ( (set1.find(i) != set1.end()) || (set2.find(i) != set2.end()) )
               texOut.item(i)=val;
          else 
               texOut.item(i)=0;
     }

     i=ext1;
     uint j;
     while (i != ext2)
     {
          std::set<uint> voisin=neigh[i];
          std::set<uint>::iterator vIt; 
          distMin=10000.0;
          for (vIt=voisin.begin(); vIt!=voisin.end(); ++vIt)
          {
               dist=distance.compute(ext2, *vIt); // APPELER ICI LA DISTANCE GEODESIQUE ENTRE *vIt et ext2
               if (dist<distMin)
               {
                    distMin=dist;
                    j=*vIt;
               }
          }
          texOut.item(j)=val;
          i=j;
     }
     std::cout << "\tReturning texOut" << std::endl;
     return(texOut);

}

}  //fin du namespace

#endif

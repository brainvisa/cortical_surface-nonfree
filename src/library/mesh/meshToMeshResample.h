#ifndef AIMS_MESH2MESH_RESAMPLE_H
#define MESH2MESH_RESAMPLE_H


#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <cortical_surface/mesh/pointDistance.h>

namespace aims
{



// When a two meshes are parameterized with surface coordinates, it is handy
// to be able to resample a texture from one to another, or to resample one mesh
// such that it has a node-to-node correspondance with the other one.
// This is the point of this class.
// Constraint : the _sx and _tx textures are the periodic ones (e.g. longitude
// for cortex, or x for sulci)

template<typename T> class Mesh2mesh
{
public:
     
     Mesh2mesh(AimsSurfaceTriangle source, AimsSurfaceTriangle target, TimeTexture<float> sx, TimeTexture<float> sy, TimeTexture<float> tx, TimeTexture<float> ty, float period) :
     _source(source), _target(target), _sx(sx), _sy(sy), _tx(tx), _ty(ty), _period(period) {checkIntegrity(); buildAlternateRep();}
     TimeTexture<T> sendTextureToTarget(TimeTexture<T> text);

     int isInTriangle(float x, float y, uint i);  // une coordonnée est-elle dans le triangle
                                                  // (de la source)

     
     
private:
     AimsSurfaceTriangle _source;
     AimsSurfaceTriangle _target;
     TimeTexture<float> _sx;
     TimeTexture<float> _sy;
     TimeTexture<float> _tx;
     TimeTexture<float> _ty;
     float _period;
     uint _ns;
     uint _nt;
     
     void checkIntegrity();
     void buildAlternateRep();
};


template<typename T> void Mesh2mesh<T>::checkIntegrity()
{
     int nsx=_sx[0].nItem(), nsy=_sy[0].nItem(), ntx=_tx[0].nItem(), nty=_ty[0].nItem();
     int ns=_source.vertex().size();
     int nt=_target.vertex().size();

     if ((nsx != nsy) || (nsx != ns) || (nsy !=ns))
     {
          std::cerr << "Mesh2mesh : number of nodes/items faulty for source" << std::endl;
          exit(EXIT_FAILURE);
     }
     else if ((ntx != nty) || (ntx != nt) || (nty !=nt))
     {
          std::cerr << "Mesh2mesh : number of nodes/items faulty for target" << std::endl;
          exit(EXIT_FAILURE);
     }
     else
     {
          _ns=ns; _nt=nt;
     }
}

template<typename T> void Mesh2mesh<T>::buildAlternateRep() // A FAIRE
{

}



template<typename T> int Mesh2mesh<T>::isInTriangle(float px, float py, uint i)
{
     AimsVector<uint, 3> triangle=_source.polygon()[i];
     uint p1=triangle[0];
     uint p2=triangle[1];
     uint p3=triangle[2];
     float p1x, p1y, p2x, p2y, p3x, p3y;
     double vp3, vp2, vp1;
     double eps=0.01;

     p1x=_sx[0].item(p1); p1y=_sy[0].item(p1);
     p2x=_sx[0].item(p2); p2y=_sy[0].item(p2);
     p3x=_sx[0].item(p3); p3y=_sy[0].item(p3);


     // The two following have to be positive

     vp1=((p3x-p1x)*(py-p1y) - (p3y-p1y)*(px-p1x))
        *((px-p1x)*(p2y-p1y) - (py-p1y)*(p2x-p1x));

     vp2=((p3x-p2x)*(py-p2y) - (p3y-p2y)*(px-p2x))
        *((px-p2x)*(p1y-p2y) - (py-p2y)*(p1x-p2x));

     // this one sort ambiguous cases (colinearities)

     vp3=((p1x-p3x)*(py-p3y) - (p1y-p3y)*(px-p3x))
        *((px-p3x)*(p2y-p3y) - (py-p3y)*(p2x-p3x));

     if ((vp1>eps) && (vp2>eps))
     {
          return(1);
     }
     else if (fabs(vp1) < eps)
     {
          if (vp2>eps)
               return(1);
          else if (fabs(vp2) < eps)
          {
               if (vp3>eps)
                    return(1);
               else return(0);
          }
          else return(0);
     }
     else if (fabs(vp2) < eps)
     {
          if (vp1 > eps)
               return(1);
          else return(0);
     }
     else return(0);
}

template<typename T> TimeTexture<T> Mesh2mesh<T>::sendTextureToTarget(TimeTexture<T> text)
{
     std::vector< AimsVector<uint,3> > poly=_source.polygon();
     uint ntr=poly.size();
     TimeTexture<T> result(1,_nt);
     uint i, j, polsrc;
     float x,y;
     int flag=0;
     float x1, x2, x3, y1, y2, y3;
     
     if (text[0].nItem() != _ns)
     {
          std::cerr << "Mesh2Mesh::sendtextureTotarget : texture does not have the same nb ofitems than the source" << std::endl;
          exit(EXIT_FAILURE);
     }

     for (i=0; i<_nt; i++)
     {
          x=_tx[0].item(i);
          y=_ty[0].item(i);
          flag=0;
          for (j=0; (j<ntr) && (flag==0); j++)
          {
               if (isInTriangle(x, y, j)==1)
               {
                    polsrc=j;
                    flag=1;
               }
          }
          if (flag==0)
          {
               result[0].item(i)=(T) -10;
          }
          else
          {
               double t1, t2, t3;
               double sum;
               double l1, l2, l3;
               x1=_sx[0].item(poly[j][0]); y1=_sy[0].item(poly[j][0]);
               x2=_sx[0].item(poly[j][1]); y2=_sy[0].item(poly[j][1]);
               x3=_sx[0].item(poly[j][2]); y3=_sy[0].item(poly[j][2]);

               l1=sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
               l2=sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1));
               l3=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

               t1=fabs((x2-x3)*(y-y3) - (y2-y3)*(x-x3))/l1;
               t2=fabs((x3-x1)*(y-y1) - (y3-y1)*(x-x1))/l2;
               t3=fabs((x1-x2)*(y-y2) - (y1-y2)*(x-x2))/l3;
               sum=t1+t2+t3; t1/=sum; t2/=sum; t3/=sum;

               result[0].item(i)=(T) t1*text[0].item(poly[j][0]) + t2*text[0].item(poly[j][1]) + t3*text[0].item(poly[j][2]);
          }
     }
     return(result);

     
}

}  //fin du namespace

#endif

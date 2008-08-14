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

class Mesh2mesh
{
public:
     
     Mesh2mesh(AimsSurfaceTriangle source, AimsSurfaceTriangle target, TimeTexture<float> sx, TimeTexture<float> sy, TimeTexture<float> tx, TimeTexture<float> ty, float period) :
     _source(source), _target(target), _sx(sx), _sy(sy), _tx(tx), _ty(ty), _period(period) {checkIntegrity(); buildTriangles();}
     TimeTexture<float> sendTextureToTarget(TimeTexture<float> text);
     int isInTriangle(float x, float y, uint i); 
     uint nearestNeighbour(float x, float y);
     AimsSurfaceTriangle remeshSourceToTarget();

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
     void buildTriangles();
     float _min(float a, float b, float c);
     float _max(float a, float b, float c);

     std::vector< AimsVector<std::pair<float, float>, 3> > _triangles;
};


//---------------------------------------------------------------------
//
//                   methods
//
//----------------------------------------------------------------------


//----------------- min and max of 3 floating point -----------------

float Mesh2mesh::_min(float a, float b, float c)
{
     if (a<=b)
          if (a<=c) return a;
          else return c;
     else 
          if (b<=c) return b;
          else return c;
}

float Mesh2mesh::_max(float a, float b, float c)
{
     if (a>=b)
          if (a>=c) return a;
          else return c;
     else 
          if (b>=c) return b;
          else return c;
}


//------------- checking dimensions of all meshes and textures ---------------

void Mesh2mesh::checkIntegrity()
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

//-------- building alternate representation of the mesh (for speed) ---------

void Mesh2mesh::buildAlternateRep() // A FAIRE
{

}


// Dealing with periodicity : all source triangles coordinates are recomputed properly

void Mesh2mesh::buildTriangles()
{
     AimsVector<uint, 3> tri;
     std::vector<AimsVector<uint, 3> > poly=_source.polygon();
     uint p1, p2, p3;
     float x1, x2, x3, y1, y2, y3;
     uint i;
     float eps=0.001;
     uint np=poly.size();
          
     for (i=0; i<np; i++)
     {
          tri=poly[i];
          p1=tri[0]; p2=tri[1]; p3=tri[2];
          x1=_sx[0].item(p1); y1=_sy[0].item(p1);
          x2=_sx[0].item(p2); y2=_sy[0].item(p2);
          x3=_sx[0].item(p3); y3=_sy[0].item(p3);
          // the following assumes that all points on the origin meridian 
          // have a value equal to 0(==px), and that no triangle 'cross'
          // this meridian (one point on each side)

          if (fabs(x1)<eps)
          {
               if ((x2 > _period - x2) || (x3 > _period - x3))
               {
                    x1=_period;
                    if (fabs(x2)<eps) x2=_period;
                    if (fabs(x3)<eps) x3=_period;
               }
          }
          else if (fabs(x2)<eps)
          {
               if ((x1>_period-x1) || (x3>_period-x3))
               {
                    x2=_period;
                    if (fabs(x3)<eps) x3=_period;
               }
          }
          else if (fabs(x3)<eps)
          {
               if ((x1>_period-x1) || (x2>_period-x2))
                    x3=_period;
          }

          AimsVector<std::pair<float, float>, 3> newTri;
          newTri[0]=std::pair<float, float>(x1, y1);
          newTri[1]=std::pair<float, float>(x2, y2);
          newTri[2]=std::pair<float, float>(x3, y3);
          _triangles.push_back(newTri);
     }

}


//------------------- is a point in a triangle ? ------------------------

int Mesh2mesh::isInTriangle(float px, float py, uint i)
{
     double p1x, p1y, p2x, p2y, p3x, p3y;
     double vp3, vp2, vp1;
     AimsVector<std::pair<float, float>, 3> newTri=_triangles[i];
     
     p1x=newTri[0].first; p1y=newTri[0].second;
     p2x=newTri[1].first; p2y=newTri[1].second;
     p3x=newTri[2].first; p3y=newTri[2].second;

     // The three following have to be positive

     vp1=((p3x-p1x)*(py-p1y) - (p3y-p1y)*(px-p1x))
        *((px-p1x)*(p2y-p1y) - (py-p1y)*(p2x-p1x));

     vp2=((p3x-p2x)*(py-p2y) - (p3y-p2y)*(px-p2x))
        *((px-p2x)*(p1y-p2y) - (py-p2y)*(p1x-p2x));

     vp3=((p1x-p3x)*(py-p3y) - (p1y-p3y)*(px-p3x))
        *((px-p3x)*(p2y-p3y) - (py-p3y)*(p2x-p3x));

     if ((vp1>=0) && (vp2>=0) && (vp3>=0))
     {
          return(1);
     }
     else return(0);
}


//------ nearestNeighbour of a coordinate pair on the source mesh -------------------------

uint Mesh2mesh::nearestNeighbour(float x, float y)
{
     uint i, nn=0;
     float dist, distMin=10000.0;
     
     for (i=0; i<_ns; i++)
     {
          dist=sqrt((x-_sx[0].item(i))*(x-_sx[0].item(i)) + (y-_sy[0].item(i))*(y-_sy[0].item(i)));
          if (dist < distMin)
          {
               nn=i;
               distMin=dist;
          }
     }
     return(nn);
}

//------- reinterpolation of a texture from the source to the target ----------------------

TimeTexture<float> Mesh2mesh::sendTextureToTarget(TimeTexture<float> text)
{
     std::vector< AimsVector<uint,3> > poly=_source.polygon();
     uint ntr_s=poly.size();
     TimeTexture<float> result(1,_nt);
     uint i, j, polsrc;
     float x,y;
     int flag=0, count=0, countZe=0;
     float x1, x2, x3, y1, y2, y3;
     
     if (text[0].nItem() != _ns)
     {
          std::cerr << "Mesh2Mesh::sendTextureTotarget : texture does not have the same nb of items than the source" << std::endl;
          exit(EXIT_FAILURE);
     }

     for (i=0; i<_nt; i++)
     {
          x=_tx[0].item(i);
          y=_ty[0].item(i);
          flag=0;
          count=0;
          for (j=0; (j<ntr_s) && (flag==0); j++)
          {
               if (isInTriangle(x, y, j)==1)
               {
                    polsrc=j;
                    flag=1;
                    count++;
               }
          }
          if (flag==0)
          {
               // temp : ce cas est particulier, parfois aucun triangle n'est trouvé
               // au voisinage des contraintes.
               // Pour l'instant on résoud en prenant le plus proche voisin, à plus
               // long terme il faut s'y prendre autrement.
               result[0].item(i)=text[0].item(nearestNeighbour(x,y));
/*               std::cout << "+" << std::endl;*/
          }
          else
          {
               double t1, t2, t3;
               float z1, z2, z3;
               AimsVector<std::pair<float, float>, 3> newTri=_triangles[polsrc];
               x1=newTri[0].first; y1=newTri[0].second;
               x2=newTri[1].first; y2=newTri[1].second;
               x3=newTri[2].first; y3=newTri[2].second;

               //-------------------------------------------------------
               // La formule ci dessous exprime l'équation du plan qui passe par les
               // trois points du triangle sous la forme (x,y,z=text(x,y) ) puis calcule
               // le z de notre point sachant qu'il appartient au plan.
               // C'est donc une interpolation linéaire.
               
               z1=text[0].item(poly[polsrc][0]);
               z2=text[0].item(poly[polsrc][1]);
               z3=text[0].item(poly[polsrc][2]);
               

               t1=(y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
               t2=(x3-x1)*(z2-z1) - (x2-x1)*(z3-z1);
               t3=(x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
               if (t3 > 0.0001)
               {
                    result[0].item(i)=((x1*t1 + y1*t2 + z1*t3 - x*t1 - y*t2)/t3);
               }
               else
               {
                    result[0].item(i)=text[0].item(nearestNeighbour(x,y));
//                // DEBUG
// 
//                     std::cout << " i="<< i << " avec x=" << x << " et y=" << y <<std::endl;
//                     std::cout << "\tx1=" << x1<< ", y1=" << y1 << ", z1=" << z1 << std::endl;
//                     std::cout << "\tx2=" << x2<< ", y2=" << y2 << ", z2=" << z2 << std::endl;
//                     std::cout << "\tx3=" << x3<< ", y3=" << y3 << ", z3=" << z3 << std::endl;
//                     countZe++;
//                // END DEBUG
               }

          }
     }
//      std::cout << "Number of Zes : " << countZe << std::endl;
//      std::cout << "Source items : " << _ns << std::endl;
//      std::cout << "Target items : " << _nt << std::endl;
     return(result);
}

//------- The source mesh is reinterpolated to the target such -------------------------------
//------- that they have a node to nodecorrespondance ----------------------------------------

AimsSurfaceTriangle Mesh2mesh::remeshSourceToTarget()
{
     TimeTexture<float> x_coord(1, _ns), y_coord(1, _ns), z_coord(1, _ns);
     TimeTexture<float> x_rs, y_rs, z_rs;
     uint i;
     std::vector<Point3df> vert_s=_source.vertex();
     std::vector<Point3df> vert_t;
     AimsSurfaceTriangle newSource;
     std::vector< AimsVector<uint,3> > poly=_target.polygon();

     for (i=0; i<_ns; i++)
     {
          x_coord[0].item(i)=vert_s[i][0];
          y_coord[0].item(i)=vert_s[i][1];
          z_coord[0].item(i)=vert_s[i][2];
     }

     x_rs=sendTextureToTarget(x_coord);
     y_rs=sendTextureToTarget(y_coord);
     z_rs=sendTextureToTarget(z_coord);

     for (i=0; i<_nt; i++)
     {
          Point3df point(x_rs[0].item(i), y_rs[0].item(i), z_rs[0].item(i));
          vert_t.push_back(point);
     }

     newSource.vertex()=vert_t;
     newSource.polygon()=poly;
     newSource.updateNormals();

     return(newSource);
}


}  //fin du namespace

#endif

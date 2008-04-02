#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <string.h>

using namespace aims;
using namespace carto;
using namespace std;




int main(int argc, const char **argv) //int argc, const char **argv)
{
     //DECLARATIONS
     std::string adressTexIn="./";
     std::string adressMesh="./";
     std::string adressTexOut="./";
     std::string adressAtlas="./";
     std::string adressIx="./";
     std::string adressIy="./";
     std::string adressAx="./";
     std::string adressAy="./";
     float px=0, py=0;

     AimsApplication     app( argc, argv, "Create an isoline mesh (tube) for a textured mesh");
     try
     {
      app.addOption( adressMesh, "-i", "input mesh");
      app.alias( "--inputMesh", "-i" );
      app.addOption( adressTexIn, "-t", "input texture (float)");
      app.alias( "--inputTex", "-t" );
      app.addOption( adressTexOut, "-o", "output texture");
      app.alias( "--outputTex", "-o" );
      app.addOption( adressAtlas, "-a", "atlas mesh");
      app.alias( "--atlasMesh", "-a" );
      app.addOption( adressIx, "-ix", "mesh x-coordinate texture");
      app.addOption( adressIy, "-iy", "mesh y-coordinate texture");
      app.addOption( adressAx, "-ax", "atlas x-coordinate texture");
      app.addOption( adressAy, "-ay", "atlas y-coordinate texture");
      app.addOption( px, "-px", "x-coord period (none=0)", 0);
      app.addOption( py, "-py", "y-coord period (none=0)", 0);
  
      app.initialize();
     }
     catch( user_interruption &)
     {
       return EXIT_FAILURE;
     }
     catch( ... )
     {
       throw;
     }

     std::cout << "Reading all mesh and textures" << endl;
     AimsSurfaceTriangle mesh, atlas;
     TimeTexture<float> texIn, texIx, texIy, texAx, texAy;
     Reader < AimsSurfaceTriangle > rm(adressMesh);
     rm.read( mesh );
     Reader < AimsSurfaceTriangle > ra(adressAtlas);
     ra.read( atlas );
     Reader < TimeTexture<float> > rt(adressTexIn);
     rt.read( texIn );
     Reader < TimeTexture<float> > rix(adressIx);
     rix.read( texIx );
     Reader < TimeTexture<float> > riy(adressIy);
     riy.read( texIy );
     Reader < TimeTexture<float> > rax(adressAx);
     rax.read( texAx );
     Reader < TimeTexture<float> > ray(adressAy);
     ray.read( texAy );

     uint na=atlas.vertex().size(), ns=mesh.vertex().size();
     TimeTexture<float> texOut(1, na), texDebug(1,na);
     float x, y, x1, y1, x2, y2, x3, y3;
     uint index1, index2, index3;
     float precisionX=10.0, precisionY=10.0;

     for (uint i=0; i<na; i++)
          texDebug[0].item(i)=0.0;

     // DABORD STRUCTURE POUR REPRESENTER LE MAILLAGE EN FONCTION DES COORDONNEES
     // CECI POUR NE PAS PARCOURIR TOUTE LA LISTE DES POINTS DU MAILLAGE A CHAQUE
     // POINT DE L'ATLAS

     std::vector< AimsVector<uint,3> >& poly=mesh.polygon();
     uint nPoly=poly.size();
     std::vector<std::set<uint> > polyVert(ns), polyVtmp(ns);
     
     std::cout << "Building alternate representation of input mesh" << endl;
     std::map<float, std::vector<std::pair<float, uint> > > mesh2;
     for (uint i=0; i<ns; i++)
     {
          x=texIx[0].item(i);
          y=texIy[0].item(i);
          mesh2[x].push_back(std::pair<float,uint>(y,i));
     }
     // alternate representation of polygons
     std::map<float, std::vector<std::pair<float, AimsVector<uint,3> > > > poly2;
     AimsVector<uint, 3> tri;
     for (uint i=0; i<nPoly; i++)
     {
          tri=poly[i];
          int j=tri[0];
          x=texIx[0].item(j);
          y=texIy[0].item(j);
          poly2[x].push_back(std::pair<float, AimsVector<uint,3> >(y, tri));
     }


    // PARCOURS DES NOEUDS DE L'ATLAS ET POUR CHACUN : INTERPOLATION

    std::cout << "Interpolating texture onto atlas" << endl;
    
    for (uint i=0; i<na; i++)
    {
          x=texAx[0].item(i);
          y=texAy[0].item(i);
          
          double dist, distMin=10000.0;
          int flag=0;
//           std::map<float, std::vector<std::pair<float, uint> > >::iterator meshIt=mesh2.begin();
//           for ( ; (meshIt!=mesh2.end()) && (flag==0) ; ++meshIt)
          std::map<float, std::vector<std::pair<float, AimsVector<uint,3> > > >::iterator polyIt=poly2.begin();
          for (; (polyIt!=poly2.end()) && (flag==0) ; ++polyIt)
          {
               x2=(*polyIt).first;
               if (x2 > (x+precisionX))
                    flag=1;
               else if (x2 >= (x-precisionX))
               {
                    std::vector<std::pair<float, AimsVector<uint,3> > >::iterator yIt=((*polyIt).second).begin();
                    for ( ; (yIt!=((*polyIt).second).end()) && (flag==0); ++yIt)
                    {
                         y2=(*yIt).first;
                         //std::cout << "\t\tDEBUG : in interval, y2=" << y2 << endl;
                         if ((y2>=(y-precisionY)) && (y2<=(y+precisionY)))
                         {
                              AimsVector<uint, 3> triangle=(*yIt).second;
                              index1=triangle[0];
                              index2=triangle[1];
                              index3=triangle[2];
                              x1=texIx[0].item(index1);y1=texIy[0].item(index1);
                              x2=texIx[0].item(index2);y2=texIy[0].item(index2);
                              x3=texIx[0].item(index3);y3=texIy[0].item(index3);
                              // gestion de la periodicite
                              if (x <= fabs(px-x))
                              {
                                   if (x1 > fabs(px-x1)) x1=x1-px;
                                   if (x2 > fabs(px-x2)) x2=x2-px;
                                   if (x3 > fabs(px-x2)) x3=x3-px;
                              }
                              else
                              {
                                   if (x1 < fabs(px-x1)) x1=px+x1;
                                   if (x2 < fabs(px-x2)) x2=px+x2;
                                   if (x3 < fabs(px-x2)) x3=px+x3;
                              }
                              if (y <= fabs(py-y))
                              {
                                   if (y1 > fabs(py-y1)) y1=y1-py;
                                   if (y2 > fabs(py-y2)) y2=y2-py;
                                   if (y3 > fabs(py-y2)) y3=y3-py;
                              }
                              else
                              {
                                   if (y1 < fabs(py-y1)) y1=py+y1;
                                   if (y2 < fabs(py-y2)) y2=py+y2;
                                   if (y3 < fabs(py-y2)) y3=py+y3;
                              }
                              float v1=((x2-x1)*(y-y1) - (y2-y1)*(x-x1))*((x-x1)*(y3-y1) - (y-y1)*(x3-x1));
                              float v2=((x1-x3)*(y-y3) - (y1-y3)*(x-x3))*((x-x3)*(y2-y3) - (y-y3)*(x2-x3));
                              if ((v1>=0) && (v2>=0))
                                   flag=1;
                         }
                    }
               }
          }
          if (flag==0)
          {
               // if we did not find the triangle, let's take the nearest neighbour
               texDebug[0].item(i)=100.0;
               cout << "+" << flush;
               std::map<float, std::vector<std::pair<float, uint> > >::iterator meshIt=mesh2.begin();
               float xbis, ybis;
               uint index;
               int flagbis=0;
               for ( ; (meshIt!=mesh2.end()) && (flagbis==0) ; ++meshIt)
               {
                    xbis=(*meshIt).first;
                    if (xbis > (x+precisionX))
                         flagbis=1;
                    else if (xbis >= (x-precisionX))
                    {
                         std::vector<std::pair<float, uint> >::iterator yIt=((*meshIt).second).begin();
                         for ( ; yIt!=((*meshIt).second).end(); ++yIt)
                         {
                              ybis=(*yIt).first;
                              //std::cout << "\t\tDEBUG : in interval, y2=" << y2 << endl;
                              if ((ybis>=(y-precisionY)) && (ybis<=(y+precisionY)))
                              {
                                   dist=(x-xbis)*(x-xbis) + (y-ybis)*(y-ybis);
                                   if (dist<distMin)
                                   {
                                        distMin=dist;
                                        index=(*yIt).second;
                                   }
                              }
                         }
                    }
               }
               texOut[0].item(i)=(float) texIn[0].item(index);
           }
           else
           {
           
           // a new interpolation scheme; should be better than before.
           // based on barycentric cordinates
               double t1, t2, t3;
               double sum;
               double l1, l2, l3;
               if (((fabs(x2-x1)<0.0001) && (fabs(x3-x1)<0.0001))
                    || ((fabs(y2-y1)<0.0001) && (fabs(y3-y1)<0.0001)))
               // three points of the triangle have the same x or y-coordinate
               // weird but it can happen on the main coordinate axis
               // (bottom ridge of a sulcus for instance)
               {
                    t1=t2=t3=1.0/3.0;
               }
               else
               {
                    l1=sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
                    l2=sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1));
                    l3=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

                    t1=fabs((x2-x3)*(y-y3) - (y2-y3)*(x-x3))/l1;
                    t2=fabs((x3-x1)*(y-y1) - (y3-y1)*(x-x1))/l2;
                    t3=fabs((x1-x2)*(y-y2) - (y1-y2)*(x-x2))/l3;
                    sum=t1+t2+t3; t1/=sum; t2/=sum; t3/=sum;
               }
               texOut[0].item(i)=(float) t1*texIn[0].item(index1) + t2*texIn[0].item(index2) + t3*texIn[0].item(index3);
           }   
    }

    std::cout << "Writing interpolated texture" << endl;
    // ECRITURE DE LA TEXTURE

    Writer< TimeTexture<float> > wt(adressTexOut);
    wt.write( texOut );
    Writer< TimeTexture<float> > wtd("debug");
    wtd.write( texDebug );
    std::cout << "Done" << endl;
    return(0);
}

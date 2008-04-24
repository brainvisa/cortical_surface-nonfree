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

float min(float a, float b, float c)
{
     if (a<=b)
          if (a<=c) return a;
          else return c;
     else 
          if (b<=c) return b;
          else return c;
}

float max(float a, float b, float c)
{
     if (a>=b)
          if (a>=c) return a;
          else return c;
     else 
          if (b>=c) return b;
          else return c;
}


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
     float px=0;

     AimsApplication     app( argc, argv, "Projet a texture from one mesh onto an atlas using spherical parameterization of both. x-coordinate is the longitude (with a period)");

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
     TimeTexture<float> texOut(1, na), texDebug(1,na), texDebug2(1, ns);
     float x=0.0, y=0.0, x1=0.0, y1=0.0, x2=0.0, y2=0.0, x3=0.0, y3=0.0, xmin, ymin, xmax, ymax;
     uint index1=0, index2=0, index3=0;
     float precisionX=10.0, precisionY=10.0;
     int countNN=0;

     for (uint i=0; i<na; i++)
          texDebug[0].item(i)=0.0;
     for (uint i=0; i<ns; i++)
          texDebug2[0].item(i)=0.0;

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
       
     
     std::cout << "Sorting triangles with peridodicity" << endl;
       
     // the following structure is far too complicated but gives access to 
     // all triangles ordered by first point coordinates and containing all 
     // indexes and coordinates corrected with peridodicity
     std::map<float, std::vector<std::pair<float, AimsVector<std::pair<uint, std::pair<float, float> >, 3> > > > poly2;
     AimsVector<uint, 3> tri;
     AimsVector<std::pair<uint, std::pair<float, float> >, 3> newtri;
     std::pair<float, float> c0, c1, c2;
     std::pair<uint, std::pair<float, float> > p0, p1, p2;
     for (uint i=0; i<nPoly; i++)
     {
          tri=poly[i];
          float x0, y0, x1, y1, x2, y2;
          x0=texIx[0].item(tri[0]); y0=texIy[0].item(tri[0]);
          x1=texIx[0].item(tri[1]); y1=texIy[0].item(tri[1]);
          x2=texIx[0].item(tri[2]); y2=texIy[0].item(tri[2]);
          
          // the following assume that all points on the origin meridian 
          // have a value equal to 0(==px), and that no triangle 'cross'
          // this meridian (one point on each side)
          
          if (fabs(x0) <0.00001) // we are on the origin meridian
          {
               if ((x1>(px-x1)) || (x2>(px-x2)))
               {
                    x0=px; 
                    if (fabs(x1) < 0.00001) x1=px;
                    if (fabs(x2) < 0.00001) x2=px;
               }
               else
               {
                    x0=0.0;
                    if (fabs(x1) < 0.00001) x1=0.0;
                    if (fabs(x2) < 0.00001) x2=0.0;
               }
          }
          else if (fabs(x1) < 0.00001)
          {
               if ((x0>(px-x0)) || (x2>(px-x2)))
               {
                    x1=px;
                    if (fabs(x2) < 0.00001) x2=px;
               }
               else 
               {
                    x1=0.0;
                    if (fabs(x2) < 0.00001) x2=0.0;
               }
          }
          else if (fabs(x2) < 0.00001)
          {
               if ((x0>(px-x0)) || (x1>(px-x1)))
                    x2=px;
               else
                    x2=0.0;
          }
          
          c0=std::pair<float, float>(x0, y0);
          p0=std::pair<uint, std::pair<float, float> >(tri[0], c0);
          c1=std::pair<float, float>(x1, y1);
          p1=std::pair<uint, std::pair<float, float> >(tri[1], c1);
          c2=std::pair<float, float>(x2, y2);
          p2=std::pair<uint, std::pair<float, float> >(tri[2], c2);
          newtri=AimsVector<std::pair<uint, std::pair<float, float> >, 3>(p0 , p1, p2);
          poly2[x0].push_back(std::pair<float, AimsVector<std::pair<uint, std::pair<float, float> >, 3> >(y0, newtri));

     }
//      std::cout << "Building alternate representation of polygons" << endl;
//      std::map<float, std::vector<std::pair<float, AimsVector<uint,3> > > > poly2;
//      AimsVector<uint, 3> tri;
//      for (uint i=0; i<nPoly; i++)
//      {
//           tri=poly[i];
//           int j=tri[0];
//           x=texIx[0].item(j);
//           y=texIy[0].item(j);
//           poly2[x].push_back(std::pair<float, AimsVector<uint,3> >(y, tri));
//      }


    // PARCOURS DES NOEUDS DE L'ATLAS ET POUR CHACUN : INTERPOLATION

    std::cout << "Interpolating texture onto atlas" << endl;

    for (uint i=0; i<na; i++)
    {
          x=texAx[0].item(i);
          y=texAy[0].item(i);
          
          double dist, distMin=10000.0;
          int flag=0;
          int found=0;
//           std::map<float, std::vector<std::pair<float, uint> > >::iterator meshIt=mesh2.begin();
//           for ( ; (meshIt!=mesh2.end()) && (flag==0) ; ++meshIt)
          std::map<float, std::vector<std::pair<float, AimsVector<std::pair<uint, std::pair<float, float> >, 3> > > >::iterator polyIt=poly2.begin();
          for (; (polyIt!=poly2.end()) && (flag==0) ; ++polyIt)
          {
               x2=(*polyIt).first;
               if (x2 > (x+precisionX))
                    flag=1;
               else if (x2 >= (x-precisionX))
               {
                    std::vector<std::pair<float, AimsVector<std::pair<uint, std::pair<float, float> >, 3> > >::iterator yIt=((*polyIt).second).begin();
                    for ( ; (yIt!=((*polyIt).second).end()) && (flag==0); ++yIt)
                    {
                         y2=(*yIt).first;
                         //std::cout << "\t\tDEBUG : in interval, y2=" << y2 << endl;
                         if ((y2>=(y-precisionY)) && (y2<=(y+precisionY)))
                         {
                              AimsVector<std::pair<uint, std::pair<float, float> >, 3> triangle=(*yIt).second;
                              index1=triangle[0].first;
                              index2=triangle[1].first;
                              index3=triangle[2].first;
                              x1=(triangle[0].second).first;y1=(triangle[0].second).second;
                              x2=(triangle[1].second).first;y2=(triangle[1].second).second;
                              x3=(triangle[2].second).first;y3=(triangle[2].second).second;
                              xmin=min(x1,x2,x3); xmax=max(x1,x2,x3);
                              ymin=min(y1,y2,y3); ymax=max(y1,y2,y3);

                              float v1=((x2-x1)*(y-y1) - (y2-y1)*(x-x1))*((x-x1)*(y3-y1) - (y-y1)*(x3-x1));
                              float v2=((x1-x2)*(y-y2) - (y1-y2)*(x-x2))*((x-x2)*(y3-y2) - (y-y2)*(x3-x2));
                              if ((v1>=0) && (v2>=0)) // CAREFUL HERE : THE STRICT > IS VERY IMPORTANT
                                                    // USING >= LEADS TO WRONG RESULTS
                              {
                                   if ((xmin<=x) && (x<=xmax) && (ymin<=y) && (y<=ymax))
                                   {
                                        found++;
                                        if ((i==5932) || (i==5142) || (i==3706) || (i==1125))
                                        {
                                             std::map<int, int> col;
                                             col[5932]=100; col[5142]=200; col[3706]=300; col[1125]=400;
                                             texDebug2[0].item(index1)=col[i];
                                             texDebug2[0].item(index2)=col[i];
                                             texDebug2[0].item(index3)=col[i];
                                        }
                                   }
                              }
                         }
                    }
               }
          }
          
          // DEBUG 
/*          found=0;*/
          //END DEBUG
          
          if (found==0)
          {
               // if we did not find the triangle, let's take the nearest neighbour
               texDebug[0].item(i)=100.0;
               countNN++;
               cout << "NN : (" << x << ", " << y << ")" << endl;
               std::map<float, std::vector<std::pair<float, uint> > >::iterator meshIt=mesh2.begin();
               float xbis, ybis;
               uint index=0;
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
               if (found > 1) 
                    cout << "+ found=" << found << endl;
/*               if (found>3) cout << "\t\t\t i=" << i << endl;*/
//                     texDebug[0].item(i)=float(found);
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
    Writer< TimeTexture<float> > wtd2("debug2");
    wtd2.write( texDebug2 );
    std::cout << "Done" << endl;
    cout << "CountNN=" << countNN << endl;
    return(0);
}

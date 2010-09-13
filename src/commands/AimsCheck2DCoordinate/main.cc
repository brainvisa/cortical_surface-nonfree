/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 *
 *  checking and fixing duplicate coordinates in the 2D surface referential
 */

#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <aims/utility/utility_g.h>
#include <aims/mesh/mesh_g.h>
#include <cortical_surface/mesh/isoLine.h>
#include <aims/io/io_g.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>
#include <aims/distancemap/meshmorphomat_d.h>
#include <aims/connectivity/meshcc.h>
#include <aims/connectivity/meshcc_d.h>
#include <cortical_surface/surfacereferential/gyri/mesh_operations.h>

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileOut, fileMesh, fileLat, fileLon;

  AimsApplication    app( argc, argv, "[UNCOMPLETE DEBUG FUNCTION] Check that a 2D surface-based coordinate system is OK (unicity)" );
  try
  {
     app.addOption( fileMesh, "-m", "input mesh" );
     app.alias( "--mesh", "-m" );
     app.addOption( fileLat, "-x", "latitude texture" );
     app.alias( "--xcoord", "-x" );
     app.addOption( fileLon, "-y", "longitude texture" );
     app.alias( "--ycoord", "-y" );
     app.addOption( fileOut, "-o", "output sign texture" );
     app.alias( "--out", "-o" );
     
     app.initialize();
     
     cout << "reading triangulation   : " << flush;
     AimsSurfaceTriangle surface;
     Reader<AimsSurfaceTriangle> triR( fileMesh );
     triR >> surface;
     cout << "done" << endl;
     
     cout << "reading textures   : " << flush;
     Reader<TimeTexture<float> > texLonR( fileLon );
     TimeTexture<float> texLon;
     texLonR >> texLon ;
     Reader<TimeTexture<float> > texLatR( fileLat );
     TimeTexture<float> texLat;
     texLatR >> texLat ;
     cout << "done " << endl;   
     
     cout << "computing neighbours  " << endl;
     vector<set<uint> >  neigh = SurfaceManip::surfaceNeighbours(surface);
     
     uint i, nVert=surface.vertex().size();
     TimeTexture<short> texOut(1, nVert);
     vector<Point3df>  vert=surface.vertex();

     cout << "computing normals" << endl;
     surface.updateNormals();
     vector<Point3df> norm=surface.normal();
     cout << "OK. Starting check" << endl;

     TimeTexture<float> latNew, lonNew;
     for (i=0; i<nVert; i++)
     {
        latNew[0].item(i)=texLat[0].item(i);
        lonNew[0].item(i)=texLon[0].item(i);
     }

     for (i=0; i< nVert; i++)
     {
          set<uint> v=neigh[i];
          set<uint>::iterator vIt=v.begin();
          uint j, jlat, jlon;
          float glat, glon, latM=-1000.0, lonM=-1000.0;
          
          for ( ; vIt!=v.end(); ++vIt)
          {
               j=*vIt;
               Point3df dist=vert[j]-vert[i];
               glat=(texLat[0].item(j)-texLat[0].item(i))/dist.norm();
               glon=(texLon[0].item(j)-texLon[0].item(i))/dist.norm();
               if (glat>latM) {latM=glat; jlat=j;}
               if (glon>lonM) {lonM=glon; jlon=j;}
          }
          Point3df gradLat, gradLon, vp, n=norm[i];
          float sign;
          gradLat=vert[jlat]-vert[i]; gradLat=gradLat/gradLat.norm();
          gradLon=vert[jlon]-vert[i]; gradLon=gradLon/gradLon.norm();
          n=n/n.norm();
          vp=vectProduct(gradLat, gradLon); 
          sign=vp.dot(n);
          if (sign>=0)
               texOut[0].item(i)=0;
          else if (sign<0)
               texOut[0].item(i)=1;
          else
               texOut[0].item(i)=0;
/*          texOut[0].item(i)=sign;*/
     }

     map<unsigned, set<pair<unsigned,float> > > poids = AimsMeshWeightFiniteElementLaplacian (surface[0], 0.98);    

     cout << "Postprocessing check" << endl;
     for (i=0; i<nVert; i++)
     {
          set<uint> v=neigh[i];
          set<uint>::iterator vIt=v.begin();
          uint j, count=0;
          if (texOut[0].item(i) == 1)
          {
               for ( ; vIt!=v.end(); ++vIt)
               {
                    j=*vIt;
                    if (texOut[0].item(i) == texOut[0].item(j))
                         count++;
               }
               if (count<=1)
               {
                    texOut[0].item(i) = 0;
               }
          }
     }
     cout << "OK. Writing texture texOut.tex" << endl;
     Writer<TimeTexture<short> > texOutW( "texOut.tex" );
     texOutW << texOut ;
     
     cout << "Closing result" << endl;

     TimeTexture<short> dil, closed;
     float sizeClosing=8.0;
     
     dil[0] = MeshDilation<short>( surface[0], texOut[0], short(0), -1 , sizeClosing, true);
     cout << "OK. Writing texture dil.tex" << endl;
     Writer<TimeTexture<short> > dilW( "dil.tex" );
     dilW << dil ;
     closed[0] = MeshErosion<short>( surface[0], dil[0], short(0), -1 , sizeClosing, true);
     cout << "OK. Writing texture closed.tex" << endl;
     Writer<TimeTexture<short> > closedW( "closed.tex" );
     closedW << closed ;
     
     cout << "Extracting connected components" << endl;
     
     TimeTexture<short> compoOut;
     compoOut[0] = AimsMeshLabelConnectedComponent<short>( surface[0], closed[0], 0, 0);

     cout << "OK. Writing texture composantes.tex" << endl;
     Writer<TimeTexture<short> > compoOutW( "composantes.tex" );
     compoOutW << compoOut ;

     std::map<short, std::vector<uint> > compos;
     short val;
     for (i=0; i<nVert; i++)
     {
          val=compoOut[0].item(i);
          if (val>0)
          {
               if (compos.find(val) == compos.end())
               {
                    compos[val]=std::vector<uint>();
                    compos[val].push_back(i);
               }
               else
               {
                    compos[val].push_back(i);
               }

          }
     }

     cout << "Found " << compos.size() << " components";
     cout << "Processing them" << endl;

     TimeTexture<short> texResult(1,nVert);
     TimeTexture<short> texBorder(1,nVert);

     for (i=0; i<nVert; i++) 
     {
        texResult[0].item(i)=0;
        texBorder[0].item(i)=0;
     }
     
     uint j,k, l, sc;
     float lat, lon, latMax, latMin, lonMax, lonMin;
     for (i=1; i<=compos.size(); i++)
     {
          cout << "Composante " << i << endl;
          lonMax=latMax=0.0; lonMin=latMin=360.0;
          std::vector<uint> compo=compos[i];
          std::vector<uint> square;
          sc=compo.size();
          cout << "\t size " << sc << endl;
          
          for (j=0; j<sc; j++)
          {
               k=compo[j];
               lat=texLat[0].item(k);
               lon=texLon[0].item(k);
               if (lat<latMin) latMin=lat;
               if (lat>latMax) latMax=lat;
               if (lon<lonMin) lonMin=lon;
               if (lon>lonMax) lonMax=lon;
          }
          cout << "bounding box : (" << latMin << ", " << lonMin << ") - (" << latMax << ", " << lonMax << ")" << endl;
          if ((lonMax-lonMin) <= (lonMin-lonMax+360.0)) // detection du sillon central ou tout est different
          {
            for (l=0; l<nVert; l++)
            {    
                lat=texLat[0].item(l);
                lon=texLon[0].item(l);
                if ( (lat<=latMax) && (lat>=latMin) && (lon<=lonMax) && (lon>=lonMin))
                {
                      texResult[0].item(l)=i;
                      square.push_back(l);
                }
            }
            uint sq=square.size();
 /**/           cout << "Looking for borders and building reparameterization vectors" << endl;
            vector<uint> inside, top, bottom, left, right;

            for (j=0; j<sq; j++)
            {
                k=square[j];
                int flagBorder=0;
                float lonOut=0.0, latOut=0.0;
                int nOut=0;
                set<uint> v=neigh[k];
                set<uint>::iterator vIt=v.begin();
                for ( ; vIt != v.end(); ++vIt)
                    if (texResult[0].item(*vIt)!=short(i))
                    {
                      //nOut++; //ici on essaye de comprendre de quel bord il s'agit
                      lonOut=texLon[0].item(*vIt); // +=
                      latOut=texLat[0].item(*vIt); // +=
                      flagBorder=1;
                    }
                if (flagBorder==0)
                {
                  inside.push_back(k);
                  texBorder[0].item(k)=5;
                }
                else
                {
                  //texResult[0].item(k)=-1;
                  //latOut/=float(nOut);
                  //lonOut/=float(nOut);
                  if (latOut<latMin) {bottom.push_back(k); texBorder[0].item(k)=1;}
                  else if (latOut>latMax) {top.push_back(k); texBorder[0].item(k)=2;}
                  else if (lonOut<lonMin) {left.push_back(k); texBorder[0].item(k)=3;}
                  else if (lonOut>lonMax) {right.push_back(k); texBorder[0].item(k)=4;}
                  else {texBorder[0].item(k)=5;}
                }
            }/**/
            vector<uint> corr;
            AimsSurfaceTriangle gyrusMesh = getGyrusMesh(surface[0], inside, corr);
            map<unsigned, set<pair<unsigned,float> > > poidsGyrus = getGyrusWeight(poids,inside,corr);

            std::vector < std::pair < std::vector <uint>, short > > constraintVert;
            std::vector < std::pair < std::vector <uint>, short > > constraintHor; 

              pair<vector<uint>, float>    con1(left, texLat[0].item(left[j]));
              constraintVert.push_back(con1);

              pair<vector<uint>, float>    con2(right, texLat[0].item(right[j]));
              constraintVert.push_back(con2);

              pair<vector<uint>, float>    con3(top, texLat[0].item(top[j]));
              constraintHor.push_back(con3);

              pair<vector<uint>, float>    con4(bottom, texLat[0].item(bottom[j]));
              constraintHor.push_back(con4);


            float criter=0.001;
            float dt=0.05;

            cout << "inside : size " << inside.size() << endl;
            cout << "top : size " << top.size() << endl;
            cout << "bottom : size " << bottom.size() << endl;
            cout << "left : size " << left.size() << endl;
            cout << "right : size " << right.size() << endl;
            cout << "Sum = " << inside.size() + top.size() + bottom.size() + left.size() + right.size() << endl;
            cout << "square : size " << square.size() << endl;     

            Texture<float> verticDiff,horizDiff;
            cout << "Starting diffusion" << endl;
            cout << "vertical" << endl;
//            verticDiff = diffusion ( poidsGyrus, gyrusMesh[0], top, bottom, constraintVert, (latMin + latMax)/2.0, inside, corr, criter, dt );
            cout << "horizontal" << endl;
//            horizDiff = diffusion ( poidsGyrus, gyrusMesh[0], left, right, constraintHor, (lonMin + lonMax)/2.0, inside, corr, criter, dt );
            cout << "Reinjecting result on coordinate field" << endl;
          /* FINIR CA ET S'OCCUPER DU SEG FAULT" */
           /* for (uint i2=0;i2<vertices[0].size();i2++)
            {
              tex.item(vertices[0][i2]) = (float) verticDiff.item(corr[vertices[0][i2]]);
                     for (uint i2=0;i2<vertices[0].size();i2++)
                        final[2].item(vertices[0][i2]) = (float) horizDiff.item(corr[vertices[0][i2]]);

            cout << "OK" << endl;   */



          }
     }

     cout << "OK. Writing texture " << fileOut << endl;
     Writer<TimeTexture<short> > texResultW( fileOut );
     texResultW << texResult ;
     cout << "OK. Writing texture borders.tex" << endl;
     Writer<TimeTexture<short> > texBorderW( "borders.tex" );
     texBorderW << texBorder ;
     return EXIT_SUCCESS;
  }
  catch( user_interruption & )
  {
  }
  catch( ... )
  {
    throw;
  }
  return EXIT_FAILURE;
}


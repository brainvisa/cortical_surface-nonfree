/*
 *  Copyright (C) 1997-2005 CEA
 *
 *  This software and supporting documentation were developed by
 *
 *      Laboratoire LSIS, �quipe LXAO,
 *      Marseille, France
 *
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 */

#include <aims/mesh/geometric.h>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>
#include <aims/math/math_g.h>
#include <aims/vector/vector.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/curv.h>
#include <aims/distancemap/meshdistance_d.h>
#include <aims/distancemap/meshmorphomat_d.h>
#include <aims/connectivity/meshcc_d.h>
#include <aims/morphology/morphology_g.h>
#include <cortical_surface/surfacereferential/shortestPath.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/getopt/getopt2.h>
#include <aims/connectivity/meshcc.h>
#include <iostream>
#include <iomanip>

using namespace aims;
using namespace carto;
using namespace std;

#define TOP2BOTTOM 0
#define BACK2FRONT 1
#define RIDGE_TOP 100
#define RIDGE_BOT 200
#define POLE_N 500
#define POLE_S 600

// SEE THE COMMENTS IN FRENCH AT THE END FOR DETAILS OF THE ALGORITHM.
// BUILDING OF THE CONSTRAINTS (POLES, TOP AND BOTTOM RIDGES) IS QUITE 
// A HEAVY TEDIOUS, BUT ROBUST, MACHINERY.

int main( int argc, const char** argv )
{
  try
    {
      string  meshfile, texFile_x, texFile_y, bottomfile, hullfile;
      string  method = "boix";
      float deltaT=0.05;
      float stop=0.000003;

/*      float tCurv=0.95;*/
      int orientation;

      AimsApplication    app( argc, argv, "Parameterize the mesh of a sulcus" );
      app.addOption( meshfile, "-i", "input mesh" );
      app.alias( "--input", "-i" );
      app.addOption( bottomfile, "-b", "sulcus bottom point image" );
      app.alias( "--bottom", "-b" );
      app.addOption( hullfile, "-t", "sulcus top point (junction with brain hull) image" );
      app.alias( "--top", "-t" );
      app.addOption( orientation, "-o", "sulcus main orientation : top->bottom = 0, back->front = 1" );
      app.alias( "--orientation", "-o" );
      app.addOption( texFile_x, "-ox", "1st output texture (depth)" );
      app.alias( "--outputx", "-ox" );
      app.addOption( texFile_y, "-oy", "2nd output texture ('along' the sulcus)" );
      app.alias( "--outputy", "-oy" );
      app.addOption( method, "-m", "curvature computation method (fem/boix/barycenter, default=boix)", "boix" );
      app.alias( "--method", "-m" );
      app.addOption( deltaT, "-d", "diffusion iteration step (default=0.05)", 0.05 );
      app.alias( "--deltaT", "-d" );
      app.addOption( stop, "-s", "laplacian variation stopping criterion (default=0.000003)", 0.000003 );
      app.alias( "--stop", "-s" );
      app.initialize();

      if ( (orientation != TOP2BOTTOM) && (orientation != BACK2FRONT) )
      {
           cerr << "-o option must be set to 0 or 1" << endl;
           return(1);
      }

      //
      // read triangulation
      //

      cout << "reading triangulation   : " << flush;
      AimsSurfaceTriangle surface;
      Reader<AimsSurfaceTriangle> triR( meshfile );
      triR >> surface;
      cout << "done" << endl;

      cout << "reading bottom image  : " << flush;
      AimsData<short> bottom;
      Reader<AimsData<short> > bottomR(bottomfile );
      bottomR >> bottom;
      cout << "done" << endl;
      cout << "reading hull image  : " << flush;
      AimsData<short> hull;
      Reader<AimsData<short> > hullR(hullfile);
      hullR >> hull;
      cout << "done" << endl;

      int x,y,z,sx,sy,sz;
      float dx, dy, dz;

      sx=bottom.dimX(); sy=bottom.dimY(); sz=bottom.dimZ();
      dx=bottom.sizeX(); dy=bottom.sizeY(); dz=bottom.sizeZ();

      AimsData<short> bottomDil(sx, sy, sz, 1, 1), hullDil(sx, sy, sz, 1, 1);
      bottomDil.setSizeXYZT(dx, dy, dz); hullDil.setSizeXYZT(dx, dy, dz);
      std::cout << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << endl;
      std::cout << "sx=" << sx << ", sy=" << sy << ", sz=" << sz << endl;

      uint i, ns=surface.vertex().size();
      TimeTexture<short> texBot(1,ns), texHull(1,ns), texPole(1,ns),
                         topDilation(1,ns), botDilation(1,ns),
                         topClosing(1,ns), botClosing(1,ns),
                         poleDilation(1, ns),
                         topLine(1,ns), botLine(1,ns);

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // computing intersection of mesh with bottom and hull image for selction of "ridges"
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      cout << "Detecting top and bottom ridges" << endl;
      
      for (i=0; i<ns; i++)
      {
          texPole[0].item(i)=short(0);
          texBot[0].item(i)=short(0);
          texHull[0].item(i)=short(0);
          Point3df vert=surface.vertex()[i];
          x=int(floor(vert[0]/dx));
          y=int(floor(vert[1]/dy));
          z=int(floor(vert[2]/dz));
          if (  (hull(x,y,z)!=0) || (hull(x,y,z+1)!=0) || (hull(x,y+1,z)!=0) || (hull(x,y+1,z+1)!=0)
             || (hull(x+1,y,z)!=0) || (hull(x+1,y,z+1)!=0) || (hull(x+1,y+1,z)!=0) || (hull(x+1,y+1,z+1)!=0) )
          {
              texHull[0].item(i)=short(RIDGE_TOP);
          }
          if (  (bottom(x,y,z)!=0) || (bottom(x,y,z+1)!=0) || (bottom(x,y+1,z)!=0) || (bottom(x,y+1,z+1)!=0)
             || (bottom(x+1,y,z)!=0) || (bottom(x+1,y,z+1)!=0) || (bottom(x+1,y+1,z)!=0) || (bottomDil(x+1,y+1,z+1)!=0) )
          {
              texBot[0].item(i)=short(RIDGE_BOT);
          }
      }

     cout << "Done... moving on to postprocessing of ridges" << endl;

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // dilation of ridges
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      topDilation[0]=MeshDilation<short>( surface[0], texHull[0], short(0), -1, 6.0, true);
      botDilation[0]=MeshDilation<short>( surface[0], texBot[0], short(0), -1, 6.0, true);
      topClosing[0]=MeshErosion<short>( surface[0], topDilation[0], short(0), -1, 6.0, true);
      botClosing[0]=MeshErosion<short>( surface[0], botDilation[0], short(0), -1, 6.0, true);


     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // skeletization/thinnning of ridges
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      vector< list<unsigned> > neighbourso(ns);
      neighbourso = AimsMeshOrderNode(surface[0]);
      topLine[0]=MeshSkeletization<short> ( surface[0], topClosing[0], short(RIDGE_TOP), short(0), neighbourso );
      botLine[0]=MeshSkeletization<short> ( surface[0], botClosing[0], short(RIDGE_BOT), short(0), neighbourso );

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // removing branches and triangles from skeletons using a graph shortest path algorithm
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     GraphPath<float> shortest;
     TimeTexture<short> topRidge(1,ns), botRidge(1,ns);
     TimeTexture<float> tmpTop(1, ns), tmpBot(1,ns), tmpRB(1,ns), tmpRT(1,ns);
          // quick convert
     for (i=0; i<ns; i++)
     {
          tmpTop[0].item(i)=float(topLine[0].item(i));
          tmpBot[0].item(i)=float(botLine[0].item(i));
     }

     // Selecting start and end point for each ridge
     // first selecting candidates than chosing them with sulci orientation

     std::vector<std::set<uint> > voisins;
     std::set<uint> topCand, botCand;
     voisins=SurfaceManip::surfaceNeighbours( surface );
     for (i=0; i<ns; i++)
     {
          if ( topLine[0].item(i) != 0 )
          {
               std::set<uint> setv=voisins[i], cand;
               std::set<uint>::iterator itV=setv.begin();
               int nbV=0;
               for ( ; itV != setv.end(); ++itV)
               {
                    if (topLine[0].item(*itV) != 0)
                    {
                         nbV++;
                         cand.insert(*itV);
                    }
               }
               if (nbV==1)
                    topCand.insert(i);
               else if (nbV==2)
               {
                    uint v1, v2;
                    itV=cand.begin(); v1=(*itV);
                    ++itV; v2=(*itV);
                    setv=voisins[v1];
                    if (setv.find(v2) != setv.end())
                    {
                         topCand.insert(i);
                    }
               }
          }
          if ( botLine[0].item(i) != 0 )
          {
               std::set<uint> setv=voisins[i], cand;
               std::set<uint>::iterator itV=setv.begin();
               int nbV=0;
               for ( ; itV != setv.end(); ++itV)
               {
                    if (botLine[0].item(*itV) != 0)
                    {
                         nbV++;
                         cand.insert(*itV);
                    }
               }
               if (nbV==1)
                    botCand.insert(i);
               else if (nbV==2)
               {
                    uint v1, v2;
                    itV=cand.begin(); v1=(*itV);
                    ++itV; v2=(*itV);
                    setv=voisins[v1];
                    if (setv.find(v2) != setv.end())
                    {
                         botCand.insert(i);
                    }
               }
          }
     }

     // the choice of the start/end points is done amongst the candidates.
     // for each pair of points the shortest path within the skezleton is
     // computed. The pair having the longest of the shortest path is the one

     std::set<uint>::iterator topIt, topIt2, botIt, botIt2;
     uint cand1, cand2, tmp1, tmp2;
     float length, lengthmax=0.0;
     uint its, itn, ibs, ibn;
     
     for (topIt=topCand.begin(); topIt!=topCand.end(); ++topIt)
     {
          tmp1=(*topIt);
          topIt2=topIt; ++topIt2;
          for ( ; topIt2!=topCand.end(); ++topIt2 )
          {
               tmp2=(*topIt2);
               length=shortest.getLongueur(tmpTop, surface, RIDGE_TOP, tmp1, tmp2);
               if (length > lengthmax)
               {
                    lengthmax=length;
                    cand1=tmp1;
                    cand2=tmp2;
               }
          }
     }

     if (orientation==TOP2BOTTOM)
     {
          if ((surface.vertex()[cand1])[2] > (surface.vertex()[cand2])[2])
          {
               itn=cand2;
               its=cand1;
          }
          else
          {
               itn=cand1;
               its=cand2;
          }
     }
     else if (orientation == BACK2FRONT)
     {
          if ((surface.vertex()[cand1])[1] > (surface.vertex()[cand2])[1])
          {
               itn=cand2;
               its=cand1;
          }
          else
          {
               itn=cand1;
               its=cand2;
          }
     }
     lengthmax=0.0;
     
     for (botIt=botCand.begin(); botIt!=botCand.end(); ++botIt)
     {
          tmp1=(*botIt);
          botIt2=botIt; ++botIt2;
          for ( ; botIt2!=botCand.end(); ++botIt2 )
          {
               tmp2=(*botIt2);
               length=shortest.getLongueur(tmpBot, surface, RIDGE_BOT, tmp1, tmp2);
               if (length > lengthmax)
               {
                    lengthmax=length;
                    cand1=tmp1;
                    cand2=tmp2;
               }
          }
     }
     if (orientation==TOP2BOTTOM)
     {
          if ((surface.vertex()[cand1])[2] > (surface.vertex()[cand2])[2])
          {
               ibn=cand2;
               ibs=cand1;
          }
          else
          {
               ibn=cand1;
               ibs=cand2;
          }
     }
     else if (orientation == BACK2FRONT)
     {
          if ((surface.vertex()[cand1])[1] > (surface.vertex()[cand2])[1])
          {
               ibn=cand2;
               ibs=cand1;
          }
          else
          {
               ibn=cand1;
               ibs=cand2;
          }
     }

     float zf;
//     , ztn=10000, zbn=10000, zts=-10000, zbs=-10000;
     float yf;
//     , ytn=10000, ybn=10000, yts=-10000, ybs=-10000;
     float xf;
     Point3df vert;


     GraphPath<float> shortest2;
          
     tmpRT=shortest2.process(tmpTop, surface, RIDGE_TOP, itn, its);
     tmpRB=shortest2.process(tmpBot, surface, RIDGE_BOT, ibn, ibs);
          
          //quick convertBack
     for (i=0; i<ns; i++)
     {
          topRidge[0].item(i)=int(floor(tmpRT[0].item(i) + 0.5));
          botRidge[0].item(i)=int(floor(tmpRB[0].item(i) + 0.5));
     }

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // computation of pole seeds
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     float xnt, xnb, xst, xsb, xn, xs;
     float ynt, ynb, yst, ysb, yn, ys;
     float znt, znb, zst, zsb, zn, zs;
     int poleS;
     int poleN;
     float dist1, dist2, distS=10000.0, distN=10000.0;
     vert=surface.vertex()[itn]; xnt=vert[0]; ynt=vert[1]; znt=vert[2];
     vert=surface.vertex()[its]; xst=vert[0]; yst=vert[1]; zst=vert[2];
     vert=surface.vertex()[ibn]; xnb=vert[0]; ynb=vert[1]; znb=vert[2];
     vert=surface.vertex()[ibs]; xsb=vert[0]; ysb=vert[1]; zsb=vert[2];
     xn=(xnt+xnb)/2.0; yn=(ynt+ynb)/2.0; zn=(znt+znb)/2.0;
     xs=(xst+xsb)/2.0; ys=(yst+ysb)/2.0; zs=(zst+zsb)/2.0;
     
     for (i=0; i<ns; i++)
     {
          vert=surface.vertex()[i];
          xf=vert[0]; yf=vert[1]; zf=vert[2];
          dist1=sqrt( (xn-xf)*(xn-xf) + (yn-yf)*(yn-yf) + (zn-zf)*(zn-zf) );
          dist2=sqrt( (xs-xf)*(xs-xf) + (ys-yf)*(ys-yf) + (zs-zf)*(zs-zf) );
          if (dist1 < distN) { poleN=i; distN=dist1; }
          if (dist2 < distS) { poleS=i; distS=dist2; }
     }
     texPole[0].item(poleN)=POLE_N;
     texPole[0].item(poleS)=POLE_S;

     poleDilation[0]=MeshDilation<short>( surface[0], texPole[0], short(0), -1, 5.0, true);

    std::set<uint> nord, sud;
    for (i=0; i<ns; i++)
    {
	if (poleDilation[0].item(i)==POLE_N)
	    nord.insert(i);
	else if (poleDilation[0].item(i)==POLE_S)
	    sud.insert(i);
    }
        
    // looking for two points in the pole dilations that are the furthest away from eachother.
    std::set<uint>::iterator nordIt, sudIt;
    uint pS, pN;
    float d, dP=0;
    Point3df vS, vN;
    for (sudIt=sud.begin(); sudIt!=sud.end(); ++sudIt)
	for (nordIt=nord.begin(); nordIt!=nord.end(); ++nordIt)
	{
	    vS=surface.vertex()[*sudIt];
	    vN=surface.vertex()[*nordIt];
	    d=sqrt( (vS[0]-vN[0])*(vS[0]-vN[0]) + 
		    (vS[1]-vN[1])*(vS[1]-vN[1]) +
		    (vS[2]-vN[2])*(vS[2]-vN[2]) ); 
	    if (d>dP) {dP=d; pS=*sudIt; pN=*nordIt;}
	}
    for (i=0; i<ns; i++)
    {
          if (i==pS)
	          poleDilation[0].item(i)=POLE_S;
	     else if (i==pN)
	          poleDilation[0].item(i)=POLE_N;
	     else
	          poleDilation[0].item(i)=0;
    }

     
     GraphPath<short> shortest3;
     TimeTexture<short> tmpNord(1, ns), tmpSud(1, ns), ptmpS(1, ns), ptmpN(1, ns), diltmpN(1, ns), diltmpS(1, ns);
     
     for (i=0; i<ns; i++)
          if (i==pS)
          {
               ptmpS[0].item(i)=1;
               ptmpN[0].item(i)=0;
          }
          else if (i==pN)
          {
               ptmpS[0].item(i)=0;
               ptmpN[0].item(i)=1;
          }
          else
          {
               ptmpS[0].item(i)=0;
               ptmpN[0].item(i)=0;
          }
     diltmpN[0]=MeshDilation<short>( surface[0], ptmpN[0], short(0), -1, 5.0, true);
     diltmpS[0]=MeshDilation<short>( surface[0], ptmpS[0], short(0), -1, 5.0, true);
      
     tmpNord=shortest3.process(diltmpN, surface, 1 , int(itn), int(pN));
     tmpSud=shortest3.process(diltmpS, surface, 1 , int(its), int(pS));
     
     for (i=0; i<ns; i++)
     {
          if ((tmpNord[0].item(i)!=0) || (tmpSud[0].item(i)!=0))
               topRidge[0].item(i)=RIDGE_TOP;
     }

     topDilation[0]=MeshDilation<short>( surface[0], topRidge[0], short(0), -1, 10.0, true);
     topClosing[0]=MeshErosion<short>( surface[0], topDilation[0], short(0), -1, 9.0, true);

//      cout << "writing texture topClosing.tex : " << flush;
//      Writer<TimeTexture<short> >  topClosingW( "topClosing.tex" );
//      topClosingW.write( topClosing );
//      cout << "done " << endl;
     
     topLine[0]=MeshSkeletization<short> ( surface[0], topClosing[0], short(RIDGE_TOP), short(0), neighbourso );
     topLine[0].item(pS)=RIDGE_TOP;
     topLine[0].item(pN)=RIDGE_TOP;
     
     list<unsigned> neighN=neighbourso[pN], neighS=neighbourso[pS];
     list<unsigned>::iterator itN=neighN.begin(), itS=neighS.begin();
     int flag=0;
     for (; itN!=neighN.end(); ++itN)
          if ((topLine[0].item(*itN) == RIDGE_TOP) && (*itN != pN))
               flag=1;
     if (flag==0)
          cout << "North pole disconnected" << endl;
     else cout << "North pole connected OK" << endl;
     flag=0;
     for (; itS!=neighS.end(); ++itS)
          if ((topLine[0].item(*itS) == RIDGE_TOP) && (*itS != pS))
               flag=1;
     if (flag==0)
          cout << "South pole disconnected" << endl;
     else cout << "South pole connected OK" << endl;
     
//      cout << "writing texture topLine.tex : " << flush;
//      Writer<TimeTexture<short> >  topLineW( "topLine.tex" );
//      topLineW.write( topLine );
//      cout << "done " << endl;

     topRidge=shortest3.process(topLine, surface, RIDGE_TOP , int(pS), int(pN));
     topRidge[0].item(pS)=0;
     topRidge[0].item(pN)=0;
  
     cout << "Generating constraints" << endl;

     TimeTexture<float> coord_x(1,ns), coord_y(1,ns);
     TimeTexture<float> init(1,ns), dist(1,ns);
     TimeTexture<short> poleCC(1,ns), topCC, botCC;


     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // Building the two sets of constraints
     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      std::set<unsigned int> pole1, pole2, topR, botR;
      std::set<unsigned int>::iterator p1It, p2It;
      float p1v=0.0, p2v=100.0;
      for (i=0; i<ns; i++)
      {
           if (poleDilation[0].item(i)==POLE_N)
           {
                dist[0].item(i)=0;
                coord_x[0].item(i)=-100.0;
                coord_y[0].item(i)=p1v;
                pole1.insert(i);
           }
           else if (poleDilation[0].item(i)==POLE_S)
           {
                dist[0].item(i)=10;
                coord_x[0].item(i)=-100.0;
                coord_y[0].item(i)=p2v;
                pole2.insert(i);
           }
           else if (topRidge[0].item(i)>0)
           {
                dist[0].item(i)=0;
                coord_x[0].item(i)=0.0;
                coord_y[0].item(i)=50.0;
                topR.insert(i);
           }
          else if (botRidge[0].item(i)>0)
          {
                dist[0].item(i)=0;
                coord_x[0].item(i)=100.0;
                coord_y[0].item(i)=50.0;
                botR.insert(i);
          }
          else
          {
                dist[0].item(i)=0;
                coord_x[0].item(i)=50.0;
                coord_y[0].item(i)=50.0;
          }
     }

//      cout << "writing texture topRidge.tex : " << flush;
//      Writer<TimeTexture<short> >  topRidgeW( "topRidge.tex" );
//      topRidgeW.write( topRidge );
//      cout << "done " << endl;
//       cout << "writing texture bottomRidge.tex : " << flush;
//       Writer<TimeTexture<short> >  botRidgeW( "bottomRidge.tex" );
//       botRidgeW.write( botRidge );
//       cout << "done " << endl;
//      cout << "writing texture poles.tex : " << flush;
//      Writer<TimeTexture<short> >  polesW( "poles.tex" );
//      polesW.write( poleDilation );
//      cout << "done " << endl;


      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Diffusion
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      // distance map for y initialisation

      init[0]=meshdistance::MeshDistance<float>(surface[0], dist[0], true);
      float distM=-10;
      for (i=0; i<ns; i++)
           if (init[0].item(i)>distM)
                distM=init[0].item(i);
      for (i=0; i<ns; i++)
           if (fabs(coord_y[0].item(i) - 50)<1)
                coord_y[0].item(i)=100.0 - (100.0*init[0].item(i)/distM);

      // Getting the list of weights for laplacian computation

      std::map<unsigned, std::set< std::pair<unsigned,float> > > weights;
      weights=AimsMeshWeightFiniteElementLaplacian( surface[0], 0.9);

      // propagating the y coordinate, pole to pole

      TimeTexture<float> laplacian(1, ns), laplacian360(1,ns);

      cout << "Diffusing y coordinate" << endl;
      int iteration=0;
      std::vector<float> laplM;
      float lMax=0, lMaxP=0;
      int flagOut=0;
      float Lmean=0, LmeanP=-10000;
      while (flagOut<1)
      {
           lMax=-10000.0;
           Lmean=0;
           laplacian[0]=AimsMeshLaplacian( coord_y[0], weights);
           for (i=0; i<ns; i++)
           {
		    coord_y[0].item(i)=coord_y[0].item(i)+deltaT*laplacian[0].item(i);
                if ((coord_y[0].item(i) != 0) && (coord_y[0].item(i) != 100) && (laplacian[0].item(i) > lMax))
                     lMax=laplacian[0].item(i);
                Lmean+=laplacian[0].item(i);

           }
           Lmean=Lmean/float(ns);

           p1It=pole1.begin(), p2It=pole2.begin();
           for ( ; p1It!=pole1.end() ; ++p1It)
                coord_y[0].item(*p1It)=p1v;
           for ( ; p2It!=pole2.end() ; ++p2It)
                coord_y[0].item(*p2It)=p2v;

           iteration++;

           if ((iteration%200) == 0)
           {
		cout << fabs(Lmean-LmeanP) << " | " << flush;
                if (fabs(Lmean - LmeanP)<=stop)
                {
                     cout << "reached at iteration " << iteration << endl;
                     flagOut++;
                }
                lMaxP=lMax;
                LmeanP=Lmean;
           }
     }

      cout << "\nDiffusion of y coordinate stopped after " << iteration << " iterations" << endl;

      // for propagation in the other direction, poles have to be virtually 'removed'
     // in the same time we identify the two sides of the origin meridian, that are concerned
     // by the periodicity of the coordinates

      TimeTexture<short> sides(1,ns), split(1,ns);
      for (i=0; i<ns; i++)
      {
          sides[0].item(i)=0;
          split[0].item(i)=0;
      }
      std::map<unsigned, std::set< std::pair<unsigned,float> > >::iterator itW=weights.begin();
      std::map<unsigned, std::set< std::pair<unsigned, float> > > newWeights;
      for ( ; itW!=weights.end(); ++itW)
      {
           i=(*itW).first;
           std::set< std::pair<unsigned,float> > setW;
           std::set< std::pair<unsigned,float> >::iterator itN=weights[i].begin();
           for ( ; itN!=weights[i].end(); ++itN)
           {
               unsigned neigh=(*itN).first;
               float poids=(*itN).second;
               if (poleDilation[0].item(neigh) > 0)
                    poids=0.0;
               setW.insert(std::pair<unsigned,float>(neigh, poids));
               if ((topRidge[0].item(i)>0) && (topRidge[0].item(neigh) == 0) && (poleDilation[0].item(neigh)==0) )
               {
                    sides[0].item(neigh)=10;
               }
           }
           newWeights[i]=setW;
      }

      // splitting 'sides' in two separate components


     split[0]=AimsMeshLabelConnectedComponent<short>( surface[0], sides[0], 9, 0);

     unsigned nbSides=AimsMeshLabelNbConnectedComponent<short>( surface[0], sides[0], 10 );
     if (nbSides < 2)
     {
          cerr << "Finding " << nbSides << " sides to the principal meridian (should be 2). Exiting ..." << endl;
          exit(1);
     }
     else if (nbSides > 2)
     {
          std::vector<int> conn;
          for (i=0; i<nbSides; i++)
               conn.push_back(0);
          cerr << "Finding " << nbSides << " sides to the principal meridian. Keeping the two biggest ones" << endl;
          for (i=0; i<ns; i++)
               if (split[0].item(i) > 0)
                    conn[split[0].item(i) - 1]++;
          int label1, label2, n1=0, n2=0;
          for (i=0; i<nbSides; i++)
          {
               if (conn[i] >= n2)
               {
                    if (conn[i] >= n1)
                    {
                         n1=conn[i];
                         label1=i+1;
                    }
                    else
                    {
                         n2=conn[i];
                         label2=i+1;
                    }
               }
          }
          cerr << "keeping labels " << label1 << " and " << label2 << endl;
          for (i=0; i<ns; i++)
          {
               if ((split[0].item(i) != label1) && (split[0].item(i) != label2))
                    sides[0].item(i)=0;
          }
          split[0]=AimsMeshLabelConnectedComponent<short>( surface[0], sides[0], 9, 0);
     }

     //propagating the x coordinate
     cout << "Diffusing x coordinate" << endl;
     iteration=0;
     lMax=0; lMaxP=0;
     Lmean=0, LmeanP=-10000;
     flagOut=0;
     std::vector<Point3df> vertV=surface[0].vertex();

     FILE *laplF;
//      laplF=fopen("laplacienX.txt", "w");

     while (flagOut<1)
     {
          lMax=-10000.0;
          Lmean=0;
          // computing the 'periodic' laplacian
          std::map<unsigned, std::set< std::pair<unsigned,float> > >::const_iterator il,el;
          int cpt=0;
          std::set< std::pair<unsigned,float> >::iterator ip,ep;
          unsigned node, voisin;
          Point3df temp, temp_current;
          float L,weight;
          for (il=newWeights.begin(), el=newWeights.end(); il!=el; ++il, cpt++)
          {
               node = il->first;
               L = 0;
               temp=vertV[cpt];

               for ( ip = (il->second).begin(), ep = (il->second).end(); ip != ep; ++ip    )
               {
                    temp_current=vertV[(ip->first)];
                    voisin = ip->first;
                    weight = ip->second;

                    // periodicity management around the meridian
                    //Side chnge with left or right hemisphere -- not managed yet
                    if(( split[0].item(node)==1)  ) // && (botR.find(voisin)!=botR.end())
                    //1 for left hemisphere, 2 for right hemisphere
                         L += weight * (- coord_x[0].item(node) + 200 );
                    else
                         if( (split[0].item(node)==2) && (botR.find(voisin)!=botR.end()) )
                         //2 for left hemisphere, 1 for righthemisphere
                              L += weight * (-coord_x[0].item(node));
                    //for all other nodes
                    else
                         L += weight * (- coord_x[0].item(node) + coord_x[0].item(voisin));
               }
               laplacian360[0].item(node) = L;
               Lmean+=L;
          }
          Lmean=Lmean/float(newWeights.size());

          for (i=0; i<ns; i++)
           {
                coord_x[0].item(i)=coord_x[0].item(i)+deltaT*laplacian360[0].item(i);
                if ((coord_x[0].item(i) != 0) && (coord_x[0].item(i) != 100) && (laplacian360[0].item(i) > lMax))
                     lMax=laplacian360[0].item(i);
           }

           std::set<unsigned>::iterator trIt=topR.begin(), brIt=botR.begin();
           for ( ; trIt!=topR.end() ; ++trIt)
                coord_x[0].item(*trIt)=0;
           for ( ; brIt!=botR.end() ; ++brIt)
                coord_x[0].item(*brIt)=100;
/*           fprintf(laplF, "%.6f\n", Lmean);*/
           iteration++;
           if ((iteration%200) == 0)
           {
                cout << fabs(Lmean-LmeanP) << " | " << flush;
                if (fabs(Lmean-LmeanP)<=stop)
                {
                     cout << "reached at iteration " << iteration << endl;
                     flagOut++;
                }
                lMaxP=lMax;
                LmeanP=Lmean;
           }
     }
/*     fclose(laplF);*/
      cout << "\nDiffusion of x coordinate stopped after " << iteration << " iterations" << endl;

      cout << "writing texture " << texFile_x << " : " << flush;
      Writer<Texture1d>  texxW( texFile_x );
      texxW.write( coord_x);
      cout << "done " << endl;
      cout << "writing texture " << texFile_y << " : " << flush;
      Writer<Texture1d>  texyW( texFile_y );
      texyW.write( coord_y);
      cout << "done " << endl;
      return( 0 );
    }

/*  catch( user_interruption & )
    {
    }*/
  catch( exception & e )
    {
      cerr << e.what() << endl;
    }
  return 1;
}

/* QUELQUES DETAILS SUR L'USINE A GAZ

Top et Bottom ridges :
----------------------

- intersection buckets/maillage
- => textures texBot et texHull à valeurs 0 ou RIDGE_TOP/RIDGE_BOT

- dilatation des deux textures => topDilation/botDilation
- fermeture des deux textures => topClosing/botClosing

- squelettisation => topLine/botLine

- enlève les branches avec un plus court chemin :
     - sélection start/end points avec le "plus long plus court chemin"
     => ibn/ibs et itn/its
     - plus court chemin itn->its ou ibn->ibs en restant sur la valeur RIDGE_TOP ou RIDGE_BOT

=> textures topRidge et botRIdge


poles : 
-------

- 2 points calculés comme le milieu de [itn, ibn] et [its, ibs] => on prend les 2 points du maillage les plus près de ca => poleN et poleS
- texture de poles à valeurs POLE_N et POLE_S sur les deux points : texPole
- dilatation de texPole => poleDilation et deux listes de points, nord et sud. 
- le couple de points de pole et sud les plus éloignés sont les poles définitifs et poleDIlation ne contient plus qu'eux. 


finalisation du méridien 0/360 :
------------------------------

- On relie le topRidge aux poles
- On ferme le résultat
- Pour chacun on refait suqlette puis le plus court chemin d'un pole à l'autre.



*/

/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 *
 *  Just my own little binary for various purposes
 */

#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <aims/utility/utility_g.h>
#include <aims/mesh/mesh_g.h>
#include <cortical_surface/mesh/isoLine.h>
#include <aims/io/io_g.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <cortical_surface/surfacereferential/gyri/mesh_operations.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>

using namespace aims;
using namespace carto;
using namespace std;



int main( int argc, const char** argv )
{
  string fileOut, fileMesh, fileGyri, fileLat, fileLon;
  float gyrus;
  float dt=0.01, delta=0.0001;

  AimsApplication    app( argc, argv, "Extract a particular gyrus from a mesh and its gyri texture" );
  app.addOption( fileMesh, "-m", "input mesh" );
  app.alias( "--mesh", "-m" );
  app.addOption( fileGyri, "-t", "gyri texture" );
  app.alias( "--tex", "-t" );
  app.addOption( gyrus, "-g", "gyrus index" );
  app.alias( "--gyrus", "-g" );
  app.addOption( fileLat, "-lat", "latitude texture");
  app.alias( "--latitude", "-lat");
  app.addOption( fileLon, "-lon", "longitude texture");
  app.alias( "--longitude", "-lon");
  app.addOption( fileOut, "-o", "output gyrus name pattern" );
  app.alias( "--out", "-o" );
  app.addOption(dt, "-dt", "diffusion time step (default=0.01)", 0.01);
  app.addOption(delta, "-delta", "diffusion stoping criterion (default=0.0001)", 0.0001);
  app.initialize();
  
  cout << "reading triangulation   : " << flush;
  AimsSurfaceTriangle surface;
  Reader<AimsSurfaceTriangle> triR( fileMesh );
  triR >> surface;
  cout << "done" << endl;

  cout << "reading gyri texture   : " << flush;
  Reader<TimeTexture<float> > texGyriR( fileGyri );
  TimeTexture<float> texG;
  texGyriR >> texG ;
  cout << "done " << endl;

  cout << "reading longitude and latitude textures   : " << flush;
  Reader<TimeTexture<float> > texLatR( fileLat );
  TimeTexture<float> texLat;
  texLatR >> texLat;
  Reader<TimeTexture<float> > texLonR( fileLon );
  TimeTexture<float> texLon;
  texLonR >> texLon;
  cout << "done " << endl;

  std::vector< Point3df > vertG;
  std::vector< uint > gyrusV;
  std::vector< AimsVector< uint, 3 > >  polyG;
  std::set<uint> cortRegion;

  std::vector<Point3df> vert=surface.vertex();
  uint nv=vert.size();
  std::vector< AimsVector< uint, 3 > > poly=surface.polygon();
  uint np=poly.size();
  uint i;
  uint count=0;
  std::map<uint, uint> trans;
  
  std::cout << "gyrus=" << gyrus << std::endl;

  cout << "Extracting nodes within " << nv << endl;
  for (i=0; i<nv; ++i)
  {
//	std::cout << texG[0].item(i) << std::endl;
    if (fabs(texG[0].item(i) - gyrus)<0.1)
    {
      vertG.push_back(vert[i]);
      trans[i]=count;
      cortRegion.insert(i);
      gyrusV.push_back(i); // this is all very redundant but does not matter so much
						  // no time to optimize, sorry.
      count++;
    }
  }

  std::cout << "count=" << count << std::endl;
  std::cout << "Reparameterizing gyrus" << std::endl;

  cout << "\t computing neighbours and borders" << endl;
  vector<set<uint> >  neigh = SurfaceManip::surfaceNeighbours(surface);
  std::map<int, set<uint> > borders;   // map(valeur derrière le bord, noeud)
  std::set<uint>::iterator gyrusIt=cortRegion.begin();

  for ( ; gyrusIt != cortRegion.end(); ++gyrusIt)
  {
	  i=*gyrusIt;

	  set<uint> vois=neigh[i];
	  set<uint>::iterator voisIt=vois.begin();
	  for (; voisIt!= vois.end(); ++voisIt)
	  {
		  if (fabs(texG[0].item(*voisIt) - gyrus) >= 0.1)
		  {
			  int val=(int) (texG[0].item(*voisIt)+0.5);
			  if (borders.find(val)!=borders.end())
				  borders[val].insert(i);
			  else
			  {
				  borders[val]=set<uint>();
				  borders[val].insert(i);
			  }
		  }
	  }
  }

  cout << "\t Sorting and merging borders" << std::endl;
  std::map<int, float> latA, latS, lonA, lonS;
  std::map<int, set<uint> >::iterator bIt;
  std::set<int> rmv;

  for (bIt=borders.begin(); bIt !=borders.end(); ++bIt)
	  if ((*bIt).second.size()<=1)
		  rmv.insert((*bIt).first);

  std::set<int>::iterator rmIt;
  for (rmIt=rmv.begin(); rmIt!=rmv.end(); rmIt++)
	  borders.erase(*rmIt);

  // computing mean and sd of lon and lat along each border
  for (bIt=borders.begin(); bIt !=borders.end(); ++bIt)
  {
	  int val=(*bIt).first;
	  if ((*bIt).second.size()>1)
	  {
		  std::set<uint>::iterator sIt=(*bIt).second.begin();
		  float mLat, mLon, sLat, sLon;
		  mLat=sLat=mLon=sLon=0.0;
		  for ( ; sIt != (*bIt).second.end(); ++sIt)
		  {
			  i=*sIt;
			  float lat=texLat[0].item(i);
			  float lon=texLon[0].item(i);
			  mLat+=lat; mLon+=lon;
		  }
		  mLat/=(float) (*bIt).second.size();
		  mLon/=(float) (*bIt).second.size();
		  for ( sIt=(*bIt).second.begin() ; sIt != (*bIt).second.end(); ++sIt)
		  {
			  i=*sIt;
			  float lat=texLat[0].item(i);
			  float lon=texLon[0].item(i);
			  sLat+=(lat-mLat)*(lat-mLat); sLon+=(lon-mLon)*(lon-mLon);
		  }
		  sLat=sqrt(sLat/((float) (*bIt).second.size()));
		  sLon=sqrt(sLon/((float) (*bIt).second.size()));
		  latA[val]=mLat; lonA[val]=mLon;
		  latS[val]=sLat; lonS[val]=sLon;
	  }
  }


  // setting borders to longitude or latitude depending on their stats
  std::map<int, float>::iterator gIt;
  set<int> borLon, borLat;

  for (gIt=latA.begin(); gIt!=latA.end(); ++gIt)
  {
	  int val=(*gIt).first;
	  if (lonS[val] < latS[val])
		  borLon.insert(val);
	  else if (latS[val] <= lonS[val])
		  borLat.insert(val);
  }

  set<int>::iterator borIt;
  std::cout << "\t equi-longitude borders" << std::endl;
  // here we merge borders which are actually the same, along longitude (then latitude below)
  for (borIt=borLon.begin(); borIt!=borLon.end(); ++borIt)
  {
	  int val=(*borIt);
	  if (borders.find(val)!=borders.end()) {
//	  std::cout << val << " : " << borders[val].size() << " points" << std::endl;
//	  std::cout << "latA="<< latA[val] << ", lonA=" << lonA[val] << std::endl;
//	  std::cout << "latS="<< latS[val] << ", lonS=" << lonS[val] << std::endl;

	  set<int>::iterator tmpIt=borIt;
	  ++tmpIt;
	  if (tmpIt!=borLon.end())
	  {
		  for ( ; tmpIt!=borLon.end(); ++tmpIt)
		  {
			  int val2=(*tmpIt);
			  if (borders.find(val2)!=borders.end()) {
			  if (fabs(lonA[val]-lonA[val2])<1.0)
			  {
//				  std::cout << "\t merging with " << val2 << std::endl;
				  set<uint>::iterator tmp2It=borders[val2].begin();
				  for ( ; tmp2It!=borders[val2].end(); ++tmp2It)
					  borders[val].insert(*tmp2It);
				  borders.erase(val2);
			  } }
		  }
	  } }
  }
  std::cout << "equi-latitude borders" << std::endl;
  // same with latitude below

   for (borIt=borLat.begin(); borIt!=borLat.end(); ++borIt)
   {
 	  int val=(*borIt);
	  if (borders.find(val)!=borders.end()) {
//	  std::cout << val << " : " << borders[val].size() << " points" << std::endl;
// 	  std::cout << "latA="<< latA[val] << ", lonA=" << lonA[val] << std::endl;
// 	  std::cout << "latS="<< latS[val] << ", lonS=" << lonS[val] << std::endl;

 	  set<int>::iterator tmpIt=borIt;
	  ++tmpIt;
	  if (tmpIt!=borLon.end())
	  {
		  for ( ; tmpIt!=borLat.end(); ++tmpIt)
		  {
			  int val2=(*tmpIt);
			  if (borders.find(val2)!=borders.end()) {
			  if (fabs(latA[val]-latA[val2])<1.0)
			  {
//				  std::cout << "\t merging with " << val2 << std::endl;
				  set<uint>::iterator tmp2It=borders[val2].begin();
				  for ( ; tmp2It!=borders[val2].end(); ++tmp2It)
					  borders[val].insert(*tmp2It);
				  borders.erase(val2);
			  } }
		  }
	  } }

   }

   vector<int> lonP, latP;
   for (bIt=borders.begin(); bIt !=borders.end(); ++bIt)
   {
	   int val=(*bIt).first;
	   if (borLon.find(val)!=borLon.end()) lonP.push_back(val);
	   else if (borLat.find(val)!=borLat.end()) latP.push_back(val);
   }
   std::cout << "\t found " << borders.size() << " borders" << std::endl;
   std::cout << "\t with " << lonP.size() << " iso-longitude and " << latP.size() << " iso-latitude" << std::endl;
   if ((borders.size()!=4) || (lonP.size()!=2) || (latP.size()!=2)) {std::cerr << "Should be 4 (2+2). This is usually due to non-unicity of coordinates in the original referential. Sorry... Stopping." << std::endl; exit(EXIT_FAILURE);}

   // building the constraints from each border for diffusion
   std::vector<uint> hautLat, basLat, hautLon, basLon;
   int h, b;
   std::set<uint>::iterator sIt;
   if (lonA[lonP[0]] < lonA[lonP[1]])
   { b=lonP[0]; h=lonP[1];}
   else
   { b=lonP[1]; h=lonP[0];}
   for (sIt=borders[b].begin(); sIt != borders[b].end(); ++sIt)
	   basLon.push_back(*sIt);
   for (sIt=borders[h].begin(); sIt != borders[h].end(); ++sIt)
   	   hautLon.push_back(*sIt);
   if (latA[latP[0]] < latA[latP[1]])
   { b=latP[0]; h=latP[1];}
   else
   { b=latP[1]; h=latP[0];}
   for (sIt=borders[b].begin(); sIt != borders[b].end(); ++sIt)
	   basLat.push_back(*sIt);
   for (sIt=borders[h].begin(); sIt != borders[h].end(); ++sIt)
   	   hautLat.push_back(*sIt);

   cout << "\t Recomputing polygons of new mesh" << endl;

    for (i=0; i<np; i++)
    {
    	if  (  (cortRegion.find(poly[i][0])	 != cortRegion.end())
    		&& (cortRegion.find(poly[i][1])	 != cortRegion.end())
    		&& (cortRegion.find(poly[i][2])	 != cortRegion.end()) )
//      if ( (fabs(texG[0].item(poly[i][0]) - gyrus)< 0.1)
//        && (fabs(texG[0].item(poly[i][1]) - gyrus)< 0.1)
//        && (fabs(texG[0].item(poly[i][2]) - gyrus)< 0.1) )
    	{
    		polyG.push_back( AimsVector< uint, 3 >(trans[poly[i][0]], trans[poly[i][1]], trans[poly[i][2]] ) );
    	}
    }

    AimsSurfaceTriangle surfaceG;
    surfaceG.vertex()=vertG;
    surfaceG.polygon()=polyG;

   //building new constraints

    cout << "\t Recomputing node indexes for constraints..." << endl;
	std::vector<uint> hautLatG, basLatG, hautLonG, basLonG;

	uint j, ng=gyrusV.size();
	//map<uint, uint> corres; //corres[oldIndex]=newIndex
    for (j=0; j<ng; ++j)
    {
    	i=gyrusV[j];
    	//corres[i]=j;
    }
    std::map<uint, uint>::iterator transIt;

    for (i=0; i<hautLat.size(); ++i) hautLatG.push_back(trans[hautLat[i]]);
    for (i=0; i<hautLon.size(); ++i) hautLonG.push_back(trans[hautLon[i]]);
    for (i=0; i<basLat.size(); ++i) basLatG.push_back(trans[basLat[i]]);
    for (i=0; i<basLon.size(); ++i) basLonG.push_back(trans[basLon[i]]);

    set<uint> bords;
    for (i=0; i<hautLatG.size(); ++i) bords.insert(hautLatG[i]);
	for (i=0; i<hautLonG.size(); ++i) bords.insert(hautLonG[i]);
	for (i=0; i<basLatG.size(); ++i) bords.insert(basLatG[i]);
	for (i=0; i<basLonG.size(); ++i) bords.insert(basLonG[i]);

    cout << "Saving gyrus patch" << flush;
    string meshName=fileOut+".mesh";
    Writer<AimsSurfaceTriangle> triW( fileOut );
    triW << surfaceG;
    cout << "done" << endl;

   //computing the laplacian weight map for the gyrus.
   std::cout << "\t Computing laplacian estimation weights" << std::endl;
   std::map<unsigned, set< pair<unsigned,float> > > weightLapl=AimsMeshWeightFiniteElementLaplacian( surface[0], 0.95 );
   std::map<unsigned, set< pair<unsigned,float> > > weightLaplG;

   std::cout << "\t Translating to gyri..." << std::endl;
   // computing lapl map for gyrus only
   std::map<unsigned, set< pair<unsigned,float> > >::iterator wLapIt;
   for (wLapIt=weightLapl.begin(); wLapIt!=weightLapl.end(); ++wLapIt)
   {
	   i=(*wLapIt).first;
	   if (cortRegion.find(i) != cortRegion.end())
	   {
		   weightLaplG[trans[i]]=set< pair<unsigned,float> >();
		   std::set<pair<unsigned, float> >::iterator pIt=(*wLapIt).second.begin();
		   float sum=0.0;
		   for (; pIt!=(*wLapIt).second.end(); ++pIt)
		   {
			   uint j=(*pIt).first;
			   if (cortRegion.find(j) != cortRegion.end())
			   {
				   sum += (*pIt).second;
			   }
		   }
		   for (pIt=(*wLapIt).second.begin() ; pIt!=(*wLapIt).second.end(); ++pIt)
		   {
			   uint j=(*pIt).first;
			   if (cortRegion.find(j) != cortRegion.end())
			   {
				   if (bords.find(trans[i]) != bords.end())
					   weightLaplG[trans[i]].insert(std::pair<unsigned,float>(trans[j], (*pIt).second/sum ));
				   else
					   weightLaplG[trans[i]].insert(std::pair<unsigned,float>(trans[j], (*pIt).second));
			   }
		   }
	   }
   }
   std::cout << "\t OK" << endl;

//   uint nv2=vertG.size();
//   TimeTexture<int> texD(1, nv2);
//   for (i=0; i<nv2; i++)
//	  texD[0].item(i)=0;
//   std::vector<uint>::iterator vIt;
//   for (vIt=basLonG.begin(); vIt!=basLonG.end(); ++vIt) texD[0].item(*vIt)=100;
//   for (vIt=hautLonG.begin(); vIt!=hautLonG.end(); ++vIt) texD[0].item(*vIt)=200;
//   for (vIt=basLatG.begin(); vIt!=basLatG.end(); ++vIt) texD[0].item(*vIt)=300;
//   for (vIt=hautLatG.begin(); vIt!=hautLatG.end(); ++vIt) texD[0].item(*vIt)=400;
//   Writer<TimeTexture<int> > texDW("debugBorders.tex");
//   texDW << texD;
   //std::cout << "Out" << std::endl;
   //exit(EXIT_FAILURE);

   std::cout << "\tStarting diffusion with delta=" << delta << " and dt=" << dt << std::endl;
   Texture<double> lonDiff, latDiff;
   std::vector< std::pair < std::vector <uint>, short > > constraints;
   std::vector<uint> dummy;
   std::cout << "\t\t Longitude..." << std::endl;
   lonDiff=diffusion(weightLaplG, surfaceG[0], hautLonG, basLonG, constraints, -1, gyrusV, dummy, delta, dt );
   std::cout << "\t\t Latitude..." << std::endl;
   latDiff=diffusion(weightLaplG, surfaceG[0], hautLatG, basLatG, constraints, -1, gyrusV, dummy, delta, dt );

   TimeTexture<float> lonG(1, ng), latG(1, ng);
   for (i=0; i< ng; ++i)
   {
	   lonG[0].item(i)=(float) lonDiff.item(i);
	   latG[0].item(i)=(float) latDiff.item(i);
   }

   std::cout << "Saving textures" << std::endl;
   string texName=fileOut+"_lon.tex";
   Writer<TimeTexture<float> > lonW(texName);
   lonW << lonG;
   texName=fileOut+"_lat.tex";
   Writer<TimeTexture<float> > latW(texName);
   latW << latG;
   cout << "Finished." << endl;

   return EXIT_SUCCESS;

  // DEBUG TEXTURE
//
//  TimeTexture<int> texD(1, nv);
//  for (i=0; i<nv; i++)
//	  texD[0].item(i)=0;
//  std::vector<uint>::iterator vIt;
//  for (vIt=basLon.begin(); vIt!=basLon.end(); ++vIt) texD[0].item(*vIt)=100;
//  for (vIt=hautLon.begin(); vIt!=hautLon.end(); ++vIt) texD[0].item(*vIt)=200;
//  for (vIt=basLat.begin(); vIt!=basLat.end(); ++vIt) texD[0].item(*vIt)=300;
//  for (vIt=hautLat.begin(); vIt!=hautLat.end(); ++vIt) texD[0].item(*vIt)=400;
//
//  Writer<TimeTexture<int> > texDW("debugBorders.tex");
//  texDW << texD;
//
//
//
//  return EXIT_SUCCESS;
}


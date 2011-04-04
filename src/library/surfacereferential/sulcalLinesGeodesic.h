/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *  	CEA/DSV/SHFJ
 *  	4 place du General Leclerc
 *  	91401 Orsay cedex
 *  	France
 */



#ifndef AIMS_SULCUS_CORTICAL_GEO_H
#define AIMS_SULCUS_CORTICAL_GEO_H

#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfacegen.h>
#include <aims/mesh/surfaceOperation.h>

#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshvoronoi.h>
#include <aims/scalespace/meshDiffuse.h>
#include <aims/distancemap/meshmorphomat.h>

#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <aims/geodesicpath/geodesic_algorithm_dijkstra.h>
#include <aims/geodesicpath/geodesic_mesh.h>
#include <aims/geodesicpath/geodesic_mesh_elements.h>

#include <iostream>
#include <fstream>

#include <queue>
#include <map>

using namespace aims;
using namespace aims::meshdistance;
using namespace std;

//namespace aims
//{
	class SulcalLinesGeodesic
	{
		public:
	
      //Parametres pour l'execution
			string _adrMesh;
			string _adrGeodesicDepth;
			string _adrBassinsDepthNorm;
			string _adrCurv;

			string _adrRootsLat;
			string _adrRootsLon;

			string _adrBassins;
			string _adrLines;

			string _adrLatGeodesicOut;
			string _adrLonGeodesicOut;
					
			int _strain;
			
			AimsSurfaceTriangle _mesh;

			TimeTexture<float> _geoDepth;

			std::vector<std::set<uint> > _neigh;

			TimeTexture<short> _texProbaPath;

			TimeTexture<short> _rootsLon;
			TimeTexture<short> _rootsLat;

			TimeTexture<float> _texCurv;
			TimeTexture<float> _texCurvSeuil;

			TimeTexture<short> _texBassins;
			TimeTexture<float> _texBassinsDepthNorm;

			TimeTexture<short> _texBassinsLat;
			TimeTexture<short> _texBassinsLon;

			geodesic::Mesh _meshSPc;
			geodesic::Mesh _meshSP;

			geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;

      vector<double> _pointsSP;
      vector<unsigned> _facesSP;

      set<int> _listIndexVertexFill;
      map<int,set<int> > _mapBassins;

      set<int> _listIndexLon;
      set<int> _listIndexLat;

      map<int, set<int> > _mapConstraintLat;
      map<int, set<int> > _mapConstraintLon;

      ofstream myHistoLat;

			//Constructor
			SulcalLinesGeodesic( string & adrMesh,string & adrCurv, string & adrGeodesicDepth, string & adrBassinsDepthNorm, string & adrRootsLon,
          string & adrRootsLat, std::string & _adrLines, string & _adrBassins, string & adrLonGeodesicOut, string & adrLatGeodesicOut,
			    int strain );
			
			~SulcalLinesGeodesic();

			//public methods
			void run();

		private :
			//private methods
			void computeGraphDijkstra (AimsSurfaceTriangle surface, int constraintType,int strain);
			double computeShortestPathSulci(unsigned source, unsigned target, vector<geodesic::SurfacePoint> & SPath, vector<int> &listIndexVertexPathSP );
			double computeDepthShortestPathSulci(unsigned source, unsigned target, vector<geodesic::SurfacePoint> & SPath, vector<int> &listIndexVertexPathSP );

			double saveHistoTemp(unsigned source, unsigned target);

	    void floodFillIter(int indexVertex, float newTextureValue, float oldTextureValue);
			void bassinsDetect();
			void bassinsDetect2();
			void bassinsDetect3();
			vector<int> maxGeodesicDistance(vector<int> points, int constraint, int* s , int* d );
			vector<int> maxGeodesicDistanceDepthStable(vector<int> points, int constraint, int* s, int *d);

	};


#endif



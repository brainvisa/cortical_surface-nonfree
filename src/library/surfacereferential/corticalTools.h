/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *  	CEA/DSV/SHFJ
 *  	4 place du General Leclerc
 *  	91401 Orsay cedex
 *  	France
 */


#ifndef AIMS_CORTICALREFERENTIAL_TOOLS_H
#define AIMS_CORTICALREFERENTIAL_TOOLS_H


#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshvoronoi.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/scalespace/meshDiffuse.h>

#include <cortical_surface/surfacereferential/corticalTools.h>

namespace aims
{
	//Tools

	//initialization of single textures
	std::pair<int,int> find_poles(TimeTexture<float> &, AimsSurfaceTriangle );
	
	float compute2D_distance(float lon1, float lon2, float lat1, float lat2);
	
	//Poles points determination (to be improved)
	void init_texture_single(TimeTexture<float> &);
	/*
    template<class T>
    Texture<float> MeshDistanceHere( const AimsSurface<3,Void> &  mesh, 
				 const Texture<T> & inittex, 
				 bool allowUnreached );
				 */
	void drawGyri(std::string & adLong, std::string & adLat, std::string & adOut, std::string & adr_cor, std::string & side);
	
	TimeTexture<float> dilate_texture(TimeTexture<float> & texture, float val, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh);
}

#endif

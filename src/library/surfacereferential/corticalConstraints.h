#ifndef AIMS_CORTICALREFERENTIAL_CONSTRAINTS_H
#define AIMS_CORTICALREFERENTIAL_CONSTRAINTS_H

#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
//#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshvoronoi.h>
//#include <aims/distancemap/projection.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/scalespace/meshDiffuse.h>

#include <cortical_surface/surfacereferential/corticalTools.h>

namespace aims
{
	
	//constraint methods
	TimeTexture<float> rescaleConstraints(TimeTexture<float> tex, std::map< int, std::map<int, std::string> > map_global);
	TimeTexture<float> rescaleLatitude(TimeTexture<float> tex);
	TimeTexture<float> rescaleLongitude(TimeTexture<float> tex);
	TimeTexture<float> constraintCleaner(TimeTexture<float> & texture, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, float contr, float curvature, float elasticity);
	std::map< int, std::map<int, std::string> > createCorrespMap( std::string & adr_cor, std::string & adr_file, std::string & side);
	TimeTexture<float> EconstraintCleaner(TimeTexture<float> & texture, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, float valeur);
	void recurs_proc(int index, TimeTexture<float> & tex, TimeTexture<float> & result, std::vector<std::set<uint> > & neigh, float & value);
	int nb_vertex(TimeTexture<float>tex, float value, int size);
	TimeTexture<float> origin_meridian(TimeTexture<float> & tex, int nord, int sud, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, TimeTexture<float> &  poles);
	void meridianLink(TimeTexture<float> & origine, TimeTexture<float> & finish, int flag, int nord, int sud, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, TimeTexture<float> & poles);
	std::vector<unsigned> buildOriginVector(TimeTexture<float> & tex, int nord, int sud, std::vector<std::set<uint> > neigh);
	void findNearNeigh(int origine, int destination, TimeTexture<float> & tex_origine, int flag, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh );
	void findNearNeighPoles(int origine, int destination, TimeTexture<float> & tex_origine, int flag, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh, TimeTexture<float> & dist_p);
	TimeTexture<float> originNeighbourgs(TimeTexture<float> originMeridian, int nord, int sud, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh , TimeTexture<float> & poles);
	TimeTexture<float> defineSides( TimeTexture<float> & sides, TimeTexture<float> & constraints, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh );
	TimeTexture<float> defineSidesPoles( TimeTexture<float> & sides, TimeTexture<float> & constraints, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh );
}

#endif

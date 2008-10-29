#ifndef AIMS_MARCHING_TRIANGLE_H
#define AIMS_MARCHING_TRIANGLE_H


#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

namespace aims
{


class IsoLine
{


public:
	
	AimsSurfaceTriangle mesh;
	TimeTexture<float> texOriginal;
	std::vector< Point3df > vertex;
	std::vector< AimsVector< uint,3> > poly;
	int value;
	float radius1, radius2;
	unsigned facets;
	bool closed, smooth;
	
	IsoLine(AimsSurfaceTriangle & mesh_read, TimeTexture<float> & texOriginal_read):mesh(mesh_read),texOriginal(texOriginal_read)
	{
		vertex = mesh.vertex();
		poly = mesh.polygon();
		radius1 = 0.2;
		radius2 = 0.2;
		facets = 6;
		closed = false;
		smooth = true;
		
	}

	AimsSurfaceTriangle makeTubes(int val);
        AimsSegments makeLine(int val);
	TimeTexture<short> setVertices();
	Point3df createNewVertex(Point3df &, Point3df &, int, int);
        void addSegment(Point3df v1, Point3df v2, AimsSegments *line);
};

}  //fin du namespace

#endif

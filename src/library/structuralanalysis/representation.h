#ifndef SURF_REPRESENTATION_H
#define SURF_REPRESENTATION_H 
#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>


using namespace aims;
using namespace carto;
using namespace std;

pair<Point2df, Point2df> getBoundingBox(set<int> &nodes_list, TimeTexture<float> &lat, TimeTexture<float> &lon);

AimsSurfaceTriangle getFlatMap(vector<vector<int> > &nodes_lists, TimeTexture<float> &lat, TimeTexture<float> &lon, TimeTexture<float> &tex);

#endif


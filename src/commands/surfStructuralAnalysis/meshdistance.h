#ifndef AIMS_RELAXMESHDISTANCE_H
#define AIMS_RELAXMESHDISTANCE_H

#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace std;
using namespace carto;

vector<map<uint,float> > CalculeCarteDistances(AimsSurfaceTriangle mesh, set<uint> nodes);
map<float, vector<pair<float, uint> > > getAlternateMesh(AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &longit);

#endif


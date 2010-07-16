#ifndef AIMS_RELAXMESHDISTANCE_H
#define AIMS_RELAXMESHDISTANCE_H

#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <cortical_surface/structuralanalysis/sites.h>

using namespace aims;
using namespace std;
using namespace carto;

//vector<map<uint,float> > CalculeCarteDistances(AimsSurfaceTriangle mesh, set<uint> nodes,float dist_thresh);
// vector<map<uint,float> > CalculeDistancesBlob(AimsSurfaceTriangle mesh, set<uint> nodes, vector<uint> &sites);
//map<float, vector<pair<float, uint> > > getAlternateMesh(AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &longit);
//map<uint,float> getDistMap( AimsSurfaceTriangle *mesh,  map<unsigned, set<unsigned> >    &neighbours,  int dep, float dist_thresh);


#endif


#ifndef AIMS_SITES_H
#define AIMS_SITES_H

using namespace aims;
using namespace carto;
using namespace std;

class Site{
  public :
    int index;
    int graph_index;
    string subject;
    int label;
    float tValue;
    float t;
    float tmin;
    float tmax;
    float trep;
    Point3df gravitycenter;
    uint node;
};

vector<Site *> ConstruireSites(Graph &primal, map<float, vector<pair<float, uint > > > &altmesh);

#endif


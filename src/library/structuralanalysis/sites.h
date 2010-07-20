#ifndef AIMS_SITES_H
#define AIMS_SITES_H
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>



class Site{
  public :
    uint index;
    int graph_index;
    std::string subject;
    int label;
    float tValue;
    float t;
    float tmin;
    float tmax;
    float trep;
    int label_occur_number;
    float significance;
    float t_rankperc;
    float sim_rankperc;
    int rank;
    Point3df gravitycenter;
    uint node;
    Point3df boundingbox_max;
    Point3df boundingbox_min;
    std::set<int> nodes_list;
};

std::vector<Site *> ConstruireSites(Graph &primal); //map<float, vector<pair<float, uint > > > &altmesh);

#endif


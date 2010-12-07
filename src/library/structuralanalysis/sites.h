#ifndef AIMS_SITES_H
#define AIMS_SITES_H

#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <graph/graph/graph.h>


class Site {
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
        Point3df boundingbox_max;
        Point3df boundingbox_min;
        std::set<int> nodes_list;
        Site() {  label_occur_number = -1; significance = -1.0; t_rankperc = -1.0; sim_rankperc = -1.0; }
};


#endif


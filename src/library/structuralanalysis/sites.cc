#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/minimization.h>

using namespace aims;
using namespace carto;
using namespace std;


std::vector<Site *> BuildSites(Graph &primal) { //map<float, vector<pair<float, uint > > > &altmesh){
  std::vector<Site *> sites;

  std::set<Vertex *>::iterator iv, jv;
  vector<float> bc1, bc2;
  float tmin_1 = -1.0,
        tmax_1 = -1.0,
        trep = -1.0,
        tvalue1 = -1.0,
        t = -1.0,
        x = -1.0,
        y = -1.0;
  int index1 = -1, rank = -1;
  string  subject1 = "",
          subject2 = "";
  int newindex = 0;
  map< float, vector< pair< float, uint > > >::iterator meshIt;
  vector< pair< float, uint > >::iterator yIt;
  vector< int > listepts;
  vector< float > bbmax, bbmin;

  for ( iv = primal.vertices().begin() ; iv != primal.vertices().end() ; ++iv ) {
    string test;
    (*iv)->getProperty( "index", index1);
    (*iv)->getProperty(  "subject", subject1);
    (*iv)->getProperty( "label", test);
    (*iv)->getProperty( "gravity_center", bc1);
    (*iv)->getProperty( "tmin", tmin_1);
    (*iv)->getProperty( "tmax", tmax_1);
    (*iv)->getProperty( "trep", trep);
    (*iv)->getProperty( "t", t);
    (*iv)->getProperty( "tValue", tvalue1);
    (*iv)->getProperty( "rank", rank);
    (*iv)->getProperty( "nodes_list", listepts);
    (*iv)->getProperty( "boundingbox_max", bbmax);
    (*iv)->getProperty( "boundingbox_min", bbmin);

    sites.push_back( new Site() );
    Site *s = sites[sites.size()-1];

    s->label = atoi(test.data());    (*iv)->setProperty("label", test);
    (*iv)->setProperty( "name", test);
    (*iv)->getProperty( "label", test);
    (*iv)->getProperty( "name", test);
    s->index = newindex++;
    s->graph_index = index1;
    s->subject = subject1;

    if ( bc1.size() == 3 ) {
        s->gravitycenter = Point3df();
        s->gravitycenter[0] = bc1[0];
        s->gravitycenter[1] = bc1[1];
        s->gravitycenter[2] = bc1[2];
        x = bc1[0];
        y = bc1[1];
    }
    s->tmin = tmin_1;
    s->tmax = tmax_1;
    s->trep = trep;
    s->rank = rank;
    if ( bbmax.size() == 3 ) {
        s->boundingbox_max[0] = bbmax[0];
        s->boundingbox_max[1] = bbmax[1];
        s->boundingbox_max[2] = bbmax[2];
    }
    if ( bbmin.size() == 3 ) {
        s->boundingbox_min[0] = bbmin[0];
        s->boundingbox_min[1] = bbmin[1];
        s->boundingbox_min[2] = bbmin[2];
    }

    s->t = t;
    s->tValue = tvalue1;
    for ( uint i = 0 ; i < listepts.size() ; i++ )
        s->nodes_list.insert( listepts[i] );


    (*iv)->setProperty( "sites_index", sites.size()-1);

  }


  return sites;
}

//##############################################################################

void SurfaceBased_StructuralAnalysis::ConvertSSBlobsToSites( std::vector<surf::ScaleSpaceBlob *> &ssblobs, std::vector<Site *> &sites ) {

    for ( uint i = 0 ; i < ssblobs.size() ; i++ ) {

        sites.push_back(new Site());
        Site *s = sites[sites.size() - 1];
        s->index = ssblobs[i]->index;

        s->subject = ssblobs[i]->subject;
        s->label = ssblobs[i]->label;
        s->tmin = ssblobs[i]->tmin;
        s->tmax = ssblobs[i]->tmax;
        s->t = ssblobs[i]->t;
        s->nodes_list = ssblobs[i]->nodes;
        set<int> sTemp(sites[i]->nodes_list);
        pair<Point2df, Point2df> bb = ssblobs[i]->get2DBoundingBox();
        s->boundingbox_min = Point3df(bb.first[0], bb.first[1], 0);
        s->boundingbox_max = Point3df(bb.second[0], bb.second[1], 0);

    }

}

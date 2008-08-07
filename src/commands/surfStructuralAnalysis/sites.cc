#include <aims/getopt/getopt2.h>
#include "sites.h"

using namespace aims;
using namespace carto;
using namespace std;


vector<Site *> ConstruireSites(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &longit) { //map<float, vector<pair<float, uint > > > &altmesh){
  vector<Site *> sites;

  std::set<Vertex *>::iterator iv, jv;
  vector<float> bc1, bc2;
  float tmin_1, tmax_1, trep, tvalue1,x,y,precisionX=10.0, precisionY=10.0;
  int index1;
  string subject1, subject2;
  int newindex=0;
  float minidis=10000.0,dis; uint mini=0;
  map<float, vector<pair<float, uint > > >::iterator meshIt;
  vector<pair<float, uint> >::iterator yIt;
  for (iv=primal.vertices().begin() ; iv!=primal.vertices().end(); ++iv){
    string test;
    (*iv)->getProperty("index", index1);
    (*iv)->getProperty( "subject", subject1);
    (*iv)->getProperty("label",test);
    (*iv)->getProperty( "subject", subject1 );
    (*iv)->getProperty( "gravity_center", bc1);
    (*iv)->getProperty( "tmin", tmin_1);
    (*iv)->getProperty( "tmax", tmax_1);
    (*iv)->getProperty( "trep", trep);
    (*iv)->getProperty( "tValue", tvalue1);
    sites.push_back(new Site());
    Site *s=sites[sites.size()-1];
    
    s->label = atoi(test.data());    (*iv)->setProperty("label", test);
    (*iv)->setProperty("name", test);
    (*iv)->getProperty("label", test);
    (*iv)->getProperty("name", test);
    s->index = newindex++;
    s->graph_index = index1;
    s->subject = subject1;
    s->gravitycenter = Point3df();
    s->gravitycenter[0] = bc1[0];
    s->gravitycenter[1] = bc1[1];
    s->gravitycenter[2] = bc1[2];
    s->tmin = tmin_1;
    s->tmax = tmax_1;
    s->trep = trep;
    s->t = tvalue1;
    s->tValue = tvalue1;
    x=bc1[0];
    y=bc1[1];
    minidis=10000.0;
    for (uint i=0;i<mesh[0].vertex().size();i++){
      dis = sqrt(pow(lat[0].item(i)-x,2)+pow(longit[0].item(i)-y,2));
      if ( dis< minidis){
        mini=i;
        minidis=dis;
      }
    }
    s->node = mini;
    (*iv)->setProperty( "sites_index", sites.size()-1);
  }

  return sites;
}

#include <aims/getopt/getopt2.h>
#include "sites.h"

using namespace aims;
using namespace carto;
using namespace std;


vector<Site *> ConstruireSites(Graph &primal) { //map<float, vector<pair<float, uint > > > &altmesh){
  vector<Site *> sites;

  std::set<Vertex *>::iterator iv, jv;
  vector<float> bc1, bc2;
  float tmin_1, tmax_1, trep, tvalue1,t1,x,y;
  int index1,rank;
  string subject1, subject2;
  int newindex=0;
  //float minidis=10000.0,dis; uint mini=0;
  map<float, vector<pair<float, uint > > >::iterator meshIt;
  vector<pair<float, uint> >::iterator yIt;
  vector<int> listepts;
  
 
  
  
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
    (*iv)->getProperty( "t", t1);
    (*iv)->getProperty( "tValue", tvalue1);
    (*iv)->getProperty( "rank", rank);
    (*iv)->getProperty( "nodes_list", listepts);

//     if (bc1[0] >91.0 && bc1[0]<150.0 && bc1[1]>257.0 && bc1[1]<316.0){
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
    s->rank = rank;
    s->t = t1;
    s->tValue = tvalue1;
    for (uint i=0;i<listepts.size();i++)
      s->nodes_list.insert(listepts[i]);
    x=bc1[0];
    y=bc1[1];
    
    
//     (*iv)->setProperty( "node", (int)mini);
    (*iv)->setProperty( "sites_index", sites.size()-1);
//     }
  }
  vector<uint> histoblobs(40),histoblobs_tvalue(40);
  for (uint i=0;i<40;i++){
    histoblobs[i] = 0;
    histoblobs_tvalue[i] = 0;
  }
  for (uint i=0;i<sites.size();i++){
    if ((sites[i]->t/1.0+5.0) < 0) histoblobs[0]++;
    else if ((sites[i]->t+5.0)/1.0 > 39) histoblobs[39]++;
    else{
    uint idx = (uint)((sites[i]->t+5.0) / 1.0);
    histoblobs[idx]++;
    }
  } 
  for (uint i=0;i<sites.size();i++){
    if ((sites[i]->tValue/1.0+5.0) < 0) histoblobs_tvalue[0]++;
    else if ((sites[i]->tValue+5.0)/1.0 > 39) histoblobs_tvalue[39]++;
    else{
      uint idx = (uint)((sites[i]->tValue+5.0) / 1.0);
      histoblobs_tvalue[idx]++;
    }
  }
  cout << "Histogramme" << endl;
  for (uint i=0;i<histoblobs.size();i++)
    cout << histoblobs[i] << " " ;
  cout << endl;
  for (uint i=0;i<histoblobs_tvalue.size();i++)
    cout << histoblobs_tvalue[i] << " " ;
  cout << endl;
  return sites;
}

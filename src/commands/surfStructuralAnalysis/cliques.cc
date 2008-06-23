#include <aims/getopt/getopt2.h>
#include "cliques.h"

using namespace aims;
using namespace carto;
using namespace std;

float Clique::ddweight, Clique::intrapsweight, Clique::simweight, Clique::lsweight, Clique::ddx2, Clique::ddx1, Clique::ddh;
void Clique::setParameters(float _ddweight, float _intrapsweight, float _simweight, float _lsweight, float _ddx2, float _ddx1, float _ddh){
  ddweight=_ddweight; intrapsweight = _intrapsweight; simweight=_simweight; lsweight=_lsweight; ddx2 =_ddx2;  ddx1 = _ddx1; ddh=_ddh;
}

void Clique::updateLabelsCount(){
  if (type == INTRAPRIMALSKETCH){
    labelscount = map<int, uint>();
    for (uint i=0;i<blobs.size();i++){
      if (labelscount.find(blobs[i]->label) == labelscount.end())
        labelscount[blobs[i]->label] = 1;
      else
        labelscount[blobs[i]->label]++;
    }
    map<int,uint>::iterator it;
    uint chksum=0;
    for (it=labelscount.begin();it!=labelscount.end();it++)
      chksum += (*it).second;
    ASSERT(chksum == blobs.size());
  }
}

vector<Clique> ConstruireCliques(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, vector<map<uint,float> > &distmap){
  
  cliquesDuSite = vector<vector<int> >(sites.size());
  vector<Clique> cliques;
  vector<set<uint> > voisins(sites.size());
  set<uint> v;
  vector<string> subjects;
  double rec;
  map<uint,float>::iterator distit;
  vector<Clique> intraps;
  for (uint i=0;i<sites.size(); i++){
    uint j=0;
    for (;j<subjects.size() && subjects[j] != sites[i]->subject;j++){}
    if (j==subjects.size()){
      subjects.push_back(sites[i]->subject);
      intraps.push_back(Clique());
      intraps[j].type = INTRAPRIMALSKETCH;
    }
    intraps[j].blobs.push_back(sites[i]);
    cliquesDuSite[sites[i]->index].push_back(j);
  }
  for (uint i=0;i<intraps.size();i++)
    cliques.push_back(intraps[i]);

  for (uint i=0;i<sites.size(); i++){
    Clique c, ls;
    c.type = DATADRIVEN; ls.type = BESTLOWERSCALE;
    cliquesDuSite[sites[i]->index].push_back(cliques.size());
    c.blobs.push_back(sites[i]);
    cliques.push_back(c);
    cliquesDuSite[sites[i]->index].push_back(cliques.size());
    ls.blobs.push_back(sites[i]);
    cliques.push_back(ls);
    for (uint j=i;j<sites.size();j++) {
      if (sites[i]->subject != sites[j]->subject) {
        distit = distmap[sites[i]->node].find(sites[j]->node);
        
        if (distit != distmap[sites[i]->node].end())
          rec = (*distit).second;
        else
          rec = 999.0;
        if ((rec < 20.0) && !((sites[j]->tmin > sites[i]->tmax) || (sites[i]->tmin > sites[j]->tmax))/* && (sites[i]->tValue * sites[j]->tValue > 40.0)*/) {
          v.insert(i);
          v.insert(j);
          voisins[i].insert(j);
          voisins[j].insert(i);
          Clique simc;
          simc.type = SIMILARITY;
          simc.rec = rec;
          cliquesDuSite[sites[i]->index].push_back(cliques.size());
          cliquesDuSite[sites[j]->index].push_back(cliques.size());
          simc.blobs.push_back(sites[i]);
          simc.blobs.push_back(sites[j]);
          cliques.push_back(simc);
        }
      }
    }
  }
  return cliques;
}

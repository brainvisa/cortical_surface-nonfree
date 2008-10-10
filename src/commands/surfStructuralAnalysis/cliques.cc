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

struct ltstr_p2d
{
  bool operator()(const AimsVector<uint,2> p1, const AimsVector<uint,2> p2)
  {
      if (p1[1]<p2[1])
        return true;
      else if (p1[1]>p2[1])
        return false;
      else
      {
        if (p1[0]<p2[0])
          return true;
        else return false;
      }
    }
};

int getNode(map<float, vector<pair<float, uint> > > &altmesh, float delta, float x, float y, TimeTexture<float> &tex1, TimeTexture<float> &tex2){
  map<float, vector<pair<float, uint> > >::iterator xit;
  float dismin=1000.0, dis;
  int mini=-1;
  for (xit=altmesh.begin();xit!=altmesh.end();xit++){
    if (xit->first>x-delta && xit->first<x+delta){
      for (uint j=0;j<xit->second.size();j++){
        if (xit->second[j].first>y-2.0 && xit->second[j].first<y+2.0){
          dis = sqrt(pow(tex1[0].item(xit->second[j].second)-x,2)+pow(tex2[0].item(xit->second[j].second)-y,2));
          if (dis<dismin){
            dismin=dis;
            mini=xit->second[j].second;
          }
        }
      }
    }
  }
  return mini;
}

vector<Clique> ConstruireCliques(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){
  uint temp=0,temp2=0,temp3=0,temp4=0;
  cliquesDuSite = vector<vector<int> >(sites.size());
  vector<Clique> cliques;
  set<uint> v;
  vector<string> subjects;
  double rec;
  cout << endl << "Construction carte de distances ..." << endl;
  set<uint> setnodes;
  float x,y;
  int mini=0,dep,arr;
  map<string,vector<map<string,uint> > > matching;
  set<uint> megasitesnodes;
  set<uint>::iterator megit;
  map<string, map<uint, map<uint,float> > > distmaptable;

  typedef map<uint, set<uint> > mapvoisins;
  map<string, mapvoisins> neighbours;
  unsigned v1b, v2, v3;

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
  for (uint s=0;s<subjects.size();s++){
    for( uint i=0; i<meshes[subjects[s]].polygon().size(); ++i ){
      v1b = meshes[subjects[s]].polygon()[i][0];
      v2 = meshes[subjects[s]].polygon()[i][1];
      v3 = meshes[subjects[s]].polygon()[i][2];

      neighbours[subjects[s]][v1b].insert( v2 );
      neighbours[subjects[s]][v2].insert( v1b );

      neighbours[subjects[s]][v1b].insert( v3 );
      neighbours[subjects[s]][v3].insert( v1b );

      neighbours[subjects[s]][v2].insert( v3 );
      neighbours[subjects[s]][v3].insert( v2 );
    }
  }


  
  for (uint i=0;i<sites.size();i++)
    for (megit=sites[i]->nodes_list.begin();megit!=sites[i]->nodes_list.end();megit++)
      megasitesnodes.insert(*megit);
    
  cout << megasitesnodes.size() << endl;

  map<string, AimsSurfaceTriangle>::iterator meshit, meshjt;
  map<string, map<float, vector<pair<float, uint> > > > altmeshes;
  for (meshit=meshes.begin();meshit!=meshes.end();meshit++){
    matching[meshit->first]=vector<map<string,uint> > (meshit->second.vertex().size());
    altmeshes[meshit->first]=getAlternateMesh(meshit->second, lats[meshit->first], lons[meshit->first]);
  }
  cout << "Construction maillages alternatifs terminée ..." << endl;
  uint cpt = megasitesnodes.size();
  megasitesnodes.clear();
  cout << "Recherche des correspondances noeud-à-noeud à travers tous les sujets..." << endl;
  for (uint i=0;i<sites.size();i++){
    for (megit=sites[i]->nodes_list.begin();megit!=sites[i]->nodes_list.end();megit++){
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << megasitesnodes.size() << "/" << cpt << flush;
      megasitesnodes.insert(*megit);
      TimeTexture<float> tex1,tex2;
      tex1=lats[sites[i]->subject]; tex2=lons[sites[i]->subject];
      for (meshjt=meshes.begin();meshjt!=meshes.end();meshjt++){

        if (matching[sites[i]->subject][*megit].find(meshjt->first)==matching[sites[i]->subject][*megit].end()){
          string atlascurr = (*meshjt).first;
          x=tex1[0].item(*megit); y=tex2[0].item(*megit);

          float delta=0.5;
          do{
            mini = getNode(altmeshes[atlascurr],delta,x,y,lats[atlascurr],lons[atlascurr]);
            delta*=2.0;
          }
          while(mini==-1);

          matching[sites[i]->subject][*megit][meshjt->first]=mini;
          matching[meshjt->first][mini][sites[i]->subject] = *megit;
        }
      }
    }
  }
  cout << "done" << endl;

  
  for (uint i=0;i<sites.size();i++){
    sites[i]->node = matching[sites[i]->subject][i][sites[0]->subject];
    setnodes.insert(sites[i]->node);
  }
//   cout << meshes[sites[0]->subject].vertex().size() << endl;
  vector<map<uint,float> > protodistmap(CalculeCarteDistances(meshes[sites[0]->subject], setnodes));
//   cout << protodistmap.size() << " " << sites[0]->node << endl;
//   cout << protodistmap[sites[0]->node].size() << endl;

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
  }
  for (uint i=0;i<sites.size();i++){
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << sites.size() << "(" << sites[i]->rank << "-" << sites[i]->subject << "-" << sites[i]->nodes_list.size() << ") "<< distmaptable.size() << flush;
    for (uint j=i;j<sites.size();j++) {
      if (sites[i]->subject != sites[j]->subject) {
        map<uint,float> dmap(protodistmap[sites[i]->node]);
        map<uint,float>::iterator distit = dmap.find(sites[j]->node);
        if (distit != dmap.end())
          rec = distit->second;
        else
          rec = 999.0;
  //         rec = sqrt(pow(mesh[0].vertex()[sites[i]->node][0]-mesh[0].vertex()[sites[j]->node][0],2)+pow(mesh[0].vertex()[sites[i]->node][1]-mesh[0].vertex()[sites[j]->node][1],2)+pow(mesh[0].vertex()[sites[i]->node][2]-mesh[0].vertex()[sites[j]->node][2],2))  // USING EUCLIDEAN DISTANCE
        if ((rec < 20.0) && !((sites[j]->tmin > sites[i]->tmax) || (sites[i]->tmin > sites[j]->tmax)) /*&& (sites[i]->tValue * sites[j]->tValue > 2.0)*/) {
          set<uint> pointsInter;
          set_intersection(sites[i]->nodes_list.begin(), sites[i]->nodes_list.end(),
                           sites[j]->nodes_list.begin(), sites[j]->nodes_list.end(),
                                               std::inserter(pointsInter, pointsInter.begin()));
          if (pointsInter.size()>0){
//             cout << "intersection" << endl;
            rec = 2*pointsInter.size()/(sites[i]->nodes_list.size()+sites[j]->nodes_list.size());
          }
          else{
//             cout << i << " " << j << ":" << endl;
            map<string, AimsSurfaceTriangle >::iterator meshit;
            set<uint> nodes_list(sites[i]->nodes_list);
            set<uint> other_list(sites[j]->nodes_list);
            set<uint>::iterator nodesit=nodes_list.begin(),otherit=other_list.begin();
            
//             map<AimsVector<uint,2>,AimsVector<float,2>,ltstr_p2d > results;
//             uint petitmini=*nodesit, grandmini=*otherit;
            float distanceMinInterblobs=10000000000.0;
            
            for (nodesit=nodes_list.begin();nodesit!=nodes_list.end();nodesit++){
              for (otherit=other_list.begin();otherit!=other_list.end();otherit++){
                
                float distance=0, ecart_type=0; uint howmany=0; vector<float> distances;
                cout << "("<< sites[i]->subject <<")\t"<< *nodesit << "\t" << "("<< sites[j]->subject <<")\t"<< *otherit << ":\t" ;
                
                for (meshit=meshes.begin();meshit!=meshes.end();meshit++){
         
                  string atlascurr = (*meshit).first;
                  dep=matching[sites[i]->subject][*nodesit][atlascurr]; //mini;
                  map<uint,float> distmap ;
//                   if (distmaptable[atlascurr].find(dep)==distmaptable[atlascurr].end())
//                     distmaptable[atlascurr][dep] = getDistMap( &meshes[atlascurr], neighbours[atlascurr],  dep);
                  
                  distmap = getDistMap( &meshes[atlascurr], neighbours[atlascurr],  dep); //distmaptable[atlascurr][dep]; //

                  arr=matching[sites[j]->subject][*otherit][atlascurr]; //mini;
                  if (distmap.find(arr)!=distmap.end()){
                    cout << "("<< atlascurr <<")"<<distmap[arr] << "-" ;
                    distance += distmap[arr];
//                   distance  +=10.0;
//                   distances.push_back(10.0);
                    distances.push_back(distmap[arr]);
                    howmany++;
                  }
                }
                
                distance /= howmany;
                ecart_type = 0.0;
                for (uint m=0;m<distances.size();m++)
                  ecart_type += (distances[m] - distance) *(distances[m] - distance) ;
                ecart_type /= distances.size();
                ecart_type = sqrt(ecart_type);
                cout << "\t=" << distance << "\t" << ecart_type << endl;
                if (distance < distanceMinInterblobs){
                  distanceMinInterblobs = distance;
//                   petitmini = *nodesit;
//                   grandmini = *otherit;
                }
              }
            }
            rec = distanceMinInterblobs;
            
  //             vector<map<uint,float> > distmap = CalculeCarteDistances(mesh, sites[i]->nodes_list);
  //             set<uint>::iterator n1it=sites[i]->nodes_list.begin(),n2it=sites[j]->nodes_list.begin();
  //             for (;n1it!=sites[i]->nodes_list.end();n1it++)
  //               for (;n2it!=sites[j]->nodes_list.end();n2it++){
  //               Point3df p1(meshes[sites[i]->subject].vertex()[*n1it]), p2(meshes[sites[j]->subject].vertex()[*n2it]);
  //                 distcurr=sqrt(pow(p1[0]-p2[0],2)+pow(p1[1]-p2[1],2)+pow(p1[2]-p2[2],2));
  //                 if (distmin<distcurr) {
  //                   distcurr=distmin;
  //                   n1min = *n1it;
  //                   n2min = *n2it;
  //                 }
  //               }
              
          }
          
          temp2++;
          Clique simc;
          simc.type = SIMILARITY;
          simc.rec = rec;
  //           cout << rec << " " ;
          cliquesDuSite[sites[i]->index].push_back(cliques.size());
          cliquesDuSite[sites[j]->index].push_back(cliques.size());
          simc.blobs.push_back(sites[i]);
          simc.blobs.push_back(sites[j]);
          cliques.push_back(simc);
  //           }
        }
      }
    }
  }
  
  cout << "TEMP :" << temp<< " " << temp2 << " " << temp3 << " " << temp4 << " "<<temp+temp4<< " ";
//   cin >> temp;

  return cliques;
}

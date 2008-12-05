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
//     cout << blobs[0]->subject << endl;
    for (it=labelscount.begin();it!=labelscount.end();it++){
//       cout << it->first << " " << it->second << "-";
      chksum += (*it).second;
    }
//     cout  << endl;
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
//   cout << ":" << delta << " " << x << " " << y;
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

// map<uint,float> test_geodmap(til::Mesh_N mesh, int dep, double distthresh){
//   std::vector< std::size_t > startPoints;
//   startPoints.push_back(dep);
//   til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
//   typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
//     
//   shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
//     
//   til::Triangle_mesh_geodesic_map<til::Mesh_N::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
//       geomap(getVertices(mesh), *pneighc, stopGhost);
//     
//     
// //     std::vector<std::size_t> startPoints(1);
// //     std::vector<double> dist(1, 0.0);
//     
//   std::vector<double> dist(( int )startPoints.size(), 0.0);
// //     startPoints[0] = vertex;
//   geomap.init(startPoints, dist);
//   geomap.process();
//   shared_ptr<til::sparse_vector<double> > sres = geomap.distanceMap(); 
//   {
//     til::sparse_vector<double>::sparse_iterator iRes = sres->sparse_begin();
//     if (false)
//     {
//       for (; iRes != sres->sparse_end(); ++iRes)
//       {
//         std::cout << iRes->first << " ";
//         std::cout << iRes->second << " ";
//       }
//       std::cout << std::endl;
//     }
//   }
//   map<uint,float> res; //(til::size(*sres));
//   {
// //     std::vector<double>::iterator iRes = res.begin();
//     til::sparse_vector<double>::const_iterator iRes2 = sres->begin();
//     uint i=0;
//     for (; iRes2 != sres->end();++iRes2,i++)
//     {
//       res[i] = *iRes2;
//     }
//   }
// //   for (std::size_t i = 0; i < til::size(res); ++i) res[i] = til::min(res[i], 100.0);
// //   
// //   {
// //     std::vector<double> res2 = res;
// //     std::sort  (res2.begin(), res2.end());
// //     std::vector<double>::iterator newend = std::unique(res2.begin(), res2.end());
// //     if (false) std::cout << "Unique : " << std::distance(res2.begin(), newend) << std::endl;
// // //     for (int i = 0; i < 10; ++i) std::cout << res[i] << " ";
// // //     std::cout << std::endl;
// //   }
//     
// 
//   return res;
// }

vector<Clique> ConstruireCliques(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){
  uint temp=0,temp2=0,temp3=0,temp4=0;
  cliquesDuSite = vector<vector<int> >(sites.size());
  vector<Clique> cliques;
  set<uint> v;
  vector<string> subjects;
  double rec;
  cout << endl << "Construction carte de distances ..." << endl;
  float dist_thresh = 50.0;
  cout << "Distance maps threshold : << " << dist_thresh <<" >>"<< endl;

  set<uint> setnodes;
  float x,y;
  int mini=0,dep,arr;
  map<string,vector<map<string,uint> > > matching;
  set<uint> megasitesnodes;
  set<uint>::iterator megit;
  map<string, map<uint, map<uint,float> > > distmaptable;
  map<string, AimsSurfaceTriangle>::iterator meshit,meshjt;
//   typedef til::Mesh_N MyMesh;
//   map<string,MyMesh> mymeshes;
//   
//   std::cout << "Adding neighbors" << std::endl;    
//   for (meshit = meshes.begin();meshit!=meshes.end();meshit++){
//     til::Mesh1 mesh0;
//     til::convert(mesh0, meshit->second);
//     mymeshes[meshit->first] = addNeighborsToMesh(mesh0);
//   }
  
  
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
    for (megit=sites[i]->nodes_list.begin();megit!=sites[i]->nodes_list.end();megit++){
//       cout << *megit << " ";
      megasitesnodes.insert(*megit);
    }
    
  cout << "megasize:" << megasitesnodes.size() << endl;

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
  cout << "depart " << setnodes.size() << endl;
//   vector<vector<double> > protodistmap(setnodes.size());
  
//   for (megit=setnodes.begin();megit!=setnodes.end();megit++){
//     cout << *megit << endl;
//     vector<double> test= test_geodmap(mymeshes[sites[0]->subject],*megit, 200.0);
//     protodistmap[*megit]=vector<double>(test);
//     
//   }
  
  vector<map<uint,float> > protodistmap(CalculeCarteDistances(meshes[sites[0]->subject], setnodes, dist_thresh));
  cout << "arrivée" << endl;
  
//   vector<double> test= test_geodmap(mymeshes[sites[0]->subject],*(setnodes.begin()), 10.0);
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
  
  float cliques_thresh=40.0;
  cout << "Similarity cliques distance limit : << " << dist_thresh << " >>" << endl;
  
  for (uint i=0;i<sites.size();i++){
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << sites.size() << "(" << sites[i]->rank << "-" << sites[i]->subject << "-" << sites[i]->nodes_list.size() << ") "<< distmaptable.size() << flush;
    for (uint j=i;j<sites.size();j++) {
      if (sites[i]->subject != sites[j]->subject) {
//         map<uint,float> dmap(protodistmap[sites[i]->node]);
//         map<uint,float>::iterator distit = dmap.find(sites[j]->node);
//         if (distit != dmap.end())
//           rec = distit->second;
//         else
//           rec = 999.0;
        rec = protodistmap[sites[i]->node][sites[j]->node];
  //         rec = sqrt(pow(mesh[0].vertex()[sites[i]->node][0]-mesh[0].vertex()[sites[j]->node][0],2)+pow(mesh[0].vertex()[sites[i]->node][1]-mesh[0].vertex()[sites[j]->node][1],2)+pow(mesh[0].vertex()[sites[i]->node][2]-mesh[0].vertex()[sites[j]->node][2],2))  // USING EUCLIDEAN DISTANCE
        if ((rec < cliques_thresh) && !((sites[j]->tmin > sites[i]->tmax) || (sites[i]->tmin > sites[j]->tmax)) /*&& (sites[i]->tValue * sites[j]->tValue > 2.0)*/) {
          set<uint> pointsInter;
          
          
          
//           (*jv)->getProperty( "index", index2 );
//           (*jv)->getProperty( "boundingbox_max", bbmax_2);
//           (*jv)->getProperty( "boundingbox_min", bbmin_2);
//           no_overlap=0;
//           if (bbmin_1[0]<=bbmin_2[0])
//             if (bbmax_1[0]<bbmin_2[0]) no_overlap=1;
//           else overlap_x= (bbmax_2[0] < bbmax_1[0] ? bbmax_2[0] : bbmax_1[0]) - bbmin_2[0] +1;
//           else
//             if (bbmax_2[0]<bbmin_1[0]) no_overlap=1;
//           else overlap_x= (bbmax_1[0] < bbmax_2[0] ? bbmax_1[0] : bbmax_2[0]) - bbmin_1[0] +1;
//           if (no_overlap==0)
//           {
//             if (bbmin_1[1]<=bbmin_2[1])
//               if (bbmax_1[1]<bbmin_2[1]) no_overlap=1;
//             else overlap_y= (bbmax_2[1] < bbmax_1[1] ? bbmax_2[1] : bbmax_1[1]) - bbmin_2[1] +1;
//             else
//               if (bbmax_2[1]<bbmin_1[1]) no_overlap=1;
//             else overlap_y= (bbmax_1[1] < bbmax_2[1] ? bbmax_1[1] : bbmax_2[1]) - bbmin_1[1] +1;
//             if (no_overlap==0)
//             {
//               if (bbmin_1[2]<=bbmin_2[2])
//                 if (bbmax_1[2]<bbmin_2[2]) no_overlap=1;
//               else overlap_z= (bbmax_2[2] < bbmax_1[2] ? bbmax_2[2] : bbmax_1[2]) - bbmin_2[2] +1;
//               else
//                 if (bbmax_2[2]<bbmin_1[2]) no_overlap=1;
//               else overlap_z= (bbmax_1[2] < bbmax_2[2] ? bbmax_1[2] : bbmax_2[2]) - bbmin_1[2] +1;
//               if (no_overlap==0)
//               {
//                 rec=overlap_x*overlap_y*overlap_z;
//                 double div=( ((bbmax_1[0]-bbmin_1[0])*(bbmax_1[1]-bbmin_1[1])*(bbmax_1[2]-bbmin_1[2]) +1)
//                       + ((bbmax_2[0]-bbmin_2[0])*(bbmax_2[1]-bbmin_2[1])*(bbmax_2[2]-bbmin_2[2]) +1) );
// 
//                 rec=2 * rec / div;
//                 (overlap->value)[std::pair<Vertex *,Vertex * >((*iv),(*jv))]=rec;
//               }
//             }
//           }
          
          
          
          
          
          
          
          
          
          
          
          
          
          
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
                  if (distmaptable[atlascurr].find(dep)==distmaptable[atlascurr].end())
                    distmaptable[atlascurr][dep] = getDistMap( &meshes[atlascurr], neighbours[atlascurr],  dep, dist_thresh);
                  
                  distmap = /*test_geodmap( mymeshes[atlascurr],   dep,100.0); //*/distmaptable[atlascurr][dep]; //

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
// 
// vector<Clique> ConstruireCliquesSimple(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons)  {
//   uint temp=0,temp2=0,temp3=0,temp4=0;
//   cliquesDuSite = vector<vector<int> >(sites.size());
//   vector<Clique> cliques;
//   set<uint> v;
//   vector<string> subjects;
//   double rec;
//   set<uint> setnodes;
//   float x,y;
//   int mini=0,dep,arr;
//   map<string,vector<map<string,uint> > > matching;
//   set<uint> megasitesnodes;
//   set<uint>::iterator megit;
//   map<string, AimsSurfaceTriangle>::iterator meshit,meshjt;
//   typedef til::Mesh_N MyMesh;
//   map<string,MyMesh> mymeshes;
//   
//   for (meshit = meshes.begin();meshit!=meshes.end();meshit++){
//     til::Mesh1 mesh0;
//     til::convert(mesh0, meshit->second);
//     mymeshes[meshit->first] = addNeighborsToMesh(mesh0);
//   }
//   
//   vector<Clique> intraps;
//   for (uint i=0;i<sites.size(); i++){
//     uint j=0;
//     for (;j<subjects.size() && subjects[j] != sites[i]->subject;j++){}
//     if (j==subjects.size()){
//       subjects.push_back(sites[i]->subject);
//       intraps.push_back(Clique());
//       intraps[j].type = INTRAPRIMALSKETCH;
//     }
//     intraps[j].blobs.push_back(sites[i]);
//     cliquesDuSite[sites[i]->index].push_back(j);
//   }
//   
//   for (uint i=0;i<sites.size();i++){
// //     cout << sites[i]->nodes_list.size() << endl;
//     for (megit=sites[i]->nodes_list.begin();megit!=sites[i]->nodes_list.end();megit++){
//     cout << *megit << " ";
//     megasitesnodes.insert(*megit);
//     }
//   }
//     
//     cout << "megasize:" << megasitesnodes.size() << " " << sites.size() << endl;
// 
//     map<string, map<float, vector<pair<float, uint> > > > altmeshes;
//     for (meshit=meshes.begin();meshit!=meshes.end();meshit++){
//       matching[meshit->first]=vector<map<string,uint> > (meshit->second.vertex().size());
//       altmeshes[meshit->first]=getAlternateMesh(meshit->second, lats[meshit->first], lons[meshit->first]);
//     }
//     cout << "Construction maillages alternatifs terminée ..." << endl;
//     
// //     uint cpt = 0;
// //     megasitesnodes.clear();
//     cout << "Recherche des correspondances noeud-à-noeud à travers tous les sujets..." << endl;
//     for (uint i=0;i<sites.size();i++){
// //       for (megit=sites[i]->nodes_list.begin();megit!=sites[i]->nodes_list.end();megit++){
//         set<uint> nodes_test;
//         nodes_test.insert(sites[i]->gravitycenter[2]);
//         megit = nodes_test.begin();
//         cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << sites.size() << flush;
// //         megasitesnodes.insert(*megit);
//         TimeTexture<float> tex1,tex2;
//         tex1=lats[sites[i]->subject]; tex2=lons[sites[i]->subject];
// //         for (meshjt=meshes.begin();meshjt!=meshes.end();meshjt++){
// //         meshjt = meshes.begin();
//           if (matching[sites[i]->subject][*megit].find(sites[0]->subject)==matching[sites[i]->subject][*megit].end()){
//             string atlascurr = sites[0]->subject;
//             x=tex1[0].item(*megit); y=tex2[0].item(*megit);
// 
//             float delta=0.5;
//             do{
//               mini = getNode(altmeshes[atlascurr],delta,x,y,lats[atlascurr],lons[atlascurr]);
//               delta*=2.0;
//             }
//             while(mini==-1);
// 
//             matching[sites[i]->subject][*megit][sites[0]->subject]=mini;
//             setnodes.insert(mini);
//             matching[sites[0]->subject][mini][sites[i]->subject] = *megit;
//           }
// //         }
// //       }
//     }
//     cout << "done" << endl;
// 
//   
//     for (uint i=0;i<sites.size();i++){
//       sites[i]->node = matching[sites[i]->subject][i][sites[0]->subject];
// //       setnodes.insert(sites[i]->node);
//     }
// //   cout << meshes[sites[0]->subject].vertex().size() << endl;
//     cout << "depart " << setnodes.size() << endl;
//     vector<vector<double> > protodistmap(setnodes.size());
//   
//     uint cpt2=0;
//     TimeTexture<float> outtest(setnodes.size(), meshes[sites[0]->subject].vertex().size());
//     for (megit=setnodes.begin();megit!=setnodes.end();megit++){
//       cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << cpt2++ << "/" << setnodes.size() << flush;
//       vector<double> test= test_geodmap(mymeshes[(meshes.begin())->first],*megit, 20.0);
//       protodistmap[*megit]=vector<double>(test);
//     
//     }
//     
//   
// //   vector<map<uint,float> > protodistmap(CalculeCarteDistances(meshes[sites[0]->subject], setnodes));
//     cout << "arrivée" << endl;
//   
// //     Texture1d test= test_geodmap(mymeshes[sites[0]->subject],*(setnodes.begin()), 10.0);
// //   cout << protodistmap.size() << " " << sites[0]->node << endl;
// //   cout << protodistmap[sites[0]->node].size() << endl;
// 
//     for (uint i=0;i<intraps.size();i++)
//       cliques.push_back(intraps[i]);
// 
//     for (uint i=0;i<sites.size(); i++){
//       Clique c, ls;
//       c.type = DATADRIVEN; ls.type = BESTLOWERSCALE;
//       cliquesDuSite[sites[i]->index].push_back(cliques.size());
//       c.blobs.push_back(sites[i]);
//       cliques.push_back(c);
//       cliquesDuSite[sites[i]->index].push_back(cliques.size());
//       ls.blobs.push_back(sites[i]);
//       cliques.push_back(ls);
//     }
//     for (uint i=0;i<sites.size();i++){
//       cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << sites.size() << "(" << sites[i]->rank << "-" << sites[i]->subject << "-" << sites[i]->nodes_list.size() << ") " << flush;
//       for (uint j=i;j<sites.size();j++) {
//         if (sites[i]->subject != sites[j]->subject) {
// //         map<uint,float> dmap(protodistmap[sites[i]->node]);
// //         map<uint,float>::iterator distit = dmap.find(sites[j]->node);
// //         if (distit != dmap.end())
// //           rec = distit->second;
// //         else
// //           rec = 999.0;
//           rec = protodistmap[sites[i]->node][sites[j]->node];
//   //         rec = sqrt(pow(mesh[0].vertex()[sites[i]->node][0]-mesh[0].vertex()[sites[j]->node][0],2)+pow(mesh[0].vertex()[sites[i]->node][1]-mesh[0].vertex()[sites[j]->node][1],2)+pow(mesh[0].vertex()[sites[i]->node][2]-mesh[0].vertex()[sites[j]->node][2],2))  // USING EUCLIDEAN DISTANCE
//           if ((rec < 20.0) && !((sites[j]->tmin > sites[i]->tmax) || (sites[i]->tmin > sites[j]->tmax)) /*&& (sites[i]->tValue * sites[j]->tValue > 2.0)*/) {
//             
//             Clique simc;
//             simc.type = SIMILARITY;
//             simc.rec = rec;
//   //           cout << rec << " " ;
//             cliquesDuSite[sites[i]->index].push_back(cliques.size());
//             cliquesDuSite[sites[j]->index].push_back(cliques.size());
//             simc.blobs.push_back(sites[i]);
//             simc.blobs.push_back(sites[j]);
//             cliques.push_back(simc);
//   //           }
//           }
//         }
//       }
//     }
//   
//     cout << "TEMP :" << temp<< " " << temp2 << " " << temp3 << " " << temp4 << " "<<temp+temp4<< " ";
// //   cin >> temp;
// 
//     return cliques;
// }

// float maxim(float


double getOverlap(Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap){

          float overlap_x,overlap_y,overlap_z,aux;
          double rec=0.0;
            if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) < 0.0001) {bbmax1[0] += 0.5; /*cout << "bbmax10+ ";*/}
            if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) < 0.0001) {bbmax1[1] += 0.5; /*cout << "bbmax11+ ";*/}
            if (sqrt(pow(bbmin2[0]-bbmax2[0],2)) < 0.0001) {bbmax2[0] += 0.5; /*cout << "bbmax20+ ";*/}
            if (sqrt(pow(bbmin2[1]-bbmax2[1],2)) < 0.0001) {bbmax2[1] += 0.5; /*cout << "bbmax21+ ";*/}
//           if (bbmin1[1]>bbmax1[1] && bbmin2[1] < bbmax2[1] ) {//alors i a bouclé autour de 360/0
            if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) >300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) <300){
//             cout << "i boucle lon " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//             ASSERT(bbmin1[1]>bbmax2[1]);
              if (360-bbmax2[1]<bbmin2[1]){
                aux = bbmax1[1];
                bbmax1[1] = bbmin1[1] + 360.0;
                bbmin1[1] = aux;
              }
              else {
                aux = bbmin1[1];
                bbmin1[1] = bbmax1[1] - 360.0;          
                bbmax1[1] = aux;
              }
          }
//           else if (bbmin1[1]<bbmax1[1] && bbmin2[1] > bbmax2[1] ) {//alors j a bouclé autour de 360/0 
          else if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) <300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) >300){  
//             cout << "j boucle lon " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//             ASSERT(bbmin2[1]>bbmax1[1]);
              if (360-bbmax1[1]<bbmin1[1]){
                aux = bbmax2[1];
                bbmax2[1] = bbmin2[1] + 360.0;
                bbmin2[1] = aux;
              }
              else {
                aux = bbmin2[1];
                bbmin2[1] = bbmax2[1] - 360.0;          
                bbmax2[1] = aux;
              }
          }
//           else if (bbmin1[1]>bbmax1[1] && bbmin2[1]>bbmax2[1] ) {//alors i&j ont bouclé
          else if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) >300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) >300){   
//               cout << "i et j bouclent lon " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
                aux = bbmin1[1];
                bbmin1[1] = bbmax1[1] - 360.0;
                bbmax1[1] = aux;
                aux = bbmin2[1];
                bbmin2[1] = bbmax2[1] - 360.0;
                bbmax2[1] = aux;
          } 
        // on s'occupe de la latitude
//           if (bbmin1[0]>bbmax1[0] && bbmin2[0] < bbmax2[0] ) {//alors i a bouclé autour de 360/0 
        if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) >150 && sqrt(pow(bbmin2[0]-bbmax2[0],2)) <150){  
//             cout << "i boucle lat" << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//             ASSERT(bbmin1[0]>bbmax2[0]);
              if (180-bbmax2[0]<bbmin2[0]){
                aux = bbmax1[0];
                bbmax1[0] = bbmin1[0] + 180.0;
                bbmin1[0] = aux;
              }
              else {
                aux = bbmin1[0];
                bbmin1[0] = bbmax1[0] - 180.0;          
                bbmax1[0] = aux;
              }
          }
//           else if (bbmin1[0]<bbmax1[0] && bbmin2[0] > bbmax2[0] ) {//alors j a bouclé autour de 360/0 
        else if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) <150 && sqrt(pow(bbmin2[0]-bbmax2[0],2)) >150){  
//             cout << "j boucle lat" << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//             ASSERT(bbmin2[0]>bbmax1[0]);
              if (180-bbmax1[0]<bbmin1[0]){
                aux = bbmax2[0];
                bbmax2[0] = bbmin2[0] + 180.0;
                bbmin2[0] = aux;
              }
              else {
                aux = bbmin2[0];
                bbmin2[0] = bbmax2[0] - 180.0;          
                bbmax2[0] = aux;
              }
            
          }
//           else if (bbmin1[0]>bbmax1[0] && bbmin2[0]>bbmax2[0] ) {//alors i&j ont bouclé
        else if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) >150 && sqrt(pow(bbmin2[0]-bbmax2[0],2)) >150){  
//               cout << "i et j bouclent lat" << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
                aux = bbmin1[0];
                bbmin1[0] = bbmax1[0] - 360.0;
                bbmax1[0] = aux;
                aux = bbmin2[0];
                bbmin2[0] = bbmax2[0] - 360.0;
                bbmax2[0] = aux;
          } 
          // prétraitements effectués on calcule le recouvrement
//           if (*no_overlap==0) cout << "rec: " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
          *no_overlap=0;
          if (bbmin1[0]<=bbmin2[0])
            if (bbmax1[0]<bbmin2[0]) *no_overlap=1;
          else overlap_x= (bbmax2[0] < bbmax1[0] ? bbmax2[0] : bbmax1[0]) - bbmin2[0] +1;
          else
            if (bbmax2[0]<bbmin1[0]) *no_overlap=1;
          else overlap_x= (bbmax1[0] < bbmax2[0] ? bbmax1[0] : bbmax2[0]) - bbmin1[0] +1;
          if (*no_overlap==0)
          {
            if (bbmin1[1]<=bbmin2[1])
              if (bbmax1[1]<bbmin2[1]) *no_overlap=1;
            else overlap_y= (bbmax2[1] < bbmax1[1] ? bbmax2[1] : bbmax1[1]) - bbmin2[1] +1;
            else
              if (bbmax2[1]<bbmin1[1]) *no_overlap=1;
            else overlap_y= (bbmax1[1] < bbmax2[1] ? bbmax1[1] : bbmax2[1]) - bbmin1[1] +1;
            if (*no_overlap==0)
            {
              rec=overlap_x*overlap_y;
              double div=( ((bbmax1[0]-bbmin1[0])*(bbmax1[1]-bbmin1[1]) +1)
                    + ((bbmax2[0]-bbmin2[0])*(bbmax2[1]-bbmin2[1]) +1) );
      
//           if (*no_overlap==0 && rec > 0.1) {cout << "rec: " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl; cout << rec << " " << div << " " << 2*rec/div << endl;}

              rec=2 * rec / div;
          //                 cout << rec << " " ;

            }
            
          }

  return rec;








}

vector<Clique> ConstruireCliquesLastChance(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons)  {
  uint temp=0,temp2=0,temp3=0,temp4=0;
  cliquesDuSite = vector<vector<int> >(sites.size());
  vector<Clique> cliques;
  set<uint> v;
  vector<string> subjects;
  double rec;
  set<uint> setnodes;
  float x,y;
  int mini=0; //,dep,arr;
  map<string,vector<map<string,uint> > > matching;
  set<uint>::iterator it;
  map<string, AimsSurfaceTriangle>::iterator meshit,meshjt;  
  float dist_thresh = 50.0;
  cout << "Distance maps threshold : << " << dist_thresh <<" >>"<< endl;


//   typedef til::Mesh_N MyMesh;
//   map<string,MyMesh> mymeshes;
// 
//   for (meshit = meshes.begin();meshit!=meshes.end();meshit++){
//     til::Mesh1 mesh0;
//     til::convert(mesh0, meshit->second);
//     mymeshes[meshit->first] = addNeighborsToMesh(mesh0);
//   }
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
  
//   for (uint i=0;i<sites.size();i++)
//     for (megit=sites[i]->nodes_list.begin();megit!=sites[i]->nodes_list.end();megit++){
// //       cout << *megit << " ";
//       megasitesnodes.insert(*megit);
//     }
  
    
  typedef map<uint, set<uint> > mapvoisins;
  map<string, mapvoisins> neighbours;
  unsigned v1b, v2, v3;
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
  
  map<string, map<float, vector<pair<float, uint> > > > altmeshes;
  for (meshit=meshes.begin();meshit!=meshes.end();meshit++){
    matching[meshit->first]=vector<map<string,uint> > (meshit->second.vertex().size());
    altmeshes[meshit->first]=getAlternateMesh(meshit->second, lats[meshit->first], lons[meshit->first]);
  }
  cout << "Construction maillages alternatifs terminée ..." << endl;
    
  cout << "Recherche des correspondances noeud-à-noeud à travers tous les sujets..." << endl;
  set<uint> sitesnodes;
  for (uint i=0;i<sites.size();i++){
    
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"  << i << "/" << sites.size() << flush;
//     x = lats[sites[i]->subject].item(i);
//     y = lons[sites[i]->subject].item(i);
    
    x = sites[i]->gravitycenter[0];
    y = sites[i]->gravitycenter[1];
//     cout << " " << x << " " << y << " " << sites[i]->gravitycenter[2] << "//" ;
//     cout << sites[i]->subject << " " << x << " " << y << endl;
//     float testdistance=0.0, testdistmin=10000000.0;
    mini=-1;
//     for (uint j=0;j<meshes[sites[0]->subject].vertex().size();j++){
//       float x0 = lats[sites[0]->subject].item(j);
//       float y0 = lons[sites[0]->subject].item(j);
//       
//       testdistance = sqrt(pow(x0-x,2)+pow(y0-y,2));
//       if (testdistance < testdistmin){
//         testdistmin = testdistance;
//         mini=j;
//       }
//     }
    float delta = 0.5;
    while (mini==-1){
    mini=getNode(altmeshes[sites[0]->subject],delta,x,y,lats[sites[0]->subject],lons[sites[0]->subject]);
    delta *= 2.0;
    }
    sites[i]->node=(uint)mini;
//     cout << mini << "(suj:"<< sites[i]->subject <<"t:" << sites[i]->t<< ")" <<endl;
    
    sitesnodes.insert(mini);
    ASSERT(mini!=-1 && sites[i]->node<meshes[sites[0]->subject].vertex().size());
//     cout << i << " " << sites[i]->node << " - " ;

  }
  map<uint,map<uint,float> > distmaps;
   cout << endl;//(sitesnodes.size(),meshes[sites[0]->subject].vertex().size());
 uint cpt=0;
 
  for (it=sitesnodes.begin();it!=sitesnodes.end();it++,cpt++){
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << cpt << "/" << sitesnodes.size() << flush;
    distmaps[*it] = getDistMap( &(meshes[sites[0]->subject]),  neighbours[sites[0]->subject], *it, 50.0); 
//     distmaps[*it] = test_geodmap(mymeshes[sites[0]->subject], *it, 40.0);
  }
    

  
  cout << endl << "done" << endl;

  
//   for (uint i=0;i<sites.size();i++){
//     sites[i]->node = matching[sites[i]->subject][i][sites[0]->subject];
//   }

  for (uint i=0;i<intraps.size();i++)
    cliques.push_back(intraps[i]);

  for (uint i=0;i<sites.size(); i++){
    Clique c;
    c.type = DATADRIVEN; /*ls.type = BESTLOWERSCALE;*/
    cliquesDuSite[sites[i]->index].push_back(cliques.size());
    c.blobs.push_back(sites[i]);
    cliques.push_back(c);
//     cliquesDuSite[sites[i]->index].push_back(cliques.size());
//     ls.blobs.push_back(sites[i]);
//     cliques.push_back(ls);
  }
  set<uint> sitesaux,sitesindexes;
  float cliques_thresh=40.0;
  cout << "Similarity cliques distance limit : << " << cliques_thresh << " >>" << endl;
  for (uint i=0;i<sites.size();i++){
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << sites.size() << "(" << sites[i]->rank << "-" << sites[i]->subject << "-" << sites[i]->nodes_list.size() << ") " << flush;
    for (uint j=i;j<sites.size();j++) {
      if (sites[i]->subject != sites[j]->subject) {

//         uint cpt0=0;
//         for (it=sitesnodes.begin();*it!=sites[i]->node;it++,cpt0++){}
        if (distmaps[sites[i]->node].find(sites[j]->node) == distmaps[sites[i]->node].end())
          rec= 600.0;
        else 
          rec = distmaps[sites[i]->node][sites[j]->node];
//         rec = distmaps[sites[i]->node][j];
        
        // DISTANCE AVEC SURFACES
//           (*jv)->getProperty( "boundingbox_max", bbmax_2);
//           (*jv)->getProperty( "boundingbox_min", bbmin_2);
          Point3df bbmax1=sites[i]->boundingbox_max, bbmax2=sites[j]->boundingbox_max;
          Point3df bbmin1=sites[i]->boundingbox_min, bbmin2=sites[j]->boundingbox_min;
         // on s'occupe de la longitude
uint no_overlap=1;
        rec = getOverlap(bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap);
//         if (rec >0.1) cout << "REC: " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//         rec = sqrt(pow(meshes[sites[0]->subject].vertex()[sites[i]->node][0]-meshes[sites[0]->subject].vertex()[sites[j]->node][0],2)+pow(meshes[sites[0]->subject].vertex()[sites[i]->node][1]-meshes[sites[0]->subject].vertex()[sites[j]->node][1],2)+pow(meshes[sites[0]->subject].vertex()[sites[i]->node][2]-meshes[sites[0]->subject].vertex()[sites[j]->node][2],2));  // USING EUCLIDEAN DISTANCE
//         cout << i<< ";" << j << ";" << sites[i]->node << "-" << sites[j]->node << "rec:" << rec << endl;
        if (rec < 0.000000001) {sitesaux.insert(sites[i]->node); sitesaux.insert(sites[j]->node); sitesindexes.insert(i); sitesindexes.insert(j);}
        
        if ((no_overlap==0 && rec < cliques_thresh) && !((sites[j]->tmin > sites[i]->tmax) || (sites[i]->tmin > sites[j]->tmax)) /*&& (sites[i]->tValue * sites[j]->tValue > 2.0)*/) {
//           if (rec < 0.2) cout << "BLIP " << endl;
          Clique simc;
          simc.type = SIMILARITY;
          simc.rec = rec;
//           cout << simc.rec << endl ;
          cliquesDuSite[sites[i]->index].push_back(cliques.size());
          cliquesDuSite[sites[j]->index].push_back(cliques.size());
          simc.blobs.push_back(sites[i]);
          simc.blobs.push_back(sites[j]);
//           cout << "(("<<simc.rec << " (lab:" << simc.blobs[0]->label << " t:"<< simc.blobs[0]->t << "(nod:" <<  simc.blobs[0]->node << " suj:"<<simc.blobs[0]->subject<<")" << "-lab:" << simc.blobs[1]->label << " t:"<< simc.blobs[1]->t << "(nod:" <<  simc.blobs[1]->node << " suj:"<<simc.blobs[1]->subject<<")" << ")=>" << simc.energie << ")) " << endl;
          cliques.push_back(simc);
  //           }
        }
      }
    }
  }
  cout << sitesaux.size() << " " << sitesindexes.size() << endl;
  cout << "TEMP :" << temp<< " " << temp2 << " " << temp3 << " " << temp4 << " "<<temp+temp4<< " ";

  return cliques;
}


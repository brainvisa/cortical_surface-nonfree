/* Copyright (c) 1995-2007 CEA
 *
 *  This software and supporting documentation were developed by
 *      CEA/DSV/SHFJ
 *      4 place du General Leclerc
 *      91401 Orsay cedex
 *      France
 *
 * This software is governed by the CeCILL license version 2 under
 * French law and abiding by the rules of distribution of free software.
 * You can  use, modify and/or redistribute the software under the
 * terms of the CeCILL license version 2 as circulated by CEA, CNRS
 * and INRIA at the following URL "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license version 2 and that you accept its terms.
 */
#include <iostream>
#include <cstdlib>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>
#include <iomanip>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <aims/graph/graphmanip.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/primalsketch/scalespace.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <aims/primalsketch/primalSketch.h>
#include <cortical_surface/structuralanalysis/representation.h>
#include <cortical_surface/structuralanalysis/cliques.h>
#include <cortical_surface/structuralanalysis/iograph.h>
#include <cortical_surface/structuralanalysis/anneal.h>
#include <cortical_surface/structuralanalysis/blobs.h>
using namespace aims;
using namespace carto;
using namespace std;


class SubjectData{
  public :
    string subject;
    AimsSurfaceTriangle mesh;
    AimsSurfaceTriangle repMesh;
    TimeTexture<float> tex;
    TimeTexture<float> lat;
    TimeTexture<float> lon;
};

//##############################################################################

vector<int> set2vector(set<int> &s){

  vector<int> v;
  set<int>::iterator it;
  for (it=s.begin();it!=s.end();it++)
    v.push_back(*it);
  return v;

}

//##############################################################################

set<int> vector2set(vector<int> &v){

  set<int> s;
  for (uint i=0;i<v.size();i++)
    s.insert(v[i]);
  return s;

}


//##############################################################################

void RecoverBlobsFromIndivGraph( Graph *graph,
                            SubjectData &subjData,
                            vector<surf::GreyLevelBlob *> &blobs,
                            vector<surf::ScaleSpaceBlob *> &ssblobs,
                            bool initNull = true){
    if (initNull){
      blobs.clear();
      ssblobs.clear();
    }

    set<Vertex *>::iterator iv;
    Edge *e;
    Vertex::iterator jv;
    Edge::iterator kv;
    string meshPath, sujet, texPath, latPath, lonPath, repMeshPath;

    graph->getProperty("mesh", meshPath);
    graph->getProperty("sujet", sujet);
    graph->getProperty("texture", texPath);
    graph->getProperty("latitude", latPath);
    graph->getProperty("longitude", lonPath);
    graph->getProperty("representation_mesh", repMeshPath);
    subjData.subject = sujet;

    cout << "Subject " << sujet << endl;
    cout << " mesh : " << meshPath << endl;
    cout << " tex : " << texPath << endl;
    cout << " lat : " << lonPath << endl << endl;
    cout << " repMesh : " << repMeshPath << endl << endl;

    cout << "Reading files..." <<  flush;
    Reader<AimsSurfaceTriangle> rdrMesh(meshPath);
    rdrMesh.read(subjData.mesh);

    Reader<AimsSurfaceTriangle> rdrRepMesh(repMeshPath);
    rdrRepMesh.read(subjData.repMesh);

    Reader<TimeTexture<float> > rdrTex(texPath);
    rdrTex.read(subjData.tex);

    Reader<TimeTexture<float> > rdrLat(latPath);
    rdrLat.read(subjData.lat);

    Reader<TimeTexture<float> > rdrLon(lonPath);
    rdrLon.read(subjData.lon);



    cout << " done" << endl;
    vector<Vertex *> listVertSSB, listVertGLB;
    map<int, set<int> > listGLBindices;
    int iNbLinks = 0;
    int iNbGLB = 0;
    int iNbSSB = 0;
    for (iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv){
      if ((*iv)->getSyntax() == "glb"){
        int index;
        float scale, t;
        vector<int> nodes_list;
        (*iv)->getProperty("scale", scale);
        (*iv)->getProperty("nodes", nodes_list);
        (*iv)->getProperty("t", t);
        blobs.push_back(new surf::GreyLevelBlob());
        surf::GreyLevelBlob *blob = blobs[blobs.size()-1];

        blob->index = iNbGLB++;
        index = blob->index;
        (*iv)->setProperty("index", (int) index );
//         (*iv)->getProperty("index", index );
//         cout << "idx: " << index << " " << flush;


        blob->nodes = vector2set(nodes_list);
        blob->scale = scale;
        blob->t = t;

      }
    }
    for (iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv){
      if ((*iv)->getSyntax() == "ssb"){
        int index;
        float tmax, tmin, t;

        (*iv)->getProperty("tmax", tmax);
        (*iv)->getProperty("tmin", tmin);
        (*iv)->getProperty("t", t);
        ssblobs.push_back(new surf::ScaleSpaceBlob());
        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];

        ssblob->index = iNbSSB++;
        index = ssblob->index;
        (*iv)->setProperty("index", (int) index );
//         cout << "idx2:" << index << flush;
        ssblob->tmax = tmax;
        ssblob->tmin = tmin;
        ssblob->t = t;
        ssblob->subject = sujet;
        if (listGLBindices.find(index) == listGLBindices.end())
          listGLBindices[index] = set<int>();

        for (jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++){
          e = *jv;
          if (e->getSyntax() == "s2g"){
            for (kv = e->begin() ; kv != e->end() ; kv++){
              if ((*kv)->getSyntax() == "ssb"){

              }
              else if ((*kv)->getSyntax() == "glb"){
                int iGLBindex;
                (*kv)->getProperty("index", iGLBindex);
//                 cout << "igl:" << iGLBindex << " " << flush;
                listGLBindices[index].insert(iGLBindex);
                iNbLinks++;
              }
            }
          }
        }

      }
    }
    cout << "Rebuilding links..." << endl;
    for (uint i = 0 ; i < ssblobs.size() ; i++){
      int index = ssblobs[i]->index;
      set<int>::iterator it;
      for (it = listGLBindices[index].begin() ; it != listGLBindices[index].end() ; it++){
        ssblobs[i]->blobs.insert(blobs[*it]);

      }
    }

    cout << iNbGLB << " blobs added" << endl;
    cout << iNbSSB << " ssblobs added " << endl;
    cout << iNbLinks << " links added" << endl;

}

//##############################################################################

double getOverlapMeasure(Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap){

  float overlap_x,overlap_y,aux;
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

    aux = bbmin1[0];
    bbmin1[0] = bbmax1[0] - 360.0;
    bbmax1[0] = aux;
    aux = bbmin2[0];
    bbmin2[0] = bbmax2[0] - 360.0;
    bbmax2[0] = aux;
  }
          // prétraitements effectués on calcule le recouvrement
  *no_overlap=0;
  if (bbmin1[0]<=bbmin2[0])
    if (bbmax1[0]<bbmin2[0]) *no_overlap=1;
  else overlap_x= (bbmax2[0] < bbmax1[0] ? bbmax2[0] : bbmax1[0]) - bbmin2[0] ;
  else
    if (bbmax2[0]<bbmin1[0]) *no_overlap=1;
  else overlap_x= (bbmax1[0] < bbmax2[0] ? bbmax1[0] : bbmax2[0]) - bbmin1[0];
  if (*no_overlap==0)
  {
    if (bbmin1[1]<=bbmin2[1])
      if (bbmax1[1]<bbmin2[1]) *no_overlap=1;
    else overlap_y= (bbmax2[1] < bbmax1[1] ? bbmax2[1] : bbmax1[1]) - bbmin2[1];
    else
      if (bbmax2[1]<bbmin1[1]) *no_overlap=1;
    else overlap_y= (bbmax1[1] < bbmax2[1] ? bbmax1[1] : bbmax2[1]) - bbmin1[1];
    if (*no_overlap==0)
    {
//       cout << "overlap_x :" << overlap_x << " overlap_y :" << overlap_y << endl;
      rec=overlap_x*overlap_y;
      double div=( ((bbmax1[0]-bbmin1[0])*(bbmax1[1]-bbmin1[1]))
            + ((bbmax2[0]-bbmin2[0])*(bbmax2[1]-bbmin2[1])));

      rec=2 * rec / div;

    }

  }

  return rec;
}


//##############################################################################

// This function takes the "ssblobs" vector and figures out which pairs of blobs
//  overlap. The resulting vector "cliques" associates to every relevant pair of
//  scale-space blobs (noted by their indices) its calculated spatial overlap.
vector<surf::SSBClique> construireCliques ( vector<surf::ScaleSpaceBlob *>   &ssblobs,
                                      vector<surf::GreyLevelBlob *>     &blobs,
                                      map<string, SubjectData> &data,
                                      vector<vector<surf::GreyLevelBlob *> > &matchingblobs){

      vector<surf::SSBClique > cliques;
      matchingblobs = vector<vector<surf::GreyLevelBlob *> > (ssblobs.size());

      set<surf::GreyLevelBlob *>::iterator itB1, itB2;
      surf::GreyLevelBlob *b1max, *b2max;

      // Start of cliques construction

      for (uint i=0 ; i < ssblobs.size() - 1 ; i++){
        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << ssblobs.size() << "(" << cliques.size() << ")" << flush;
        for (uint j=i+1 ; j < ssblobs.size() ; j++){

          // For every single pair of scale-space blobs, computes a maximal overlap
          //   between every possible pair of grey-level blobs.

          // We consider only pairs of scale-space blobs from different subjects.
          if (ssblobs[i]->subject != ssblobs[j]->subject){

            float overmax=-1.0;

            for (itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() ; itB1++){
              for (itB2 = ssblobs[j]->blobs.begin() ; itB2 != ssblobs[j]->blobs.end() ; itB2++){

                // For every possible pair of grey-level blobs between these two scale-
                //   space blobs, we figure out their possible spatial overlap.

//                 vector<int> listNodesB1(set2vector((*itB1)->nodes_list)),
//                             listNodesB2(set2vector((*itB2)->nodes_list));

                pair<Point2df,Point2df> bbi = getBoundingBox((*itB1)->nodes, data[ssblobs[i]->subject].lat, data[ssblobs[i]->subject].lon),
                                        bbj = getBoundingBox((*itB2)->nodes, data[ssblobs[j]->subject].lat, data[ssblobs[j]->subject].lon);

                Point3df bbmin1 (bbi.first[0], bbi.first[1], 0.0),
                        bbmax1 (bbi.second[0], bbi.second[1], 0.0),
                        bbmin2 (bbj.first[0], bbj.first[1], 0.0),
                        bbmax2 (bbj.second[0], bbj.second[1], 0.0) ;

                uint no_overlap = 2;
                double overlap = getOverlapMeasure( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );


                if (no_overlap == 0 ){

                  // If the current pair's overlap is maximal, then the glb indices are stored.

//                   cout << "bbi("<< (*itB1)->nodes_list.size() << "):" << bbi.first[0] << "-" << bbi.first[1] << " " <<
//                       bbi.second[0] << "-" << bbi.second[1] << " " <<
//                       "bbj("<< (*itB2)->nodes_list.size() << "):" << bbj.first[0] << " " << bbj.first[1] << " " <<
//                       bbj.second[0] << " " << bbj.second[1] << " over:" << overlap << endl;
//                   cout << (*itB1)->scale << " " << (*itB2)->scale << endl;

                  if (overlap > overmax) {
                    overmax = overlap;
                    b1max = *itB1;
                    b2max = *itB2;
                  }
                }

              }
            }


            // Here all the possible glb pairs have been processed for the two current ssb

            if (overmax > 0.10 &&
                !((ssblobs[j]->tmin > ssblobs[i]->tmax) || (ssblobs[i]->tmin > ssblobs[j]->tmax))){

              // If the two scale-space blobs have at least one pair of grey-level
              //   overlapping (bounding-boxes) (+ scales overlapping), then a clique
              // is created between these two ssb and the max-overlapping pair of glb
              // is stored in "matchingblobs".



              cliques.push_back(surf::SSBClique(ssblobs[i], ssblobs[j], overmax));
//               assert(!( res.first[0] == res.first[1] && ssblobs[i]->subject == ssblobs[j]->subject));
              matchingblobs[i].push_back(b1max);
              matchingblobs[j].push_back(b2max);
//               cout << "max (" << ssblobs[i]->index <<","<< ssblobs[j]->index << ") between:" << b1max->index << " "
//                    << b2max->index << " overmax:" << overmax << endl;
//               cout << "scales: " << b1max->scale << " " << b1max->scale << endl;

            }
          }

          // The next pair of scale-space blobs will now be processed.
        }
      }
      cout << ssblobs.size() << "/" << ssblobs.size() << "(" << cliques.size() << ")" << endl;
      // Construction of a representation blob for each scale-space blob
      for (uint i = 0 ; i < ssblobs.size() ; i++){

        // For every scale-space blob, we create a representation blob
        //   from the set of grey-level blobs found to be max-matching
        //   with some others (from other scale-space blobs)
        set<uint>::iterator it;

        if (matchingblobs[i].size()!=0)
          cout << i << ":";

//         for (it = matchingblobs[i].begin() ; it != matchingblobs[i].end() ; it++){
        for (uint j = 0 ; j < matchingblobs[i].size() ; j++){
          set<int> blobNodes(matchingblobs[i][j]->nodes);
          ssblobs[i]->nodes.insert(blobNodes.begin(), blobNodes.end());
          cout << ssblobs[i]->nodes.size() << " " << flush;
        }

        if (matchingblobs[i].size()!=0)
          cout << endl ;

      }
//       float test;
//       cin >> test;
      return cliques;
}




//##############################################################################

void ConstruireGrapheGroupe(Graph *graph,
//                       vector<Blob *> &blobs,
                            vector<surf::ScaleSpaceBlob *> &ssblobs,
                            vector<surf::SSBClique> &cliques,
                            string texPaths,
                            string meshPaths,
                            string latPaths,
                            string lonPaths,
                            string indivGraphPaths,
                            string repMeshPaths,
                            vector<string> &listSujets){

  cerr << "Construction du Graphe" << endl;
  vector<float> resolution,bbmin2D,bbmax2D;
  vector<int> bbmin, bbmax;
  resolution.push_back(1.0); resolution.push_back(1.0); resolution.push_back(1.0);
  //   bbmin.push_back(mesh[0].minimum()[0]-1); bbmin.push_back(mesh[0].minimum()[1]-1); bbmin.push_back(mesh[0].minimum()[2]-1);
  //   bbmax.push_back(mesh[0].maximum()[0]+1); bbmax.push_back(mesh[0].maximum()[1]+1); bbmax.push_back(mesh[0].maximum()[2]+1);
  bbmin.push_back(-10); bbmin.push_back(-10); bbmin.push_back(-10);
  bbmax.push_back(10); bbmax.push_back(10); bbmax.push_back(10);
  graph->setProperty( "filename_base", "*");

  graph->setProperty("voxel_size", resolution);
  graph->setProperty("boundingbox_min", bbmin);
  graph->setProperty("boundingbox_max", bbmax);
  graph->setProperty("meshes", splitGraphFile(meshPaths));
  graph->setProperty("sujets", listSujets);
  graph->setProperty("textures", splitGraphFile(texPaths));
  graph->setProperty("latitudes", splitGraphFile(latPaths));
  graph->setProperty("longitudes", splitGraphFile(lonPaths));
  graph->setProperty("indiv_graphs", splitGraphFile(indivGraphPaths));
  graph->setProperty("representation_meshes", splitGraphFile(repMeshPaths));

//   AimsSurfaceTriangle *objects = new AimsSurfaceTriangle();
//   *objects = getBlobsSphericalMeshes( ssblobs, repMesh[repMesh.size()-1], lat[0], lon[0], nodes_lists);

  Vertex *vert;
  carto::rc_ptr<AimsSurfaceTriangle> ptr;
  aims::GraphManip manip;
  vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( ssblobs.size() );


  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {

    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index

    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("ssb");
//     int index1 = i;
    ssblobs[i]->index = i;
    vert->setProperty("index", i);
    vert->setProperty("label", "0");
    vert->setProperty("t", ssblobs[i]->t);
    vert->setProperty( "subject", ssblobs[i]->subject);
    vert->setProperty( "tmin", ssblobs[i]->tmin);
    vert->setProperty( "tmax", ssblobs[i]->tmax);
    vert->setProperty( "tValue", 100.0);
    vert->setProperty("nodes", set2vector(ssblobs[i]->nodes));


//     // We associate the proper mesh patch from "objects" to the vertex
//     ptr=carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
//     (*ptr)[0]=(*objects)[i];
//     manip.storeAims(*graph, vert, "glb", ptr);
//     vert->setProperty("glb_label",i);

    listVertSSB[ i ] = vert;
  }

  // CONSTRUIRE LES ARETES
  cout << "Construction cliques" << endl;
  for (uint i = 0 ; i < cliques.size() ; i++){

    // For every clique, we get the two corresponding vertices from listVertices
    Vertex *v1, *v2;
    int aux1, aux2;
    for (aux1 = 0 ;
         aux1 < (int) ssblobs.size() &&
             !(cliques[i].ssb1->index == ssblobs[aux1]->index &&
                 cliques[i].ssb1->subject == ssblobs[aux1]->subject) ;
         aux1++){  }
    for (aux2 = 0 ;
         aux2 < (int) ssblobs.size() &&
             !(cliques[i].ssb2->index == ssblobs[aux2]->index &&
             cliques[i].ssb2->subject == ssblobs[aux2]->subject) ;
         aux2++){  }
    assert(aux1 != ssblobs.size() && aux2 != ssblobs.size());
    v1 = listVertSSB[aux1];
    v2 = listVertSSB[aux2];

    Edge *edge= graph->addEdge(v1, v2, "b2b");

    edge->setProperty("similarity", cliques[i].similarity);

  }

}

//##############################################################################

void recoverSubjectBlobs ( string sujet,
                           vector<surf::GreyLevelBlob *> &blobs,
                           vector<surf::ScaleSpaceBlob *> &ssblobs,
                           vector<surf::GreyLevelBlob *> &subjBlobs,
                           vector<surf::ScaleSpaceBlob *> &subjSsblobs){
  subjBlobs.clear();
  subjSsblobs.clear();
//   for (uint i = 0 ; i < blobs.size() ; i++){
//     blobs[i]->index = i;
//   }
  set<surf::GreyLevelBlob *>::iterator it;
  for (uint i = 0 ; i < ssblobs.size() ; i++){
    ssblobs[i]->index = i;
    if (ssblobs[i]->subject == sujet){
      subjSsblobs.push_back(ssblobs[i]);
      for (it = ssblobs[i]->blobs.begin() ; it != ssblobs[i]->blobs.end() ; it++){
        subjBlobs.push_back(*it);
      }
    }
  }

}


//##############################################################################


int main( int argc, const char **argv ){
  try {

    string indivGraphPaths = "",
           groupGraphPath = "",
           meshPaths = "",
           texPaths = "",
           latPaths = "",
           lonPaths = "",
           repMeshPaths = "",
           sujets = "";

    AimsApplication app( argc, argv, "surfMesh2Graph" );
    app.addOption( meshPaths, "-m", "mesh");
    app.addOption( texPaths, "-t", "texture");
    app.addOption( indivGraphPaths, "-g", "indiv graphs");
    app.addOption( groupGraphPath, "-G", "group graph");
    app.addOption( sujets, "-s", "sujet");
    app.addOption( latPaths, "--lat", "latitude");
    app.addOption( lonPaths, "--lon", "longitude");
    app.addOption( repMeshPaths, "--repM", "repMesh",1);
    app.initialize();








        vector<string> listSujets = splitGraphFile(sujets);
        vector<string> listGraphPaths = splitGraphFile(indivGraphPaths);
        map<string, SubjectData> data;
        vector<surf::GreyLevelBlob *> blobs;
        vector<surf::ScaleSpaceBlob *> ssblobs;
        // Processing every subject...
        for (uint i = 0 ; i < listSujets.size() ; i++) {

          string sujet = listSujets[i];
          cout << "sujet : " << sujet << endl;
          Reader<Graph> rdrIndivGraph(listGraphPaths[i]);
          Graph *graph = new Graph("BlobsArg");
          rdrIndivGraph.read(*graph);
          pair<string, SubjectData> subjData;
          subjData.first = sujet;
          RecoverBlobsFromIndivGraph(graph, subjData.second, blobs, ssblobs, false);
          data.insert(subjData);

        }

        // Building the cliques...
        vector<surf::SSBClique> cliques;
        vector<vector<surf::GreyLevelBlob *> > matchingblobs;
        cliques = construireCliques(ssblobs, blobs, data, matchingblobs);
        cout << endl << cliques.size() << " cliques de similarité " << endl;

        // Building the group graph...
        Graph *graph = new Graph("BlobsArg");
        ConstruireGrapheGroupe(graph, /*blobs,*/ ssblobs, cliques, texPaths, meshPaths, latPaths, lonPaths, indivGraphPaths, repMeshPaths, listSujets);

        Writer<Graph> wtrGraph(groupGraphPath);
        wtrGraph.write(*graph);
    return EXIT_SUCCESS;
  }
  catch( carto::user_interruption & )
  {
  }
  catch( exception & e )
  {
    cerr << e.what() << endl;
  }
}


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
#include "blobs.h"




using namespace aims;
using namespace carto;
using namespace std;


//##############################################################################

// Function that builds a collection of Blob and SSBlob objects from a previously
//   computed Primal Sketch
void construireBlobs(PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch, 
                     vector<Blob *> &blobs, vector<SSBlob *> &ssblobs){
  
//     // Inititalization of the results vectors "blobs" and "ssblobs"
//     blobs = vector<Blob*>();
//     ssblobs = vector<SSBlob*>();
    
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*> listBlobs 
         = sketch.BlobSet();
    uint iBlob=0, iSSBlob=0;
    
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*>::iterator itSSB;
    list<GreyLevelBlob<SiteType<AimsSurface<3, Void> >::type > *>::iterator itGLB;
    
    ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *ssb;
    set< SiteType<AimsSurface<3, Void> >::type, 
      ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> >::iterator itPoints;
    
    
    for (itSSB=listBlobs.begin();itSSB!=listBlobs.end();itSSB++){
      
      // For each scale-space blob, we create a SSBlob in "ssblobs" containing 
      //    various Blob objects (being themselves contained in a general resulting 
      // "blobs" vector).
      ssb = *itSSB;
      ssblobs.push_back(new SSBlob());
      SSBlob *ssblob = ssblobs[ssblobs.size()-1];
      ssblob->index = iSSBlob;
      ssblob->subject = sketch.Subject();
      ssblob->tmin = 999.0;
      ssblob->tmax = -999.0;
      
      for (itGLB = ssb->glBlobs.begin(); itGLB != ssb->glBlobs.end(); itGLB++){
        
        // For each grey-level blob, we create a Blob
        blobs.push_back(new Blob());
        Blob *blob = blobs[blobs.size()-1];
        
        // Each Blob has a specific index iBlob, and a SSBlob has an iSSBlob
        blob->index = iBlob++;
        blob->parent = iSSBlob;
//         blob->subject = sketch.Subject();
        blob->t = (*itGLB)->measurements.t;
        blob->scale = (*itGLB)->GetScale();
        
        // The Blob's nodeslist contains its corresponding nodes indices on the 
        //    mesh it was extracted from.
        set<SiteType<AimsSurface<3, Void> >::type, 
             ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> > listePoints 
                 = (*itGLB)->GetListePoints();
        for (itPoints = listePoints.begin() ; itPoints != listePoints.end() ; itPoints++)
          (blob->nodes_list).insert((*itPoints).second);
        
        ssblob->blobs.insert(blob);
        
        if (blob->scale < ssblob->tmin)
          ssblob->tmin=blob->scale;
        if (blob->scale > ssblob->tmax)
          ssblob->tmax=blob->scale;
      }
      
      ssblob->t = ssb->GetMeasurements().t;
//       cout << ssblob->t << " " << ssb->GetMeasurements().tValue << " " << ssblob->blobs.size() << endl;
//       cout << ssblob->tmin << "/" << ssb->ScaleMin() << " " << ssblob->tmax << "/" << ssb->ScaleMax() << endl;
      
      
      iSSBlob++;
    }
    
    cout << " iBlob : " << iBlob << 
            " iSSblob : " << iSSBlob <<
            " blobs.size : " << blobs.size() << 
            " ssblobs.size : " << ssblobs.size() << endl;
}

//##############################################################################

// Function that creates a Scale Space from a texture, a mesh and two coordinates textures
ScaleSpace<AimsSurface<3, Void>, Texture<float> > getScaleSpace(
                                  TimeTexture<float> &laTexture, 
                                  AimsSurfaceTriangle &laMesh, 
                                  TimeTexture<float> &lat, 
                                  TimeTexture<float> &longit){
                                    
      float moy=0.0;
      for (uint i=0;i<laTexture.nItem();i++)
        if (laTexture.item(i) == laTexture.item(i))
          moy += laTexture.item(i);
         
      moy /= laTexture.nItem();
      for (uint i=0;i<laTexture.nItem();i++)
        if (laTexture.item(i) != laTexture.item(i))
          laTexture.item(i) = moy;

      cout << "Size=" << laMesh[0].vertex().size() << endl;
      cout << "Ssmoother creation" << endl;

      FiniteElementSmoother<3, float> *smooth;
      smooth=new FiniteElementSmoother<3, float>(0.01, &(laMesh[0]));
      
      ScaleSpace< AimsSurface<3, Void>, Texture<float> > scale_space(
            &(laMesh[0]), &(laTexture[0]), smooth);
            
      cout << "Scale-space creation" << endl;
      
      vector<Point3df> *coordinates;
      coordinates=new vector<Point3df>();
      
      for (uint i=0;i<lat[0].nItem();i++)
        (*coordinates).push_back(Point3df(lat[0].item(i), longit[0].item(i),i));
      
      scale_space.PutCoordinates(coordinates);
      scale_space.GenerateDefaultScaleSpace(8.0);
      return scale_space;
}


//##############################################################################

vector<set<uint> > getTriangles(AimsSurfaceTriangle &mesh){
  vector<set<uint> > triangles(mesh[0].vertex().size());
  uint p1,p2,p3;
  for (uint i=0;i<mesh[0].polygon().size();i++){
     p1=mesh[0].polygon()[i][0];
     p2=mesh[0].polygon()[i][1];
     p3=mesh[0].polygon()[i][2];
     triangles[p1].insert(i);
     triangles[p2].insert(i);
     triangles[p3].insert(i);
  }
  return triangles;
}


//##############################################################################

int find(const vector<int> &v, int item){
  int c=-1;
  for (uint i=0; i<v.size() && c==-1 ; i++){
    if (item == v[i]) c=(int)i;
  }
  return c;
}

//##############################################################################

// Function that extracts mesh patches from a "mesh", being given a "blobs" vector,
//   and returning a collection of "objects" plus a vector of "nodes_lists".
AimsSurfaceTriangle getBlobsMeshes( vector<Blob *> &blobs, 
                                    AimsSurfaceTriangle &mesh, 
                                    vector<vector<int> > &nodes_lists){
  AimsSurfaceTriangle objects;
  uint p1,p2,p3;
  nodes_lists=vector<vector<int> >(blobs.size());

  set<uint>::iterator it;

  for (uint i=0;i<blobs.size();i++){
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " << objects.size() << flush ;
    set<uint> tri,comp;
    vector<uint> corres;
    for (uint j=0;j<mesh[0].polygon().size();j++){
      
      p1=mesh[0].polygon()[j][0];
      p2=mesh[0].polygon()[j][1];
      p3=mesh[0].polygon()[j][2];
    
      if ( blobs[i]->nodes_list.find(p1)!=blobs[i]->nodes_list.end() &&     
              blobs[i]->nodes_list.find(p2)!=blobs[i]->nodes_list.end() &&    
              blobs[i]->nodes_list.find(p3)!=blobs[i]->nodes_list.end() )
        tri.insert(j);
      
    }

    for (it=tri.begin();it!=tri.end();it++){
      p1=mesh[0].polygon()[*it][0];
      p2=mesh[0].polygon()[*it][1];
      p3=mesh[0].polygon()[*it][2];
      comp.insert(p1); comp.insert(p2); comp.insert(p3);
    }
    corres=vector<uint>(mesh[0].vertex().size());
    for (it=comp.begin();it!=comp.end();it++){
      assert(*it<corres.size());
      assert(*it<mesh[0].vertex().size());
      assert(i<nodes_lists.size());
      (objects)[i].vertex().push_back(mesh[0].vertex()[*it]);
      corres[*it]=(objects)[i].vertex().size()-1;
      nodes_lists[i].push_back(*it);
    }

    for (it=tri.begin();it!=tri.end();it++){
      p1=mesh[0].polygon()[*it][0];
      p2=mesh[0].polygon()[*it][1];
      p3=mesh[0].polygon()[*it][2];
      (objects)[i].polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
    }
    
  }
  cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << blobs.size() << endl;
  
  return objects;
}

//##############################################################################

// That function takes a vector of SSBlob and build mesh patches corresponding to
//   previously computed representation blobs
AimsSurfaceTriangle getBlobsMeshes( vector<SSBlob *> &blobs, 
                                    AimsSurfaceTriangle &mesh, 
                                    vector<vector<int> > &nodes_lists){
                                      
                                      cout << "mesh.vertex:" << mesh[0].vertex().size() << endl;
                                      cout << "mesh.polygon:" << mesh[0].polygon().size() << endl;
                                      AimsSurfaceTriangle objects;
                                      uint p1,p2,p3;
                                      nodes_lists = vector<vector<int> >(blobs.size());

                                      set<uint>::iterator it;

                                      for (uint i = 0 ; i < blobs.size() ; i++) {
    
                                        cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " <<
                                                objects.size() << flush ;
                                        set<uint> tri,comp;
                                        vector<uint> corres;
                                        cout << endl << blobs[i]->representation.size() << endl;
                                        for (uint j=0;j<mesh[0].polygon().size();j++){
      
                                          p1=mesh[0].polygon()[j][0];
                                          p2=mesh[0].polygon()[j][1];
                                          p3=mesh[0].polygon()[j][2];
    
                                          if ( blobs[i]->representation.find(p1) != blobs[i]->representation.end()
                                            && blobs[i]->representation.find(p2) != blobs[i]->representation.end()
                                            && blobs[i]->representation.find(p3) != blobs[i]->representation.end()
                                             )
                                            tri.insert(j);
      
                                        }
                                        
                                        cout << "t " << tri.size() << " " << flush;
                                        
                                        for (it = tri.begin() ; it != tri.end() ; it++){
                                          p1=mesh[0].polygon()[*it][0];
                                          p2=mesh[0].polygon()[*it][1];
                                          p3=mesh[0].polygon()[*it][2];
                                          comp.insert(p1); comp.insert(p2); comp.insert(p3);
                                        }
                                        
                                        corres = vector<uint>( mesh[0].vertex().size() );
                                        
                                        for (it = comp.begin() ; it != comp.end() ; it++){
                                          
                                          assert( *it < corres.size() );
                                          assert( *it < mesh[0].vertex().size() );
                                          assert( i < nodes_lists.size() );
                                          
                                          (objects)[i].vertex().push_back( mesh[0].vertex()[*it] );
                                          corres[*it] = (objects)[i].vertex().size() - 1;
                                          nodes_lists[i].push_back( *it );
                                          
                                        }

                                        for (it = tri.begin() ; it != tri.end() ; it++){
                                          
                                          p1 = mesh[0].polygon()[*it][0];
                                          p2 = mesh[0].polygon()[*it][1];
                                          p3 = mesh[0].polygon()[*it][2];
                                          (objects)[i].polygon().push_back( AimsVector<uint,3>(corres[p1], corres[p2], corres[p3]) );
                                          
                                        }
    
                                      }
                                      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << blobs.size() << endl;
  
                                      return objects;
                                    }


//##############################################################################

// This function extracts normal mesh patches corresponding to the different label
//   components existing in a texture. It returns a vector of nodes lists.                                   
AimsSurfaceTriangle getObjects( TimeTexture<short> &tex, 
                                AimsSurfaceTriangle &mesh, 
                                vector<vector<int> > &nodes_lists){
    
    int labelmax=0;
    for (uint i=0;i<tex[0].nItem();i++){
      if (tex[0].item(i)>labelmax)
        labelmax=tex[0].item(i);
    }
    AimsSurfaceTriangle objects;
    nodes_lists=vector<vector<int> >(labelmax+1);
    vector<set<uint> > triangles(labelmax+1);
    
    set<uint>::iterator it;
    set<uint> tri,comp;
    uint p1,p2,p3;  int L1,L2,L3;
    
    for (uint i=0;i<mesh[0].polygon().size();i++){
      p1=mesh[0].polygon()[i][0];
      p2=mesh[0].polygon()[i][1];
      p3=mesh[0].polygon()[i][2];

      L1=tex[0].item(p1); L2=tex[0].item(p2); L3=tex[0].item(p3);
      
      if (L1==L2 || L1==L3)
        triangles[L1].insert(i);
      else if (L2==L3)
        triangles[L2].insert(i);
    }
    vector<uint> corres;
    for (uint i=0;i<triangles.size();i++){
      tri = triangles[i];

      comp.clear();
      for (it=tri.begin();it!=tri.end();it++){
        p1=mesh[0].polygon()[*it][0];
        p2=mesh[0].polygon()[*it][1];
        p3=mesh[0].polygon()[*it][2];
        comp.insert(p1); comp.insert(p2); comp.insert(p3);
      }
      corres=vector<uint>(mesh[0].vertex().size());
      for (it=comp.begin();it!=comp.end();it++){
        assert(*it<corres.size());
        assert(*it<mesh[0].vertex().size());
        assert(i<nodes_lists.size());
        (objects)[i].vertex().push_back(mesh[0].vertex()[*it]);
        corres[*it]=(objects)[i].vertex().size()-1;
        nodes_lists[i].push_back(*it);
      }

      for (it=tri.begin();it!=tri.end();it++){
        p1=mesh[0].polygon()[*it][0];
        p2=mesh[0].polygon()[*it][1];
        p3=mesh[0].polygon()[*it][2];
        (objects)[i].polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
      }
    }
    return objects;
}

//##############################################################################

// This function builds different spheres around the barycenters of every blob described
//  as a nodes list
AimsSurfaceTriangle getBarycenters ( AimsSurfaceTriangle  &mesh, 
                                     vector<vector<int> > &nodes_lists, 
                                     float                radius ) {
    AimsSurfaceTriangle objects;
    uint jmin;
    for (uint i=0;i<nodes_lists.size();i++){
      if (nodes_lists[i].size()!=0){
        jmin = nodes_lists[i][nodes_lists[i].size()/2];
        AimsSurfaceTriangle *msh;
        msh = SurfaceGenerator::sphere( mesh[0].vertex()[jmin], radius, 10 );
        objects[i]= (*msh)[0];
      }
    }
    return objects;
}

//##############################################################################

vector<int> set2vector(set<uint> &s){
  
  vector<int> v;
  set<uint>::iterator it;
  for (it=s.begin();it!=s.end();it++)
    v.push_back(*it);     
  return v;  
  
}

//##############################################################################

// This function takes the "ssblobs" vector and figures out which pairs of blobs 
//  overlap. The resulting vector "cliques" associates to every relevant pair of
//  scale-space blobs (noted by their indices) its calculated spatial overlap.
vector<pair<Point2d, float> > construireCliques ( vector<SSBlob *>   &ssblobs,
                                                  vector<Blob *>     &blobs, 
                                                  map<string, TimeTexture<float> > &latitudes,
                                                  map<string, TimeTexture<float> > &longitudes, 
                                                  vector<set<uint> > &matchingblobs){
    
      vector<pair<Point2d, float> > cliques;
      matchingblobs = vector<set<uint> > (ssblobs.size());
      
      set<Blob *>::iterator itB1, itB2;
      Blob *b1max, *b2max;
      
      // Start of cliques construction
            
      for (uint i=0 ; i < ssblobs.size() - 1 ; i++){
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
                
                vector<int> listNodesB1(set2vector((*itB1)->nodes_list)),
                            listNodesB2(set2vector((*itB2)->nodes_list));
                            
                pair<Point2df,Point2df> bbi = getBoundingBox(listNodesB1, latitudes[ssblobs[i]->subject], longitudes[ssblobs[i]->subject]),
                                        bbj = getBoundingBox(listNodesB2, latitudes[ssblobs[j]->subject], longitudes[ssblobs[j]->subject]);
                                                        
                Point3df bbmin1 (bbi.first[0], bbi.first[1], 0.0), 
                        bbmax1 (bbi.second[0], bbi.second[1], 0.0), 
                        bbmin2 (bbj.first[0], bbj.first[1], 0.0),
                        bbmax2 (bbj.second[0], bbj.second[1], 0.0) ;
                        
                uint no_overlap = 2;              
                double overlap = getOverlap( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
  
                
                if (no_overlap == 0 ){
                  
                  // If the current pair's overlap is maximal, then the glb indices are stored.
                  
  //                 cout << "bbi:" << bbi.first[0] << "-" << bbi.first[1] << " " <<
  //                     bbi.second[0] << "-" << bbi.second[1] << " " <<
  //                     "bbj:" << bbj.first[0] << " " << bbj.first[1] << " " <<
  //                     bbj.second[0] << " " << bbj.second[1] << " over:" << overlap << endl;
  //                 cout << (*itB1)->scale << " " << (*itB2)->scale << endl;
                  
                  if (overlap > overmax) {
                    overmax = overlap;
                    b1max = *itB1;
                    b2max = *itB2;
                  }
                }
            
              }
            }
          
          
            // Here all the possible glb pairs have been processed for the two current ssb
            
            if (overmax > 0.0 && 
                !((ssblobs[j]->tmin > ssblobs[i]->tmax) || (ssblobs[i]->tmin > ssblobs[j]->tmax))){
              
              // If the two scale-space blobs have at least one pair of grey-level
              //   overlapping (bounding-boxes) (+ scales overlapping), then a clique 
              // is created between these two ssb and the max-overlapping pair of glb
              // is stored in "matchingblobs".
              
              pair<Point2d, float> res;
              res.first = Point2d(i,j);
              res.second = overmax;
              cliques.push_back(res);
              matchingblobs[i].insert(b1max->index);
              matchingblobs[j].insert(b2max->index);            
              cout << "max (" << i <<","<<j<< ") between:" << b1max->index << " " 
                   << b2max->index << " overmax:" << overmax << endl;
              cout << "scales: " << b1max->scale << " " << b1max->scale << endl;
  
            }
          }
          
          // The next pair of scale-space blobs will now be processed.
        }
      }
      
      // Construction of a representation blob for each scale-space blob
      for (uint i = 0 ; i < ssblobs.size() ; i++){
        
        // For every scale-space blob, we create a representation blob
        //   from the set of grey-level blobs found to be max-matching 
        //   with some others (from other scale-space blobs)
        set<uint>::iterator it;
        
        if (matchingblobs[i].size()!=0) 
          cout << i << ":";
        
        for (it = matchingblobs[i].begin() ; it != matchingblobs[i].end() ; it++){
          
          set<uint> blobNodes(blobs[*it]->nodes_list);
          ssblobs[i]->representation.insert(blobNodes.begin(), blobNodes.end());
          cout << ssblobs[i]->representation.size() << " " << flush;
        }
        
        if (matchingblobs[i].size()!=0) 
          cout << endl ;
        
      }
      float test;
      cin >> test;
      return cliques;
}


//##############################################################################

void construireGraph (  Graph                         *graph,
                        vector<SSBlob *>              &ssblobs, 
                        vector<pair<Point2d, float> > &cliques,
                        AimsSurfaceTriangle           *objects ){

  
  Vertex *vert;
  carto::rc_ptr<AimsSurfaceTriangle> ptr;
  aims::GraphManip manip;
  vector<Vertex *> listVertices;
    
      
  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("blob");
    vert->setProperty("index", i);
    vert->setProperty("name", i);
    vert->setProperty("label", "0");
    vert->setProperty("t", 100.0);
    vert->setProperty("rank", i);
    vert->setProperty( "subject", ssblobs[i]->subject);
    vert->setProperty( "tmin", 1);
    vert->setProperty( "tmax", 4);
    vert->setProperty( "trep", 2);
    vert->setProperty( "depth", 100.0);
    vert->setProperty( "tValue", 100.0);
    vert->setProperty("nodes_list", set2vector(ssblobs[i]->representation));
      
    
    // We associate the proper mesh patch from "objects" to the vertex
    ptr=carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
    (*ptr)[0]=(*objects)[i];
    manip.storeAims(*graph, vert, "blob", ptr);
    vert->setProperty("blob_label",i);
    
    listVertices.push_back(vert);
  }

  // CONSTRUIRE LES ARETES
  cout << "Construction cliques" << endl;
  for (uint i = 0 ; i < cliques.size() ; i++){
    
    // For every clique, we get the two corresponding vertices from listVertices
    Vertex *v1, *v2;
    v1 = listVertices[cliques[i].first[0]];
    v2 = listVertices[cliques[i].first[1]];
    int index1, index2;
    v1->getProperty("index", index1);
    v2->getProperty("index", index2);
    cout << index1 << " " << index2 << endl;
    
    
    Edge *edge= graph->addEdge(v1,v2,"b2b");
    edge->setProperty("blob_first", cliques[i].first[0]);
    edge->setProperty("blob_second", cliques[i].first[1]);
    edge->setProperty("similarity", cliques[i].second);

  }    
}


//##############################################################################

// int main( int argc, const char **argv ){
//   try {
//   
//     int mode=0;
//     string outpath = "", 
//            meshPath, 
//            texPath, 
//            latpath="", 
//            lonpath="", 
//            flatpath="", 
//            sujet;
// 
//     AimsApplication app( argc, argv, "surfLabelsTex2Graph" );
//     app.addOption( meshPath, "-m", "mesh");
//     app.addOption( texPath, "-t", "texture");
//     app.addOption( outpath, "-o", "output file");
//     app.addOption( mode, "-M", "mode (0: normal - 1:barycenters (provide the lat/lon textures)",1);
//     app.addOption( sujet,"-s", "sujet");
//     app.addOption( latpath, "--lat", "latitude");
//     app.addOption( lonpath, "--lon", "longitude");
//     app.addOption( flatpath, "--flat", "flat",1);
//     app.initialize();
//     assert(latpath!="");
//     assert(lonpath!="");
//     
//     // Read files (the mesh, the tex and the coordinates)
//     Reader<AimsSurfaceTriangle> rdrMesh(meshPath);
//     Reader<TimeTexture<float> > rdrTex(texPath),
//                                 rdrLat(latpath),
//                                 rdrLon(lonpath);
//     AimsSurfaceTriangle mesh;
//     TimeTexture<float> tex, lat, lon;
//     rdrMesh.read(mesh);
//     rdrTex.read(tex);
//     rdrLat.read(lat);
//     rdrLon.read(lon);
//     
//     // Computes a radius for the "barycenters" sphere-based representation mode
//     mesh[0].setMini(); mesh[0].setMaxi();
//     cerr << mesh[0].minimum()[0] << " " << mesh[0].minimum()[1] << " " << mesh[0].minimum()[2] << ";" << mesh[0].maximum()[0] << " " << mesh[0].maximum()[1] << " " << mesh[0].maximum()[2] << endl;    
//     float dist = sqrt(pow(mesh[0].minimum()[0]-mesh[0].maximum()[0],2)+pow(mesh[0].minimum()[1]-mesh[0].maximum()[1],2)+pow(mesh[0].minimum()[2]-mesh[0].maximum()[2],2));
//     float radius = dist / 300.0;
//     
//     // Definition of blobs, nodes lists and objects vectors..
//     vector<vector<int> > nodes_lists;
//     AimsSurfaceTriangle *objects;
//     
//     vector<SSBlob *> ssblobs;
//     vector<Blob *> blobs;
//     vector<pair<Point2d, float> > cliques;
//     vector<set<uint> > matchingblobs;
// 
// 
//     if (mode == 1){
//       
//       // Barycenters mode
//       objects = new AimsSurfaceTriangle(getBarycenters(mesh,nodes_lists,radius));
//     }
//     else if (mode == 2){
//       
//       // Construction of a primal-sketch
//       ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss(getScaleSpace(tex,mesh,lat,lon));
// 
//       PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch(sujet, &ss, SURFACE);
// 
//       // Launching the computation of the PS (tmin, tmax, statfile, intersection_criterium)
//       sketch.ComputePrimalSketch(1.0, 8.0, "", 10);
// 
//       cout << "CONSTRUIRE BLOBS" << endl;
//       construireBlobs(sketch, blobs, ssblobs);
// 
//       cliques = construireCliques(ssblobs, blobs, lat, lon, matchingblobs);
//       
//       
// 
//       
//       cout << "FIN CONSTRUIRE BLOBS" << endl;
//       
//       // Extracting mesh patches before building the Aims graph
//       objects = new AimsSurfaceTriangle(getBlobsMeshes(ssblobs,mesh,nodes_lists));
// 
//     // End of mode 2
//     }
//     
// //     if (flatpath != ""){
// //       cerr << "écriture flat mesh:" << flatpath << endl;
// //       TimeTexture<float> texflat;
// //       AimsSurfaceTriangle flat(getFlatMap(nodes_lists,lat,lon,texflat));
// //       cerr << flat[0].vertex().size() << " !=! " << texflat[0].nItem() << endl;
// //       Writer<AimsSurfaceTriangle> wflat(flatpath);
// //       wflat.write(flat);
// 
// //     }
// //     
//     cerr << "construction graphe" << endl;
//     Graph graph("BlobsArg");
//     vector<float> resolution,bbmin2D,bbmax2D;
//     vector<int> bbmin, bbmax;
//     resolution.push_back(1.0); resolution.push_back(1.0); resolution.push_back(1.0);
//     bbmin.push_back(mesh[0].minimum()[0]-1); bbmin.push_back(mesh[0].minimum()[1]-1); bbmin.push_back(mesh[0].minimum()[2]-1);
//     bbmax.push_back(mesh[0].maximum()[0]+1); bbmax.push_back(mesh[0].maximum()[1]+1); bbmax.push_back(mesh[0].maximum()[2]+1);
//     graph.setProperty( "filename_base", "*");
// 
//     graph.setProperty("voxel_size", resolution);
//     graph.setProperty("boundingbox_min", bbmin);
//     graph.setProperty("boundingbox_max", bbmax);
//     graph.setProperty("mesh", meshPath);
//     graph.setProperty("sujet", sujet);
//     graph.setProperty("texture", texPath);
//     
//     
//     construireGraph(&graph, ssblobs, cliques, objects);
//     
// 
// 
// 
//     cerr << "graph.order:" << graph.order() << endl;
//     cerr << "graph.edgesSize:" << graph.edgesSize() << endl;
// 
//     Writer<Graph> graphWtr(outpath);
//     graphWtr.write(graph);
//     
//     return EXIT_SUCCESS;
//   }
//   catch( carto::user_interruption & )
//   {
//   }
//   catch( exception & e )
//   {
//     cerr << e.what() << endl;
//   }
// }


//##############################################################################

void setupData( map<string, TimeTexture<float> > &textures,
                map<string, AimsSurfaceTriangle > &meshes,
                map<string, TimeTexture<float> > &latitudes,
                map<string, TimeTexture<float> > &longitudes,
                string texPaths,
                string meshPaths,
                string latPaths,
                string lonPaths,
                vector<string> &listSujets){

  vector<string> paths;

  // Let's split the textures pathes and then read/store every texture
  paths = splitGraphFile(texPaths);
  assert( paths.size() == listSujets.size() );
  for ( uint i = 0 ; i < listSujets.size() ; i++ ) {
    Reader<TimeTexture<float> > texRdr ( paths[i] ) ;
    pair<string, TimeTexture<float> > tex;
    tex.first = listSujets[i];
    texRdr.read(tex.second);
    textures.insert(tex);
  }

  // Let's do the same for the meshes
  paths = splitGraphFile(meshPaths);
  assert( paths.size() == listSujets.size() );
  for ( uint i = 0 ; i < listSujets.size() ; i++ ) {
    Reader<AimsSurfaceTriangle> meshRdr (paths[i] ) ;
    pair<string, AimsSurfaceTriangle > mesh;
    mesh.first = listSujets[i];
    meshRdr.read(mesh.second);
    meshes.insert(mesh);
  }

  // Let's do it for the latitudes
  paths = splitGraphFile(latPaths);
  assert( paths.size() == listSujets.size() );
  for ( uint i = 0 ; i < listSujets.size() ; i++ ) {
    Reader<TimeTexture<float> > latRdr ( paths[i] );
    pair<string, TimeTexture<float> > lat;
    lat.first = listSujets[i];
    latRdr.read(lat.second);
    latitudes.insert(lat);
  }

  // Let's do it for the longitudes
  paths = splitGraphFile(lonPaths);
  assert( paths.size() == listSujets.size() );
  for ( uint i = 0 ; i < listSujets.size() ; i++ ) {
    Reader<TimeTexture<float> > lonRdr ( paths[i] );
    pair<string, TimeTexture<float> > lon;
    lon.first = listSujets[i];
    lonRdr.read(lon.second);
    longitudes.insert(lon);
  }
  
}

//##############################################################################

void ConstruireGraphe(Graph *graph,
                      vector<Blob *> &blobs,
                      vector<SSBlob *> &ssblobs,
                      vector<pair<Point2d, float> > &cliques,
                      vector<set<uint> > &matchingblobs,
                      string texPaths,
                      string meshPaths,
                      string latPaths,
                      string lonPaths,
                      vector<string> &listSujets){

  graph = new Graph("BlobsArg");

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

  Vertex *vert;
  carto::rc_ptr<AimsSurfaceTriangle> ptr;
  aims::GraphManip manip;
  vector<Vertex *> listVertices;
    
      
  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("ssb");
    vert->setProperty("index", i);
    vert->setProperty("name", i);
    vert->setProperty("label", "0");
    vert->setProperty("t", ssblobs[i]->t);
    vert->setProperty("rank", i);
    vert->setProperty( "subject", ssblobs[i]->subject);
    vert->setProperty( "tmin", ssblobs[i]->tmin);
    vert->setProperty( "tmax", ssblobs[i]->tmax);
    vert->setProperty( "tValue", 100.0);
    vert->setProperty("nodes_list", set2vector(ssblobs[i]->representation));
      
    
    // We associate the proper mesh patch from "objects" to the vertex
//     ptr=carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
//     (*ptr)[0]=(*objects)[i];
//     manip.storeAims(*graph, vert, "blob", ptr);
//     vert->setProperty("blob_label",i);
    
    listVertices.push_back(vert);
  }

  // CONSTRUIRE LES ARETES
  cout << "Construction cliques" << endl;
  for (uint i = 0 ; i < cliques.size() ; i++){
    
    // For every clique, we get the two corresponding vertices from listVertices
    Vertex *v1, *v2;
    v1 = listVertices[cliques[i].first[0]];
    v2 = listVertices[cliques[i].first[1]];
    int index1, index2;
    v1->getProperty("index", index1);
    v2->getProperty("index", index2);
    cout << index1 << " " << index2 << endl;
    
    
    Edge *edge= graph->addEdge(v1,v2,"b2b");
    edge->setProperty("blob_first", cliques[i].first[0]);
    edge->setProperty("blob_second", cliques[i].first[1]);
    edge->setProperty("similarity", cliques[i].second);

  }


  // AJOUTER LES GREY LEVEL BLOBS
  for (int i = 0 ; i < (int) blobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("glb");
    vert->setProperty("index", blobs[i]->index);
    vert->setProperty("t", blobs[i]->t);
    vert->setProperty( "scale", blobs[i]->scale);

  }

  // AJOUTER LES RELATIONS SSB -> GLB


}


//##############################################################################

int main( int argc, const char **argv ){
  try {

    int mode=0;
    string outPaths = "", outPath = "",
           meshPaths = "", meshPath = "",
           texPaths = "", texPath = "",
           latPaths = "", latPath = "",
           lonPaths = "", lonPath = "",
           flatPaths = "", flatPath = "",
           sujets = "";

    AimsApplication app( argc, argv, "surfLabelsTex2Graph" );
    app.addOption( meshPaths, "-m", "mesh");
    app.addOption( texPaths, "-t", "texture");
    app.addOption( outPaths, "-o", "output file");
    app.addOption( mode, "-M", "mode (0: normal - 1:barycenters (provide the lat/lon textures)",1);
    app.addOption( sujets, "-s", "sujet");
    app.addOption( latPaths, "--lat", "latitude");
    app.addOption( lonPaths, "--lon", "longitude");
    app.addOption( flatPaths, "--flat", "flat",1);
    app.initialize();

    map<string, TimeTexture<float> > textures;
    map<string, AimsSurfaceTriangle > meshes;
    map<string, TimeTexture<float> > latitudes;
    map<string, TimeTexture<float> > longitudes;
    vector<string> listSujets = splitGraphFile(sujets);
    cerr << "split string sujets -> " << listSujets.size() << " sujets" << endl;
    
    setupData(textures, meshes, latitudes, longitudes, texPaths, meshPaths, latPaths, lonPaths, listSujets);

    vector<Blob *> blobs;
    vector<SSBlob *> ssblobs;
    vector<pair<Point2d, float> > cliques;
    vector<set<uint> > matchingblobs;

    
    for (uint i = 0 ; i < listSujets.size() ; i++) {
      string sujet = listSujets[i];
      cout << "sujet : " << sujet << endl;;
      // Construction of a primal-sketch
      ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss(
          getScaleSpace( textures[sujet], meshes[sujet], latitudes[sujet], longitudes[sujet]) );

      PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch(sujet, &ss, SURFACE);

      // Launching the computation of the PS (tmin, tmax, statfile, intersection_criterium)
      sketch.ComputePrimalSketch(1.0, 8.0, "", 10);

      cout << "CONSTRUIRE BLOBS" << endl;
      construireBlobs(sketch, blobs, ssblobs);      

      cout << "blobs.size() = " << blobs.size() << endl;
      cout << "ssblobs.size() = " << ssblobs.size() << endl;
    }

    cliques = construireCliques(ssblobs, blobs, latitudes, longitudes, matchingblobs);

    Graph *graph;
    ConstruireGraphe(graph, blobs, ssblobs, cliques, matchingblobs, texPaths, meshPaths, latPaths, lonPaths, listSujets);

    
    Writer<Graph> wtrGraph("/volatile/operto/graph.arg");
    wtrGraph.write(*graph);
    // On vient d'écrire le graphe pour l'analyse (à créer : une fonction qui lit
    //   ce graphe et qui récupère les sites et cliques)
        
    // Eventuellement une fonction pour séparer le graphe en sous-graphes correspondants
    //   aux différents sujets
     
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

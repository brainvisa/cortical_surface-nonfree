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
#include "blobs.h"




using namespace aims;
using namespace carto;
using namespace std;


//##############################################################################

// Function that builds a collection of Blob and SSBlob objects from a previously
//   computed Primal Sketch
void construireBlobs(PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch, 
                     vector<Blob *> &blobs, vector<SSBlob *> &ssblobs, bool initNull = true){
  
//     // Inititalization of the results vectors "blobs" and "ssblobs"
    if (initNull){
      blobs.clear();
      ssblobs.clear();
    }
    
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*> listBlobs 
         = sketch.BlobSet();
    uint iBlob=blobs.size(), iSSBlob=ssblobs.size();
    
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*>::iterator itSSB;
    list<GreyLevelBlob<SiteType<AimsSurface<3, Void> >::type > *>::iterator itGLB;
    
    ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *ssb;
    set< SiteType<AimsSurface<3, Void> >::type, 
      ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> >::iterator itPoints;
    
    
    for (itSSB = listBlobs.begin() ; itSSB != listBlobs.end() ; itSSB++){
      
      // For each scale-space blob, we create a SSBlob in "ssblobs" containing 
      //    various Blob objects (being themselves contained in a general resulting 
      // "blobs" vector).
      ssb = *itSSB;
      ssblobs.push_back(new SSBlob());
      SSBlob *ssblob = ssblobs[ssblobs.size() - 1];
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
          (blob->nodes_set).insert((*itPoints).second);
        
        ssblob->blobs.insert(blob);
        
        if (blob->scale < ssblob->tmin)
          ssblob->tmin=blob->scale;
        if (blob->scale > ssblob->tmax)
          ssblob->tmax=blob->scale;
      }
      
      ssblob->t = ssb->GetMeasurements().t;
      
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
    
      if ( blobs[i]->nodes_set.find(p1)!=blobs[i]->nodes_set.end() &&
              blobs[i]->nodes_set.find(p2)!=blobs[i]->nodes_set.end() &&
              blobs[i]->nodes_set.find(p3)!=blobs[i]->nodes_set.end() )
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
vector<SSBClique> construireCliques ( vector<SSBlob *>   &ssblobs,
                                                  vector<Blob *>     &blobs, 
                                                  map<string, SubjectData> &data,
                                                  vector<vector<Blob *> > &matchingblobs){
    
      vector<SSBClique > cliques;
      matchingblobs = vector<vector<Blob *> > (ssblobs.size());
      
      set<Blob *>::iterator itB1, itB2;
      Blob *b1max, *b2max;
      
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

                pair<Point2df,Point2df> bbi = getBoundingBox((*itB1)->nodes_list, data[ssblobs[i]->subject].lat, data[ssblobs[i]->subject].lon),
                                        bbj = getBoundingBox((*itB2)->nodes_list, data[ssblobs[j]->subject].lat, data[ssblobs[j]->subject].lon);

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

              
              
              cliques.push_back(SSBClique(ssblobs[i], ssblobs[j], overmax));
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
          set<int> blobNodes(matchingblobs[i][j]->nodes_set);
          ssblobs[i]->representation.insert(blobNodes.begin(), blobNodes.end());
          cout << ssblobs[i]->representation.size() << " " << flush;
        }
        
        if (matchingblobs[i].size()!=0) 
          cout << endl ;
        
      }
//       float test;
//       cin >> test;
      return cliques;
}




//##############################################################################

// int main( int argc, const char **argv ){
//     // Computes a radius for the "barycenters" sphere-based representation mode
//     mesh[0].setMini(); mesh[0].setMaxi();
//     cerr << mesh[0].minimum()[0] << " " << mesh[0].minimum()[1] << " " << mesh[0].minimum()[2] << ";" << mesh[0].maximum()[0] << " " << mesh[0].maximum()[1] << " " << mesh[0].maximum()[2] << endl;    
//     float dist = sqrt(pow(mesh[0].minimum()[0]-mesh[0].maximum()[0],2)+pow(mesh[0].minimum()[1]-mesh[0].maximum()[1],2)+pow(mesh[0].minimum()[2]-mesh[0].maximum()[2],2));
//     float radius = dist / 300.0;
//     
// 
//     if (mode == 1){
//       
//       // Barycenters mode
//       objects = new AimsSurfaceTriangle(getBarycenters(mesh,nodes_lists,radius));
//     }

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

// }


//##############################################################################

void setupData( map<string, SubjectData> &data,
                vector<string> &listTexPaths,
                vector<string> &listMeshPaths,
                vector<string> &listLatPaths,
                vector<string> &listLonPaths,
                vector<string> &listSujets){
  
  assert( listLatPaths.size() == listSujets.size() );
  assert( listLonPaths.size() == listSujets.size() );
  assert( listMeshPaths.size() == listSujets.size() );
  assert( listTexPaths.size() == listSujets.size() );
  
  for ( uint i = 0 ; i < listSujets.size() ; i++ ) {
    
    pair<string, SubjectData > pSubjData;
//     SubjectData subjData;
    pSubjData.first = listSujets[i];
//     pSubjData.second = subjData;
    pSubjData.second.subject = listSujets[i];
    
    Reader<TimeTexture<float> > texRdr ( listTexPaths[i] ) ;
    texRdr.read(pSubjData.second.tex);
    
    Reader<AimsSurfaceTriangle> meshRdr (listMeshPaths[i] ) ;
    meshRdr.read(pSubjData.second.mesh);
    
    Reader<TimeTexture<float> > latRdr ( listLatPaths[i] );
    latRdr.read(pSubjData.second.lat);
  
    Reader<TimeTexture<float> > lonRdr ( listLonPaths[i] );
    lonRdr.read(pSubjData.second.lon);
    
    // Checking the data
    cout << " subject : " << pSubjData.second.subject << endl;
    cout << "  texture : " << pSubjData.second.tex[0].nItem() << " values" << endl;
    cout << "  mesh : " << pSubjData.second.mesh[0].vertex().size() << " nodes" << endl;
    cout << "  lat : " << pSubjData.second.lat[0].nItem() << " values" << endl;
    cout << "  lon : " << pSubjData.second.lon[0].nItem() << " values" << endl;
    
    data.insert(pSubjData);
    
  }
  
}

//##############################################################################

void ConstruireGraphe(Graph *graph,
//                       vector<Blob *> &blobs,
                      vector<SSBlob *> &ssblobs,
                      vector<SSBClique> &cliques,
                      string texPaths,
                      string meshPaths,
                      string latPaths,
                      string lonPaths,
                      string indivGraphPaths,
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

  Vertex *vert;
  carto::rc_ptr<AimsSurfaceTriangle> ptr;
//   aims::GraphManip manip;
  vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( ssblobs.size() );
  
      
  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("ssb");
    int index1=ssblobs[i]->index;
    vert->setProperty("index", index1);
    vert->setProperty("label", "0");
    vert->setProperty("t", ssblobs[i]->t);
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
    v1 = listVertSSB[aux1];
    v2 = listVertSSB[aux2];
    
    Edge *edge= graph->addEdge(v1, v2, "b2b");


//     edge->setProperty("blob_first", cliques[i].first[0]);
//     edge->setProperty("blob_second", cliques[i].first[1]);
    edge->setProperty("similarity", cliques[i].similarity);

  }


//   // AJOUTER LES GREY LEVEL BLOBS
//   for (int i = 0 ; i < (int) blobs.size() ; i++) {
//         
//     // For every scale-space blob, we create a vertex in the Aims graph : we define
//     //   its properties and store a link between the created vertex and the blob index
//     
//     cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
//     vert = graph->addVertex("glb");
//     vert->setProperty("index", blobs[i]->index);
//     vert->setProperty("t", blobs[i]->t);
//     vert->setProperty( "scale", blobs[i]->scale);
//     vert->setProperty( "nodes_list", blobs[i]->nodes_list);
//     listVertGLB.push_back(vert);
//     
//   }
// 
//   // AJOUTER LES RELATIONS SSB -> GLB
//   
//   for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
//     
//     set<Blob *>::iterator itB1;
//     set<Blob *> &listGLB = ssblobs[i]->blobs;
//     for (itB1 = listGLB.begin(); itB1 != listGLB.end() ; itB1++) {
//       
//       Vertex *v1, *v2;
//       
//       v1 = listVertSSB[ssblobs[i]->index];
//       v2 = listVertGLB[(*itB1)->index];
//       Edge *edge= graph->addEdge(v1,v2,"s2g");
// //       edge->setProperty("ssb_index", ssblobs[i]->index);
// //       edge->setProperty("glb_index", (*itB1)->index);
//       
//     }
//   }



}

//##############################################################################

// Creates an Aims Graph for only ONE subject with scale-space blobs, grey-level 
//  blobs and links between both types

void ConstruireIndividualGraph( Graph *graph,
                                vector<Blob *> &blobs,
                                vector<SSBlob *> &ssblobs,
                                string meshPath,
                                string texPath,
                                string latPath,
                                string lonPath,
                                string sujet){
  // From the two vectors, we build a graph containing the ssb, the glb and
  //  the links between ssb and glb                                 

  cerr << "Graph construction..." << endl;
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
  graph->setProperty("mesh", meshPath);
  graph->setProperty("sujet", sujet);
  graph->setProperty("texture", texPath);
  graph->setProperty("latitude", latPath);
  graph->setProperty("longitude", lonPath);

  Vertex *vert;
  carto::rc_ptr<AimsSurfaceTriangle> ptr;
//   aims::GraphManip manip;
  vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( blobs.size() );
    
      
  // Let's add the scale-space blobs
  
  cout << "Adding scale-space blobs..." << endl;
  
  int iNbSSB = 0;
  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("ssb");
    iNbSSB++;
    
    vert->setProperty("index", ssblobs[i]->index);
    vert->setProperty("label", "0");
    vert->setProperty("t", ssblobs[i]->t);
    vert->setProperty( "subject", ssblobs[i]->subject);
    vert->setProperty( "tmin", ssblobs[i]->tmin);
    vert->setProperty( "tmax", ssblobs[i]->tmax);
    vert->setProperty( "tValue", 100.0);

      
     
    // We associate the proper mesh patch from "objects" to the vertex
    // ptr=carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
    // (*ptr)[0]=(*objects)[i];
    // manip.storeAims(*graph, vert, "blob", ptr);
    // vert->setProperty("blob_label",i);
        
    assert(ssblobs[i]->index < ssblobs.size());
    listVertSSB[ ssblobs[i]->index ] = vert;
  }
  cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbSSB << " blobs added... done" << endl; 
    
  // Let's add the grey-level blobs
  
  cout << "Adding grey-level blobs..." << endl;
  
  int iNbGLB = 0;
  for (int i = 0 ; i < (int) blobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("glb");
    iNbGLB++;
    vert->setProperty("index", blobs[i]->index);
    vert->setProperty("t", blobs[i]->t);
    vert->setProperty( "scale", blobs[i]->scale);
    vert->setProperty( "nodes_list", blobs[i]->nodes_list);

    assert( blobs[i]->index < blobs.size() );
    listVertGLB[ blobs[i]->index ] = vert;
    
  }
  cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbGLB << " blobs added... done" << endl; 
   
  // Let's add the links between scale-space and grey-level blobs
  
  cout << "Adding links between both types..." << endl;
  
  uint iNbLinks=0;
  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
    
    set<Blob *>::iterator itB1;
    set<Blob *> &listGLB = ssblobs[i]->blobs;
    for (itB1 = listGLB.begin(); itB1 != listGLB.end() ; itB1++) {
      
      Vertex *v1, *v2;
      
      v1 = listVertSSB[ssblobs[i]->index];
      v2 = listVertGLB[(*itB1)->index];
      graph->addEdge(v1,v2,"s2g");
      iNbLinks++;
//       edge->setProperty("ssb_index", ssblobs[i]->index);
//       edge->setProperty("glb_index", (*itB1)->index);
      
    }
  }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << endl; 

}

//##############################################################################

void RecoverBlobsFromIndivGraph( Graph *graph,
                            SubjectData &subjData,
                            vector<Blob *> &blobs,
                            vector<SSBlob *> &ssblobs,
                            bool initNull = true){
    if (initNull){
      blobs.clear();
      ssblobs.clear();                               
    }
                             
    set<Vertex *>::iterator iv;
    Edge *e;
    Vertex::iterator jv;
    Edge::iterator kv;
    string meshPath, sujet, texPath, latPath, lonPath;
    
    graph->getProperty("mesh", meshPath);
    graph->getProperty("sujet", sujet);
    graph->getProperty("texture", texPath);
    graph->getProperty("latitude", latPath);
    graph->getProperty("longitude", lonPath);
    
    subjData.subject = sujet;
    
    cout << "Subject " << sujet << endl;
    cout << " mesh : " << meshPath << endl;
    cout << " tex : " << texPath << endl;
    cout << " lat : " << lonPath << endl << endl;
    
    cout << "Reading files..." <<  flush;
    Reader<AimsSurfaceTriangle> rdrMesh(meshPath);
    rdrMesh.read(subjData.mesh);
    
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
    for (iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv){
      if ((*iv)->getSyntax() == "glb"){
        int index;
        float scale, t;
        vector<int> nodes_list;
        (*iv)->getProperty("index", index);
        (*iv)->getProperty("scale", scale);
        (*iv)->getProperty("nodes_list", nodes_list);
        (*iv)->getProperty("t", t);
        blobs.push_back(new Blob());
        Blob *blob = blobs[blobs.size()-1];
        blob->index = index;
        blob->nodes_list = nodes_list;
        blob->nodes_set = vector2set(nodes_list);
        blob->scale = scale;
        blob->t = t;
        
      }
      else if ((*iv)->getSyntax() == "ssb"){
        int index;
        float tmax, tmin, t;
        (*iv)->getProperty("index", index);
        (*iv)->getProperty("tmax", tmax);
        (*iv)->getProperty("tmin", tmin);
        (*iv)->getProperty("t", t);
        ssblobs.push_back(new SSBlob());
        SSBlob *ssblob = ssblobs[ssblobs.size()-1];
        ssblob->index = index;
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
    

      
            
    
    

    
    cout << blobs.size() << " blobs added" << endl;
    cout << ssblobs.size() << " ssblobs added " << endl;
    cout << iNbLinks << " links added" << endl;



}

//##############################################################################

// void RecoverBlobsFromGroupGraph( Graph *graph,
//                             vector<SubjectData> &data,
//                             vector<Blob *> &blobs,
//                             vector<SSBlob *> &ssblobs,
//                             bool initNull = true){
//     if (initNull){
//       blobs.clear();
//       ssblobs.clear();                               
//     }
//                              
//     set<Vertex *>::iterator iv;
//     Edge *e;
//     Vertex::iterator jv;
//     Edge::iterator kv;
//     vector<string> listMeshPaths, listSujets, listTexPaths, listLatPaths, listLonPaths;
//     
//     graph->getProperty("mesh", listMeshPaths);
//     graph->getProperty("sujet", listSujets);
//     graph->getProperty("texture", listTexPaths);
//     graph->getProperty("latitude", listLatPaths);
//     graph->getProperty("longitude", listLonPaths);
// 
//     for (uint i = 0 ; i < listSujets.size() ; i++){
//       data.push_back(SubjectData());
//       SubjectData &subjData = data[data.size() - 1];
//       subjData.subject = listSujets[i];
//       
//       cout << "Subject " << listSujets[i] << endl;
//       cout << " mesh : " << listMeshPaths[i] << endl;
//       cout << " tex : " << listTexPaths[i] << endl;
//       cout << " lat : " << listLonPaths[i] << endl << endl;
//       
//       cout << "Reading files..." <<  flush;
//       Reader<AimsSurfaceTriangle> rdrMesh(listMeshPaths[i]);
//       rdrMesh.read(subjData.mesh);
//       
//       Reader<TimeTexture<float> > rdrTex(listTexPaths[i]);
//       rdrTex.read(subjData.tex);
//       
//       Reader<TimeTexture<float> > rdrLat(listLatPaths[i]);
//       rdrLat.read(subjData.lat);
//       
//       Reader<TimeTexture<float> > rdrLon(listLonPaths[i]);
//       rdrLon.read(subjData.lon);
//       
//       
//       cout << " done" << endl;
//     }
// //     vector<Vertex *> listVertSSB, listVertGLB;
// //     map<int, set<int> > listGLBindices;
// //     int iNbLinks = 0;
// //     for (iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv){
// //       if ((*iv)->getSyntax() == "glb"){
// //         int index;
// //         float scale, t;
// //         vector<int> nodes_list;
// //         (*iv)->getProperty("index", index);
// //         (*iv)->getProperty("scale", scale);
// //         (*iv)->getProperty("nodes_list", nodes_list);
// //         (*iv)->getProperty("t", t);
// //         blobs.push_back(new Blob());
// //         Blob *blob = blobs[blobs.size()-1];
// //         blob->index = index;
// //         blob->nodes_list = nodes_list;
// //         blob->nodes_set = vector2set(nodes_list);
// //         blob->scale = scale;
// //         blob->t = t;
// //         
// //       }
// //       else if ((*iv)->getSyntax() == "ssb"){
// //         int index;
// //         float tmax, tmin, t;
// //         (*iv)->getProperty("index", index);
// //         (*iv)->getProperty("tmax", tmax);
// //         (*iv)->getProperty("tmin", tmin);
// //         (*iv)->getProperty("t", t);
// //         ssblobs.push_back(new SSBlob());
// //         SSBlob *ssblob = ssblobs[ssblobs.size()-1];
// //         ssblob->index = index;
// //         ssblob->tmax = tmax;
// //         ssblob->tmin = tmin;
// //         ssblob->t = t;                
// //         ssblob->subject = sujet;
// //         if (listGLBindices.find(index) == listGLBindices.end())
// //           listGLBindices[index] = set<int>();
// //             
// //         for (jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++){
// //           e = *jv;
// //           if (e->getSyntax() == "s2g"){
// //             for (kv = e->begin() ; kv != e->end() ; kv++){
// //               if ((*kv)->getSyntax() == "ssb"){
// // 
// //               }
// //               else if ((*kv)->getSyntax() == "glb"){
// //                 int iGLBindex;
// //                 (*kv)->getProperty("index", iGLBindex);
// //                 listGLBindices[index].insert(iGLBindex);
// //                 iNbLinks++;
// //               }
// //             }
// //           }
// //         }
// // 
// //       }
// //     }
// //     cout << "Rebuilding links..." << endl;
// //     for (uint i = 0 ; i < ssblobs.size() ; i++){
// //       int index = ssblobs[i]->index;
// //       set<int>::iterator it;
// //       for (it = listGLBindices[index].begin() ; it != listGLBindices[index].end() ; it++){
// //         ssblobs[i]->blobs.insert(blobs[*it]);
// //         
// //       }
// //     }
// //     
// // 
// //       
// //             
// //     
// //     
// // 
// //     
// //     cout << blobs.size() << " blobs added" << endl;
// //     cout << ssblobs.size() << " ssblobs added " << endl;
// //     cout << iNbLinks << " links added" << endl;
// 
// 
// 
// }


//##############################################################################

// void getSubjectBlobs( vector<Blob *> &blobs,
//                       vector<SSBlob *> &ssblobs,
//                       string sujet,
//                       vector<Blob *> &subjBlobs,
//                       vector<SSBlob *> &subjSsblobs){
//   
// 
// 
// }

vector<string> getVectorStringFromGraph(Graph &graph, string graph_property){
  vector<string> v;
  if( graph.hasProperty( graph_property ) )  {
    Object slist = graph.getProperty( graph_property ); // note the different getProperty() method
    cout << "node with 'sujets' property:\n";
    Object oit = slist->objectIterator();  // iterator on the list
    while( oit->isValid() ) {
      Object s = oit->currentValue(); // the list element, type Object
      string ss = s->getString(); // extract as std::string or convert to string      
      cout << ss << ", ";
      v.push_back(ss);
      oit->next();
    }
    cout << endl;    
  }     
  
  return v;
  
}

//##############################################################################

void recoverGroupData ( Graph &graph,
                        map<string, SubjectData> &data){
  vector<string> listSujets, listGraphPaths, listMeshPaths, listTexPaths, listLatPaths, listLonPaths;;
  listSujets = getVectorStringFromGraph(graph, "sujets");
  listGraphPaths = getVectorStringFromGraph(graph, "indiv_graphs");
  listTexPaths = getVectorStringFromGraph(graph, "textures");
  listLatPaths = getVectorStringFromGraph(graph, "latitudes");
  listLonPaths = getVectorStringFromGraph(graph, "longitudes");
  listMeshPaths = getVectorStringFromGraph(graph, "meshes");
    
  
  for (uint i = 0 ; i < listSujets.size() ; i++){
    pair<string, SubjectData> pSubjData;
    pSubjData.first = listSujets[i];
    Reader<TimeTexture<float> > texRdr ( listTexPaths[i] ) ;
    texRdr.read(pSubjData.second.tex);
    
    Reader<AimsSurfaceTriangle> meshRdr (listMeshPaths[i] ) ;
    meshRdr.read(pSubjData.second.mesh);
    
    Reader<TimeTexture<float> > latRdr ( listLatPaths[i] );
    latRdr.read(pSubjData.second.lat);
  
    Reader<TimeTexture<float> > lonRdr ( listLonPaths[i] );
    lonRdr.read(pSubjData.second.lon);
    
    // Checking the data
    cout << " subject : " << pSubjData.second.subject << endl;
    cout << "  texture : " << pSubjData.second.tex[0].nItem() << " values" << endl;
    cout << "  mesh : " << pSubjData.second.mesh[0].vertex().size() << " nodes" << endl;
    cout << "  lat : " << pSubjData.second.lat[0].nItem() << " values" << endl;
    cout << "  lon : " << pSubjData.second.lon[0].nItem() << " values" << endl;
    
    data.insert(pSubjData); 
  }                          
                          
}

//##############################################################################

void readGroupGraph ( Graph &graph,
                      map<string, SubjectData> &data,
                      vector<SSBlob *> &ssblobs, 
                      vector<SSBClique> &cliques,
                      vector<Vertex *> &listVertex){

  cout << "Recovering the data..." << endl;
  recoverGroupData(graph, data);
    
  set<Vertex *>::iterator iv;
  // Recovering the scale-space blobs
  string sujet;
  int index;
  int newindex=0;
  cout << " Recovering the scale-space blobs..." << endl;
  for (iv = graph.vertices().begin() ; iv != graph.vertices().end(); ++iv){
    string label;
    
    float t, tmin, tmax, tvalue;
    vector<int> representation;

    if ((*iv)->getSyntax() == "ssb"){
      (*iv)->getProperty("index", index);
      (*iv)->getProperty( "subject", sujet);
      (*iv)->getProperty("label", label);      
      (*iv)->getProperty( "tmin", tmin);
      (*iv)->getProperty( "tmax", tmax);
      (*iv)->getProperty( "t", t);
//       (*iv)->getProperty( "tValue", tvalue);
//       (*iv)->getProperty( "rank", rank);
//       (*iv)->setProperty("label", label);
      (*iv)->getProperty( "nodes_list", representation);      

      ssblobs.push_back(new SSBlob());
      SSBlob *s=ssblobs[ssblobs.size()-1];
      listVertex.push_back(*iv);


//       (*iv)->getProperty("label", label);

      s->index = newindex++;
      s->label = atoi(label.data());    

      s->graph_index = index;
      s->subject = sujet;
      s->tmin = tmin;
      s->tmax = tmax;
      s->t = t;
      for (uint i=0;i<representation.size();i++)
        s->representation.insert(representation[i]);      
      (*iv)->setProperty( "sites_index", (int)(ssblobs.size()-1));
      index =0;
      (*iv)->getProperty( "sites_index", index);
//       cout << "index:" << index << " " << flush;
    }

  }


  cout << "ssblobs.size() :"<< ssblobs.size() << endl;

  
  // Recovering the links...
  Edge *e;
  Vertex::iterator jv;
  Edge::iterator kv;
  
  cout << " Recovering the similarity cliques..." << endl;
  for (iv = graph.vertices().begin() ; iv!=graph.vertices().end(); ++iv){
    if ((*iv)->getSyntax() == "ssb"){
//         (*iv)->getProperty( "sites_index", index);
//         cout << "idx:" << index << " " << flush;
        (*iv)->getProperty("index", index);
        (*iv)->getProperty( "subject", sujet);
        for (jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++){
        e = *jv;
        if (e->getSyntax() == "b2b"){
          float similarity;
          e->getProperty("similarity", similarity);
          for (kv = e->begin() ; kv != e->end() ; kv++){
            if ((*kv)->getSyntax() == "ssb"){
              int indexB2;
              string sujetB2;
              (*kv)->getProperty("index", indexB2);
              (*kv)->getProperty("subject", sujetB2);
              if (!(indexB2 == index && sujetB2 == sujet)){
                int blobs_index1, blobs_index2;
                (*iv)->getProperty("sites_index", blobs_index1);
                (*kv)->getProperty("sites_index", blobs_index2);
//                 cout << cliques.size() << " " << ssblobs.size() << " " << blobs_index1 << " " << blobs_index2 << flush;
                if (blobs_index1 < blobs_index2){
                  cliques.push_back(SSBClique(ssblobs[blobs_index1], ssblobs[blobs_index2], similarity));
                }
                
              }

            }
          }
        }
      }
    }
  }

  cout << "ssbcliques.size() :"<< cliques.size() << endl << endl;


}

//##############################################################################

void convertSSBlobsToSites(vector<SSBlob *> &ssblobs, vector<Site *> &sites){
//    int newindex= 0;
  for (uint i = 0 ; i < ssblobs.size() ; i++){
//     if (ssblobs[i]->t>2.0){
     sites.push_back(new Site());
     Site *s = sites[sites.size() - 1];
     s->index = ssblobs[i]->index;
     s->graph_index = ssblobs[i]->index;
     s->subject = ssblobs[i]->subject;
     s->label = ssblobs[i]->label;
     s->tmin = ssblobs[i]->tmin;
     s->tmax = ssblobs[i]->tmax;
     s->t = ssblobs[i]->t;
     s->nodes_list = ssblobs[i]->representation;
//     }
  }    
    
}

//##############################################################################

void getCliquesFromSSBCliques ( vector<SSBClique> &ssbcliques, 
                                vector<Site *> &sites,
                                vector<Clique> &cliques,
                                vector<vector<int> > &cliquesDuSite){
  
 cliquesDuSite = vector<vector<int> >(sites.size()); 
 
 for (uint i = 0 ; i < ssbcliques.size() ; i++){
   Clique simc;
   simc.type = SIMILARITY;
   simc.rec = ssbcliques[i].similarity;
   SSBlob *ssb1, *ssb2; 
   int iSSB1, iSSB2;   
   ssb1 = ssbcliques[i].ssb1;
   ssb2 = ssbcliques[i].ssb2;
   iSSB1 = ssb1->index;
   iSSB2 = ssb2->index;

   cliquesDuSite[sites[iSSB1]->index].push_back(i);
   cliquesDuSite[sites[iSSB2]->index].push_back(i);
   
   simc.blobs.push_back(sites[iSSB1]);
   simc.blobs.push_back(sites[iSSB2]);
   cliques.push_back(simc);
   
 }
  for (uint i=0; i<sites.size();i++){
  
    for (uint n=0;n<cliquesDuSite[i].size();n++){
        uint aux = cliquesDuSite[i][n];
        if (cliques[aux].type == SIMILARITY){
          if ( cliques[aux].blobs[0]->index == (uint) i ) {}
            else if (cliques[aux].blobs[1]->index ==(uint) i ) {}
            else {
                cout << i << " " << aux << " " << cliques[aux].type << " " << cliques[aux].blobs.size() << " " << cliques[aux].blobs[0]->index << " " << cliques[aux].blobs[1]->index << endl;
                ASSERT(false);
              }          
        }
    }
  
  }
  cout << "ok" << endl;
  
}

//##############################################################################

void computeSitesBoundingBoxes(vector<Site *> &sites, map<string, SubjectData> &data){

  for (uint i = 0 ; i < sites.size() ; i++) {
    vector<int> vTemp(set2vector(sites[i]->nodes_list));
    pair<Point2df, Point2df> bb = getBoundingBox(vTemp, data[sites[i]->subject].lat, data[sites[i]->subject].lon);
    sites[i]->boundingbox_min = Point3df(bb.first[0], bb.first[1], 0);
    sites[i]->boundingbox_max = Point3df(bb.second[0], bb.second[1], 0);
  }

}

//##############################################################################

void FromRawTexturesToIndividualGraphsViaPrimalSketches ( string sujets,
                                                          string indivGraphPaths,
                                                          string meshPaths,
                                                          string texPaths,
                                                          string latPaths,
                                                          string lonPaths){
  
  map<string, SubjectData> data;
  vector<string> listSujets = splitGraphFile(sujets);
  vector<string> listGraphPaths = splitGraphFile(indivGraphPaths);
  vector<string> listMeshPaths = splitGraphFile(meshPaths);
  vector<string> listTexPaths = splitGraphFile(texPaths);
  vector<string> listLatPaths = splitGraphFile(latPaths);
  vector<string> listLonPaths = splitGraphFile(lonPaths);
      
  cerr << "  split string sujets -> " << listSujets.size() << " sujets" << endl << endl;

  cout << "Reading the data..." << endl;
  setupData(data, listTexPaths, listMeshPaths, listLatPaths, listLonPaths, listSujets);
  cout << "done" << endl<< endl;
  
  
  vector<Blob *> blobs;
  vector<SSBlob *> ssblobs;

  
  // Processing every subject...
  for (uint i = 0 ; i < listSujets.size() ; i++) {
    
    string sujet = listSujets[i];
    
    // Checking the data
    cout << " subject : " << sujet << endl;
    cout << "  texture : " << data[sujet].tex[0].nItem() << " values" << endl;
    cout << "  mesh : " << data[sujet].mesh[0].vertex().size() << " nodes" << endl;
    cout << "  lat : " << data[sujet].lat[0].nItem() << " values" << endl;
    cout << "  lon : " << data[sujet].lon[0].nItem() << " values" << endl;        
    
    // Generating a scale-space
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss(
        getScaleSpace( data[sujet].tex, data[sujet].mesh, data[sujet].lat, data[sujet].lon) );

    // Constructing a primal-sketch        
    PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch(sujet, &ss, SURFACE);

    // Launching the computation of the PS (tmin, tmax, statfile, intersection_criterium)
    sketch.ComputePrimalSketch(1.0, 8.0, "", 10);

    // Getting the blobs from PS structure
    cout << "Blobs vectors construction..." << endl;
    construireBlobs(sketch, blobs, ssblobs);      
    cout << "blobs.size() = " << blobs.size() << endl;
    cout << "ssblobs.size() = " << ssblobs.size() << endl;
    
    // Converting the blobs into an Aims Graph
    Graph *tmpGraph = new Graph("BlobsArg");
    ConstruireIndividualGraph(tmpGraph, blobs, ssblobs, listMeshPaths[i], listTexPaths[i], listLatPaths[i], listLonPaths[i], sujet);
    
    // Storing the graph on the hard disk
    cout << "Writing graph .. " << listGraphPaths[i] << endl;
    Writer<Graph> wtrGraph(listGraphPaths[i]);
    wtrGraph.write(*tmpGraph); 
    cout << "done" << endl << endl;
    
  // Next subject...
  }
  
  // Deleting the blobs vectors
  for (uint i = 0 ; i < blobs.size() ; i++)
    delete(blobs[i]);
  for (uint i = 0 ; i < ssblobs.size() ; i++)
    delete(ssblobs[i]);
    
}


//##############################################################################

int main( int argc, const char **argv ){
  try {

    int mode=0;
    string indivGraphPaths = "",
           groupGraphPath = "",
           meshPaths = "",
           texPaths = "",
           latPaths = "",
           lonPaths = "",
           flatPaths = "",
           sujets = "";

    AimsApplication app( argc, argv, "surfLabelsTex2Graph" );
    app.addOption( meshPaths, "-m", "mesh");
    app.addOption( texPaths, "-t", "texture");
    app.addOption( indivGraphPaths, "-g", "indiv graphs");
    app.addOption( groupGraphPath, "-G", "group graph");
    app.addOption( mode, "-M", "mode : 0 (compute the primal sketches) - 1 (take previously computed primal sketches)",1);
    app.addOption( sujets, "-s", "sujet");
    app.addOption( latPaths, "--lat", "latitude");
    app.addOption( lonPaths, "--lon", "longitude");
    app.addOption( flatPaths, "--flat", "flat",1);
    app.initialize();


    if (mode == 0){
        FromRawTexturesToIndividualGraphsViaPrimalSketches
            ( sujets, indivGraphPaths, meshPaths, texPaths, latPaths, lonPaths);
      
    }
    else if (mode == 1){

        vector<string> listSujets = splitGraphFile(sujets);
        vector<string> listGraphPaths = splitGraphFile(indivGraphPaths);
        map<string, SubjectData> data;
        vector<Blob *> blobs;
        vector<SSBlob *> ssblobs;
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
        vector<SSBClique> cliques;
        vector<vector<Blob *> > matchingblobs;
        cliques = construireCliques(ssblobs, blobs, data, matchingblobs);
        cout << endl << cliques.size() << " cliques de similarité " << endl;
    
        // Building the group graph...
        Graph *graph = new Graph("BlobsArg");
        ConstruireGraphe(graph, /*blobs,*/ ssblobs, cliques, texPaths, meshPaths, latPaths, lonPaths, indivGraphPaths, listSujets);
        
        Writer<Graph> wtrGraph(groupGraphPath);
        wtrGraph.write(*graph);
    }
    else if (mode == 2){
        cout << "Reading the group graph..." << endl;
        Graph graph; 
        Reader<Graph> rdrGroupGraph(groupGraphPath);
        rdrGroupGraph.read(graph);
        cout << "Recovering the sites and cliques..." << endl;
        // Getting the sites and cliques
        vector<SSBlob *> ssblobs;
        vector<SSBClique> ssbcliques;
        vector<Vertex *> listVertex;
        map<string, SubjectData> data;
        
        readGroupGraph(graph, data, ssblobs, ssbcliques, listVertex);
        
        // faire l'analyse = étiquetter
        Anneal swc;
//         vector<Site *> sites;
//         vector<Clique> cliques;
//         vector<vector<int> > cliquesDuSite;
        
        // Sauvegarder les labels dans les graphes individuels
        //  sur le modèle de comment on faisait avant
        
        convertSSBlobsToSites(ssblobs, swc.sites);
        cout << swc.sites.size() << " sites generated" << endl;
        
        getCliquesFromSSBCliques(ssbcliques, swc.sites, swc.cliques, swc.cliquesDuSite);
        
        ConstruireCliquesDataDriven(swc.sites, swc.cliquesDuSite, swc.cliques);        
        ConstruireCliquesIntraPS(swc.sites, swc.cliquesDuSite, swc.cliques);
        cout << swc.cliques.size() << " cliques created" << endl;

        
        set<string> subjects;

        cout << endl << "  done" << endl;
        for (uint i=0;i< swc.sites.size();i++)
          subjects.insert(swc.sites[i]->subject);
        swc.nbsujets = subjects.size();
        
        uint nb_cl_sim=0, nb_cl_dd=0, nb_cl_intraps=0, nb_cl_lower=0;
        for (uint i=0;i<swc.cliques.size();i++){
          if (swc.cliques[i].type == SIMILARITY) nb_cl_sim++;
          else if (swc.cliques[i].type == DATADRIVEN) nb_cl_dd++;
          else if (swc.cliques[i].type == BESTLOWERSCALE) nb_cl_lower++;
          else if (swc.cliques[i].type == INTRAPRIMALSKETCH) nb_cl_intraps++;
        }
        cout << " done (" << nb_cl_sim << " cliques de similarité ; " << nb_cl_dd << " cliques datadriven ; " << nb_cl_lower << " cliques lower ; " << nb_cl_intraps << " cliques intraps ; " << swc.cliques.size() << " cliques en tout)" << endl;
        computeSitesBoundingBoxes(swc.sites, data);
        swc.prepareLabelsZones();
        
        
  
        float _ddweight=0.8, _intrapsweight = 4.0, _simweight=1.0, _lsweight=1.0, _ddx1 = 8.0, _ddx2 = 4.0, _ddh=0.0001;
        cout << _ddweight << "-" << _intrapsweight << "-" << _simweight << "-" << _lsweight << "-" << _ddx2 << "-" << _ddx1 << "-" << _ddh << endl;
        swc.setModelParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh);
        swc.run = 1;
        swc.Initialization();
      
        swc.Run(1);
        swc.SummaryLabels();
 
    }
      
      
  

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

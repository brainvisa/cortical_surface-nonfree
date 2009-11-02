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

#include "blobs.h"




using namespace aims;
using namespace carto;
using namespace std;

vector<Blob *> construireBlobs(PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch){
    vector<Blob *> blobs;
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*> blobList=sketch.BlobSet();
    uint test=0,test2=0;
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*>::iterator blobIt;
    list<GreyLevelBlob<SiteType<AimsSurface<3, Void> >::type > *>::iterator glbit;
    ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *ssb;
    set<SiteType<AimsSurface<3, Void> >::type,ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> >::iterator ptsit;
    
    
    for (blobIt=blobList.begin();blobIt!=blobList.end();blobIt++){
      ssb = *blobIt;
      for (glbit = ssb->glBlobs.begin(); glbit != ssb->glBlobs.end(); glbit++){
        blobs.push_back(new Blob());
        Blob *blob=blobs[blobs.size()-1];
        blob->index = test++;
        blob->parent = test2;
        blob->subject = sketch.Subject();
        blob->t = 0.0;
        blob->scale = (*glbit)->GetScale();
        set<SiteType<AimsSurface<3, Void> >::type,ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> > listePoints;
        listePoints = (*glbit)->GetListePoints();
        for (ptsit=listePoints.begin();ptsit!=listePoints.end();ptsit++){
          (blob->nodes_list).insert((*ptsit).second);
        }
        blobs.push_back(blob);
        //#########################################################
        // CLASSES BLOB ET SSBLOB
        // FONCTION QUI GENERE DES VECTEURS DE BLOB et de SSBLOB
        // ATTENTION LES GLB APPARTENANT AU MEME SSB DOIVENT ETRE RELIES D'UNE FACON OU UNE 
        //    AUTRE
        // FONCTION QUI EXTRAIT LES MAILLAGES QUI CONVIENNENT DU VECTEUR DE BLOBS A PARTIR
        //    DU MAILLAGE
        
      }
      test2++;
    }
    cout << "nb de blobs : " << test << " / glb : " <<blobs.size()<< endl;
    return blobs;
}




ScaleSpace<AimsSurface<3, Void>, Texture<float> > getScaleSpace(TimeTexture<float> &laTexture, AimsSurfaceTriangle &laMesh, TimeTexture<float> &lat, TimeTexture<float> &longit){
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
      ScaleSpace<AimsSurface<3, Void>, Texture<float> > scale_space(&(laMesh[0]), &(laTexture[0]), smooth);
      cout << "Scale-space creation" << endl;
      vector<Point3df> *coordinates;
      coordinates=new vector<Point3df>();
      for (uint i=0;i<lat[0].nItem();i++){
        (*coordinates).push_back(Point3df(lat[0].item(i), longit[0].item(i),i));
      }
      scale_space.PutCoordinates(coordinates);

      scale_space.GenerateDefaultScaleSpace(8.0);
      return scale_space;
     
}




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

int find(const vector<int> &v, int item){
  int c=-1;
  for (uint i=0; i<v.size() && c==-1 ; i++){
    if (item == v[i]) c=(int)i;
  }
  return c;
}

AimsSurfaceTriangle getBlobsMeshes(vector<Blob *> &blobs, AimsSurfaceTriangle &mesh, vector<vector<int> > &nodes_lists){
  AimsSurfaceTriangle objects;
  uint p1,p2,p3;
  nodes_lists=vector<vector<int> >(blobs.size());

  set<uint>::iterator it;

  for (uint i=0;i<blobs.size();i++){
    cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " << objects.size() << flush ;
    set<uint> tri,comp;
    vector<uint> corres;
    for (uint j=0;j<mesh[0].polygon().size();j++){
      p1=mesh[0].polygon()[j][0];
      p2=mesh[0].polygon()[j][1];
      p3=mesh[0].polygon()[j][2];
    
      if (blobs[i]->nodes_list.find(p1)!=blobs[i]->nodes_list.end() && blobs[i]->nodes_list.find(p2)!=blobs[i]->nodes_list.end() &&     blobs[i]->nodes_list.find(p3)!=blobs[i]->nodes_list.end())
        tri.insert(j);
    }
//     if (tri.size() == 0) {countzero++; cout << endl; }

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

AimsSurfaceTriangle getObjects(TimeTexture<short> &tex, AimsSurfaceTriangle &mesh, vector<vector<int> > &nodes_lists){
    
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


AimsSurfaceTriangle getBarycenters(AimsSurfaceTriangle &mesh,  vector<vector<int> > &nodes_lists, float radius){
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



int main( int argc, const char **argv ){
  try {
  
    int mode=0;
    string outpath = "", meshPath, texPath, latpath="", lonpath="", flatpath="", sujet;

    AimsApplication app( argc, argv, "surfLabelsTex2Graph" );
    app.addOption( meshPath, "-m", "mesh");
    app.addOption( texPath, "-t", "texture");
    app.addOption( outpath, "-o", "output file");
    app.addOption( mode, "-M", "mode (0: normal - 1:barycenters (provide the lat/lon textures)",1);
    app.addOption( sujet,"-s", "sujet");
    app.addOption( latpath, "--lat", "latitude");
    app.addOption( lonpath, "--lon", "longitude");
    app.addOption( flatpath, "--flat", "flat",1);
    app.initialize();
    Reader<AimsSurfaceTriangle> meshRdr(meshPath);
    Reader<TimeTexture<float> > texRdr(texPath);
    AimsSurfaceTriangle mesh;
    meshRdr.read(mesh);
    TimeTexture<float> intex;
    texRdr.read(intex);
    TimeTexture<float> latex(1,intex[0].nItem());
    TimeTexture<short> tex(1,intex[0].nItem());
    for (uint i=0;i<intex[0].nItem();i++){
      tex[0].item(i) = (short) intex[0].item(i); //+1;
      latex[0].item(i) = (float) intex[0].item(i); //+1;
    }
    
    mesh[0].setMini(); mesh[0].setMaxi();
    cerr << mesh[0].minimum()[0] << " " << mesh[0].minimum()[1] << " " << mesh[0].minimum()[2] << ";" << mesh[0].maximum()[0] << " " << mesh[0].maximum()[1] << " " << mesh[0].maximum()[2] << endl;
    
    float dist = sqrt(pow(mesh[0].minimum()[0]-mesh[0].maximum()[0],2)+pow(mesh[0].minimum()[1]-mesh[0].maximum()[1],2)+pow(mesh[0].minimum()[2]-mesh[0].maximum()[2],2));
    float radius = dist / 300.0;
    
    vector<vector<int> > nodes_lists;
    AimsSurfaceTriangle *objects;
    objects = new AimsSurfaceTriangle(getObjects(tex,mesh,nodes_lists));
    TimeTexture<float> lat,lon;
    assert(latpath!="");
    assert(lonpath!="");
    Reader<TimeTexture<float> > rlat(latpath),rlon(lonpath);
    rlat.read(lat);
    rlon.read(lon);
    
    
    vector<SSBlob> ssblobs;
    vector<Blob *> blobs;
    
    if (mode == 1){
      objects = new AimsSurfaceTriangle(getBarycenters(mesh,nodes_lists,radius));
    }
    else if (mode == 2){
      ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss(getScaleSpace(latex,mesh,lat,lon));
    
      PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch(sujet, &ss, SURFACE);

      sketch.ComputePrimalSketch(1.0, 8.0, "", 10);
      cout << "CONSTRUIRE BLOBS" << endl;
      blobs = construireBlobs(sketch);
      cout << "FIN CONSTRUIRE BLOBS" << endl;
      objects = new AimsSurfaceTriangle(getBlobsMeshes(blobs,mesh,nodes_lists));

    }
    
    if (flatpath != ""){
      cerr << "écriture flat mesh:" << flatpath << endl;
      TimeTexture<float> texflat;
      AimsSurfaceTriangle flat(getFlatMap(nodes_lists,lat,lon,texflat));
      cerr << flat[0].vertex().size() << " !=! " << texflat[0].nItem() << endl;
      Writer<AimsSurfaceTriangle> wflat(flatpath);
      wflat.write(flat);
//       Writer<TimeTexture<float> > wtex("/volatile/operto/test.tex");
//       wtex.write(texflat);
    }
    
    cerr << "construction graphe" << endl;
    
    Graph graph("BlobsArg");
    vector<float> resolution,bbmin2D,bbmax2D;
    vector<int> bbmin, bbmax;
    resolution.push_back(1.0); resolution.push_back(1.0); resolution.push_back(1.0); 
    bbmin.push_back(mesh[0].minimum()[0]-1); bbmin.push_back(mesh[0].minimum()[1]-1); bbmin.push_back(mesh[0].minimum()[2]-1); 
    bbmax.push_back(mesh[0].maximum()[0]+1); bbmax.push_back(mesh[0].maximum()[1]+1); bbmax.push_back(mesh[0].maximum()[2]+1); 
    graph.setProperty( "filename_base", "*");

    graph.setProperty("voxel_size", resolution);
    graph.setProperty("boundingbox_min", bbmin);
    graph.setProperty("boundingbox_max", bbmax);
    graph.setProperty("mesh", meshPath);
    graph.setProperty("sujet", sujet);
    graph.setProperty("texture", texPath);
    
    Vertex *vert;
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
//     for (int i=0;i<(int)blobs.size();i++){
    for (int i=0;i<(int)(*objects).size();i++){
      if ((*objects)[i].vertex().size()!=0){
        pair<Point2df, Point2df> bb(getBoundingBox(nodes_lists[i],lat,lon));
        float area = (bb.second[0]-bb.first[0])*(bb.second[1]-bb.first[1]);
        if(area<1000.0){
          cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph.order() << flush ;
          vert = graph.addVertex("blob");
          vert->setProperty("index", i);
          vert->setProperty("name", i);
          vert->setProperty("label", "0");
          vert->setProperty("t", 100.0);
          vert->setProperty("rank", i);
          vert->setProperty( "subject", sujet);
          vert->setProperty( "tmin", 1);
          vert->setProperty( "tmax", 4);
          vert->setProperty( "trep", 2);
          vert->setProperty( "depth", 100.0);
          vert->setProperty( "tValue", 100.0);
          vert->setProperty("nodes_list", nodes_lists[i]);
          
          bbmin2D.clear(); bbmax2D.clear();
          bbmin2D.push_back(bb.first[0]);
          bbmin2D.push_back(bb.first[1]);
          bbmin2D.push_back(-1);
          bbmax2D.push_back(bb.second[0]);
          bbmax2D.push_back(bb.second[1]);
          bbmax2D.push_back(-1);
          vert->setProperty( "gravity_center", bbmin2D);
          vert->setProperty("boundingbox_min", bbmin2D);
          vert->setProperty("boundingbox_max", bbmax2D);
          ptr=carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
          (*ptr)[0]=(*objects)[i];
          manip.storeAims(graph, vert, "blob", ptr);
          vert->setProperty("blob_label",i);
        }
      }
    }

    cerr << "graph.order:" << graph.order() << endl;

    Writer<Graph> graphWtr(outpath);
    graphWtr.write(graph);
    
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


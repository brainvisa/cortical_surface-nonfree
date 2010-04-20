
#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/representation.h>
#include <cortical_surface/structuralanalysis/blobs.h>
#include <aims/mesh/surfacegen.h>

using namespace aims;
using namespace carto;
using namespace std;


//##############################################################################

// Function that extracts mesh patches from a "mesh", being given a "blobs" vector.
//   On a sphere version

AimsSurfaceTriangle getBlobsSphericalMeshes ( vector<surf::GreyLevelBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     Texture<float> &lat,
                                     Texture<float> &lon,
                                     vector<set<int> > &nodes_lists){
  AimsSurfaceTriangle objects;
  nodes_lists=vector<set<int> >(blobs.size());

  for (uint i = 0 ; i < blobs.size() ; i++){

    cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " << objects.size() << flush ;
    objects[i] = blobs[i]->getAimsPatchOnASphere(mesh, lat, lon, nodes_lists[i]);

  }

  cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << blobs.size() << endl;
  return objects;
}

AimsSurfaceTriangle getBlobsSphericalMeshes ( vector<surf::ScaleSpaceBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     Texture<float> &lat,
                                     Texture<float> &lon,
                                     vector<set<int> > &nodes_lists){
  AimsSurfaceTriangle objects;
  nodes_lists=vector<set<int> >(blobs.size());

  for (uint i = 0 ; i < blobs.size() ; i++){

    cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " << objects.size() << flush ;
    objects[i] = blobs[i]->getAimsPatchOnASphere(mesh, lat, lon, nodes_lists[i]);

  }

  cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << blobs.size() << endl;
  return objects;
}


//##############################################################################
// Function that extracts mesh patches from a "mesh", being given a "blobs" vector.
//   All mesh patches are flat, belong to a same plane, and their z varies with
//   their scale
AimsSurfaceTriangle getBlobs2DMeshes ( vector<surf::GreyLevelBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     Texture<float> &lat,
                                     Texture<float> &lon,
                                     vector<set<int> > &nodes_lists){
  AimsSurfaceTriangle objects;
  nodes_lists=vector<set<int> >(blobs.size());


  for (uint i = 0 ; i < blobs.size() ; i++){

    cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " << objects.size() << flush ;
    objects[i] = blobs[i]->getAimsPatchOnAPlane(mesh, lat, lon, nodes_lists[i]);


  }

  cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << blobs.size() << endl;
  return objects;
}


//##############################################################################

// Function that extracts mesh patches from a "mesh", being given a "blobs" vector,
//   and returning a collection of "objects" plus a vector of "nodes_lists".
AimsSurfaceTriangle getBlobsMeshes ( vector<surf::GreyLevelBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     vector<set<int> > &nodes_lists){
  AimsSurfaceTriangle objects;
  nodes_lists=vector<set<int> >(blobs.size());


  for (uint i = 0 ; i < blobs.size() ; i++){

    cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " << objects.size() << flush ;
    objects[i] = blobs[i]->getAimsMeshPatch(mesh, nodes_lists[i]);


  }

  cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << blobs.size() << endl;
  return objects;
}

//##############################################################################


// That function takes a vector of SSBlob and build mesh patches corresponding to
//   previously computed representation blobs
AimsSurfaceTriangle getBlobsMeshes ( vector<surf::ScaleSpaceBlob *> &blobs,
                                     AimsSurface<3, Void> &mesh,
                                     vector<vector<int> > &nodes_lists){

    cout << "mesh.vertex:" << mesh.vertex().size() << endl;
    cout << "mesh.polygon:" << mesh.polygon().size() << endl;
    AimsSurfaceTriangle objects;
    uint p1,p2,p3;
    nodes_lists = vector<vector<int> >(blobs.size());

    set<uint>::iterator it;

    for (uint i = 0 ; i < blobs.size() ; i++) {

      cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << " " <<
              objects.size() << flush ;
      set<uint> tri,comp;
      vector<uint> corres;
      cout << endl << blobs[i]->nodes.size() << endl;
      for (uint j = 0 ; j < mesh.polygon().size() ; j++){

        p1=mesh.polygon()[j][0];
        p2=mesh.polygon()[j][1];
        p3=mesh.polygon()[j][2];

        if ( blobs[i]->nodes.find(p1) != blobs[i]->nodes.end()
          && blobs[i]->nodes.find(p2) != blobs[i]->nodes.end()
          && blobs[i]->nodes.find(p3) != blobs[i]->nodes.end()
            )
          tri.insert(j);

      }

      cout << "t " << tri.size() << " " << flush;

      for (it = tri.begin() ; it != tri.end() ; it++){
        p1=mesh.polygon()[*it][0];
        p2=mesh.polygon()[*it][1];
        p3=mesh.polygon()[*it][2];
        comp.insert(p1); comp.insert(p2); comp.insert(p3);
      }

      corres = vector<uint>( mesh.vertex().size() );

      for (it = comp.begin() ; it != comp.end() ; it++){

        assert( *it < corres.size() );
        assert( *it < mesh.vertex().size() );
        assert( i < nodes_lists.size() );

        (objects)[i].vertex().push_back( mesh.vertex()[*it] );
        corres[*it] = (objects)[i].vertex().size() - 1;
        nodes_lists[i].push_back( *it );

      }

      for (it = tri.begin() ; it != tri.end() ; it++){

        p1 = mesh.polygon()[*it][0];
        p2 = mesh.polygon()[*it][1];
        p3 = mesh.polygon()[*it][2];
        (objects)[i].polygon().push_back( AimsVector<uint,3>(corres[p1], corres[p2], corres[p3]) );

      }

    }
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << blobs.size() << endl;

    return objects;
  }



//##############################################################################





AimsSurfaceTriangle getFlatMap ( vector<set<int> > &nodes_lists,
                                 TimeTexture<float> &lat,
                                 TimeTexture<float> &lon,
                                 TimeTexture<float> &tex){
  AimsSurfaceTriangle objects;
  for (uint i=0;i<nodes_lists.size();i++){
    if (nodes_lists[i].size()!=0){
      pair<Point2df,Point2df> bb(getBoundingBox(nodes_lists[i],lat,lon));
      assert(bb.first[0]<=bb.second[0] || !(cout << bb.first[0] << " /\\" << bb.second[0] << endl));
      assert(bb.first[1]<=bb.second[1]|| !(cout << bb.first[1] << " /\\" << bb.second[1] << endl));
      float area = (bb.second[0]-bb.first[0])*(bb.second[1]-bb.first[1]);
      if(area<1000.0){
        tex[0].push_back(area);
        tex[0].push_back(area);
        tex[0].push_back(area);
        tex[0].push_back(area);
//         cerr << tex[0].nItem() << " " << flush;
        objects[0].vertex().push_back(Point3df(bb.first[0],bb.first[1],0.001));
        objects[0].vertex().push_back(Point3df(bb.first[0],bb.second[1],0.002));
        objects[0].vertex().push_back(Point3df(bb.second[0],bb.second[1],0.003));
        objects[0].vertex().push_back(Point3df(bb.second[0],bb.first[1],0.0005));
//         cerr << objects[0].vertex().size() << endl;
        objects[0].polygon().push_back(AimsVector<uint,3>(objects[0].vertex().size()-4,objects[0].vertex().size()-3,objects[0].vertex().size()-2));
        objects[0].polygon().push_back(AimsVector<uint,3>(objects[0].vertex().size()-2,objects[0].vertex().size()-1,objects[0].vertex().size()-4));
      }
    }
  }

  return objects;
}

//##############################################################################

// This function extracts normal mesh patches corresponding to the different label
//   components existing in a texture. It returns a vector of nodes lists.
AimsSurfaceTriangle getLabelObjectsOnASphere( TimeTexture<short> &tex,
                                AimsSurface<3,Void> &mesh,
                                Texture<float> &lat,
                                Texture<float> &lon,
                                vector<set<int> > &nodes_lists){

    int labelmax=0;
    for (uint i=0;i<tex[0].nItem();i++){
      if (tex[0].item(i)>labelmax)
        labelmax=tex[0].item(i);
    }
    AimsSurfaceTriangle objects;
    nodes_lists=vector<set<int> >(labelmax+1);
    vector<set<uint> > triangles(labelmax+1);

    set<uint>::iterator it;
    set<uint> tri,comp;
    uint p1,p2,p3;  int L1,L2,L3;
    cout << "TEST" << endl;
    for (uint i=0;i<mesh.polygon().size();i++){
      p1=mesh.polygon()[i][0];
      p2=mesh.polygon()[i][1];
      p3=mesh.polygon()[i][2];

      L1=tex[0].item(p1); L2=tex[0].item(p2); L3=tex[0].item(p3);

      if (L1==L2 || L1==L3)
        triangles[L1].insert(i);
      else if (L2==L3)
        triangles[L2].insert(i);
    }
    vector<uint> corres;

    cout << "TRI:" << triangles.size() << endl;

    for (uint i=0;i<triangles.size();i++){
      tri = triangles[i];

      comp.clear();
      for (it=tri.begin();it!=tri.end();it++){
        p1=mesh.polygon()[*it][0];
        p2=mesh.polygon()[*it][1];
        p3=mesh.polygon()[*it][2];
        comp.insert(p1); comp.insert(p2); comp.insert(p3);
      }
      corres=vector<uint>(mesh.vertex().size());
      for (it=comp.begin();it!=comp.end();it++){
        assert(*it<corres.size());
        assert(*it<mesh.vertex().size());
        assert(i<nodes_lists.size());
        (objects)[i].vertex().push_back(Point3df ( 1.0 * cos((lat.item(*it)-90.)/180.0*3.1415957) * cos(lon.item(*it)/180.0*3.1415957),
                    1.0 * cos((lat.item(*it)-90.)/180.0*3.1415957) * sin(lon.item(*it)/180.0*3.1415957),
                    1.0 * sin((lat.item(*it)-90.)/180.0*3.1415957) ));
        corres[*it]=(objects)[i].vertex().size()-1;
        nodes_lists[i].insert(*it);

      }
      cout << (objects)[i].vertex().size() << endl;
      for (it=tri.begin();it!=tri.end();it++){
        p1=mesh.polygon()[*it][0];
        p2=mesh.polygon()[*it][1];
        p3=mesh.polygon()[*it][2];
        (objects)[i].polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
      }
    }
    return objects;
}

//##############################################################################


AimsSurfaceTriangle getLabelObjectsOnAMesh( TimeTexture<short> &tex,
                                AimsSurface<3,Void> &mesh,
                                vector<set<int> > &nodes_lists){

    int labelmax=0;
    for (uint i=0;i<tex[0].nItem();i++){
      if (tex[0].item(i)>labelmax)
        labelmax=tex[0].item(i);
    }
    AimsSurfaceTriangle objects;
    nodes_lists=vector<set<int> >(labelmax+1);
    vector<set<uint> > triangles(labelmax+1);

    set<uint>::iterator it;
    set<uint> tri,comp;
    uint p1,p2,p3;  int L1,L2,L3;

    for (uint i=0;i<mesh.polygon().size();i++){
      p1=mesh.polygon()[i][0];
      p2=mesh.polygon()[i][1];
      p3=mesh.polygon()[i][2];

      L1=tex[0].item(p1); L2=tex[0].item(p2); L3=tex[0].item(p3);

      if (L1==L2 || L1==L3)
        triangles[L1].insert(i);
      else if (L2==L3)
        triangles[L2].insert(i);
    }
    vector<uint> corres;

    cout << "TRI:" << triangles.size() << endl;

    for (uint i=0;i<triangles.size();i++){
      tri = triangles[i];

      comp.clear();
      for (it=tri.begin();it!=tri.end();it++){
        p1=mesh.polygon()[*it][0];
        p2=mesh.polygon()[*it][1];
        p3=mesh.polygon()[*it][2];
        comp.insert(p1); comp.insert(p2); comp.insert(p3);
      }
      corres=vector<uint>(mesh.vertex().size());
      for (it=comp.begin();it!=comp.end();it++){
        assert(*it<corres.size());
        assert(*it<mesh.vertex().size());
        assert(i<nodes_lists.size());
        (objects)[i].vertex().push_back(mesh.vertex()[*it]);
        corres[*it]=(objects)[i].vertex().size()-1;
        nodes_lists[i].insert(*it);

      }
      cout << (objects)[i].vertex().size() << endl;
      for (it=tri.begin();it!=tri.end();it++){
        p1=mesh.polygon()[*it][0];
        p2=mesh.polygon()[*it][1];
        p3=mesh.polygon()[*it][2];
        (objects)[i].polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
      }
    }
    return objects;
}

//##############################################################################

Point3df getBlobBarycenterOnASphere( surf::GreyLevelBlob *glb ) {
    float latMoy = 0.0, lonMoy = 0.0;
    set<int>::iterator it;

    for ( it = glb->nodes.begin() ; it != glb->nodes.end() ; it ++ ) {

        latMoy += glb->coordinates[*it][0];
        lonMoy += glb->coordinates[*it][1];
    }
    ASSERT(glb->nodes.size() == glb->coordinates.size());

    latMoy /= glb->nodes.size();
    lonMoy /= glb->nodes.size();
    return Point3df( log(glb->scale) * cos((latMoy-90.)/180.0*3.1415957) * cos(lonMoy/180.0*3.1415957),
                log(glb->scale) * cos((latMoy-90.)/180.0*3.1415957) * sin(lonMoy/180.0*3.1415957),
                log(glb->scale) * sin((latMoy-90.)/180.0*3.1415957) );

}


float compareBlobsScales(const surf::GreyLevelBlob *s1, const surf::GreyLevelBlob *s2){
    ASSERT(s1->scale != s2->scale);
    return (s1->scale - s2->scale);
}

struct ltstr
{
  bool operator()(const surf::GreyLevelBlob * s1, const surf::GreyLevelBlob * s2) const
  {
    return compareBlobsScales(s1, s2) < 0.0;
  }
};





AimsSurfaceTriangle getLinkMesh ( surf::ScaleSpaceBlob *ssb, vector< pair<uint, uint> > &blobsIndices ) {
    AimsSurfaceTriangle mesh, *cyl;
    set<surf::GreyLevelBlob *>::iterator itB;
    set<surf::GreyLevelBlob *> &unsortedListGLB = ssb->blobs;
    set<surf::GreyLevelBlob *, ltstr> listGLB;
    set<surf::GreyLevelBlob *, ltstr>::iterator itB1, itB2;
    for ( itB = unsortedListGLB.begin() ; itB != unsortedListGLB.end() ; itB ++ )
        listGLB.insert(*itB);
    
    itB1 = listGLB.begin();
    itB2 = itB1;
    if ( itB2 != listGLB.end() )
        itB2++;
    else
        ASSERT(false);
    uint i = 0;
    while ( itB2 != listGLB.end() ) {
        surf::GreyLevelBlob *glb1, *glb2;
        glb1 = (*itB1);
        glb2 = (*itB2);
        Point3df    p1 ( getBlobBarycenterOnASphere( glb1 ) ),
                    p2 ( getBlobBarycenterOnASphere( glb2 ) ) ;
        cout << glb1->index << ":" << p1[0] << " " << p1[1] << " " << p1[2] << ";" << glb2->index << ":" << p2[0] << " " << p2[1] << " " << p2[2] << endl;
        cyl = SurfaceGenerator::cylinder(p1, p2, 0.001, 0.001, 10, true, true);
        pair<uint, uint> ind ( glb1->index, glb2->index );
        mesh[i] = (*cyl)[0];
        blobsIndices.push_back( ind );
        itB1++, itB2++, i++;
    }

    return mesh;
}

AimsSurfaceTriangle getG2GRelationsMeshes ( vector<surf::ScaleSpaceBlob *> &ssblobs,
                                            vector< pair<uint, uint> > &blobsIndices ) {
    
    AimsSurfaceTriangle meshes;
    uint j = 0;

    for ( uint i = 0 ; i < ssblobs.size() ; i ++ ) {
        AimsSurfaceTriangle linkMesh ;
        
        linkMesh = getLinkMesh( ssblobs[i], blobsIndices );
        cout << linkMesh.size() << " " << endl;
        for ( uint k = 0 ; k < linkMesh.size() ; k++, j++ ) 
            meshes[j] = linkMesh[k];
        j = j + ssblobs[i]->blobs.size();
    }
    ASSERT( meshes.size() == blobsIndices.size() );
    cout << "G2G : " << meshes.size() << flush;
    return meshes;
}

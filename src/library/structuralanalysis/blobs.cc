#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/blobs.h>
#include <time.h>

using namespace aims;
using namespace carto;
using namespace std;

//##############################################################################

Point3df Point3dfOnSphere ( int i, float radius,
                            Texture<float> &lat,
                            Texture<float> &lon ){

    return  Point3df ( log(radius) * cos((lat.item(i)-90.)/180.0*3.1415957) * cos(lon.item(i)/180.0*3.1415957),
                    log(radius) * cos((lat.item(i)-90.)/180.0*3.1415957) * sin(lon.item(i)/180.0*3.1415957),
                    log(radius) * sin((lat.item(i)-90.)/180.0*3.1415957) );    //(float)(rand()/RAND_MAX) * 0.001 ));
}

Point3df Point3dfOnMesh ( int i, AimsSurface<3,Void> &mesh ) {

    return mesh.vertex()[i];

}

Point3df Point3dfOnPlane ( int i, float height,
                            Texture<float> &lat,
                            Texture<float> &lon ){

    return Point3df(lat.item(i), lon.item(i), height * 10.0 + (float)(rand()/RAND_MAX) * 0.001 );

}


//##############################################################################

AimsSurface<3, Void> surf::Blob::getAimsPatchOnAPlane ( AimsSurface<3, Void> &mesh,
                                              Texture<float> &lat,
                                              Texture<float> &lon,
                                              float height,
                                              set<int> &nodes_list ) {
  srand( (unsigned)time( NULL ) );

  AimsSurface<3, Void> patch;
  uint p1,p2,p3;

  set<uint>::iterator it;
  set<uint> tri, comp;
  vector<uint> corres;

  for (uint j = 0 ; j < mesh.polygon().size() ; j++){

    p1=mesh.polygon()[j][0];
    p2=mesh.polygon()[j][1];
    p3=mesh.polygon()[j][2];

    if ( nodes.find(p1) != nodes.end() &&
            nodes.find(p2) != nodes.end() &&
            nodes.find(p3) != nodes.end() )
      tri.insert(j);
  }

  for (it=tri.begin();it!=tri.end();it++){
    p1=mesh.polygon()[*it][0];
    p2=mesh.polygon()[*it][1];
    p3=mesh.polygon()[*it][2];
    comp.insert(p1); comp.insert(p2); comp.insert(p3);
  }

  corres = vector<uint>(mesh.vertex().size());

  for (it = comp.begin() ; it != comp.end() ; it++){
    assert(*it<corres.size());
    assert(*it<mesh.vertex().size());

    patch.vertex().push_back( Point3dfOnPlane(*it, height, lat, lon) );

    corres[*it] = patch.vertex().size()-1;
    nodes_list.insert(*it);
  }

  for (it=tri.begin();it!=tri.end();it++){
    p1=mesh.polygon()[*it][0];
    p2=mesh.polygon()[*it][1];
    p3=mesh.polygon()[*it][2];
    patch.polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
  }
  return patch;

}

//##############################################################################

AimsSurface<3, Void> surf::Blob::getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh,
                                              Texture<float> &lat,
                                              Texture<float> &lon,
                                              float radius,
                                              set<int> &nodes_list ) {
  AimsSurface<3, Void> patch;
  uint p1,p2,p3;

  set<uint>::iterator it;
  set<uint> tri, comp;
  vector<uint> corres;

  for (uint j = 0 ; j < mesh.polygon().size() ; j++){

    p1=mesh.polygon()[j][0];
    p2=mesh.polygon()[j][1];
    p3=mesh.polygon()[j][2];

    if ( nodes.find(p1) != nodes.end() &&
            nodes.find(p2) != nodes.end() &&
            nodes.find(p3) != nodes.end() )
      tri.insert(j);
  }

  for (it=tri.begin();it!=tri.end();it++){
    p1=mesh.polygon()[*it][0];
    p2=mesh.polygon()[*it][1];
    p3=mesh.polygon()[*it][2];
    comp.insert(p1); comp.insert(p2); comp.insert(p3);
  }

  corres = vector<uint>(mesh.vertex().size());

  for (it = comp.begin() ; it != comp.end() ; it++){
    assert(*it<corres.size());
    assert(*it<mesh.vertex().size());

    patch.vertex().push_back( Point3dfOnSphere(*it, radius, lat, lon) );

    corres[*it] = patch.vertex().size()-1;
    nodes_list.insert(*it);
  }

  for (it=tri.begin();it!=tri.end();it++){
    p1=mesh.polygon()[*it][0];
    p2=mesh.polygon()[*it][1];
    p3=mesh.polygon()[*it][2];
    patch.polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
  }
  return patch;

}

//##############################################################################


AimsSurface<3, Void> surf::Blob::getAimsMeshPatch ( AimsSurface<3, Void> &mesh,
                                                    set<int> &nodes_list ){

  AimsSurface<3, Void> patch;
  uint p1,p2,p3;

  set<uint>::iterator it;
  set<uint> tri, comp;
  vector<uint> corres;

  for (uint j = 0 ; j < mesh.polygon().size() ; j++){

    p1=mesh.polygon()[j][0];
    p2=mesh.polygon()[j][1];
    p3=mesh.polygon()[j][2];

    if ( nodes.find(p1) != nodes.end() &&
            nodes.find(p2) != nodes.end() &&
            nodes.find(p3) != nodes.end() )
      tri.insert(j);
  }

  for (it=tri.begin();it!=tri.end();it++){
    p1=mesh.polygon()[*it][0];
    p2=mesh.polygon()[*it][1];
    p3=mesh.polygon()[*it][2];
    comp.insert(p1); comp.insert(p2); comp.insert(p3);
  }

  corres = vector<uint>(mesh.vertex().size());

  for (it = comp.begin() ; it != comp.end() ; it++){
    assert(*it<corres.size());
    assert(*it<mesh.vertex().size());

    patch.vertex().push_back( Point3dfOnMesh(*it, mesh) );

    corres[*it] = patch.vertex().size()-1;
    nodes_list.insert(*it);
  }

  for (it=tri.begin();it!=tri.end();it++){
    p1=mesh.polygon()[*it][0];
    p2=mesh.polygon()[*it][1];
    p3=mesh.polygon()[*it][2];
    patch.polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
  }
  return patch;
}

//##############################################################################

AimsSurface<3, Void> surf::GreyLevelBlob::getAimsMeshPatch ( AimsSurface<3, Void> &mesh,
                                                              set<int> &nodes_list ){

  return surf::Blob::getAimsMeshPatch(mesh, nodes_list);

}

AimsSurface<3, Void> surf::GreyLevelBlob::getAimsPatchOnAPlane ( AimsSurface<3, Void> &mesh,
                                                                  Texture<float> &lat,
                                                                  Texture<float> &lon,
                                                                  set<int> &nodes_list ){

  return surf::Blob::getAimsPatchOnAPlane(mesh, lat, lon, scale, nodes_list);

}

AimsSurface<3, Void> surf::GreyLevelBlob::getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh,
                                                              Texture<float> &lat,
                                                              Texture<float> &lon,
                                                              set<int> &nodes_list ){

  return surf::Blob::getAimsPatchOnASphere(mesh, lat, lon, scale, nodes_list);

}

//##############################################################################

AimsSurface<3, Void> surf::ScaleSpaceBlob::getAimsMeshPatch ( AimsSurface<3, Void> &mesh,
                                                              set<int> &nodes_list ){

  return surf::Blob::getAimsMeshPatch(mesh, nodes_list);

}

AimsSurface<3, Void> surf::ScaleSpaceBlob::getAimsPatchOnAPlane ( AimsSurface<3, Void> &mesh,
                                                                  Texture<float> &lat,
                                                                  Texture<float> &lon,
                                                                  set<int> &nodes_list ){

  return surf::Blob::getAimsPatchOnAPlane(mesh, lat, lon, tmin, nodes_list);

}

AimsSurface<3, Void> surf::ScaleSpaceBlob::getAimsPatchOnASphere ( AimsSurface<3, Void> &mesh,
                                                              Texture<float> &lat,
                                                              Texture<float> &lon,
                                                              set<int> &nodes_list ){

  return surf::Blob::getAimsPatchOnASphere(mesh, lat, lon, tmin, nodes_list);

}

//##############################################################################


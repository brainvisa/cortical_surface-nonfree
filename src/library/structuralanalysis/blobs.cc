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

double getOverlapMeasure( Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap ){

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

void filteringBlobs (  vector<surf::ScaleSpaceBlob *> & ssblobs,
                       vector<surf::GreyLevelBlob *> &filteredBlobs,
                       vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
                       Point3df bbmin2,
                       Point3df bbmax2  ){
                           
    // Filtering according to positions

    for (uint i = 0 ; i < ssblobs.size() ; i++ ){
        string subject;
        subject = ssblobs[i]->subject;
        bool firstGLB = false;

        surf::ScaleSpaceBlob *ssb;
        ssb = new surf::ScaleSpaceBlob();
        ssb->index = ssblobs[i]->index;
        ssb->t = ssblobs[i]->t;
        ssb->label = ssblobs[i]->label;
        ssb->subject = ssblobs[i]->subject;
        ssb->tmin = ssblobs[i]->tmin;
        ssb->tmax = ssblobs[i]->tmax;
 
        set<surf::GreyLevelBlob *>::iterator itB1;
        for (itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() ; itB1++){

            pair<Point2df,Point2df> bbi = getBoundingBox( (*itB1)->nodes, (*itB1)->coordinates );

            Point3df bbmin1 (bbi.first[0], bbi.first[1], 0.0),
                bbmax1 (bbi.second[0], bbi.second[1], 0.0);
//                 bbmin2 (105.0, 256.0, 0.0), // SUPERIOR TEMPORAL SULCUS
//                 bbmax2 (107.0, 306.0, 0.0) ;  // STS

            uint no_overlap = 2;
            double overlap = getOverlapMeasure( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
            if (no_overlap == 0){
                surf::GreyLevelBlob *glb;
                glb = new surf::GreyLevelBlob();
                glb->t = (*itB1)->t;
                glb->scale = (*itB1)->scale;
                glb->boundingbox_max = (*itB1)->boundingbox_max;
                glb->boundingbox_min = (*itB1)->boundingbox_min;
                glb->nodes = (*itB1)->nodes;
                glb->coordinates = (*itB1)->coordinates;

                firstGLB = true;

                glb->ssb_parent = ssb;
                ssb->blobs.insert(glb);
                filteredBlobs.push_back(glb);
            }

        }
        if (firstGLB)
            filteredSsblobs.push_back(ssb);
        else
            delete(ssb);
    }

    // Now that the blobs are filtered, we add the correct bifurcations
    set< uint > filteredIndices;
    for ( uint i = 0 ; i < filteredSsblobs.size() ; i ++ ) {
        filteredIndices.insert ( filteredSsblobs[i]->index ) ;
    }
    for ( uint i = 0 ; i < filteredSsblobs.size() ; i ++ ) {
        set<surf::ScaleSpaceBlob *> auxTop( filteredSsblobs[i]->topBlobs );
        set<surf::ScaleSpaceBlob *>::iterator it;
        filteredSsblobs[i]->topBlobs.clear();
        for ( it = auxTop.begin() ; it != auxTop.end() ; it ++ ) {
            if (filteredIndices.find((*it)->index) != filteredIndices.end())
               filteredSsblobs[i]->topBlobs.insert(*it);
        }
        set<surf::ScaleSpaceBlob *> auxBot( filteredSsblobs[i]->bottomBlobs );
        filteredSsblobs[i]->bottomBlobs.clear();
        for ( it = auxBot.begin() ; it != auxBot.end() ; it ++ ) {
            if (filteredIndices.find((*it)->index) != filteredIndices.end())
               filteredSsblobs[i]->bottomBlobs.insert(*it);
        }
    }

}

// #############################################################################

pair<Point2df, Point2df> getBoundingBox ( set<int> &nodes_list,
                                          map<int, float> &lat,
                                          map<int, float> &lon ){
  Point2df bbmin, bbmax;
  bbmin[0] = 181.0;
  bbmin[1] = 361.0;
  bbmax[0] = -1.0;
  bbmax[1] = -1.0;
  set<int>::iterator it;
  pair<Point2df, Point2df> bb;
  for (it = nodes_list.begin() ; it != nodes_list.end() ; it ++){
    if (lat[*it] < bbmin[0])
      bbmin[0] = lat[*it];
    if (lon[*it] < bbmin[1])
      bbmin[1] = lon[*it];

    if (lat[*it] > bbmax[0])
      bbmax[0] = lat[*it];
    if (lon[*it] > bbmax[1])
      bbmax[1] = lon[*it];
  }

  if ( bbmax[1] > 300.0 && bbmin[1] < 60.0 ) {
    for ( uint i = 0 ; i < nodes_list.size() ; i++ ) {
      if ( lon[*it] >300.0 && lon[*it] < bbmax[1])
        bbmax[1] = lon[*it];
      if (lon[*it] < 60.0 && lon[*it] > bbmin[1])
        bbmin[1] = lon[*it];
    }
  }

  bb.first = bbmin;
  bb.second = bbmax;
  return bb;
}


pair<Point2df, Point2df> getBoundingBox ( set<int> &nodes_list,
                                          TimeTexture<float> &lat,
                                          TimeTexture<float> &lon ){
    map<int, float> mapLat, mapLon;
    set<int>::iterator it;
    for (it = nodes_list.begin() ; it != nodes_list.end() ; it ++){
        mapLat[*it] = lat[0].item(*it);
        mapLon[*it] = lon[0].item(*it);
    }
    return getBoundingBox( nodes_list, mapLat, mapLon );
}

pair<Point2df, Point2df> getBoundingBox ( set<int> &nodes_list,
                                          map<int, vector<float> > &coordinates ){
    map<int, float> mapLat, mapLon;
    set<int>::iterator it;
    for (it = nodes_list.begin() ; it != nodes_list.end() ; it ++){
        mapLat[*it] = coordinates[*it][0];
        mapLon[*it] = coordinates[*it][1];
    }
    return getBoundingBox( nodes_list, mapLat, mapLon );
}

//##############################################################################


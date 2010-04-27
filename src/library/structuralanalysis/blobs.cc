#include <cstdlib>
#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/blobs.h>
#include <time.h>
#include <aims/mesh/surfacegen.h>


using namespace aims;
using namespace carto;
using namespace std;


float compareBlobsScales(const surf::GreyLevelBlob *s1, const surf::GreyLevelBlob *s2){
    ASSERT(s1->scale != s2->scale);
    return (s1->scale - s2->scale);
}


//##############################################################################

Point3df Point3dfOnSphere ( int i, float radius,
                            Texture<float> &lat,
                            Texture<float> &lon ){

    return  Point3df ( log(radius) * cos((lat.item(i)-90.)/180.0*3.1415957) * cos(lon.item(i)/180.0*3.1415957),
                    log(radius) * cos((lat.item(i)-90.)/180.0*3.1415957) * sin(lon.item(i)/180.0*3.1415957),
                    log(radius) * sin((lat.item(i)-90.)/180.0*3.1415957) );    //(float)(rand()/RAND_MAX) * 0.001 ));
}

Point3df Point3dfOnMesh ( int i, AimsSurface<3,Void> &mesh, float radius = 1.0) {

    return Point3df( log(radius) * mesh.vertex()[i][0],
                     log(radius) * mesh.vertex()[i][1],
                     log(radius) * mesh.vertex()[i][2] );

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
    ASSERT ( mesh.vertex().size() != 0 ) ;
    if ( mesh.vertex().size() > 1 ) {
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
    }
    else {
        AimsSurfaceTriangle *sphere;
        Point3df p1 ( mesh.vertex()[0] );
        for ( uint i = 0 ; i < 2 ; i ++ )
            p1[i] *= height;
        sphere = SurfaceGenerator::sphere(p1, 0.02, 10, true);
        patch = (*sphere)[0];
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
    ASSERT ( mesh.vertex().size() != 0 ) ;
    if ( mesh.vertex().size() > 1 ) {
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
    }
    else {
        AimsSurfaceTriangle *sphere;
        Point3df p1 ( mesh.vertex()[0] );
        for ( uint i = 0 ; i < 2 ; i ++ )
            p1[i] *= radius;
        sphere = SurfaceGenerator::sphere(p1, 0.02, 10);
        patch = (*sphere)[0];
    }
    return patch;

}

//##############################################################################


AimsSurface<3, Void> surf::Blob::getAimsMeshPatch ( AimsSurface<3, Void> &mesh,
                                                    set<int> &nodes_list,
                                                    float radius){
    
    AimsSurface<3, Void> patch;
    ASSERT ( mesh.vertex().size() != 0 ) ;
        
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
    if ( tri.size() >= 1 ) {

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

            patch.vertex().push_back( Point3dfOnMesh(*it, mesh, radius) );

            corres[*it] = patch.vertex().size()-1;
            nodes_list.insert(*it);
        }

        for (it=tri.begin();it!=tri.end();it++){
            p1=mesh.polygon()[*it][0];
            p2=mesh.polygon()[*it][1];
            p3=mesh.polygon()[*it][2];
            patch.polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
        }
    }
    else {
        AimsSurfaceTriangle *sphere;
        
        Point3df p1 ( mesh.vertex()[*(nodes.begin())] );
        for ( uint i = 0 ; i < 3 ; i ++ )
            p1[i] *= log( radius );
        sphere = SurfaceGenerator::sphere(p1, 0.2, 10);
        patch = (*sphere)[0];

    }
    return patch;
}


//##############################################################################

AimsSurface<3, Void> surf::GreyLevelBlob::getAimsMeshPatch ( AimsSurface<3, Void> &mesh,
                                                              set<int> &nodes_list ){

  return surf::Blob::getAimsMeshPatch(mesh, nodes_list, scale + 1.0);

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

int getEcartMaxIndice( set<float> &longitudes ) {

    uint i = 0, imax = -1 ;
    set<float>::iterator itf, itf2;

    itf = longitudes.begin();
    itf2 = longitudes.end();
    itf2--;
    
    if ( *itf2 - *itf > 300.0 ) {
        imax = 0;
        itf2 = itf;
        itf2 ++;
        float ecart = *itf2 - *itf;
        for ( ; itf2 != longitudes.end() ; ) {
            if ( *itf2 - *itf > ecart ) {
                imax = i;
                ecart = *itf2 - *itf;
            }
            i++, itf++, itf2++;
        }
    }
    return imax;

}

//##############################################################################

Point3df surf::GreyLevelBlob::getBlobBarycenterOnASphere( void ) {
    surf::GreyLevelBlob *glb;
    glb = this;
    float latMoy = 0.0, lonMoy = 0.0;
    ASSERT(glb->nodes.size() == glb->coordinates.size());
    
    set<float> latitudes, longitudes;
    
    set<int>::iterator it;
    set<float>::iterator itf, itf2;
    ASSERT( glb->nodes.size() >= 1 );
    it = glb->nodes.begin();
    
    if ( glb->nodes.size() == 1 ) {
        latMoy = glb->coordinates[*it][0];
        lonMoy = glb->coordinates[*it][1];
    }
    else {
        
        for ( it = glb->nodes.begin() ; it != glb->nodes.end() ; it ++ ) {
            latMoy += glb->coordinates[*it][0];
            lonMoy += glb->coordinates[*it][1];
            latitudes.insert( glb->coordinates[*it][0] );
            longitudes.insert( glb->coordinates[*it][1] );
        }
        uint imax = getEcartMaxIndice( longitudes );
        lonMoy += ( imax + 1) * 360.0;

        latMoy /= glb->nodes.size();
        lonMoy /= glb->nodes.size();
    }
    return Point3df( log(glb->scale + 1) * cos((latMoy-90.)/180.0*3.1415957) * cos(lonMoy/180.0*3.1415957),
                log(glb->scale + 1) * cos((latMoy-90.)/180.0*3.1415957) * sin(lonMoy/180.0*3.1415957),
                log(glb->scale + 1) * sin((latMoy-90.)/180.0*3.1415957) );
}

Point3df surf::GreyLevelBlob::getBlobBarycenter( void ) {
    surf::GreyLevelBlob *glb;
    glb = this;
    float xMoy = 0.0, yMoy = 0.0, zMoy = 0.0;
    ASSERT(glb->nodes.size() == glb->raw_coordinates.size());    
    
    set<float> latitudes, longitudes;
    
    set<int>::iterator it;
    set<float>::iterator itf;
    ASSERT( glb->nodes.size() >= 1 );
            
    for ( it = glb->nodes.begin() ; it != glb->nodes.end() ; it ++ ) {
        ASSERT( glb->raw_coordinates[*it].size() == 3 );
        xMoy += glb->raw_coordinates[*it][0];
        yMoy += glb->raw_coordinates[*it][1];
        zMoy += glb->raw_coordinates[*it][2];
    }

    xMoy /= glb->nodes.size();
    yMoy /= glb->nodes.size();
    zMoy /= glb->nodes.size();
    return Point3df( log(glb->scale + 1) * xMoy,
                     log(glb->scale + 1) * yMoy,
                     log(glb->scale + 1) * zMoy );
}

// PAS VERIFIEE
Point3df surf::GreyLevelBlob::getBlobBarycenterOnAPlane( void ) {
    surf::GreyLevelBlob *glb;
    glb = this;
    float latMoy = 0.0, lonMoy = 0.0;
    ASSERT(glb->nodes.size() == glb->coordinates.size());
    
    set<float> latitudes, longitudes;
    
    set<int>::iterator it;
    set<float>::iterator itf, itf2;
    ASSERT( glb->nodes.size() >= 1 );
    it = glb->nodes.begin();
    
    if ( glb->nodes.size() == 1 ) {
        latMoy = glb->coordinates[*it][0];
        lonMoy = glb->coordinates[*it][1];
    }
    else {
        
        for ( it = glb->nodes.begin() ; it != glb->nodes.end() ; it ++ ) {
            latMoy += glb->coordinates[*it][0];
            lonMoy += glb->coordinates[*it][1];
            latitudes.insert( glb->coordinates[*it][0] );
            longitudes.insert( glb->coordinates[*it][1] );
        }
        uint imax = getEcartMaxIndice( longitudes );
        lonMoy += ( imax + 1) * 360.0;
        
        
        latMoy /= glb->nodes.size();
        lonMoy /= glb->nodes.size();
    }
    return Point3df( log(glb->scale + 1) * latMoy,
                     log(glb->scale + 1) * lonMoy,
                     log(glb->scale + 1) * 1.0 );
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



double getOverlapMeasure( Point3df bbmin1,
                          Point3df bbmax1,
                          Point3df bbmin2,
                          Point3df bbmax2,
                          uint *no_overlap ){

    float overlap_x, overlap_y, aux;
    double rec = 0.0;

    if ( sqrt(pow( bbmin1[0] - bbmax1[0], 2) ) < 0.0001 )  bbmax1[0] += 0.5; 
    if ( sqrt(pow( bbmin1[1] - bbmax1[1], 2) ) < 0.0001 )  bbmax1[1] += 0.5; 
    if ( sqrt(pow( bbmin2[0] - bbmax2[0], 2) ) < 0.0001 )  bbmax2[0] += 0.5; 
    if ( sqrt(pow( bbmin2[1] - bbmax2[1], 2) ) < 0.0001 )  bbmax2[1] += 0.5; 
  
    if (sqrt(pow(bbmin1[1] - bbmax1[1], 2)) > 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) < 300) {

        if ( 360 - bbmax2[1] < bbmin2[1] ) {
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

    else if (sqrt(pow( bbmin1[1] - bbmax1[1], 2)) < 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) > 300) {

        if ( 360 - bbmax1[1] < bbmin1[1] ) {
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
    else if (sqrt(pow( bbmin1[1] - bbmax1[1], 2)) > 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) > 300) {

        aux = bbmin1[1];
        bbmin1[1] = bbmax1[1] - 360.0;
        bbmax1[1] = aux;
        aux = bbmin2[1];
        bbmin2[1] = bbmax2[1] - 360.0;
        bbmax2[1] = aux;
    }
  
    // ON S'OCCUPE DE LA LATITUDE
    if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) > 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) < 150) {

        if ( 180 - bbmax2[0] < bbmin2[0] ) {
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
    else if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) < 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) > 150){

        if ( 180 - bbmax1[0] < bbmin1[0] ) {
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

    else if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) > 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) > 150){

        aux = bbmin1[0];
        bbmin1[0] = bbmax1[0] - 360.0;
        bbmax1[0] = aux;
        aux = bbmin2[0];
        bbmin2[0] = bbmax2[0] - 360.0;
        bbmax2[0] = aux;
    }
    
    // PRÉTRAITEMENTS EFFECTUÉS ON CALCULE LE RECOUVREMENT

    *no_overlap = 0;
    if ( bbmin1[0] <= bbmin2[0] )
        if ( bbmax1[0] < bbmin2[0] )
            *no_overlap = 1;
        else
            overlap_x = ( bbmax2[0] < bbmax1[0] ? bbmax2[0] : bbmax1[0] ) - bbmin2[0];
    else
        if ( bbmax2[0] < bbmin1[0] )
            *no_overlap = 1;
        else
            overlap_x = ( bbmax1[0] < bbmax2[0] ? bbmax1[0] : bbmax2[0] ) - bbmin1[0];

    if ( *no_overlap == 0 ) {

        if ( bbmin1[1] <= bbmin2[1] )
            if ( bbmax1[1] < bbmin2[1] )
                *no_overlap = 1;
            else
                overlap_y = ( bbmax2[1] < bbmax1[1] ? bbmax2[1] : bbmax1[1] ) - bbmin2[1];
        else
            if ( bbmax2[1] < bbmin1[1] )
                *no_overlap = 1;
            else
                overlap_y = ( bbmax1[1] < bbmax2[1] ? bbmax1[1] : bbmax2[1] ) - bbmin1[1];
        if ( *no_overlap == 0 ) {
            rec = overlap_x * overlap_y;
            double div=  ( bbmax1[0] - bbmin1[0] ) * ( bbmax1[1] - bbmin1[1] )
                + ( bbmax2[0] - bbmin2[0] ) * ( bbmax2[1] - bbmin2[1] ) ;

            rec = 2 * rec / div;
        }
    }

    return rec;
}





//##############################################################################

void filteringBlobs (  vector<surf::ScaleSpaceBlob *> & ssblobs,
                       vector<surf::GreyLevelBlob *> &filteredBlobs,
                       vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
                       Point2df bbmin,
                       Point2df bbmax  ){
                           
    // Filtering according to positions
    Point3df bbmin2( bbmin[0], bbmin[1], 0.0 ), bbmax2( bbmax[0], bbmax[1], 0.0 );
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
                glb->raw_coordinates = (*itB1)->raw_coordinates;
                

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

void filteringBlobs (  vector<surf::ScaleSpaceBlob *> & ssblobs,
                       vector<surf::GreyLevelBlob *> &filteredBlobs,
                       vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
                       set<int> &nodes  ) {
    
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
            
            set<int>::iterator it;
            set<int> intersection;
            for ( it = (*itB1)->nodes.begin() ; it != (*itB1)->nodes.end() ; it++ )
                if ( nodes.find(*it) != nodes.end() )
                    intersection.insert(*it);
            
            if ( intersection.size() != 0 ){
                surf::GreyLevelBlob *glb;
                glb = new surf::GreyLevelBlob();
                glb->t = (*itB1)->t;
                glb->scale = (*itB1)->scale;
                glb->boundingbox_max = (*itB1)->boundingbox_max;
                glb->boundingbox_min = (*itB1)->boundingbox_min;
                glb->nodes = (*itB1)->nodes;
                glb->coordinates = (*itB1)->coordinates;
                glb->raw_coordinates = (*itB1)->raw_coordinates;
                glb->index = (*itB1)->index;
                
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


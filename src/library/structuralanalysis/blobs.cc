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


Point3df Point3dfOnSphere ( float radius,
                            float lat,
                            float lon ){

    return  Point3df ( log(radius) * cos((lat-90.)/180.0*3.1415957) * cos(lon/180.0*3.1415957),
                    log(radius) * cos((lat-90.)/180.0*3.1415957) * sin(lon/180.0*3.1415957),
                    log(radius) * sin((lat-90.)/180.0*3.1415957) );
}

Point3df Point3dfOnMesh ( vector<float> &coordinates, float radius = 1.0) {

    return Point3df( log(radius) * coordinates[0],
                     log(radius) * coordinates[1],
                     log(radius) * coordinates[2] );
}

Point3df Point3dfOnPlane (  float height,
                            float lat,
                            float lon ){
    return Point3df(lat, lon, height );
}


void surf::Blob::getAimsMesh ( AimsSurface<3, Void> &inMesh) {

    this->mesh = AimsSurface<3,Void>();
	vector<int> gyrusVertices, corres;
	set<int>::iterator it;
		// on extrait un gyrus, le maillage a moins de vertex que l'hemisphere, du coup on cree un vecteur qui renseigne
		// sur les correspondances entre points homologues..s
	set<uint> gyrusSet;
	corres = *(new vector<int>(inMesh.vertex().size()));
	uint i = 0;
	for ( it = nodes.begin() ; it != nodes.end() ; it ++ ) {
	  gyrusSet.insert(*it);
	  this->mesh.vertex().push_back(inMesh.vertex()[*it]);
	  corres[*it] = i++;
	}
	for ( uint i = 0 ; i < inMesh.polygon().size() ; i++ ) {
	  if (gyrusSet.find(inMesh.polygon()[i][0])!=gyrusSet.end() &&
			gyrusSet.find(inMesh.polygon()[i][1])!=gyrusSet.end() &&
			gyrusSet.find(inMesh.polygon()[i][2])!=gyrusSet.end()){
	          this->mesh.polygon().push_back(AimsVector<uint,3>(corres[inMesh.polygon()[i][0]],
			   corres[inMesh.polygon()[i][1]],
			   corres[inMesh.polygon()[i][2]]));
			}

	}
}

//void surf::GreyLevelBlob::getAimsMesh ( AimsSurface<3, Void> &mesh ){
//    surf::Blob::getAimsMesh ( mesh );
//    assert( this->mesh.vertex().size() > 0 );
//    if ( this->mesh.polygon().size() == 0 ) {
////        cout << " sphere " << endl;
//        AimsSurfaceTriangle *sph;
//        sph = SurfaceGenerator::sphere(this->mesh.vertex()[0], 0.9, 10);
//        this->mesh = (*sph)[0];
//    }
//}

//void surf::ScaleSpaceBlob::getAimsMesh ( AimsSurface<3, Void> &inMesh ){
//    this->mesh = AimsSurface<3,Void>();
//    vector<int> gyrusVertices, corres;
//    set<int>::iterator it;
//
//        // on extrait un gyrus, le maillage a moins de vertex que l'hemisphere, du coup on cree un vecteur qui renseigne
//        // sur les correspondances entre points homologues..s
//    set<uint> gyrusSet;
//    corres = *(new vector<int>(inMesh.vertex().size()));
//    uint i = 0;
//    assert(blobs.size() == 1);
//    set<int> nodes2 ( (*(blobs.begin()))->nodes );
//    for ( it = nodes2.begin() ; it != nodes2.end() ; it ++ ) {
//      gyrusSet.insert(*it);
//      this->mesh.vertex().push_back(inMesh.vertex()[*it]);
//      corres[*it] = i++;
//    }
//    for ( uint i = 0 ; i < inMesh.polygon().size() ; i++ ) {
//      if (gyrusSet.find(inMesh.polygon()[i][0])!=gyrusSet.end() &&
//            gyrusSet.find(inMesh.polygon()[i][1])!=gyrusSet.end() &&
//            gyrusSet.find(inMesh.polygon()[i][2])!=gyrusSet.end()){
//              this->mesh.polygon().push_back(AimsVector<uint,3>(corres[inMesh.polygon()[i][0]],
//               corres[inMesh.polygon()[i][1]],
//               corres[inMesh.polygon()[i][2]]));
//            }
//
//    }
//    assert( this->mesh.vertex().size() > 0 );
//    if ( this->mesh.polygon().size() == 0 ) {
////        cout << " sphere " << endl;
//        AimsSurfaceTriangle *sph;
//        sph = SurfaceGenerator::sphere(this->mesh.vertex()[0], 0.9, 10);
//        this->mesh = (*sph)[0];
//    }
//}

//##############################################################################


void surf::Blob::getAimsEllipsoid ( float abscissa, float depth, float height, float area ) {
    AimsSurfaceTriangle *ellipse;

    Point3df p1(0.0, 0.0, 0.0);

    ellipse = SurfaceGenerator::ellipse(p1, height, area, 10);
    for ( uint i = 0 ; i < (*ellipse)[0].vertex().size() ; i ++ ) {
        (*ellipse)[0].vertex()[i][0] += abscissa*10.0;
        (*ellipse)[0].vertex()[i][1] += log(height)*100.0;
        (*ellipse)[0].vertex()[i][2] += depth*10.0;
    }
    this->mesh = (*ellipse)[0];
}

void surf::GreyLevelBlob::getAimsEllipsoid ( void ) {

    set<int>::iterator it;
    float moyX = 0.0, moyY = 0.0;
    int nbX = 0, nbY = 0;

    for ( it = nodes.begin() ; it != nodes.end() ; it ++ ) {

        ASSERT ( coordinates[*it].size() > 0);
        if ( coordinates[*it][0] >= 0.0 ) {
            moyX += coordinates[*it][0];
            nbX++;
        }
        if ( coordinates[*it].size() == 2 && coordinates[*it][1] >= 0.0 ) {
            moyY += coordinates[*it][1];
            nbY++;
        }
    }

    it = nodes.begin();
    if (nbX > 0)
        moyX /= nbX;
    else
        moyX = 0.0;
    if (nbY > 0)
        moyY /= nbY;
    else
        moyY = 0.0;



    surf::Blob::getAimsEllipsoid ( moyX, moyY, scale, 2.0 );
}


void surf::Blob::moveMeshToSphericalAtlas ( float radius ) {

    if ( mesh.vertex().size() > 0 ) {

        if ( mesh.polygon().size() > 0 ) {
            // If at least one triangle exists in the blob mesh

            for ( uint i = 0 ; i < mesh.vertex().size() ; i++ )
                mesh.vertex()[i] = Point3dfOnSphere( radius, coordinates[i][0], coordinates[i][1]);

        }
        else {
            // If the blob mesh has no triangle, the representing mesh is a sphere
            //  centered at the first blob node
            AimsSurfaceTriangle *sphere;
            Point3df p ( Point3dfOnSphere( radius, coordinates[0][0], coordinates[0][1]) );
//            for ( uint i = 0 ; i < 3 ; i ++ )
//                p[i] *= log( radius );
            sphere = SurfaceGenerator::sphere(p, 0.2, 10);
            mesh = (*sphere)[0];

        }

    }
}

void surf::GreyLevelBlob::moveMeshToSphericalAtlas ( ) {
    return surf::Blob::moveMeshToSphericalAtlas ( scale + 1.0 );
}

void surf::Blob::moveMeshToPlaneAtlas ( float height ) {

    if ( mesh.vertex().size() > 0 ) {

        if ( mesh.polygon().size() > 0 ) {
            // If at least one triangle exists in the blob mesh

            for ( uint i = 0 ; i < mesh.vertex().size() ; i++ )
                mesh.vertex()[i] = Point3dfOnPlane( log(height)*100.0, coordinates[i][0], coordinates[i][1]);

        }
        else {
            // If the blob mesh has no triangle, the representing mesh is a sphere
            //  centered at the first blob node
            AimsSurfaceTriangle *sphere;
            Point3df p ( Point3dfOnPlane( log(height)*100.0, coordinates[0][0], coordinates[0][1]) );
            sphere = SurfaceGenerator::sphere(p, 0.2, 10);
            mesh = (*sphere)[0];

        }

    }

}

void surf::GreyLevelBlob::moveMeshToPlaneAtlas ( ) {
    return surf::Blob::moveMeshToPlaneAtlas ( scale + 1.0 );
}
//##############################################################################

void surf::Blob::getAimsSphereAtMaxNode (  AimsSurface<3, Void> &mesh, Texture<float> &tex ) {
    set<int>::iterator it;
    int maxim_node = getMaximumNode(tex);
    AimsSurfaceTriangle *sph;
//    sph = SurfaceGenerator::sphere(mesh.vertex()[maxim_node], 0.9, 10);
    Point3df p(raw_coordinates[maxim_node][0],raw_coordinates[maxim_node][1],raw_coordinates[maxim_node][2]); 
    sph = SurfaceGenerator::sphere(p, 0.9, 10);
    this->mesh = (*sph)[0];

}



//##############################################################################
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

Point3df surf::Blob::getBlobBarycenterOnASphere( void ) {
    surf::Blob *glb;
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
    return Point3df(  cos((latMoy-90.)/180.0*3.1415957) * cos(lonMoy/180.0*3.1415957),
                 cos((latMoy-90.)/180.0*3.1415957) * sin(lonMoy/180.0*3.1415957),
                 sin((latMoy-90.)/180.0*3.1415957) );
}

Point3df surf::Blob::getBlobBarycenter( void ) {
    surf::Blob *glb;
    glb = this;
    float xMoy = 0.0, yMoy = 0.0, zMoy = 0.0;
    ASSERT(glb->nodes.size() == glb->raw_coordinates.size());

    set<int>::iterator it;
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
    return Point3df(  xMoy, yMoy, zMoy );
}


// PAS VERIFIEE

Point3df surf::Blob::getBlobBarycenterOnAPlane( void ) {
    surf::Blob *glb;
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
    return Point3df(  latMoy, lonMoy, 1.0 );
}

Point3df surf::Blob::getBlobBarycenterFromMesh( void ) {
    assert(mesh.vertex().size() > 0);
    Point3df bc(0.0, 0.0, 0.0);
    for ( uint i = 0 ; i < mesh.vertex().size() ; i++ )
        bc += mesh.vertex()[i];
    bc /= mesh.vertex().size();
    return bc;
}

//##############################################################################



double getOverlapMeasure( Point2df bbmin1,
                          Point2df bbmax1,
                          Point2df bbmin2,
                          Point2df bbmax2,
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



bool isInside2DBox( Point2df p1, Point2df bbmin, Point2df bbmax) {
    uint no_overlap = 2;

    Point2df bbmin1 (p1[0] - 0.0001, p1[1] - 0.0001),
        bbmax1 (p1[0] + 0.0001, p1[1] + 0.0001),
        bbmin2 (bbmin[0], bbmin[1]),
        bbmax2 (bbmax[0], bbmax[1]);

    getOverlapMeasure( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
    if (no_overlap == 0)
        return true;
    else if (no_overlap == 1)
        return false;
    assert(false);
    return false;
}


void filteringBlobs (  vector<surf::ScaleSpaceBlob *> & ssblobs,
                       vector<surf::GreyLevelBlob *> &filteredBlobs,
                       vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
                       set<int> &nodes ){

    for (uint i = 0 ; i < ssblobs.size() ; i++ )
        ssblobs [i] -> index = i;

    set< uint > filteredIndices;

    // Filtering according to positions
    for (uint i = 0 ; i < ssblobs.size() ; i++ ){
        string subject;
        subject = ssblobs[i]->subject;
        bool firstGLB = false;

        surf::ScaleSpaceBlob *ssb;
        ssb = new surf::ScaleSpaceBlob( ssblobs[i] );
        ssb->blobs.clear();
        ssb->topBlobs.clear();
        ssb->bottomBlobs.clear();

        set<surf::GreyLevelBlob *>::iterator itB1;
        for (itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() && !firstGLB ; itB1++){
            uint no_overlap = 2;

            set<int>::iterator it;
            set<int> intersection;
            for ( it = (*itB1)->nodes.begin() ; it != (*itB1)->nodes.end() ; it++ )
                if ( nodes.find(*it) != nodes.end() )
                    intersection.insert(*it);

            if ( intersection.size() != 0 )
                no_overlap = 0;
            else
                no_overlap = 1;
//             }
            if (no_overlap == 0)
                firstGLB = true;
        }

        if (firstGLB) {
            for ( itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() ; itB1++ ) {

                surf::GreyLevelBlob *glb;
                glb = new surf::GreyLevelBlob( *itB1 );
                ASSERT(glb->nodes.size() == glb->raw_coordinates.size());
                glb->ssb_parent = ssb;
                ssb->blobs.insert(glb);
                filteredBlobs.push_back(glb);
            }
            filteredSsblobs.push_back(ssb);
            filteredIndices.insert ( ssb->index );
        }
        else
            delete(ssb);
    }
    cerr << filteredBlobs.size() << " filtered blobs - " << filteredSsblobs.size() << " filtered ssblobs" << endl;


    // Now that the blobs are filtered, we add the correct bifurcations

    for ( uint i = 0 ; i < filteredSsblobs.size() ; i ++ ) {

        set<surf::ScaleSpaceBlob *> auxTop = ssblobs[filteredSsblobs[i]->index]->topBlobs;
        set<surf::ScaleSpaceBlob *>::iterator it;
        for ( it = auxTop.begin() ; it != auxTop.end() ; it ++ ) {
            if (filteredIndices.find((*it)->index) != filteredIndices.end()) {
                uint i1 = 0;
                for ( ; i1 < filteredSsblobs.size() && filteredSsblobs[i1]->index != (*it)->index ; i1 ++ ) {}
                filteredSsblobs[i]->topBlobs.insert( filteredSsblobs[i1] );
            }
        }

        set<surf::ScaleSpaceBlob *> auxBot = ssblobs[filteredSsblobs[i]->index]->bottomBlobs;
        for ( it = auxBot.begin() ; it != auxBot.end() ; it ++ ) {
            if (filteredIndices.find((*it)->index) != filteredIndices.end()){
                uint i1 = 0;
                for ( ; i1 < filteredSsblobs.size() && filteredSsblobs[i1]->index != (*it)->index ; i1 ++ ) {}
                filteredSsblobs[i]->bottomBlobs.insert(filteredSsblobs[i1]);
            }
        }

    }

}

// #############################################################################


pair<Point2df, Point2df> surf::Blob::get2DBoundingBox ( void ){
  Point2df bbmin, bbmax;
  bbmin[0] = 181.0;
  bbmin[1] = 361.0;
  bbmax[0] = -1.0;
  bbmax[1] = -1.0;
  set<int>::iterator it;
  pair<Point2df, Point2df> bb;
  for (it = nodes.begin() ; it != nodes.end() ; it ++){
    if (coordinates[*it][0] < bbmin[0])
      bbmin[0] = coordinates[*it][0];
    if (coordinates[*it][1] < bbmin[1])
      bbmin[1] = coordinates[*it][1];

    if (coordinates[*it][0] > bbmax[0])
      bbmax[0] = coordinates[*it][0];
    if (coordinates[*it][1] > bbmax[1])
      bbmax[1] = coordinates[*it][1];
  }

  if ( bbmax[1] > 300.0 && bbmin[1] < 60.0 ) {
    for ( uint i = 0 ; i < nodes.size() ; i++ ) {
      if ( coordinates[*it][1] >300.0 && coordinates[*it][1] < bbmax[1])
        bbmax[1] = coordinates[*it][1];
      if (coordinates[*it][1] < 60.0 && coordinates[*it][1] > bbmin[1])
        bbmin[1] = coordinates[*it][1];
    }
  }

  bb.first = bbmin;
  bb.second = bbmax;
  return bb;
}


pair<Point2df, Point2df> surf::GreyLevelBlob::get2DBoundingBox ( void ) {
    return surf::Blob::get2DBoundingBox ( );
}

pair<Point2df, Point2df> surf::ScaleSpaceBlob::get2DBoundingBox ( void ) {
    return surf::Blob::get2DBoundingBox ( );
}



//##############################################################################

int surf::Blob::getMaximumNode( Texture<float> &tex ) {
	assert(nodes.size() > 0);
    int maximum = *(nodes.begin());
    set<int>::iterator it;
    for ( it = nodes.begin() ; it != nodes.end() ; it ++ ) {
        if ( tex.item(*it) > tex.item(maximum) )
            maximum = *it;
    }
    return maximum;
}


// PREVIOUS VERSION THE NEW ONE HAS BEEN COPIED FROM CORTICAL/SURF_REF/GYRI/MESH_OPERATIONS.CC

//void surf::Blob::getAimsMesh ( AimsSurface<3, Void> &mesh ) {
//
//    ASSERT ( mesh.vertex().size() != 0 ) ;
//
//    uint p1,p2,p3;
//
//    set<uint>::iterator it;
//    set<uint> tri, comp;
//    vector<uint> corres;
//
//    for (uint j = 0 ; j < mesh.polygon().size() ; j++){
//
//        p1=mesh.polygon()[j][0];
//        p2=mesh.polygon()[j][1];
//        p3=mesh.polygon()[j][2];
//
//        if ( nodes.find(p1) != nodes.end() &&
//            nodes.find(p2) != nodes.end() &&
//            nodes.find(p3) != nodes.end() )
//        tri.insert(j);
//    }
//
//    for (it=tri.begin();it!=tri.end();it++){
//        p1=mesh.polygon()[*it][0];
//        p2=mesh.polygon()[*it][1];
//        p3=mesh.polygon()[*it][2];
//        comp.insert(p1); comp.insert(p2); comp.insert(p3);
//    }
//
//    corres = vector<uint>(mesh.vertex().size());
//
//    for (it = comp.begin() ; it != comp.end() ; it++){
//        assert(*it<corres.size());
//        assert(*it<mesh.vertex().size());
//        this->mesh.vertex().push_back( Point3dfOnMesh( raw_coordinates[*it] ) );
//
//        corres[*it] = this->mesh.vertex().size()-1;
//    }
//
//    for (it=tri.begin();it!=tri.end();it++){
//        p1=mesh.polygon()[*it][0];
//        p2=mesh.polygon()[*it][1];
//        p3=mesh.polygon()[*it][2];
//        this->mesh.polygon().push_back(AimsVector<uint,3>(corres[p1],corres[p2],corres[p3]));
//    }
//
//}

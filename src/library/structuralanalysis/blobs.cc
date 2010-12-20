#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/blobs.h>
#include <time.h>
#include <aims/mesh/surfacegen.h>


using namespace aims;
using namespace carto;
using namespace std;


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

//##############################################################################

void surf::Blob::getAimsEllipsoid ( float abscissa, float height, float depth, float area ) {
    AimsSurfaceTriangle *ellipse;

    Point3df p1(0.0, 0.0, 0.0);

    ellipse = SurfaceGenerator::ellipse(p1,  depth / 5.0, area, 10);
    for ( uint i = 0 ; i < (*ellipse)[0].vertex().size() ; i ++ ) {
        (*ellipse)[0].vertex()[i][0] += abscissa;
        (*ellipse)[0].vertex()[i][1] += depth * 5.0;
        (*ellipse)[0].vertex()[i][2] += height;
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

    surf::Blob::getAimsEllipsoid ( moyX*15.0, 0.0 /*moyY*10.0*/, log(scale)*10.0, 2.0 );
}

void surf::GreyLevelBlob::getAimsEllipsoidAtMaxNode (  Texture<float> &tex ) {
    int maxim_node = getMaximumNode(tex);
    surf::Blob::getAimsEllipsoid ( coordinates[maxim_node][0]*15.0, 0.0 /*moyY*10.0*/, log(scale)*10.0, 2.0 );
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

void surf::Blob::getAimsSphereAtMaxNode (  Texture<float> &tex, float radius ) {
    set<int>::iterator it;
    assert(tex.nItem() > 0);
    int maxim_node = getMaximumNode(tex);
    AimsSurfaceTriangle *sph;
//    sph = SurfaceGenerator::sphere(mesh.vertex()[maxim_node], 0.9, 10);
    assert( nodes.size() == raw_coordinates.size() );
    Point3df p( raw_coordinates[maxim_node][0], raw_coordinates[maxim_node][1], raw_coordinates[maxim_node][2] );
    sph = SurfaceGenerator::sphere( p, radius, 10 );
    this->mesh = (*sph)[0];

}

void surf::Blob::getNodesFromBlob ( surf::Blob *b ){
	std::set<int>::iterator it;
	nodes = set<int>( b->nodes );
	for (it = nodes.begin() ; it != nodes.end() ; it++ ){
		raw_coordinates[*it] = std::vector<float>(b->raw_coordinates[*it]);
		coordinates[*it] = std::vector<float>(b->coordinates[*it]);
	}
	assert( nodes.size() == raw_coordinates.size() || raw_coordinates.size() == 0) ;
	assert( nodes.size() == coordinates.size() || coordinates.size() == 0);
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



// #############################################################################


pair<Point2df, Point2df> surf::Blob::get2DBoundingBox ( void ){
    Point2df bbmin, bbmax;
    bbmin[0] = 181.0;
    bbmin[1] = 361.0;
    bbmax[0] = -1.0;
    bbmax[1] = -1.0;
    std::set<int>::iterator it;
    std::pair<Point2df, Point2df> bb;

    for ( it = nodes.begin() ; it != nodes.end() ; it ++ ) {
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
        for ( it = nodes.begin() ; it != nodes.end() ; it ++ ) {
            if ( coordinates[*it][1] > 300.0 && coordinates[*it][1] < bbmax[1] )
                bbmax[1] = coordinates[*it][1];
            if (coordinates[*it][1] < 60.0 && coordinates[*it][1] > bbmin[1] )
                bbmin[1] = coordinates[*it][1];
        }
    }

    bb.first = bbmin;
    bb.second = bbmax;
    return bb;
}


std::pair<Point2df, Point2df> surf::GreyLevelBlob::get2DBoundingBox ( void ) {
    return surf::Blob::get2DBoundingBox ( );
}

std::pair<Point2df, Point2df> surf::ScaleSpaceBlob::get2DBoundingBox ( void ) {
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

//##############################################################################

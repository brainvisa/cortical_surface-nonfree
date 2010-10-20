
#include <aims/getopt/getopt2.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <cortical_surface/surfacereferential/gyri/mesh_operations.h>
#include <aims/distancemap/meshdistance.h>
#include <cortical_surface/structuralanalysis/region.h>

using namespace aims;
using namespace carto;
using namespace std;

void surf::Region::getGlobalFromLocalNodes ( surf::Blob *blob ) {
    set<int>::iterator it;
    set<int> nodes_aux (blob->nodes);
    map<int, vector<float> > coordinates_aux (blob->coordinates);
    map<int, vector<float> > raw_coordinates_aux (blob->raw_coordinates);

    blob->nodes.clear();
    blob->coordinates.clear();
    blob->raw_coordinates.clear();
    for ( it =  nodes_aux.begin() ; it != nodes_aux.end() ; it ++ ) {
        blob->nodes.insert( nodes[*it] );
        blob->coordinates[ nodes[*it]] = vector<float>(coordinates_aux[*it]);
        blob->raw_coordinates[ nodes[*it]] = vector<float>(raw_coordinates_aux[*it]);
    }
}

surf::Region::Region( SubjectData *subject, std::vector<uint> &gyrusVertices  ) {

    this->subject = subject;

    nodes = std::vector<uint>( gyrusVertices );

//    vector<uint> corres;
    std::map<unsigned, std::set< std::pair<unsigned,float> > > globalWeightLapl = AimsMeshWeightFiniteElementLaplacian ( *(subject->mesh), 0.98 );
	// Extracting the gyrus mesh
	AimsSurfaceTriangle gyrusMesh ( getGyrusMesh ( *(subject->mesh), gyrusVertices, corres ) ) ;

	regionMesh = AimsSurface<3,Void>( gyrusMesh[0] );
	weightLapl = getGyrusWeight( globalWeightLapl, gyrusVertices, corres);




//              regionMesh.vertex() = gyrusMesh[0].vertex();
//              regionMesh.polygon() = gyrusMesh[0].polygon();

//              // Writing the gyrus mesh on the disk
//              Writer<AimsSurfaceTriangle > w1(tex2ps->_meshPath + "Gyrus.mesh");
//              w1.write(gyrusMesh);

//              // Flattening the gyrus mesh and writing it on the disk
//              AimsSurfaceTriangle flatMesh;
//              TimeTexture<float> param(3, lat[0].nItem());
//              param[0] = tex[0];
//              param[1] = lat[0];
//              param[2] = lon[0];
//              flatMesh[0] = getFlatMesh( gyrusMesh[0], gyrusVertices, corres, param );
//              Writer<AimsSurfaceTriangle> wFlat(tex2ps->_meshPath + "GyrusFlat.mesh" );
//              wFlat.write(flatMesh);

                // Extracting the region associated to the gyrus from the input texture
//                regionTex = Texture<float> ( gyrusVertices.size());
//                regionLat = Texture<float> ( gyrusVertices.size());
//                regionLon = Texture<float> ( gyrusVertices.size());
//
//                for ( uint i = 0 ; i < gyrusVertices.size() ; i++ ){
//                    uint corr = gyrusVertices[i];
//                    regionTex.item(i) = tex->item(corr);
//                    regionLat.item(i) = lat->item(corr);
//                    regionLon.item(i) = lon->item(corr);
//                }
//              Writer<TimeTexture<float> > w(tex2ps->_texPath + "Gyrus.tex");
//              w.write(gyrusTex);

//              // Storing the original data in SubjectData original
//              globalMesh = AimsSurfaceTriangle (mesh);
//              globalTex = TimeTexture<float> (tex);
//              globalLat = TimeTexture<float> (lat);
//              globalLon = TimeTexture<float> (lon);

//              // Replacing the data by the gyrus-sized
//              mesh = AimsSurfaceTriangle (gyrusMesh);
//              tex = TimeTexture<float> (gyrusTex);
//              lat = TimeTexture<float> (gyrusLat);
//              lon = TimeTexture<float> (gyrusLon);
                cout << " ░ gyrus : " << gyrusVertices.size() << " nodes ░" << endl;


}

Texture<float> surf::Region::getGlobalTexture( float background_value ) {
                Texture<float> outtex ( nodes.size() );
                for ( uint i = 0 ; i < outtex.nItem() ; i ++ )
                    outtex.item(i) = background_value;
                for ( uint i = 0 ; i < nodes.size() ; i ++ )
                    outtex.item( nodes[i] ) = regionTex.item(i);

                return outtex;
}

Texture<float> surf::Region::getLocalFromGlobalTexture ( Texture<float> &corticalTex ) {
    Texture<float> outtex( nodes.size() );

    for ( uint i = 0 ; i < outtex.nItem() ; i ++ )
        outtex.item(i) = corticalTex.item( nodes[i] );

    return outtex;
}

Texture<float> surf::Region::getGlobalFromLocalTexture ( Texture<float> &gyrusTex, float background_value ) {
    Texture<float> outtex ( subject->mesh->vertex().size() );

    for ( uint i = 0 ; i < outtex.nItem() ; i ++ )
        outtex.item(i) = background_value;
    for ( uint i = 0 ; i < gyrusTex.nItem() ; i ++ )
        outtex.item( nodes[i] ) = gyrusTex.item(i);

    return outtex;
}

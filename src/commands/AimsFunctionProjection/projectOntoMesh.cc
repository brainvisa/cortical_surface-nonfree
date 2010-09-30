#include "projectOntoMesh.h"



std::vector< AimsData<float> > load_kernel ( std::string path ) {
    aims::Reader< AimsData<float> > rK(path);
    AimsData<float> data;
    std::vector< AimsData<float> > kernel;
    rK.read(data);
    int size = data.dimX();
    float vsizeX = data.sizeX();
    float vsizeY = data.sizeY();
    float vsizeZ = data.sizeZ();
    for ( int i = 0 ; i < data.dimT() ; i++ ) {
        AimsData<float> convol ( size, size, size, 1 );
        convol.setSizeXYZT ( vsizeX, vsizeY, vsizeZ, 1.0 );
      
        for ( int x = 0 ; x < size ; x++ )
            for ( int y = 0 ; y < size ; y++ )
                for ( int z = 0 ; z < size ; z++ ) {
                    convol ( x, y, z, 0 ) = data ( x, y, z, i );
//                sum_weight += convol(x,y,z,0);
                }
//       for (int x=0;x<size;x++)
//          for (int y=0;y<size;y++)
//             for (int z=0;z<size;z++){
//                convol(x,y,z,0) /= sum_weight;
//                //sum += convol(x,y,z,0);
//             }
        kernel.push_back ( convol );
    }
    std::cout << kernel.size() << std::endl;
    return kernel;
}



Texture<float> deconvolve ( AimsData<float> inFuncData, 
                            std::vector< AimsData<float> > & kernel, 
                            AimsSurfaceTriangle mesh ) {

    Texture<float> tex;
    int size = kernel[0].dimX();
    float vsizeX = kernel[0].sizeX();
    float vsizeY = kernel[0].sizeY();
    float vsizeZ = kernel[0].sizeZ();
    for ( uint i = 0 ; i < mesh[0].vertex().size() ; i++ ) {
        tex.push_back(0.0);
    }
    for ( uint i = 0 ; i < kernel.size() ; i++ ) {
        Point3df p ( mesh[0].vertex()[i] );
        Point3d nearest_voxel( (int) ( ( p[0] + ( vsizeX/2.0 ) ) / vsizeX ),
			    (int) ( ( p[1] + ( vsizeY/2.0 ) ) / vsizeY ),
			    (int) ( ( p[2] + ( vsizeZ/2.0 ) ) / vsizeZ) );
        for ( int x = -size/2, x0 = 0 ; x <= size/2 ; x++, x0++ ) {
            for ( int y = -size/2, y0 = 0 ; y <= size/2 ; y++, y0++ ) {
                for ( int z = -size/2, z0 = 0 ; z <= size/2 ; z++, z0++ ) {
                    Point3d vxl ( nearest_voxel );
                    vxl += Point3d ( x, y, z );
                    tex.item(i) += kernel[i] ( x0, y0, z0, 0 ) * inFuncData ( vxl[0], vxl[1], vxl[2], 0 );
                }
            }
        }
    }
    return tex;
}

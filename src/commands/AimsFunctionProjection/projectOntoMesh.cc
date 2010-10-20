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



TimeTexture<float> deconvolve ( AimsData<float> inFuncData, 
                            const std::vector< AimsData<float> > & kernel, 
                            AimsSurfaceTriangle mesh,
			    bool verbose ) {

    TimeTexture<float> tex( inFuncData.dimT(), mesh[0].vertex().size() );
    int size = kernel[0].dimX();
    float vsizeX = kernel[0].sizeX();
    float vsizeY = kernel[0].sizeY();
    float vsizeZ = kernel[0].sizeZ();
    
    for ( uint timepoint = 0 ; timepoint < inFuncData.dimT() ; timepoint++ ) {
	std::cout << timepoint << "/" << inFuncData.dimT() << std::endl;
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
			assert( x0 < size && y0 < size && z0 < size );
			//if (vxl[0]<0) {
			//    std::cout << vxl << std::endl;
			//}
			
			
			int Zoffset =  0; // - (int)(26.0/vsizeZ) ;

<<<<<<< .mine
			//assert( vxl[0] >= 0 );
			//assert( vxl[1] >= 0 );
			//assert( vxl[2] + Zoffset >= 0 );
			//assert( vxl[0] < inFuncData.dimX() );
			//assert( vxl[1] < inFuncData.dimY() );
			//if (vxl[2] + Zoffset >= inFuncData.dimZ() )
			//    std::cout << vxl[2] << " " << inFuncData.dimZ() << std::endl;
			//assert( vxl[2] + Zoffset < inFuncData.dimZ() );
			if ( vxl[0] >= 0 && 
				vxl[1] >= 0 &&
				vxl[2] >= 0 &&
				vxl[2] + Zoffset >= 0 &&
				vxl[0] < inFuncData.dimX() &&
				vxl[1] < inFuncData.dimY() &&
				vxl[2] + Zoffset < inFuncData.dimZ() ) {
			    
			    tex[timepoint].item(i) += kernel[i] ( x0, y0, z0, 0 ) * inFuncData ( vxl[0], vxl[1], vxl[2] + Zoffset, timepoint );
			}
			else {
			    if ( verbose )
				std::cout << "Problem accessing the volume data (" << vxl << ")" << std::endl;
			}
		    }
		}
	    }
	}		
=======
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
                    assert( x0 < size && y0 < size && z0 < size );
                    if (vxl[0]<0) {
                        std::cout << vxl << std::endl;
                    }
                    
                    
                    int Zoffset =  0; 

                    // introduction of an offset in Z to cope with partial brain acquisitions
                                        
                    if ( vxl[0] >= 0 && 
                            vxl[1] >= 0 &&
                            vxl[2] >= 0 &&
                            vxl[2] + Zoffset >= 0 &&
                            vxl[0] < inFuncData.dimX() &&
                            vxl[1] < inFuncData.dimY() &&
                            vxl[2] + Zoffset < inFuncData.dimZ() ) {
                        
                        tex.item(i) += kernel[i] ( x0, y0, z0, 0 ) * inFuncData ( vxl[0], vxl[1], vxl[2] + Zoffset, 0 );
                    }
                    else {
                        std::cout << "Problem accessing the volume data" << std::endl;
                    }
                }
            }
        }
>>>>>>> .r39575
    }
    return tex;
}

#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include "createKernels.h"
#include "projectOntoMesh.h"

using namespace aims;
using namespace carto;

int main ( int argc, const char **argv ) {
    try
    {
        std::string outpath, 
                    meshpath, 
                    datapath, 
                    datapath1,
                    mixtype;
        int operation, 
            size = 7, 
            type = 0, 
            index = -1;
        float vsizeX = 3.0,
              vsizeY = 3.0,
              vsizeZ = 3.0, 
              geod_decay = 5.0, 
              norm_decay = 2.0;
        AimsApplication app ( argc, argv, "AimsFunctionProjection : first computes \
            anatomically-informed kernels from one anatomy and uses them to project\
             some functional data onto a cortical mesh" );
        app.addOption( operation, "-op", "0 : computes convolution kernels from one anatomy ; 1 : projects functional volumes onto the surface (using kernels)");
        app.addOption( meshpath , "-m", "Grey/white matter mesh (.mesh)" );
        app.addOption( datapath , "-d", "Convolution kernels (.ima) to be used for projection (-op=1)", "" );
        app.addOption( datapath1 , "-d1", "4D functional volume (.ima/.img) to project onto the mesh (-op=1)", "" );
        app.addOption( size, "-i", "Size of computed kernels (integer)", 7);
        app.addOption( vsizeX, "-vx", "X-resolution of kernels (float)", 3.0);
        app.addOption( vsizeY, "-vy", "Y-resolution of kernels (float)", 3.0);
        app.addOption( vsizeZ, "-vz", "Z-resolution of kernels (float)", 3.0);
        app.addOption( geod_decay, "-g", "Geodesic decay (in mm;default = 5.0) ", 5.0);
        app.addOption( norm_decay, "-n", "Normal decay (in mm;default = 2.0)", 2.0);
        app.addOption( type , "-t", "For computing convolution kernels (-op=0), selects the cortical thickness evaluation method : 0 for 3mm constant (only for now) ", 1); //; 1 for cortical mask-based evaluation ; 2 if using thickness texture\n  For projecting volumes onto the surface (-op=1), sets the output type : 0 for single volume projection ; 1 for several volumes -> several textures ; 2 for several volumes -> one timetexture", 0 );
        app.addOption( outpath , "-o", "Output file : convolution kernels (-op=0) or projection texture (-op=1)" );
        app.addOption( index, "-I", "[DEBUG] Index of a precise kernel to be computed", 1);
        app.addOption( mixtype , "-mt", "Mix type : label (majority) or level (default), used only in projection mode (-op=1)", true );
        app.initialize();

        if ( operation == 0 ) {
            // calcul des noyaux de convolution
            // besoin : un meshpath, un entier, un flottant, un outpath
            AimsSurfaceTriangle mesh;
            Reader<AimsSurfaceTriangle> r1 ( meshpath );
            r1.read ( mesh );
            getMaskPath ( datapath );
            get_kernelindex ( index );
           
            Point3df vsize ( vsizeX, vsizeY, vsizeZ );
            AimsData<float> kernel ( fast_marching_kernels ( meshpath, size, vsize, geod_decay, norm_decay) );
//            float mini = 5.0, maxi = -5.0;
//            for ( uint x = 0 ; x < size ; x++ )
//                for ( uint y = 0 ; y < size ; y++ )
//                    for ( uint z = 0 ; z < size ; z++ ) {
//                        if ( kernel ( x, y, z, 3114 ) < mini )
//                            mini = kernel ( x, y, z, 3114 );
//                        if ( kernel ( x, y, z, 3114 ) > maxi )
//                            maxi = kernel ( x, y, z, 3114 );
//                    }
//            std::cout << std::endl << " warning : MINI= " << mini << " MAXI=" << maxi << " time=" << time << " node=" << kernel(3,3,3,3114) << std::endl;
//            std::cout << "milieu node avant= " << kernel(3,3,3,3110) << std::endl;
//            for (uint i= 0;i<3115;i++) 
//                std::cout << kernel(3,3,3,i) << " " << std::flush;
            
            Writer<AimsData<float> > w ( outpath );
            w.write ( kernel );
                       
        }
        else if ( operation == 1 ) {
            // cration de la texture de projection  partir de l'anatomie, des noyaux prcalculs et in fine d'un volume fonctionnel
            // besoin : un meshpath, un datapath1  projeter, un datapath (kernels), un outpath (texture)
            AimsSurfaceTriangle mesh;
            Reader < AimsSurfaceTriangle > r1 ( meshpath );
            r1.read ( mesh );
            std::vector< AimsData<float> > kernel = load_kernel( datapath );
         
            Reader < Volume<float> > r ( datapath1 );
            VolumeRef<float> funcdata;
            funcdata = r.read();

            std::cout << "kernel resolution : " << kernel[0].sizeX() << ";" << kernel[0].sizeY() << ";" << kernel[0].sizeZ() << std::endl;
            std::cout << "kernel size : " << kernel[0].dimX() << std::endl;
            std::cout << "volume resolution : " << funcdata->getVoxelSize()[0] << ";" << funcdata->getVoxelSize()[1] << ";" << funcdata->getVoxelSize()[2] << std::endl;
	    std::cout << "number of scans : " << funcdata.getSizeT() << std::endl;
	    
            TimeTexture<float> tex;

            if( mixtype == "label" )
              tex = deconvolve<Label>( funcdata, kernel, mesh );
            else
              tex = deconvolve<Level>( funcdata, kernel, mesh );
            //TimeTexture<float> ttex;
            //ttex[0] = tex;
            Writer < TimeTexture<float> > w ( outpath );
            w.write ( tex );
        }

    }
    catch ( user_interruption & )
    {
        return EXIT_FAILURE;
    }
    catch ( std::exception & e )
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
   
    return EXIT_SUCCESS;
}



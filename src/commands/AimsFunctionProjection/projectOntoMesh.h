#ifndef AIMS_FUNCTION_PROJECTION_H
#define AIMS_FUNCTION_PROJECTION_H

#include <aims/io/reader.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>

std::vector< AimsData<float> > load_kernel( std::string path );

TimeTexture<float> deconvolve ( AimsData<float> inFuncData, 
                            const std::vector< AimsData<float> > & kernel, 
                            AimsSurfaceTriangle mesh,
			    bool verbose = false );
 
#endif


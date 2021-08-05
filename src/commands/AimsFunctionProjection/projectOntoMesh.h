#ifndef AIMS_FUNCTION_PROJECTION_H
#define AIMS_FUNCTION_PROJECTION_H

#include <aims/io/reader.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <cartodata/volume/volume.h>
#include <aims/data/data.h>


enum MixType
{
  Level,
  Label
};


std::vector< AimsData<float> > load_kernel( std::string path );

template <MixType MT>
TimeTexture<float> deconvolve ( const carto::VolumeRef<float> inFuncData,
                                const std::vector< AimsData<float> > & kernel,
                                const AimsSurfaceTriangle mesh,
                                bool verbose = false );
 
#endif


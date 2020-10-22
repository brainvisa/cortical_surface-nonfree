#include "projectOntoMesh.h"
#include <cartodata/volume/volumeutil.h>

using namespace carto;
using namespace std;


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


template <MixType MT>
float store_convolve( vector<TimeTexture<float> > & tex, uint timepoint,
                      float w, const VolumeRef<float> & inFuncData,
                      Point3d pos,
                      uint i );


template <>
float store_convolve<Level>( vector<TimeTexture<float> > & tex, uint timepoint,
                             float w, const VolumeRef<float> & inFuncData,
                             Point3d pos, uint i )
{
  tex[0][timepoint].item(i) += w * inFuncData ( pos[0], pos[1], pos[2],
                                                timepoint );
}


template <>
float store_convolve<Label>( vector<TimeTexture<float> > & tex, uint timepoint,
                             float w, const VolumeRef<float> & inFuncData,
                             Point3d pos, uint i )
{
  int label = int( inFuncData ( pos[0], pos[1], pos[2], timepoint ) );
  tex[label][timepoint].item(i) += w;
//   if( w != 0. )
//     cout << label << ", " << i << ", " << tex[label][timepoint].item(i) << ": " << inFuncData ( pos[0], pos[1], pos[2], timepoint ) << " " << pos << endl;
}


template <MixType MT>
void alloc_tex( vector<TimeTexture<float> > & tex,
                const VolumeRef<float> & inFuncData,
                const AimsSurfaceTriangle mesh );


template <>
void alloc_tex<Level>( vector<TimeTexture<float> > & tex,
                       const VolumeRef<float> & inFuncData,
                       const AimsSurfaceTriangle mesh )
{
  tex.push_back( TimeTexture<float>() );
  tex[0] = TimeTexture<float>( inFuncData.getSizeT(),
                               mesh.begin()->second.vertex().size(), 0. );
}


template <>
void alloc_tex<Label>( vector<TimeTexture<float> > & tex,
                       const VolumeRef<float> & inFuncData,
                       const AimsSurfaceTriangle mesh )
{
  unsigned max_label = unsigned( max( inFuncData ) );
  tex.reserve( max_label +  1 );

  for( unsigned l=0; l<=max_label; ++l )
  {
    tex.push_back( TimeTexture<float>() );
    tex[l] = TimeTexture<float>( inFuncData.getSizeT(),
                                 mesh.begin()->second.vertex().size(), 0. );
  }
}


template <MixType MT>
void mix_results( vector<TimeTexture<float> > & tex );


template <>
void mix_results<Level>( vector<TimeTexture<float> > & tex )
{
  // nothing to do
}


template <>
void mix_results<Label>( vector<TimeTexture<float> > & tex )
{
  unsigned i, l, n = tex.size(), nv = tex[0][0].nItem(), t, nt = tex[0].size();
  float m;
  unsigned am;

  for( t=0; t<nt; ++t )
    for( i=1; i<nv; ++i )
    {
      m = 0.;
      am = 0.;
      for( l=0; l<n; ++l )
        if( tex[l][t][i] > m )
        {
          m = tex[l][t][i];
          am = l;
        }
//       cout << "vertex " << t << ", " << i << ": " << am << ", " << m << endl;
      tex[0][t][i] = am;
    }
}


template <MixType MT>
TimeTexture<float> deconvolve ( const VolumeRef<float> inFuncData,
                                const std::vector< AimsData<float> > & kernel,
                                const AimsSurfaceTriangle mesh,
                                bool verbose )
{
  vector<TimeTexture<float> > tex;
  alloc_tex<MT>( tex, inFuncData, mesh );
  int size = kernel[0].dimX();
  float vsizeX = kernel[0].sizeX();
  float vsizeY = kernel[0].sizeY();
  float vsizeZ = kernel[0].sizeZ();

  unsigned dimX = inFuncData.getSizeX();
  unsigned dimY = inFuncData.getSizeY();
  unsigned dimZ = inFuncData.getSizeZ();
  unsigned dimT = inFuncData.getSizeT();

  for ( uint timepoint = 0 ; timepoint < dimT ; timepoint++ )
  {
    std::cout << timepoint << "/" << dimT << std::endl;
    for ( uint i = 0 ; i < kernel.size() ; i++ )
    {
      Point3df p( mesh.vertex()[i] );
      Point3d nearest_voxel( (int) ( ( p[0] + ( vsizeX/2.0 ) ) / vsizeX ),
                            (int) ( ( p[1] + ( vsizeY/2.0 ) ) / vsizeY ),
                            (int) ( ( p[2] + ( vsizeZ/2.0 ) ) / vsizeZ) );
      for ( int x = -size/2, x0 = 0 ; x <= size/2 ; x++, x0++ )
      {
        for ( int y = -size/2, y0 = 0 ; y <= size/2 ; y++, y0++ )
        {
          for ( int z = -size/2, z0 = 0 ; z <= size/2 ; z++, z0++ )
          {
            Point3d vxl ( nearest_voxel );
            vxl += Point3d ( x, y, z );
            assert( x0 < size && y0 < size && z0 < size );

            int Zoffset =  0; // - (int)(26.0/vsizeZ) ;
            vxl[2] += Zoffset;
//             cout << x0 << ", " << y0 << ", " << z0 << ": " << kernel[i]( x0, y0, z0, 0 ) << ", " << vxl << endl;

            if ( vxl[0] >= 0 &&
                vxl[1] >= 0 &&
                vxl[2] >= 0 &&
                vxl[2] >= 0 &&
                vxl[0] < dimX &&
                vxl[1] < dimY &&
                vxl[2] < dimZ )
            {
              store_convolve<MT>( tex, timepoint, kernel[i]( x0, y0, z0, 0 ),
                                  inFuncData, vxl, i );
            }
            else
            {
              if ( verbose )
                cout << "Problem accessing the volume data (" << vxl
                     << ")" << std::endl;
            }
          }
        }
      }
    }
  }

  mix_results<MT>( tex );

  return tex[0];
}


template TimeTexture<float> deconvolve<Level>(
  const VolumeRef<float>, const vector< AimsData<float> > &,
  const AimsSurfaceTriangle, bool );
template TimeTexture<float> deconvolve<Label>(
  const VolumeRef<float>, const vector< AimsData<float> > &,
  const AimsSurfaceTriangle, bool );



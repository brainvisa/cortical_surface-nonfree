#ifndef AIMS_CREATE_KERNELS_H
#define AIMS_CREATE_KERNELS_H

#include <aims/data/data.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

template<class T>
      Texture<float> distance_map ( const AimsSurface<3,Void> &  mesh,
                                    const Texture<T> & inittex,
                                    bool allowUnreached,
                                    float limit);

inline float calcule_distance ( const Point3df p, 
                                const Point3df t );

inline float calcule_distance ( const Point3df p, 
                                const Point3d t );

inline float geod_weight_function ( float d, 
                                    float d0 );

inline float cortical_distance ( Point3d nv3, 
                                 Point3df vsize, 
                                 Point3df v, 
                                 Point3df n, 
                                 float width, 
                                 float seuil);

void getMaskPath ( std::string path );

float cortical_distance_via_tex ( Point3df v, short vertex, Point3d nv3, Point3df vsize, Point3df n );

inline std::vector<uint> nearest_vertices ( Point3df pf, 
                                            AimsSurfaceTriangle &mesh, 
                                            float rayon);

inline std::pair<int,float> plus_proche_point ( Point3df p, 
                                                AimsSurfaceTriangle &mesh);

struct ltstr;

void compute_kernel ( AimsData<float> &kernels, 
                      uint time, 
                      AimsSurfaceTriangle &mesh, 
                      uint node, 
                      const std::vector< std::map <uint, float> > &voisins2, 
                      AimsData<long> &vertex,  
                      AimsData<short> &classe, 
                      float &geod_decay, 
                      float &norm_decay, 
                      Point3df &vsize, 
                      uint size);

void get_kernelindex( int index );

AimsData<float> fast_marching_kernels ( std::string meshpath, 
                                        int size, 
                                        Point3df vsize, 
                                        float geod_decay, 
                                        float norm_decay );

#endif


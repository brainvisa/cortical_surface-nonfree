#include <cstdlib>
#include <float.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <cortical_surface/mesh/pointDistance.h>

using namespace aims;

float MeshPointDistance::compute(uint p1, uint p2)
{
     uint i;
     std::multimap<float,uint>              front1, front2;
     std::multimap<float,uint>              *cfront = &front1, *nfront = &front2, *tmpf;
     std::multimap<float,uint>::iterator    iv, fv;
     std::set<uint>                         neigh;
     std::set<uint>::iterator               in, fn;
     float                                  d, d2, l;
     Point3df                               pos;
     float                                  dist=0;

     front1.insert( std::pair<float,uint>( 0, p1 ) );

     if (p1 == p2)
          return(0.0);

     std::map<uint,float> distmap;
     std::map<uint,float>::iterator distit;
     distmap[p1]=0.0;
     int taille=0, ns=_mesh.vertex().size();

     while(  taille < ns )
     {
/*          std::cerr << "Taille : " << taille << std::endl;*/
          nfront->clear();
          neigh.clear();
/*          std::cerr << "\tFirst loop" << std::endl;*/
          
          for( iv=cfront->begin(), fv=cfront->end(); iv!=fv; ++iv )
          {
               i = (*iv).second;
               d = (*iv).first;
               for( in=_neigh[i].begin(), fn=_neigh[i].end(); in!=fn; ++in )
               {
                    distit = distmap.find(*in);
                    if (distit == distmap.end()) distmap[*in] = FLT_MAX;
                    d2 = distmap[*in];
                    pos = _mesh.vertex()[i] - _mesh.vertex()[*in];
                    l = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
                    if( d2 > d + l )
                    {
                         distmap[*in] = d+l;
                         neigh.insert( *in );
                         dist=d+l;
                    }
               }
          }
//           std::cerr << "\tSecond loop" << std::endl;

          for( in=neigh.begin(), fn=neigh.end(); in!=fn; ++in )
          {
               if ( *in == p2)
               {
/*                    std::cerr << "\tFound it" << std::endl;*/
                    front1.clear();
                    front2.clear();
//                     (*cfront).clear();
//                     (*nfront).clear();
//                     (*tmpf).clear();
//                     neigh.clear();
                    return distmap[*in];
               }
               else
                    nfront->insert( std::pair<float,unsigned>( distmap[*in], *in ) );
          }
                 
/*          std::cerr << "\tOK" << std::endl;*/
          tmpf = cfront;
          cfront = nfront;
          nfront = tmpf;
          taille=distmap.size();
     }

//      std::cerr << "Pb !!! MeshPointDistance::run(uint p1, uint p2) !!!" << std::endl;
//      std::cerr << "Points : " << p1 << " and " << p2 << ". Second point not found" << std::endl;
     exit(EXIT_FAILURE);
}


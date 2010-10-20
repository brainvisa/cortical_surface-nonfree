#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <float.h>
#include <cortical_surface/structuralanalysis/meshdistance.h>

using namespace aims;
using namespace std;
using namespace carto;

std::map<uint,float> getDistMap( AimsSurfaceTriangle *mesh,  map<unsigned, set<unsigned> >    &neighbours,  int dep, float dist_thresh){

  unsigned i;
  multimap<float,unsigned>    front1, front2;
  multimap<float,unsigned>    *cfront = &front1, *nfront = &front2, *tmpf;
  multimap<float,unsigned>::iterator    iv, fv;
  set<unsigned>                neigh;
  set<unsigned>::iterator        in, fn;
  float                    d, d2, l;
  Point3df                pos;
  float dist=0;
 
  front1.insert( pair<float,unsigned>( 0, dep ) );

  map<uint,float> distmap;
  map<uint,float>::iterator distit;
  distmap[dep]=0.0;
  if (dep != -1){
    while(  dist<=dist_thresh )
    {
//       cout << "TEST" << endl;
      nfront->clear();
      neigh.clear();
 
      for( iv=cfront->begin(), fv=cfront->end(); iv!=fv; ++iv )
      {
        i = (*iv).second;
        d = (*iv).first;
        for( in=neighbours[i].begin(), fn=neighbours[i].end(); in!=fn; ++in )
        {
          distit = distmap.find(*in);
          if (distit == distmap.end()) distmap[*in] = FLT_MAX;
          d2 = distmap[*in];
          pos = (*mesh).vertex()[i] - (*mesh).vertex()[*in];
          l = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
          if( d2 > d + l )
          {
            distmap[*in] = d+l;
            neigh.insert( *in );
            dist=d+l;
          }
        }
      }

      for( in=neigh.begin(), fn=neigh.end(); in!=fn; ++in )
        nfront->insert( pair<float,unsigned>( distmap[*in], *in ) );
 
      tmpf = cfront;
      cfront = nfront;
      nfront = tmpf;
    }
     

    front1.clear();
    front2.clear();
    (*cfront).clear();
    (*nfront).clear();
    (*tmpf).clear();
    neigh.clear();
//   tex.item(dep) = FLT_MAX;
  }
  return( distmap );
}

// vector<map<uint,float> > CalculeDistancesBlob(AimsSurfaceTriangle mesh, vector<uint> &sites){
//   // SURFACIQUE DISTMAP
//   map<unsigned, set<unsigned> >    neighbours;
//   unsigned v1b, v2, v3;
//  
//   for( uint i=0; i<mesh.polygon().size(); ++i ){
//     v1b = mesh.polygon()[i][0];
//     v2 = mesh.polygon()[i][1];
//     v3 = mesh.polygon()[i][2];
// 
//     neighbours[v1b].insert( v2 );
//     neighbours[v2].insert( v1b );
// 
//     neighbours[v1b].insert( v3 );
//     neighbours[v3].insert( v1b );
// 
//     neighbours[v2].insert( v3 );
//     neighbours[v3].insert( v2 );
//   }
//   
//   vector<map<uint,float> > distmap;
//   map<uint,float>::iterator it;
//   uint ns=mesh[0].vertex().size();
//   
//   for (uint i=0;i<ns;i++)
//     distmap.push_back(map<uint,float>());
//   for (uint i=0;i<ns;i++){
//     cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << ns << flush;
//     if (nodes.find(i) != nodes.end()){
//       map<uint,float> dmap(getDistMap(&mesh,&neighbours,i));
//       for (it=dmap.begin();it!=dmap.end();it++){
//         distmap[i][it->first]=it->second;
//         distmap[it->first][i]=it->second;
//       }
//     }
//   }
// 
//   return distmap;
// }

std::vector< std::map<uint,float> > CalculeCarteDistances(AimsSurfaceTriangle mesh, set<uint> nodes, float dist_thresh){
  // SURFACIQUE DISTMAP
  map<unsigned, set<unsigned> >    neighbours;
  unsigned v1b, v2, v3;
 
  for( uint i=0; i<mesh.polygon().size(); ++i ){
    v1b = mesh.polygon()[i][0];
    v2 = mesh.polygon()[i][1];
    v3 = mesh.polygon()[i][2];

    neighbours[v1b].insert( v2 );
    neighbours[v2].insert( v1b );

    neighbours[v1b].insert( v3 );
    neighbours[v3].insert( v1b );

    neighbours[v2].insert( v3 );
    neighbours[v3].insert( v2 );
  }
  
  vector<map<uint,float> > distmap;
  map<uint,float>::iterator it;
  uint ns=mesh[0].vertex().size();
  
  for (uint i=0;i<ns;i++)
    distmap.push_back(map<uint,float>());
  for (uint i=0;i<ns;i++){
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << ns << flush;
    if (nodes.find(i) != nodes.end()){
      map<uint,float> dmap(getDistMap(&mesh,neighbours,i,dist_thresh));
      for (it=dmap.begin();it!=dmap.end();it++){
        distmap[i][it->first]=it->second;
        distmap[it->first][i]=it->second;
      }
    }
  }

  return distmap;
}

std::map<float, std::vector<std::pair<float, uint> > > getAlternateMesh ( AimsSurfaceTriangle &mesh, 
                                                                          TimeTexture<float> &lat, 
                                                                          TimeTexture<float> &longit ) {
    std::set<uint> v;
    // LECTURE DU MAILLAGE "ATLAS" ET CREATION DE LA STRUCTURE ALTERNATIVE
    uint ns = mesh.vertex().size();
    float x, y;
    std::vector< std::set<uint> > polyVert(ns), 
                                polyVtmp(ns);
    std::map<float, std::vector<std::pair<float, uint> > > mesh2;
    for ( uint imesh = 0 ; imesh < ns ; imesh++ ) {
        x = lat[0].item ( imesh );
        y = longit[0].item ( imesh );
        mesh2[x].push_back ( std::pair<float,uint>(y,imesh) );
    }
    return mesh2;
}


std::map<uint, float> LocalMeshDistanceMap ( AimsSurface<3,Void> *mesh,
                                             const std::vector< std::set<unsigned> >    &neighbours,
                                             uint max_node,
                                             float ind ) {
    unsigned i;
    std::multimap<float,unsigned> front1, front2;
    std::multimap<float,unsigned> *cfront = &front1, *nfront = &front2, *tmpf;
    std::multimap<float,unsigned>::iterator    iv, fv;
    std::set<unsigned> neigh_local;
    std::set<unsigned>::const_iterator in, fn;
    float d, d2, l;
    Point3df pos;
    float dist = 0;
    std::map< uint, float > distances;
    distances[max_node] = 0.0;
    std::map<uint,float>::iterator it;
    for ( it = distances.begin() ; it != distances.end() ; it++ )
        if ( it->second == 0.0 )
            front1.insert( std::pair<float,unsigned>( 0, (*it).first ) );

    while( dist <= ind ) {
        nfront->clear();
        neigh_local.clear();

        for( iv = cfront->begin(), fv = cfront->end() ; iv != fv ; ++iv ) {

            i = (*iv).second;
            d = (*iv).first;
            for( in = neighbours[i].begin(), fn = neighbours[i].end() ; in != fn ; ++in ) {

                //d2 = tex.item( *in );
                if ( distances.find(*in) == distances.end() )
                    distances[*in] = FLT_MAX;
                d2 = distances[*in];

                pos = mesh->vertex()[i] - mesh->vertex()[*in];
                l = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );

                if( d2 > d + l ) {
                    //tex.item( *in ) = d + l;
                    distances[*in] = d+l;
                    neigh_local.insert( *in );

                    dist = distances [ *in ];
                }

            }

        }

        for ( in = neigh_local.begin(), fn = neigh_local.end() ; in != fn ; ++in )
            nfront->insert( std::pair<float,unsigned>( distances[*in], *in ) );

        tmpf = cfront;
        cfront = nfront;
        nfront = tmpf;
    }

    front1.clear();
    front2.clear();
    (*cfront).clear();
    (*nfront).clear();
    (*tmpf).clear();

    return ( distances );
}

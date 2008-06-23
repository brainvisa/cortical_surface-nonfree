#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <float.h>

using namespace aims;
using namespace std;
using namespace carto;

map<uint,float> getDistMap( AimsSurfaceTriangle *mesh,  map<unsigned, set<unsigned> >    *neighbours,  int dep){

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
    while(  dist<=20.0 )
    {
//       cout << "TEST" << endl;
      nfront->clear();
      neigh.clear();
 
      for( iv=cfront->begin(), fv=cfront->end(); iv!=fv; ++iv )
      {
        i = (*iv).second;
        d = (*iv).first;
        for( in=(*neighbours)[i].begin(), fn=(*neighbours)[i].end(); in!=fn; ++in )
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


vector<map<uint,float> > CalculeCarteDistances(AimsSurfaceTriangle mesh, set<uint> nodes){
  // SURFACIQUE DISTMAP
  
  map<unsigned, set<unsigned> >    neighbours;
  unsigned v1b, v2, v3;
 
  for( uint i=0; i<mesh.polygon().size(); ++i )
  {
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
  
  uint ns=mesh.vertex().size();
  for (uint i=0;i<ns;i++){
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << ns << flush;
    if (nodes.find(i) != nodes.end())
      distmap.push_back(getDistMap(&mesh,&neighbours,i));
    else
      distmap.push_back(getDistMap(&mesh,&neighbours,-1));
  }

  return distmap;
}

map<float, vector<pair<float, uint> > > getAlternateMesh(AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &longit){
  set<uint> v;
  // LECTURE DU MAILLAGE "ATLAS" ET CREATION DE LA STRUCTURE ALTERNATIVE
  uint ns=mesh.vertex().size();
  float x, y;
  vector<set<uint> > polyVert(ns), polyVtmp(ns);
  map<float, vector<pair<float, uint> > > mesh2;
  for (uint imesh=0;imesh<ns;imesh++) {
    
    x=lat[0].item(imesh);
    y=longit[0].item(imesh);
    mesh2[x].push_back(pair<float,uint>(y,imesh));
  }
  return mesh2;
}


/* Copyright (c) 1995-2007 CEA
 *
 *  This software and supporting documentation were developed by
 *      CEA/DSV/SHFJ
 *      4 place du General Leclerc
 *      91401 Orsay cedex
 *      France
 *
 * This software is governed by the CeCILL license version 2 under 
 * French law and abiding by the rules of distribution of free software.
 * You can  use, modify and/or redistribute the software under the 
 * terms of the CeCILL license version 2 as circulated by CEA, CNRS
 * and INRIA at the following URL "http://www.cecill.info". 
 *  
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 * 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license version 2 and that you accept its terms.
 */
#include <iostream>
#include <cstdlib>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>
#include <iomanip>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <aims/math/random.h>
#include <math.h>
#include <aims/distancemap/meshvoronoi.h>







using namespace aims;
using namespace carto;
using namespace std;
using namespace aims::meshdistance;

Texture<float> MeshDistance_adapt_tex(AimsSurface<3,Void> &mesh, const Texture<float> & inittex, bool allowUnreached, float ind )
{
  Texture<float> tex;
  TimeTexture<float> result_lim(1,mesh.vertex().size());
  std::vector< AimsVector<uint,3> > poly= mesh.polygon();
  unsigned                i, n = mesh.vertex().size();

  ASSERT( inittex.nItem() == n );
  tex.reserve( n );

 // neighbours map

  allowUnreached=true;

  std::map<unsigned, std::set<unsigned> >    neighbours;
  unsigned v1, v2, v3;

  for( i=0; i<poly.size(); ++i )
  {
    v1 = poly[i][0];
    v2 = poly[i][1];
    v3 = poly[i][2];
    if(inittex.item(v1)!=MESHDISTANCE_FORBIDDEN
       && inittex.item(v2)!=MESHDISTANCE_FORBIDDEN)
    {
      neighbours[v1].insert( v2 );
      neighbours[v2].insert( v1 );
    }
    if(inittex.item(v1)!=MESHDISTANCE_FORBIDDEN
       && inittex.item(v3)!=MESHDISTANCE_FORBIDDEN)
    {
      neighbours[v1].insert( v3 );
      neighbours[v3].insert( v1 );
    }
    if(inittex.item(v2)!=MESHDISTANCE_FORBIDDEN
       && inittex.item(v3)!=MESHDISTANCE_FORBIDDEN)         {
      neighbours[v2].insert( v3 );
      neighbours[v3].insert( v2 );
       }
  }

 // init texture

  for( i=0; i<n; ++i )
  {
    if( inittex.item(i) == 0 )
      tex.push_back( FLT_MAX );
    else if( inittex.item(i) == MESHDISTANCE_FORBIDDEN )
      tex.push_back( MESHDISTANCE_FORBIDDEN );
    else
      tex.push_back( 0 );
  }

  std::multimap<float,unsigned>    front1, front2;
  std::multimap<float,unsigned>    *cfront = &front1, *nfront = &front2, *tmpf;
  std::multimap<float,unsigned>::iterator    iv, fv;
  std::set<unsigned>                neigh_local;
  std::set<unsigned>::iterator        in, fn;
  float                    d, d2, l;
  Point3df                pos;
  float dist=0;

  for( i=0; i<n; ++i )
    if( tex.item(i) == 0 )
      front1.insert( std::pair<float,unsigned>( 0, i ) );

  while( dist<=ind )
  {
    nfront->clear();
    neigh_local.clear();

    for( iv=cfront->begin(), fv=cfront->end(); iv!=fv; ++iv )
    {
      i = (*iv).second;
      d = (*iv).first;
      for( in=neighbours[i].begin(), fn=neighbours[i].end(); in!=fn; ++in )
      {
        d2 = tex.item( *in );
        pos = mesh.vertex()[i] - mesh.vertex()[*in];
        l = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
        if( d2 > d + l )
        {
          tex.item( *in ) = d + l;
//                     result_lim[0].item( *in )=tex.item( *in );
          neigh_local.insert( *in );
          dist=tex.item( *in );
        }
      }
    }

    for( in=neigh_local.begin(), fn=neigh_local.end(); in!=fn; ++in )
      nfront->insert( std::pair<float,unsigned>( tex.item( *in ), *in ) );

    tmpf = cfront;
    cfront = nfront;
    nfront = tmpf;
  }
  neighbours.clear();
  front1.clear();
  front2.clear();
  (*cfront).clear();
  (*nfront).clear();
  (*tmpf).clear();
  return( tex );
} 

void graphe(vector<float> tab){
  
  vector<uint> votes;
  for (uint i=0;i<100;i++)
    votes.push_back(0);
  float tmin=1000000000000.0,tmax=-1000000000000.0;
  
  for (uint i=0; i<tab.size() ; i++)  {
    float valeur = tab[i];
    if (valeur<tmin) tmin=valeur;
    if (valeur>tmax) tmax=valeur;
  }
  printf("min: %f - max: %f\n", tmin, tmax);
  for (uint i=0; i<tab.size() ; i++)  {
    float valeur = tab[i];
//     printf("%d ", (uint) (((valeur - tmin) / (tmax-tmin))*votes.size()));
    votes[(uint) (((valeur - tmin) / (tmax-tmin))*votes.size())]++;
  }
  printf("\n");
  for (uint i=0;i<100;i++){
    int aux=0;
//     for (uint j=i;j<100;j++){
      aux += votes[i];
      
//     }
    printf("%.3f %d\n",tmin+i*(tmax-tmin)/100.0, aux);
  }
  printf("\n");
  
  
}

Texture<float> createSynthData(AimsSurfaceTriangle mesh, vector<pair<uint, float> > &sites, float smooth, float intensitynoise, float backgroundnoise, float simlocationnoise){
  cout << "test" << endl;
  vector<set<uint> >  voisins(SurfaceManip::surfaceNeighbours(mesh));
  Texture<float> result(mesh[0].vertex().size());
  float random;
  float distance = 1000.0;
  cout << "test" << endl;
  FiniteElementSmoother<3, float> *smoother;
  smoother=new FiniteElementSmoother<3, float>(0.05, &(mesh[0]));
  cout << "test" << endl;
    // bruit
  for (uint i=0;i<result.nItem();i++){
//     random= fabs(GaussianRandom(0.0,1.0))*0.0; 
    random =(float)GaussianRandom((float)backgroundnoise, (float)0.05);//s + ((float)UniformRandom() * 5.0 - 7.0) ;

//     cout << random << " ";
    result.item(i) = random;
  }
  cout << "testa" << endl;
  cout << smooth << endl;
//   lissage
  result=smoother->doSmoothing(result, ((int)((smooth)/0.05))*0.05);

//   for (uint i=0;i<result.nItem();i++)
//       result.item(i) += ((float)UniformRandom()*3.0 - 1.5);
//   for (uint i=0;i<result.nItem();i++){
//     random= ((float)UniformRandom() * backgroundnoise - backgroundnoise/2.0);//* 15.0 - 7.0) ;
//     result.item(i) += random;
//   }
//   result=smoother->doSmoothing(result, ((int)(smooth/0.05))*0.05);
// 
//   
  // texture où on crée les diracs, on les lisse, et qu'on fusionne ensuite à result
  Texture<float> tex(result.nItem());
  for (uint j=0;j<result.nItem();j++)
    tex.item(j)=0.0;
  
  // lissage de l'activation
  double smooth1 = smooth;
  
  Texture<float> tex1a,aux;   
  for (uint k=0;k<result.nItem();k++){
    aux.push_back(0.0);
  }
  uint backnode=0;
  

  // création des pics au niveau des sites
  for (uint i=0;i<sites.size();i++){
    cout << "site " << i  << endl;
    uint node=sites[i].first; 

    // creation de texture de distance avant seuillage
    aux.item(backnode)= 0.0;
    aux.item(node)=100.0;
    backnode=node;

    tex1a = MeshDistance_adapt_tex(mesh[0], aux, false,30.0);

    vector<uint> vois1;
    set<uint> vois2;
    for (uint k=0;k<result.nItem();k++){
      if (tex1a.item(k) < simlocationnoise+1 && tex1a.item(k)> simlocationnoise-0.01) vois1.push_back(k);
      if (tex1a.item(k) < 7.0) vois2.insert(k);
    }


    // introduction du bruit de similarité dans la position
    uint randomnode = (uint)(UniformRandom() * vois1.size());

    
    cout << "bruit de la position : " << randomnode << "/" << vois1.size() << "-(" << mesh[0].vertex()[vois1[randomnode]][0] << ";" << mesh[0].vertex()[vois1[randomnode]][1] << ";" << mesh[0].vertex()[vois1[randomnode]][2] << "/" << mesh[0].vertex()[sites[i].first][0] << ";" << mesh[0].vertex()[sites[i].first][1] << ";" << mesh[0].vertex()[sites[i].first][2] <<  ") " << sqrt(pow(mesh[0].vertex()[vois1[randomnode]][0]-mesh[0].vertex()[sites[i].first][0],2)+pow(mesh[0].vertex()[vois1[randomnode]][1]-mesh[0].vertex()[sites[i].first][1],2)+pow(mesh[0].vertex()[vois1[randomnode]][2]-mesh[0].vertex()[sites[i].first][2],2)) << endl;
    sites[i].first = vois1[randomnode];
     
    // on reparcourt le maillage pour trouver le noeud correspondant après changement
    node=vois1[randomnode];

    aux.item(backnode)= 0.0;
    aux.item(node)=100.0;
    backnode=node;
    tex1a = MeshDistance_adapt_tex(mesh[0], aux, false,30.0);

    vois1.clear();
    vois2.clear();
    for (uint k=0;k<result.nItem();k++){
      if (tex1a.item(k) < simlocationnoise) vois1.push_back(k);
      if (tex1a.item(k) < 10.0)  vois2.insert(k); 
    }
    
    // bruit de la tvalue max
    random = (UniformRandom() * intensitynoise) - intensitynoise/2.0;
    
//     cout << node << endl;
    cout << "tvalue2 : " << (sites[i].second) << endl;
    cout << "bruit de la tvalue : " << random << endl;
    // bruit de l'étendue
    set<uint>::iterator it,it2;
    
    cout << (sites[i].second + random)*sqrt(smooth1)/2.0*log(2.) << endl;
    for (it=vois2.begin();it!=vois2.end();it++){
      
      tex.item(*it) += (sites[i].second + random)*sqrt(8.0)/2.0*log(2.);
    }
    cout << "===" << endl;

  }
    // lissage des diracs
  tex=smoother->doSmoothing(tex, ((int)(200.0/0.05))*0.05);
  
  //fusion des deux textures
  for (uint j=0;j<result.nItem();j++)
//     if (result.item(j) < tex.item(j) && tex.item(j) > 1.0) result.item(j) = 0.3*result.item(j) + 0.7* tex.item(j);
    result.item(j) += tex.item(j);
  
  
  // introduction du bruit de similarité dans l'étendue (dans la valeur du lissage)
//   random = (UniformRandom() * simareanoise);

  return result;
}



// Texture<float> controldistance(AimsSurfaceTriangle mesh, pair<Point3df,float> site, float distance, float smooth, float intensitynoise, float backgroundnoise, float simlocationnoise, float simareanoise){
//   vector<pair<Point3df, float> > sites;
//   // creation de texture de distance avant seuillage
// 
//   Texture<float> tex1a,aux;   
//   for (uint k=0;k<mesh[0].vertex().size();k++){
//     aux.push_back(0.0);
//   }
//   uint node=0; float dist=1000.0;
//   for (uint j=0;j<mesh[0].vertex().size();j++){
//     float temp = sqrt(pow(mesh[0].vertex()[j][0]-site.first[0],2) + pow(mesh[0].vertex()[j][1]-site.first[1],2) +pow(mesh[0].vertex()[j][2]-site.first[2],2));
//     if (temp < dist) {
//       dist = temp;
//       node = j;
//     }
//   }
//   aux.item(node)= 100.0;
//   tex1a = MeshDistance(mesh[0], aux, false);         
//   vector<uint> vois1;
//   for (uint k=0;k<mesh[0].vertex().size();k++){
//     if (tex1a.item(k) < distance+0.05 && tex1a.item(k) > distance) vois1.push_back(k);
//   }
//   uint randomnode = (uint)(UniformRandom() * vois1.size());
//    
//   sites.push_back(site);
//   sites.push_back(pair<Point3df,float>(mesh[0].vertex()[vois1[randomnode]],10.0));
//   return createSynthData(mesh,sites,smooth,intensitynoise,backgroundnoise,simlocationnoise,simareanoise);
// }




int main( int argc, const char **argv )
{
  try {
  
    std::vector< std::size_t > startPoints;
    vector<std::size_t> intensities;

    std::string outpath = "";
    
   
    float noise;
    float location;
    float intensity;
    float smooth;
    Reader<AimsSurfaceTriangle> rmesh;
    AimsApplication app( argc, argv, "surfTexActivationSimulation" );
    app.addOption( rmesh, "-m", "mesh for testGeomap ");
    app.addOption( outpath, "-o", "output File");
    app.addOption( noise , "-n", "bruit de fond", true );
    app.addOption( location , "-l", "bruit de position", true );
    app.addOption( intensity , "-i", "bruit d'intensité", true );
    app.addOption( smooth , "-s", "lissage", true);
    app.addOptionSeries(startPoints, "-F", "vertex index list",false);
    app.addOptionSeries( intensities, "-I", "intensities", false);

//     AimsSurfaceTriangle *mesh;
//     mesh = SurfaceGenerator::sphere(Point3df(0.0,0.0,0.0) , 70.0, 22000 );
//     Writer<AimsSurfaceTriangle> w("/home/grg/data/simulations/sphere.mesh");
//     w.write(*mesh);
//     TimeTexture<float> lat(1,mesh[0].vertex().size()), lon(1,mesh[0].vertex().size());
//     for (uint i=0;i<mesh[0].vertex().size();i++){
//       Point3df p(mesh[0].vertex()[i]);
//       lat[0].item(i)= atan(p[1]/p[0])*57.2957795;
//       lon[0].item(i)= acos(p[2]/70.0)*57.2957795;
//     }
//     Writer<TimeTexture<float> > w1("/home/grg/data/simulations/latitude2.tex");
//     w1.write(lat);
//     Writer<TimeTexture<float> > w2("/home/grg/data/simulations/longitude2.tex");
//     w2.write(lon);

    app.initialize();
    cout << "SIMULATION" << endl;
    if ( startPoints.empty() )
      throw runtime_error( "bad startPoints input" );
   
    int startPoints_nb = ( int )startPoints.size();
    
    AimsSurfaceTriangle msh;
    rmesh.read(msh);
    
    for ( int i = 0; i < startPoints_nb; i++ )
    {
      std::cout << "input Vertex [" << i << "]: " << startPoints[i] << std::endl;
    }
    
    Writer<TimeTexture<float> > writer(outpath);
    TimeTexture<float> synthtex;
    vector<pair<uint,float> > sites; 
    
    for (uint i=0;i<startPoints.size();i++){
      pair<uint,float> newfocus;
      newfocus.first = startPoints[i];
      newfocus.second = intensities[i];
      sites.push_back(newfocus);
    }
        
      
      // float smooth, float intensitynoise, float backgroundnoise, float simlocationnoise, float simareanoise
    synthtex[0] = createSynthData(msh, sites, smooth, intensity, noise, location);
    writer.write(synthtex);
    return EXIT_SUCCESS;
  }
  catch( carto::user_interruption & )
  {
  // Exceptions thrown by command line parser (already handled, simply exit)
  }
  catch( exception & e )
  {
    cerr << e.what() << endl;
  }
}


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

#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <aims/math/random.h>
#include <aims/distancemap/meshdistance_d.h>






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

Texture<float> createSynthData(AimsSurfaceTriangle mesh, vector<pair<Point3df, float> > &sites, float smooth, float intensitynoise, float backgroundnoise, float simlocationnoise){
  vector<set<uint> >  voisins(SurfaceManip::surfaceNeighbours(mesh));
  Texture<float> result(mesh[0].vertex().size());
  float random;
  float distance = 1000.0;
  FiniteElementSmoother<3, float> *smoother;
  smoother=new FiniteElementSmoother<3, float>(0.05, &(mesh[0]));
  
    // bruit
  for (uint i=0;i<result.nItem();i++){
//     random= fabs(GaussianRandom(0.0,1.0))*0.0; //(float)GaussianRandom((float)backgroundnoise, (float)7.0);//s + ((float)UniformRandom() * 5.0 - 7.0) ;
// //     cout << random << " ";
    result.item(i) = 0.0; // random;
  }

  
  // lissage
//   result=smoother->doSmoothing(result, 0, ((int)(smooth/0.05))*0.05);

  
  for (uint i=0;i<result.nItem();i++){
    random= ((float)UniformRandom() * backgroundnoise - backgroundnoise/2.0);//* 15.0 - 7.0) ;
    result.item(i) += random;
  }
  result=smoother->doSmoothing(result, 0, ((int)(smooth/0.05))*0.05);

  
  // texture où on crée les diracs, on les lisse, et qu'on fusionne ensuite à result
  Texture<float> tex(result.nItem());
  for (uint j=0;j<result.nItem();j++)
    tex.item(j)=0.0;
  
  // lissage de l'activation
  double smooth1 = 16.0;
  
  Texture<float> tex1a,aux;   
  for (uint k=0;k<result.nItem();k++){
    aux.push_back(0.0);
  }
  uint backnode=0;
  
  // création des pics au niveau des sites
  for (uint i=0;i<sites.size();i++){
    
    uint node=0; distance =1000.0;
    for (uint j=0;j<result.nItem();j++){
      float temp = sqrt(pow(mesh[0].vertex()[j][0]-sites[i].first[0],2) + pow(mesh[0].vertex()[j][1]-sites[i].first[1],2) +pow(mesh[0].vertex()[j][2]-sites[i].first[2],2));
      if (temp < distance) {
        distance = temp;
        node = j;
      }
    }
    
    // creation de texture de distance avant seuillage
    aux.item(backnode)= 0.0;
    aux.item(node)=100.0;
    backnode=node;
    tex1a = MeshDistance_adapt_tex(mesh[0], aux, false,10.0);         


    vector<uint> vois1;
    set<uint> vois2;
    for (uint k=0;k<result.nItem();k++){
      if (tex1a.item(k) < simlocationnoise) vois1.push_back(k);
      if (tex1a.item(k) < 7.0) vois2.insert(k);
    }
    
    // introduction du bruit de similarité dans la position
    uint randomnode = (uint)(UniformRandom() * vois1.size());
    cout << "bruit de la position : " << randomnode << "/" << vois1.size() << "-(" << mesh[0].vertex()[vois1[randomnode]][0] << ";" << mesh[0].vertex()[vois1[randomnode]][1] << ";" << mesh[0].vertex()[vois1[randomnode]][2] << "/" << sites[i].first[0] << ";" << sites[i].first[1] << ";" << sites[i].first[2] <<  ") " << sqrt(pow(mesh[0].vertex()[vois1[randomnode]][0]-sites[i].first[0],2)+pow(mesh[0].vertex()[vois1[randomnode]][1]-sites[i].first[1],2)+pow(mesh[0].vertex()[vois1[randomnode]][2]-sites[i].first[2],2)) << endl;
    sites[i].first = mesh[0].vertex()[vois1[randomnode]];
    
    // on reparcourt le maillage pour trouver le noeud correspondant après changement
    node=0; distance = 1000.0;
    for (uint j=0;j<result.nItem();j++){
      float temp = sqrt(pow(mesh[0].vertex()[j][0]-sites[i].first[0],2) + pow(mesh[0].vertex()[j][1]-sites[i].first[1],2) +pow(mesh[0].vertex()[j][2]-sites[i].first[2],2));
      if (temp < distance) {
        distance = temp;
        node = j;
      }
    }

    aux.item(backnode)= 0.0;
    aux.item(node)=100.0;
    backnode=node;
    tex1a = MeshDistance_adapt_tex(mesh[0], aux, false,10.0);


    vois1.clear();
    vois2.clear();
    for (uint k=0;k<result.nItem();k++){
      if (tex1a.item(k) < simlocationnoise) vois1.push_back(k);
      if (tex1a.item(k) < 7.0) vois2.insert(k);
    }
    
    
    // bruit de la tvalue max
    random = (UniformRandom() * intensitynoise) - intensitynoise/2.0;

    cout << "tvalue2 : " << (sites[i].second) << endl;
    cout << "bruit de la tvalue : " << random << endl;
    // bruit de l'étendue
    set<uint>::iterator it,it2;
//     vois2.insert(node);
//     for (it=voisins[node].begin();it!=voisins[node].end();it++){
//       vois2.insert(*it);
//       for (it2=voisins[*it].begin();it2!=voisins[*it].end();it2++){
//         vois2.insert(*it2);
//       }
//     }
//     for (uint ii=0;ii<vois1.size();ii++)
//       vois2.insert(vois1[ii]);
    
    cout << (sites[i].second + random)*sqrt(smooth1)/2.0*log(2.) << endl;
    for (it=vois2.begin();it!=vois2.end();it++){
      
      tex.item(*it) += (sites[i].second + random)*sqrt(smooth1)/2.0*log(2.);
    }
    cout << "===" << endl;

  }
  
    // lissage des diracs
  tex=smoother->doSmoothing(tex, 0, ((int)(smooth1/0.05))*0.05);

  
  //fusion des deux textures
  for (uint j=0;j<result.nItem();j++)
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


int main(int argc, const char **argv){
   try
   {
      string outpath, meshpath, meshpath2, datapath, datapath1;
      float noise;
      float location;
      float intensity;
      float smooth;
      float f1x,f1y,f1z,f1t,f2x,f2y,f2z,f2t,f3x,f3y,f3z,f3t,f4x,f4y,f4z,f4t,f5x,f5y,f5z,f5t;
//       int operation, size=7, type = 0, index=-1;
//       float vsizeX=3.0,vsizeY=3.0,vsizeZ=3.0, geod_decay=5.0, norm_decay=2.0;
      AimsApplication app( argc, argv, "AimsFunctionProjection : first computes anatomically-informed kernels from one anatomy and uses them to project some functional data onto a cortical mesh" );
//       app.addOption( operation, "-op", "0 : computes convolution kernels from one anatomy ; 1 : projects functional volumes onto the surface (using kernels)");
      app.addOption( meshpath , "-m", "Grey/white matter mesh (.mesh)" );
      app.addOption( noise , "-n", "bruit de fond", 0.1 );
      app.addOption( location , "-l", "bruit de position", 0.1 );
      app.addOption( intensity , "-i", "bruit d'intensité", 0.1 );
      app.addOption( smooth , "-s", "lissage", 10.0);
      app.addOption( f1x , "-f1x", "lissage", 10.0);
      app.addOption( f1y , "-f1y", "lissage", 10.0);
      app.addOption( f1z , "-f1z", "lissage", 10.0);
      app.addOption( f1t , "-f1t", "lissage", 10.0);
      app.addOption( f2x , "-f2x", "lissage", 10.0);
      app.addOption( f2y , "-f2y", "lissage", 10.0);
      app.addOption( f2z , "-f2z", "lissage", 10.0);
      app.addOption( f2t , "-f2t", "lissage", 10.0);
      app.addOption( f3x , "-f3x", "lissage", 10.0);
      app.addOption( f3y , "-f3y", "lissage", 10.0);
      app.addOption( f3z , "-f3z", "lissage", 10.0);
      app.addOption( f3t , "-f3t", "lissage", 10.0);
      app.addOption( f4x , "-f4x", "lissage", 10.0);
      app.addOption( f4y , "-f4y", "lissage", 10.0);
      app.addOption( f4z , "-f4z", "lissage", 10.0);
      app.addOption( f4t , "-f4t", "lissage", 10.0);
      app.addOption( f5x , "-f5x", "lissage", 10.0);
      app.addOption( f5y , "-f5y", "lissage", 10.0);
      app.addOption( f5z , "-f5z", "lissage", 10.0);
      app.addOption( f5t , "-f5t", "lissage", 10.0);



//       app.addOption( type , "-t", "For computing convolution kernels (-op=0), selects the cortical thickness evaluation method : 0 for 3mm constant (only for now) ", 1); //; 1 for cortical mask-based evaluation ; 2 if using thickness texture\n  For projecting volumes onto the surface (-op=1), sets the output type : 0 for single volume projection ; 1 for several volumes -> several textures ; 2 for several volumes -> one timetexture", 0 );
      app.addOption( outpath , "-o", "Output file : convolution kernels (-op=0) or projection texture (-op=1)" );
//       app.addOption( index, "-I", "[DEBUG] Index of a precise kernel to be computed", 1);
      app.initialize();
      AimsSurfaceTriangle msh;
      AimsSurfaceTriangle *msh2;
      msh2 = SurfaceGenerator::sphere( Point3df(0.0,0.0,0.0), 1.0, 100);
//       Writer<AimsSurfaceTriangle> writer2("/home/grg/data/synth/sphere.mesh");
//       writer2.write(*msh2);
      Reader<AimsSurfaceTriangle> r(meshpath);
      r.read(msh);

      Writer<TimeTexture<float> > writer(outpath);
      TimeTexture<float> result;
      vector<pair<Point3df,float> > sites; //126.764, 23.888, 46.989
//           sites.push_back( pair<Point3df,float>(Point3df(126.764, 23.888, 46.989), 5.0));
      sites.push_back(pair<Point3df,float>(Point3df(f1x, f1y, f1z), f1t));
      sites.push_back(pair<Point3df,float>(Point3df(f2x, f2y, f2z), f2t));
      sites.push_back(pair<Point3df,float>(Point3df(f3x, f3y, f3z), f3t));
      sites.push_back(pair<Point3df,float>(Point3df(f4x, f4y, f4z), f4t));
      sites.push_back(pair<Point3df,float>(Point3df(f5x, f5y, f5z), f5t));
//       sites.push_back(pair<Point3df,float>(Point3df(113.669, 69.4213, 29.1237), 5.0));
// 
//       sites.push_back(pair<Point3df,float>(Point3df(141.079, 81.536, 58.6003), 5.0));
// 
//       sites.push_back(pair<Point3df,float>(Point3df(129.447, 131.799, 40.0678), 5.0));
// 
//       sites.push_back(pair<Point3df,float>(Point3df(123.486, 49.0954, 62.2587), 5.0));
      
      



      // float smooth, float intensitynoise, float backgroundnoise, float simlocationnoise, float simareanoise
      result[0] = createSynthData(msh, sites, smooth, intensity, noise, location);
      writer.write(result);

   }
   catch( user_interruption & )
   {
      return EXIT_FAILURE;
   }
   catch( exception & e )
   {
      cerr << e.what() << endl;
      return EXIT_FAILURE;
   }
   
   return EXIT_SUCCESS;
}



#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include "iograph.h"
#include "icm.h"
#include "anneal.h"
#include "cluster.h"

#include <aims/math/random.h>

using namespace aims;
using namespace std;
using namespace carto;

int main(int argc, const char **argv){
  string  graphFile, output, atlaspath, recuitpath = "/home/grg/recuit.txt", energypath = "/home/grg/energy.txt";
  Graph primal;
  float _ddweight=0.7, _intrapsweight = 2.0, _simweight=1.4, _lsweight=0.01, _ddx2 = 14.0, _ddx1 = 6.0, _ddh=0.0001;

  int run=1;
  AimsApplication     app( argc, argv, "Initialize");
  app.addOption( graphFile, "-p", "PS graph");
  app.alias( "--primal", "-p" );
  app.addOption( output, "-o", "PS graph", "");
  app.alias( "--output", "-o" );
  app.addOption( atlaspath, "-m" , "Chemin de base des données d'atlas");
  app.addOption(run,"--run","run","");
  app.addOption(_ddweight, "--ddw", "ddweight", 1.0);
  app.addOption(_intrapsweight, "--ipsw", "intrapsweight",  1.0);
  app.addOption(_simweight, "--simw", "simweight",  1.0);
  app.addOption(_lsweight, "--lsw", "lsweight",  1.0);
  app.addOption(_ddx1, "--ddx1", "ddx1",  1.0);
  app.addOption(_ddx2, "--ddx2", "ddx2",  1.0);
  app.addOption(_ddh, "--ddh", "ddh",  1.0);
  app.addOption(energypath, "--energypath", "energypath",  1.0);
  app.addOption(recuitpath, "--recuitpath", "recuitpath",  1.0);
  app.initialize();

  LireGraphes(graphFile,primal);
  set<string> sujets(RecupererSujets(primal));
  map<string, AimsSurfaceTriangle> meshes;
  map<string, TimeTexture<float> > lats, lons;
  set<string>::iterator it=sujets.begin();


  for (;it!=sujets.end();it++){
    string meshpath, latpath, lonpath;
    if (atlaspath.find("nmr_surface")!=string::npos){
    meshpath = atlaspath + *it + "/mesh/" + *it + "_Lwhite.mesh";
    latpath = atlaspath + *it + "/surface/" + *it + "_L_lat.tex";
    lonpath = atlaspath + *it + "/surface/" + *it + "_L_lon.tex";
    }
    else if (atlaspath.find("simulations")!=string::npos){
    meshpath = atlaspath + "sphere.mesh";
    latpath = atlaspath + "latitude.tex";
    lonpath = atlaspath + "longitude.tex";
    }
    else ASSERT(false);
    cout << "chargement" << endl;
    Reader<AimsSurfaceTriangle> rmesh(meshpath);
    Reader<TimeTexture<float> > rlat(latpath);
    Reader<TimeTexture<float> > rlon(lonpath);
    AimsSurfaceTriangle newmesh;
    rmesh.read(newmesh);

    meshes[*it] = AimsSurfaceTriangle(newmesh);

    TimeTexture<float> lat, lon;
    rlat.read(lat);
    rlon.read(lon);
    
    lats[*it]=TimeTexture<float>(lat);
    lons[*it]=TimeTexture<float>(lon);
  }
//   SWC swc(primal, mesh, lat, lon);
  Anneal swc(primal, meshes, lats,lons);
  
  swc.recuitpath = recuitpath;
  swc.energypath = energypath;
  
  cout << _ddweight << "-" << _intrapsweight << "-" << _simweight << "-" << _lsweight << "-" << _ddx2 << "-" << _ddx1 << "-" << _ddh << endl;
  swc.setModelParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh);
   
  swc.Run();
  swc.SummaryLabels();
  swc.StoreToGraph(primal);
  SauvegarderGraphes(primal, graphFile, output);

  return(0);

}


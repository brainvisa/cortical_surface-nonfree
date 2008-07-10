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
  string  graphFile, output, meshpath, latpath, lonpath;
  Graph primal;
  float _ddweight=0.7, _intrapsweight = 2.0, _simweight=1.0, _lsweight=0.01, _ddx2 = 5.0, _ddx1 = 1.0, _ddh=0.0001;

  int run=1;
  AimsApplication     app( argc, argv, "Initialize");
  app.addOption( graphFile, "-p", "PS graph");
  app.alias( "--primal", "-p" );
  app.addOption( output, "-o", "PS graph", "");
  app.alias( "--output", "-o" );
  app.addOption( meshpath, "-m" , "Maillage sur lequel on calcule les distances intersujets");
  app.addOption( latpath , "-lat", "Texture latitudes de ce maillage");
  app.addOption( lonpath , "-lon", "Texture longitudes de ce maillage");
  app.addOption(run,"--run","run","");
  app.addOption(_ddweight, "--ddw", "ddweight", 1.0);
  app.addOption(_intrapsweight, "--ipsw", "intrapsweight",  1.0);
  app.addOption(_simweight, "--simw", "simweight",  1.0);
  app.addOption(_lsweight, "--lsw", "lsweight",  1.0);
  app.addOption(_ddx1, "--ddx1", "ddx1",  1.0);
  app.addOption(_ddx2, "--ddx2", "ddx2",  1.0);
  app.addOption(_ddh, "--ddh", "ddh",  1.0);
  app.initialize();

  LireGraphes(graphFile,output,primal);

  AimsSurfaceTriangle mesh;
  Reader<AimsSurfaceTriangle> r(meshpath);
  r.read(mesh);
  TimeTexture<float> lat, lon;
  Reader<TimeTexture<float> > rlat(latpath);
  rlat.read(lat);
  Reader<TimeTexture<float> > rlongit(lonpath);
  rlongit.read(lon);
  SWC swc(primal, mesh, lat, lon);
  cout << _ddweight << "-" << _intrapsweight << "-" << _simweight << "-" << _lsweight << "-" << _ddx2 << "-" << _ddx1 << "-" << _ddh << endl;
  swc.setModelParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh);
  
  swc.Run();
  swc.SummaryLabels();
  SauvegarderGraphes(primal, graphFile, output);

  return(0);

}


#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <cortical_surface/structuralanalysis/iograph.h>
#include <cortical_surface/structuralanalysis/icm.h>
#include <cortical_surface/structuralanalysis/old_anneal.h>
#include <cortical_surface/structuralanalysis/anneal.h>
#include <cortical_surface/structuralanalysis/cluster.h>
// #include "validation.h"

#include <aims/math/random.h>

using namespace aims;
using namespace std;
using namespace carto;

int main(int argc, const char **argv){
  string  graphFile, output, atlaspath= "", recuitpath="", energypath="";
  int verbose=1;
  Graph primal;
  uint save=1;
  float _ddweight=0.8, _intrapsweight = 4.0, _simweight=1.0, _lsweight=1.0, _ddx1 = 3.125, _ddx2 = 4.50, _ddh=0.0001;
  _ddx1 = 8.0;
  _ddx2 = 4.0;


  int run=1;
  AimsApplication     app( argc, argv, "Initialize");
  app.addOption( graphFile, "-p", "PS graph");
  app.alias( "--primal", "-p" );
  app.addOption( output, "-o", "PS graph", "");
  app.alias( "--output", "-o" );
  app.addOption(run,"--run","run",0);
  app.addOption(save,"--save","save",0);
  app.addOption(verbose,"--verbose","verbose", 1.0);
  app.addOption(_ddweight, "--ddw", "ddweight", 1.0);
  app.addOption(_intrapsweight, "--ipsw", "intrapsweight",  1.0);
  app.addOption(_simweight, "--simw", "simweight",  1.0);
  app.addOption(_lsweight, "--lsw", "lsweight",  1.0);
  app.addOption(_ddx1, "--ddx1", "ddx1",  1.0);
  app.addOption(_ddx2, "--ddx2", "ddx2",  1.0);
  app.addOption(_ddh, "--ddh", "ddh",  1.0);
  app.addOption(energypath, "--energypath", "energypath",  0);
  app.addOption(recuitpath, "--recuitpath", "recuitpath",  0);
  app.initialize();

  LireGraphes(graphFile,primal);
//   set<string> sujets(RecupererSujets(primal));
  Anneal swc(primal);
  
  swc.recuitpath = recuitpath;
  swc.energypath = energypath;
  
  cout << _ddweight << "-" << _intrapsweight << "-" << _simweight << "-" << _lsweight << "-" << _ddx2 << "-" << _ddx1 << "-" << _ddh << endl;
  swc.setModelParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh);
  swc.run = run;
  swc.Initialization();

  swc.Run(verbose);
  swc.SummaryLabels();
  swc.StoreToGraph(primal);
  if (save) SauvegarderGraphes(primal, graphFile, output);
  
  return(0);

}


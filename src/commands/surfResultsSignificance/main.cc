#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

#include <cortical_surface/structuralanalysis/validation.h>
#include <cortical_surface/structuralanalysis/anneal.h>
#include <cortical_surface/structuralanalysis/iograph.h>

#include <aims/math/random.h>

using namespace aims;
using namespace std;
using namespace carto;

int main(int argc, const char **argv){
  string  graphFile, output, atlaspath= "", recuitpath="", energypath="";
  int verbose=1;
  Graph primal;
  float _ddweight=0.8, _intrapsweight = 4.0, _simweight=1.0, _lsweight=1.0, _ddx1 = 3.125, _ddx2 = 4.50, _ddh=0.0001;
  _ddx1 = 8.0;
  _ddx2 = 4.0;


  AimsApplication     app( argc, argv, "Initialize");
  app.addOption( graphFile, "-p", "PS graph");
  app.alias( "--primal", "-p" );
  app.addOption( output, "-o", "PS graph", "");
  app.alias( "--output", "-o" );
  app.addOption(verbose,"--verbose","verbose", 1.0);
  app.addOption(_ddweight, "--ddw", "ddweight", 1.0);
  app.addOption(_intrapsweight, "--ipsw", "intrapsweight",  1.0);
  app.addOption(_simweight, "--simw", "simweight",  1.0);
  app.addOption(_lsweight, "--lsw", "lsweight",  1.0);
  app.addOption(_ddx1, "--ddx1", "ddx1",  1.0);
  app.addOption(_ddx2, "--ddx2", "ddx2",  1.0);
  app.addOption(_ddh, "--ddh", "ddh",  1.0);
  app.initialize();

  LireGraphes(graphFile,primal);
  set<string> sujets(RecupererSujets(primal));
  Anneal swc(primal);
  swc.recuitpath = recuitpath;
  swc.energypath = energypath;

  swc.setModelParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh);
  cout << "Estimating the results significance : " << endl;

  swc.Initialization();
  cout << "Estimating the results significance : " << endl;
  StructuralAnalysis_Validation valid(&swc);
  valid.ValidAround();
  swc.StoreToGraph(primal);
  SauvegarderGraphes(primal, graphFile, output);
  return(0);

}


#include <cstdlib>
#include <aims/io/writer.h>
#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>
#include "gyri_parameterization.h"

using namespace aims;
using namespace std;
using namespace carto;

int main( int argc, const char** argv )
{
  try
    {
      string meshfile; 
      string intexfile;
      string outtexfile;
      string gyrifile;
      string diffmodfile;
      string spmtexfile;
      uint constraint_method;
      float criter;
      float dt; 
   //
   // Parser of options
   //
      AimsApplication app( argc, argv, "Paramétrisation de la surface corticale à partir d'une parcellisation en gyri");
      app.addOption( meshfile, "-i", "inputMesh");
      app.alias( "--inMesh", "-i" );
      app.addOption( intexfile, "-t", "gyral parcellation inputTexture");
      app.alias( "--inTex", "-t" );   
      app.addOption( gyrifile, "-g", "gyri to texture correspondance table");
      app.alias( "--gyri", "-g" );
      app.addOption( diffmodfile, "-m", "text file model of the relations between gyri and the system constraints");
      app.alias( "--model", "-m" );
      app.addOption( spmtexfile, "-s", "activation texture for surface-based mapping ( default = \"\" )", "");
      app.alias( "--spmtex", "-s" );
      app.addOption( constraint_method, "-c", "constraint method (0 = w/o constraints, 1 = iso-extracted constraints, 2 = mixed constraints, 3 = functional constraints iso-extracted mode, 4 = functional constraints mixed mode, 5 = spot constraints)", 0);
      app.alias( "--constraint", "-c" );
      app.addOption( criter, "-C", "Limit for the diffusion process ( default = 1e-4 )", 1e-4);
      app.alias( "--criter", "-C" );
      app.addOption( dt, "-d", "Iterative step ( default = 0.05 )", 0.05);
      app.alias( "--dt", "-d" );
      app.addOption( outtexfile, "-o", "outputTexture");
      app.alias( "--outTex", "-o" );
      app.initialize();      
      
      TimeTexture<float> outTex(GyriParamTexture(meshfile, intexfile, gyrifile, diffmodfile, spmtexfile, constraint_method, criter, dt));

      stringstream sstr;  string s;   sstr << constraint_method;   sstr >> s;
      outtexfile = outtexfile + "_"+s+".tex";
      Writer<TimeTexture<float> > w(outtexfile);
      w.write(outTex);

      return 0;  
    }
  catch( user_interruption & )  {  }
  catch( exception & e )
    {
      cerr << e.what() << endl;
      return EXIT_FAILURE;
    }
  return( EXIT_SUCCESS );
}


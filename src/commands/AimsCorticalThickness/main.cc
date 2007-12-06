#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <cortical_surface/mesh/median.h>

using namespace aims;
using namespace carto;
using namespace std;


int main(int argc, const char **argv){
   try
   {
      string intmeshpath, extmeshpath, neighbourspath, neighOutPath, outtexpath,neighbourspath2,medianmeshpath;
      uint direction=0;
      AimsApplication app( argc, argv, "AimsMeshMedianSurface" );
      app.addOption( intmeshpath , "-i", "White Hemisphere Mesh");
      app.addOption( extmeshpath , "-e", "Grey/CSF Hemisphere Mesh");
      app.addOption( neighbourspath , "-v", "Neighbours table file (from \"arrival\" mesh)", "" );
      app.addOption( direction, "-d", "Process direction (int -> ext || ext -> int)");
      app.addOption( neighOutPath, "-vout", "Output file in case the neighbours table would not be available and should be computed (6th order by default)", "" );
      app.addOption( outtexpath , "-o", "Output width texture", "");
      app.addOption( neighbourspath2 , "-v2", " [DEBUG ONLY] Neighbours table file (from \"other arrival\" mesh)", "" );
      app.addOption( medianmeshpath, "-m", "If provided, saves the created median mesh used for thickness computation", "");
      
      app.initialize();
      
      if (neighbourspath2 == ""){
         vector<set<uint> > voisins;
         AimsSurfaceTriangle intmesh, extmesh;
         Reader<AimsSurfaceTriangle> ri(intmeshpath);
         ri.read(intmesh);
         Reader<AimsSurfaceTriangle> re(extmeshpath);
         re.read(extmesh);
         if (neighbourspath != "")
            voisins = readVoisinsFromDisk(neighbourspath);
         else
         {
            voisins = compute_neighbours_order(extmesh, 4);
            if (neighOutPath != "")
               writeVoisinsToDisk(neighOutPath,voisins);
         }
         
         pair<AimsSurfaceTriangle,TimeTexture<float> > output(build_median_surface(intmesh,extmesh,voisins,direction));
         Writer<TimeTexture<float> > w(outtexpath);
         w.write(output.second);
         if (medianmeshpath != ""){
            Writer<AimsSurfaceTriangle> w1(medianmeshpath);
            w1.write(output.first);
         }
      }
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



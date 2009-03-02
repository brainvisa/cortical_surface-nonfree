#include <cstdlib>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <cortical_surface/mesh/median.h>

using namespace aims;
using namespace carto;
using namespace std;

pair<AimsSurfaceTriangle,AimsSurfaceTriangle> build_time_median_surface(AimsSurfaceTriangle &intmesh, AimsSurfaceTriangle &extmesh, vector<set<uint> > &voisins, vector<set<uint> > &voisins2, uint iter){
   AimsSurfaceTriangle outmeshInt,outmeshExt;
   pair<AimsSurfaceTriangle,TimeTexture<float> > output(build_median_surface(intmesh,extmesh,voisins,0));
   outmeshInt[0] = output.first[0];
   output= build_median_surface(extmesh,intmesh,voisins2,1);
   outmeshExt[0] = output.first[0];

   for (uint i=0;i<iter;i++){
      cout << i <<"th iteration" << endl;
                        
      intmesh[0] = outmeshInt[i];
      extmesh[0] = output.first[0];
            
      pair<AimsSurfaceTriangle,TimeTexture<float> > output(build_median_surface(intmesh,extmesh,voisins,0));
      outmeshInt[i+1] = output.first[0];
            
      output = build_median_surface(extmesh,intmesh,voisins2,1);
      outmeshExt[i+1] = output.first[0];            
   }
   pair<AimsSurfaceTriangle,AimsSurfaceTriangle> result(outmeshInt, outmeshExt);
   return result;
}

int main(int argc, const char **argv){
   try
   {
      string intmeshpath, extmeshpath, neighbourspath, neighOutPath, outmeshpath,neighbourspath2;
      uint direction=0;
      AimsApplication app( argc, argv, "AimsMeshMedianSurface" );
      app.addOption( intmeshpath , "-i", "White Hemisphere Mesh");
      app.addOption( extmeshpath , "-e", "Grey/CSF Hemisphere Mesh");
      app.addOption( neighbourspath , "-v", "Neighbours table file (from \"arrival\" mesh)", "" );
      app.addOption( direction, "-d", "Process direction (int -> ext || ext -> int)");
      app.addOption( neighOutPath, "-vout", "Output file in case the neighbours table would not be available and should be computed (8th order by default)", "" );
      app.addOption( outmeshpath , "-o", "Output mesh", "");
      app.addOption( neighbourspath2 , "-v2", " [DEBUG ONLY] Neighbours table file (from \"other arrival\" mesh)", "" );
      
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
         Writer<AimsSurfaceTriangle> w(outmeshpath);
         w.write(output.first);
      }
      else {
         vector<set<uint> > voisinsInt,voisinsExt;
         voisinsInt = readVoisinsFromDisk(neighbourspath);
         voisinsExt = readVoisinsFromDisk(neighbourspath2);
         AimsSurfaceTriangle intmesh, extmesh;
         Reader<AimsSurfaceTriangle> ri(intmeshpath);
         ri.read(intmesh);
         Reader<AimsSurfaceTriangle> re(extmeshpath);
         re.read(extmesh);

         pair<AimsSurfaceTriangle,AimsSurfaceTriangle> result(build_time_median_surface(intmesh,extmesh,voisinsInt,voisinsExt, 16));
         
         Writer<AimsSurfaceTriangle> w(outmeshpath + "_int.mesh");
         w.write(result.first);
         Writer<AimsSurfaceTriangle> w2(outmeshpath + "_ext.mesh");
         w2.write(result.second);
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



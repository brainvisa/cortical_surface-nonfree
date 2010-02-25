#include <cstdlib>
#include <aims/distancemap/meshdistance_d.h>
#include <aims/mesh/texture.h>
#include <aims/io/io_g.h>
#include "vector_operations.h"

using namespace std;
using namespace aims;
using namespace aims::meshdistance;

AimsSurfaceTriangle getGyrusMesh(AimsSurface<3, Void> &inMesh, const vector<uint> &gyrusVertices, vector<uint> &corres){
   AimsSurfaceTriangle gyrusMesh;
        // on extrait un gyrus, le maillage a moins de vertex que l'h�misph�re, du coup on cr�e un vecteur qui renseigne
        // sur les correspondances entre points homologues..s
   set<uint> gyrusSet;
   corres = *(new vector<uint>(inMesh.vertex().size()));
   for (uint i=0;i<gyrusVertices.size();i++){
      gyrusSet.insert(gyrusVertices[i]);
      gyrusMesh[0].vertex().push_back(inMesh.vertex()[gyrusVertices[i]]);
      corres[gyrusVertices[i]] = i;
   }
   for (uint i=0;i<inMesh.polygon().size();i++){
      if (gyrusSet.find(inMesh.polygon()[i][0])!=gyrusSet.end() &&
            gyrusSet.find(inMesh.polygon()[i][1])!=gyrusSet.end() &&
            gyrusSet.find(inMesh.polygon()[i][2])!=gyrusSet.end()){
               gyrusMesh[0].polygon().push_back(AimsVector<uint,3>(corres[inMesh.polygon()[i][0]],
               corres[inMesh.polygon()[i][1]],
               corres[inMesh.polygon()[i][2]]));
            }
         
   }
   return gyrusMesh;
}

map<unsigned, set<pair<unsigned,float> > > getGyrusWeight (map<unsigned, set<pair<unsigned,float> > > &poids, vector<uint> &gyrusVertices, vector<uint> &corres){
   map<unsigned, set<pair<unsigned,float> > > poidsGyrus;
   map<unsigned, set<pair<unsigned,float> > >::iterator mapIt;
   set<pair<unsigned,float> >::iterator setIt;
   set<uint> gyrusSet;
   for (uint i=0;i<gyrusVertices.size();i++)
      gyrusSet.insert(gyrusVertices[i]);
   for (uint i=0;i<gyrusVertices.size();i++){
      set<pair<unsigned,float> > aux(poids[gyrusVertices[i]]);
      set<pair<unsigned,float> > aux2;
      for (setIt = aux.begin();setIt != aux.end();setIt++){
         if (gyrusSet.find(setIt->first)!=gyrusSet.end()) aux2.insert(*(new pair<unsigned,float>(corres[setIt->first],setIt->second)));
      }
      poidsGyrus[corres[gyrusVertices[i]]] = aux2 ;
   }
   return poidsGyrus;
}

AimsSurface<3,Void> getFlatMesh(const AimsSurface<3,Void> &mesh, const vector<uint> &vertices, const vector<uint> &corres,
      TimeTexture<float> &paramTex){
   
   set<uint> gyrusSet;
   for (uint i=0;i<vertices.size();i++)
      gyrusSet.insert(vertices[i]);
   AimsSurface<3,Void> outMesh(mesh);
   TimeTexture<short> outTex(1,mesh.vertex().size());
   
   for (uint i=0;i<paramTex[0].nItem();i++)
      if (gyrusSet.find(i)!=gyrusSet.end()) {
         outMesh.vertex()[corres[i]] = *new AimsVector<float,3>(paramTex[1].item(i)+(i*0.01)/paramTex[0].nItem(), paramTex[2].item(i)+(i*0.01)/paramTex[0].nItem(), (i*0.01)/paramTex[0].nItem());
      }
   
   return outMesh;
}

Texture<float> AimsMeshLaplacian( const Texture<float> &inittex,
                      const map<unsigned, set< pair<unsigned,float> > > &lapl)
{
  unsigned              neigh,node, n =inittex.nItem();
  Texture<float>           tex(n);
  map<unsigned, set< pair<unsigned,float> > >::const_iterator il,el;
  set< pair<unsigned,float> >::iterator      ip,ep;
  float                                        L;
  float weight;
  ASSERT ( lapl.size() == n);
  
  for (il = lapl.begin(), el = lapl.end(); il != el; ++il)
    {
      node = il->first;
      L = 0;

      //Pondered sum on the neighbour of the node 
      for ( ip = (il->second).begin(), ep = (il->second).end(); ip != ep; ++ip    ) 
   {
     neigh = ip->first;
     weight = ip->second;
     L += weight * (- inittex.item(node) + inittex.item(neigh));
   }
      tex.item(node) = L;
    }

  return(tex);
}

Texture<float> diffusion(map<unsigned, set<pair<unsigned,float> > > &poidsGyrus,
      AimsSurface<3,Void> &mesh_base, vector<uint> &haut, vector<uint> &bas, const vector<pair<uint, float> > &constraints
      , short init, vector<uint> &gyrusvertices, vector<uint> &corr, const float criter, const float dt){
         
   Texture<float> tex1,lapl,tex2;
   Texture<float> tex1a, aux;
   float maxi=0.0;
   
   for (uint i=0;i<mesh_base.vertex().size();i++){
      aux.push_back(0.0);
      tex1.push_back(0.0);
      tex2.push_back(0.0);
      }
      
   for (uint i=0;i<haut.size();i++){
      aux.item(haut[i])= 100.0;
      }
   switch(init){
      case -1:
         tex1a = MeshDistance(mesh_base, aux, false);         
         for (uint i=0;i<mesh_base.vertex().size();i++)
            maxi = max(maxi, tex1a.item(i));
         printf("maxi : %f\n",maxi);
         for (uint i=0;i<mesh_base.vertex().size();i++){
            tex1.item(i) = 100-tex1a.item(i)/maxi*100;
            tex2.item(i) = tex1.item(i);
         }      
      break;
      default:
         for (uint i=0;i<mesh_base.vertex().size();i++){
            tex1.item(i) = (float)init;
            tex2.item(i) = (float)init;
         }
   }
      

      
   for (uint i=0;i<haut.size();i++)
      tex1.item(haut[i])= 100.0;            
   for (uint i=0;i<bas.size();i++)
      tex1.item(bas[i])= 0.0;        

   for (uint i=0;i<constraints.size();i++)
         tex1.item(constraints[i].first) = constraints[i].second;

   /*Texture<double> test;
   for (uint i=0;i<corr.size();i++)
      test.push_back(0);          
   writeTexture(test, "/opt/goperto/tempTex.tex");
   writeTexture(test, "/opt/goperto/laplTex.tex");*/

   
   printf("Diffusion en cours (seuil=0.001) : \n");
   fflush(NULL);

   float stop=9.0, moy=0.0, stop2=9.0, moy2=0.0;
   printf("%.9f", stop);
//   FILE *f1= fopen ("/home/grg/lapl.txt","a");
//   fprintf(f1,"============================================================================================\n");  
   

   for (long int n=1;stop2>criter;n++){
      stop=0.0;
      stop2=0.0;
      moy=0.0;
      moy2=0.0;
      lapl = AimsMeshLaplacian( tex1, poidsGyrus );
      for (uint i=0;i<mesh_base.vertex().size();i++){
         tex2.item(i) = dt*lapl.item(i) + tex1.item(i);
      }
      for (uint i=0;i<haut.size();i++)
         tex2.item(haut[i])= 100.0;
      for (uint i=0;i<bas.size();i++)
         tex2.item(bas[i])= 0.0;
      for (uint i=0;i<constraints.size();i++)
            tex2.item(constraints[i].first) = constraints[i].second;
      for (uint i=0;i<lapl.nItem();i++){
         moy +=lapl.item(i);
         moy2 += abs(lapl.item(i));
         stop = max(stop,abs(lapl.item(i)));
         stop2 = max(stop2, abs(tex2.item(i) - tex1.item(i)));
         tex1.item(i) = tex2.item(i);
         }
      moy = abs(moy/lapl.nItem());
      moy2 = abs(moy2/lapl.nItem());
      if (n%1000){
         cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << f2str((float)n).substr(0,7) << ":"<< f2str(stop).substr(0,15) << "  " << f2str(stop2).substr(0,15) << "   ";
         fflush(NULL);
      }
      /*if (n%7000==0) {
         addToTexture(tex1, "/opt/goperto/tempTex.tex", gyrusvertices, corr);
         //addToTexture(lapl, "/opt/goperto/laplTex.tex", gyrusvertices, corr);
      }*/
      if (n%100==0) {
//         fprintf(f1, "%.7f\t%.7f\t%.7f\t%.7f\n", stop, stop2, moy, moy2);
      }
      
   }
//   fclose(f1);
   printf("\n");

   return tex2;
}


#include <stdio.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surface.h>
#include "gyri_operations.h"
#include "vertices_operations.h"
#include "vector_operations.h"

using namespace std;
using namespace aims;


bool verifMeshTexMatch(const AimsSurface<3,Void> &inMesh, const Texture<short> &inTex){
   return (inMesh.vertex().size() == inTex.nItem());
}

bool verifTableTexMatch(const Texture<short> &inTex, const vector<uint> &corres){
   short mini=500, maxi=-500;
   for (uint i=0;i<inTex.nItem();i++){
      maxi = max(maxi, inTex.item(i));
      mini = min(mini, inTex.item(i));
   }
   //nb de gyri = maxi - mini +1
   ASSERT(mini != 500 && maxi != -500);
   
   return (mini==0 && (short)corres.size() == maxi - mini + 1);      
}

bool verifVoisins(const vector<set<uint> > &voisins, const Texture<short> &inTex, const vector<uint> &corres){
   uint test = 0;
   for (uint j=0;j<corres.size();j++)
      if (getNeighbours(j,voisins,inTex).size()<=2) {
         cerr << "Attention : le gyrus" << j << " n'a pas plus de 2 voisins.\n" << endl;
         test = 1;
      }   
   return (test == 0);
}

bool verifIntersectionsModel(const vector<set<uint> > &voisins, const Texture<short> &inTex, const vector<uint> &corres){
   uint test = 0;
   for (uint i=0;i<corres.size();i++)
      for (uint j=i+1;j<corres.size();j++)
         if (getIntersection(i,j, voisins, inTex).size() !=0)
            if (getNeighbouringLabels(i,j, voisins, inTex).size()!=2) {
               cerr << "Attention : les gyri " << i << " et " << j << " n'a pas plus de 2 voisins.\n" << endl;
               test=1;
            }
   return (test == 0);
}

bool verifGyriExistence(short g, const Point3d &gHaut, const Point3d &gBas, const vector<uint> &corres){
   uint test=0;
   if (find(corres, g)==-1){
      cerr << "Attention : le gyrus " << g << " référencé dans le modèle de diffusion n'existe pas dans la parcellisation.\n" << endl;
      test=1;
   }
   for (uint i=0;i<3;i++)
      if (find(corres, gHaut[i])==-1) {
         cerr << "Attention : le gyrus " << gHaut[i] << " référencé dans le modèle de diffusion n'existe pas dans la parcellisation.\n" << endl;
         test=1;
      }
      else if (find(corres, gBas[i])==-1){
         cerr << "Attention : le gyrus " << gBas[i] << " référencé dans le modèle de diffusion n'existe pas dans la parcellisation.\n" << endl;
         test = 1;
      }
   return (test==0);
}


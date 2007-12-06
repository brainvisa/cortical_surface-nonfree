#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include "vector_operations.h"
#include "model_operations.h"
#include "vertices_operations.h"

using namespace std;
using namespace aims;

pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr,
      const vector<set<uint> > &voisins, const Texture<short> &inTex){

      
   vector<vector<short> > constraints(getGyrusConstraints(gyrus, constrMod));
   pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > result;
   
   for (uint i=0;i<constraints.size();i++){
      vector<uint> inters(getRealIntersection(gyrus,constraints[i][0],constraints[i][1], voisins, inTex));
      if (inters.size()!=0){
         
         if (find(hautBas.first, inters[0])!=-1 || find(hautBas.second, inters[0])!=-1)
            result.second.push_back(*new pair<vector<uint>,short>(getCorresVector(inters,corr), constraints[i][2]));
         else if (find(gaucheDroite.first, inters[0])!=-1 || find(gaucheDroite.second, inters[0])!=-1)
            result.first.push_back(*new pair<vector<uint>,short>(getCorresVector(inters,corr), constraints[i][2]));
         else ASSERT(false);
         
      }
      else printf("Attention : the constraint (%d %d %d) has found no corresponding vertex.\n", gyrus, constraints[i][0], constraints[i][1]);
   }
   return result;
}



pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &vertTex, const Texture<double> &horizTex){

      
   vector<vector<short> > constraints(getGyrusConstraints(gyrus, constrMod));   
   pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > result;
   for (uint i=0;i<constraints.size();i++){      
      vector<uint> inters(getRealIntersection(gyrus,constraints[i][0],constraints[i][1], voisins, inTex));
      if (inters.size()!=0){        
         if (find(hautBas.first, inters[0])!=-1 || find(hautBas.second, inters[0])!=-1) {
            vector<uint> isoline(isolineExtraction(horizTex.item(inters[0]),gyrusMesh,horizTex));
            result.second.push_back(*new pair<vector<uint>,short>(isoline, constraints[i][2]));
         }
         else if (find(gaucheDroite.first, inters[0])!=-1 || find(gaucheDroite.second, inters[0])!=-1){
            vector<uint> isoline(isolineExtraction(vertTex.item(corr[inters[0]]),gyrusMesh,vertTex));
            result.first.push_back(*new pair<vector<uint>,short>(isoline, constraints[i][2]));
         }
         else ASSERT(false);
      }
      else printf("Attention : the constraint (%d %d %d) has found no corresponding vertex.\n", gyrus, constraints[i][0], constraints[i][1]);
   }
   return result;
}

void init_constvector(vector<uint> &haut, const vector<uint> &prohaut, const Texture<double> &horizTex){
   uint maxi,mini;
   maxi=0; mini=0;
   for (uint i=0;i<prohaut.size();i++){
      if (horizTex.item(prohaut[i])>horizTex.item(prohaut[maxi])) maxi = i;
      if (horizTex.item(prohaut[i])<horizTex.item(prohaut[mini])) mini = i;
      haut[i]=102;
   }
   haut[maxi] = 100;
   haut[mini] = 0;
}

double getRatioValue(uint valX, uint val1, uint val2, double valA, double valB){   
   double result =(valB + (double)(valX-val2)*(valA-valB)/(double)(val1-val2));
   return result;
}

void update_result(vector<pair<vector<uint>,short> > &result, const vector<uint> &haut,
      const vector<uint> &bas, const vector<uint> &prohaut, const vector<uint> &probas, const vector<uint> &corr,
      AimsSurface<3,Void> &flatMesh, const Texture<double> &horizTex){
   uint i=0;

   while (i<haut.size()){
      while ((haut[i]>=100 || haut[i]==0) && i<haut.size()) i++;
      if (i<haut.size()){
         uint up=find(bas,100), down=find(bas,0);
         for (uint j=0;j<bas.size() && !(bas[up]==haut[i] && bas[down]==haut[i]);j++){
            if (bas[j]==haut[i]) up = down = j;
            else if (bas[j]<bas[up] && bas[j]>haut[i]) up = j;
            else if (bas[j]>bas[down] && bas[j]<haut[i]) down = j;
         }

         double value = getRatioValue(haut[i], bas[up], bas[down], horizTex[corr[probas[up]]],
                   horizTex[corr[probas[down]]]);

         vector<uint> bas2(getCorresVector(probas,corr));

         uint target = getNearestPoint(bas2, horizTex, value);

         result.push_back(*new pair<vector<uint>,short>(lineExtraction(corr[prohaut[i]], target, flatMesh), haut[i]));
         i++;
      }
   }
}


pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, const Texture<double> &vertTex, const Texture<double> &horizTex){

   pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > result;
   
   vector<vector<short> > constraints(getGyrusConstraints(gyrus, constrMod));
   vector<uint> haut(hautBas.first.size());
   vector<uint> bas(hautBas.second.size());
   vector<uint> gauche(gaucheDroite.first.size());
   vector<uint> droite(gaucheDroite.second.size());   
   init_constvector(haut,getCorresVector(hautBas.first,corr),horizTex);
   init_constvector(bas,getCorresVector(hautBas.second,corr),horizTex);
   init_constvector(gauche,getCorresVector(gaucheDroite.first,corr),vertTex);
   init_constvector(droite,getCorresVector(gaucheDroite.second,corr),vertTex);
   
   for (uint i=0;i<constraints.size();i++){
      vector<uint> inters(getRealIntersection(gyrus,constraints[i][0],constraints[i][1], voisins, inTex));
      if (inters.size()!=0){
         int index;
         index = find(hautBas.first, inters[0]);
         ASSERT(index != (int)hautBas.first.size()-1 && index != 0);
         if (index!=-1) {haut[index] = constraints[i][2]; /*printf("haut\n");*/}
         else {
            index = find(hautBas.second, inters[0]);
            ASSERT(index != (int)hautBas.second.size()-1 && index != 0);
            if (index!=-1) {bas[index] = constraints[i][2]; /*printf("bas\n");*/}
            else {
               index = find(gaucheDroite.first, inters[0]);
               ASSERT(index != (int)gaucheDroite.first.size()-1 && index != 0);
               if (index!=-1) {gauche[index] = constraints[i][2]; /*printf("gauche\n");*/}
               else {
                  index = find(gaucheDroite.second, inters[0]);
                  ASSERT(index != (int)gaucheDroite.second.size()-1 && index != 0);
                  if (index!=-1) {droite[index] = constraints[i][2]; /*printf("droite\n");*/}
                  else {
                     ASSERT(false);
                  }
               }
            }
         }
      }
      else printf("Attention : the constraint (%d %d %d) has found no corresponding vertex.\n", gyrus, constraints[i][0], constraints[i][1]);
   }

   update_result(result.second, haut, bas, hautBas.first, hautBas.second, corr, flatMesh, horizTex);
   
   update_result(result.second, bas, haut, hautBas.second, hautBas.first, corr, flatMesh, horizTex);
   
   update_result(result.first, gauche, droite, gaucheDroite.first, gaucheDroite.second, corr, flatMesh, vertTex);
   
   update_result(result.first, droite, gauche, gaucheDroite.second, gaucheDroite.first, corr, flatMesh, vertTex);
   
   return result;
}

pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, const Texture<double> &vertTex, const Texture<double> &horizTex,
      const Texture<float> &spmTex){

      pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > result;
      result = getConstraints(gyrus, constrMod, hautBas, gaucheDroite, corr, voisins, inTex, flatMesh, vertTex, horizTex);
      //pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > result;

      vector<uint> haut(hautBas.first.size());
      vector<uint> bas(hautBas.second.size());
      vector<uint> gauche(gaucheDroite.first.size());
      vector<uint> droite(gaucheDroite.second.size());
      
      init_constvector(haut,getCorresVector(hautBas.first,corr),horizTex);
      init_constvector(bas,getCorresVector(hautBas.second,corr),horizTex);
      init_constvector(gauche,getCorresVector(gaucheDroite.first,corr),vertTex);
      init_constvector(droite,getCorresVector(gaucheDroite.second,corr),vertTex);
      
      int maxi=-1;
      for (uint i=0;i<spmTex.nItem() && maxi==-1;i++)
         if (inTex.item(i) == gyrus)
            maxi = i;
      for (uint i=0;i<spmTex.nItem();i++)
         if (inTex.item(i) == gyrus)
            if (spmTex.item(i)>spmTex.item(maxi)) maxi = i;
      vector<uint> inters;
      printf("MAXI : %d %.3f\n", maxi, vertTex.item(corr[maxi]));
      inters.push_back(find(corr,getNearestPoint(getCorresVector(gaucheDroite.first, corr), vertTex, vertTex.item(corr[maxi]))));
      int index;
      index = find(hautBas.first, inters[0]);
      if (index != (int)hautBas.first.size()-1 && index != 0){
         if (index!=-1) {haut[index] = 58; printf("haut\n");}
         else {
            index = find(hautBas.second, inters[0]);
            if (index != (int)hautBas.second.size()-1 && index != 0){
               if (index!=-1) {bas[index] = 58; printf("bas\n");}
               else {
                  index = find(gaucheDroite.first, inters[0]);
                  if (index != (int)gaucheDroite.first.size()-1 && index != 0){
                     if (index!=-1) {gauche[index] = 58; printf("gauche\n");}
                     else {
                        index = find(gaucheDroite.second, inters[0]);
                        if (index != (int)gaucheDroite.second.size()-1 && index != 0){
                           if (index!=-1) {droite[index] = 58; printf("droite\n");}
                           else {
                              ASSERT(false);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      update_result(result.second, haut, bas, hautBas.first, hautBas.second, corr, flatMesh, horizTex);
      
      update_result(result.second, bas, haut, hautBas.second, hautBas.first, corr, flatMesh, horizTex);
      
      update_result(result.first, gauche, droite, gaucheDroite.first, gaucheDroite.second, corr, flatMesh, vertTex);
      
      update_result(result.first, droite, gauche, gaucheDroite.second, gaucheDroite.first, corr, flatMesh, vertTex);
            
      /*for (uint i=0;i<constraints.first.size();i++)
         result.first.push_back(constraints.first[i]);
      for (uint i=0;i<constraints.second.size();i++)
         result.second.push_back(constraints.second[i]);
      */
      return result;

}

pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &vertTex,
      const Texture<double> &horizTex, const Texture<float> &spmTex){
      
      pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > result;
      result = getConstraints(gyrus, constrMod, hautBas, gaucheDroite, corr, voisins, inTex, flatMesh, vertTex, horizTex);
      
      //pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > result;
      int maxi=-1;
      for (uint i=0;i<spmTex.nItem() && maxi==-1;i++)
         if (inTex.item(i) == gyrus)
            maxi = i;
      
      for (uint i=0;i<spmTex.nItem();i++)
         if (inTex.item(i) == gyrus)
            if (spmTex.item(i)>spmTex.item(maxi)) maxi = i;
      vector<uint> isoline(isolineExtraction(vertTex.item(corr[maxi]),gyrusMesh,vertTex));
      result.first.push_back(*new pair<vector<uint>,short>(isoline, 57));
      
      return result;

}

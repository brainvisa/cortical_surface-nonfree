#include <cstdlib>
#include <math.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include "vector_operations.h"

using namespace std;
using namespace aims;


vector<vector<uint> > sortVertices(uint gyruslabel, const vector<set<uint> > &voisins, const Texture<short> &inTex) {

   // méthode qui pour un gyrus donné renvoit les indices des points du gyrus, des points du bord du gyrus, des points du bord étant voisins de deux gyri distincts

   uint test;
   set<uint>::const_iterator it;
   
   vector<vector<uint> > result;
   vector<uint> gyrusVertices;
   vector<uint> borderVertices;
   vector<uint> hotVertices;

   for (uint j=0;j<inTex.nItem();j++){
      if (!voisins[j].empty()) {
         test=0;
         short ref = inTex.item(j);
         short ref2 = 0;
         if ((uint) ref==gyruslabel){
            gyrusVertices.push_back(j);
            for (it = voisins[j].begin(); it != voisins[j].end() && test !=2 ; it++){
               if (test ==0 && inTex.item(*it) != ref){
                     test=1;
                     ref2 = inTex.item(*it);
                     borderVertices.push_back(j);
                  }
               else if (test==1 && inTex.item(*it)!=ref2 && inTex.item(*it)!=ref){
                  test = 2;
                  hotVertices.push_back(j);
               }
            }
         }
         else {}
      }
      else {}
   }
   result.push_back(gyrusVertices);
   result.push_back(borderVertices);
   result.push_back(hotVertices);
  
   return result;  
}



vector<uint> getIntersection(short gyruslabel, short voisin1, short voisin2, const vector<set<uint> > &voisins,
      const Texture<short> &inTex){
   
   // cette méthode renvoit les points d'une couleur donnée ayant pour voisins des points de deux couleurs données

   vector<uint> result;
   set<uint>::iterator it;
   short test;
   ASSERT(voisin1!=voisin2 && voisin1 != gyruslabel && gyruslabel != voisin2);
   for (uint i=0;i<voisins.size();i++){
      if (inTex.item(i) == gyruslabel){
         test=-1;
         for (it=voisins[i].begin();it!=voisins[i].end() && test!=-2;it++){
            if (test==-1 && (inTex.item(*it)==voisin1 || inTex.item(*it)==voisin2)){
               test = inTex.item(*it);
               //it=voisins[i].begin();
            }
            else if (test!=-1 && (inTex.item(*it)==voisin1 || inTex.item(*it)==voisin2) && test!=inTex.item(*it)){
               test=-2;
            }
         }
         if (test==-2){
            result.push_back(i);
         }
      }
   }
   return result;
}

vector<uint> getRealIntersection(short gyruslabel, short voisin1, short voisin2, const vector<set<uint> > &voisins,
      const Texture<short> &inTex){
   
   // cette méthode renvoit les points d'une couleur donnée ayant pour voisins des points de deux couleurs données

   set<uint> result;
   set<uint>::iterator it,it2,it3;
   
   ASSERT(voisin1!=voisin2 && voisin1 != gyruslabel && gyruslabel != voisin2);
   vector<uint> candidats(getIntersection(gyruslabel, voisin1, voisin2, voisins, inTex));
   for (uint i=0 ; i<candidats.size();i++){
      for (it = voisins[candidats[i]].begin();it!=voisins[candidats[i]].end(); it++){
         if (inTex.item(*it)==voisin1)
            for (it2=voisins[*it].begin();it2!=voisins[*it].end();it2++)
               if (inTex.item(*it2)==voisin2 && voisins[candidats[i]].find(*it2) != voisins[candidats[i]].end())
                  result.insert(candidats[i]);
      }
   }
   return setToVector(result);
}

vector<uint> getIntersection(short gyruslabel, short voisin1, const vector<set<uint> > &voisins, const Texture<short> &inTex){
    // cette méthode renvoit les points d'une couleur donnée ayant pour voisins des points d'une couleur donnée

   vector<uint> result;
   set<uint>::iterator it;
   uint test;
   ASSERT(voisin1 != gyruslabel);
   for (uint i=0;i<voisins.size();i++){
      if (inTex.item(i) == gyruslabel){
         test=0;
         for (it=voisins[i].begin();it!=voisins[i].end() && test != 1 ;it++){
            if (inTex.item(*it) == voisin1)  test=1;
         }
         if (test==1){
            result.push_back(i);
         }
      }
   }
   return result;
}

vector<uint> getThirdPoints(uint a, uint b, const vector<set<uint> > &voisins){

   // méthode qui pour deux indices de points donnés renvoit la liste de points voisins de l'un et de l'autre simultanément

   vector<uint> result;
   ASSERT(a!=b);
   set<uint>::iterator it;
   for (it=voisins[a].begin();it!=voisins[a].end();it++){
      if (voisins[b].find(*it)!=voisins[b].end()) result.push_back(*it);      
   }         
   return result;
}


vector<uint> getBorderNeighbours(uint v, const vector<uint> &borderVertices, const vector<set<uint> > &voisins){

   // à un point en entrée la méthode renvoit ses voisins contenus dans le bord

   set<uint>::iterator it;
   set<uint> result;
   ASSERT(v<voisins.size());
   for (it = voisins[v].begin();it != voisins[v].end();it++)
      if (find(borderVertices,*it)!=-1) result.insert(*it);
   return setToVector(result);
}

double getdistance(Point3df &a, Point3df &b){
   
   return ( sqrt((double) ((a[0]-b[0])*(a[0]-b[0]) +  (a[1]-b[1])*(a[1]-b[1])  +  (a[2]-b[2])*(a[2]-b[2]) ) )) ;

}

bool distanceCompare(Point3df &a, Point3df &b, Point3df &c){
   double d1 = getdistance(a,c);
   double d2 = getdistance(b,c);
   return (d1<d2);
}


vector<uint> isolineExtraction(double value, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &inTex){

   // méthode qui renvoit une liste d'indice de points correspondant 

   vector<set<uint> > voisins(SurfaceManip::surfaceNeighbours(gyrusMesh));
   //set<uint> result;
   vector<uint> result,aux;
   int index=-1;
   set<uint>::iterator it;

   for (uint i=0;i<inTex.nItem();i++)
      if (inTex.item(i)<=value) aux.push_back(i);

   for (uint i=0;i<aux.size();i++){
      index = -1;
      for (it = voisins[aux[i]].begin();it!=voisins[aux[i]].end() && index ==-1;it++)
         if (find(aux,*it)==-1) index = 0;
      if (index == 0) result.push_back(aux[i]);
   }
   
   
   return result;

}


vector<uint> lineExtraction2(uint start, uint end, const pair<vector<uint>, vector<uint> > &hautBas, AimsSurface<3,Void> &gyrusSurf){
   vector<uint> result;
   AimsSurfaceTriangle gyrusMesh;
   gyrusMesh[0] = gyrusSurf;
   vector<set<uint> > voisins(SurfaceManip::surfaceNeighbours(gyrusMesh));
   uint curr = hautBas.first[start];
   set<uint>::iterator it;
   int best;
   result.push_back(curr);
   while (curr != hautBas.second[end] && best != -1){      
      best = *voisins[curr].begin();
      //best = -1;
      //for (it = voisins[curr].begin() ; it != voisins[curr].end() && best == -1 ; it ++)
      //   if (find(result,*it)==-1) best = *it;
      if (best != -1){
         for (it = voisins[curr].begin() ; it != voisins[curr].end() ; it ++)
            if (find(result,*it)==-1 && getdistance(gyrusMesh[0].vertex()[*it], gyrusMesh[0].vertex()[hautBas.second[end]])< getdistance(gyrusMesh[0].vertex()[best], gyrusMesh[0].vertex()[hautBas.second[end]])) best = *it;
         result.push_back(best);
         curr = best;
      }  
   }
   return result;
}

Point3df cross(Point3df p1, Point3df p2){
   return *new Point3df(p1[1]*p2[2]-p1[2]*p2[1],p1[0]*p2[2]-p1[2]*p2[0],p1[0]*p2[1]-p1[1]*p2[0]);
}


vector<uint> lineExtraction(uint start, uint end, const pair<vector<uint>, vector<uint> > &hautBas, AimsSurface<3,Void> &gyrusSurf){
   vector<uint> result;
   AimsSurfaceTriangle gyrusMesh;
   gyrusMesh[0] = gyrusSurf;
   vector<set<uint> > voisins(SurfaceManip::surfaceNeighbours(gyrusMesh));
   set<uint>::iterator it;
   vector<uint> aux;
   int index;
   uint curr = hautBas.first[start];

   result.push_back(curr);
   Point3df vecteur(gyrusSurf.vertex()[hautBas.second[end]][0]-gyrusSurf.vertex()[hautBas.first[start]][0],
      gyrusSurf.vertex()[hautBas.second[end]][1]-gyrusSurf.vertex()[hautBas.first[start]][1],
      gyrusSurf.vertex()[hautBas.second[end]][2]-gyrusSurf.vertex()[hautBas.first[start]][2]);

   for (uint i=0;i<gyrusSurf.vertex().size();i++){
      if (cross(*new Point3df(gyrusSurf.vertex()[i][0]-gyrusSurf.vertex()[hautBas.first[start]][0],
      gyrusSurf.vertex()[i][1]-gyrusSurf.vertex()[hautBas.first[start]][1],
      gyrusSurf.vertex()[i][2]-gyrusSurf.vertex()[hautBas.first[start]][2]),vecteur)[2]>0.0) aux.push_back(i);
   }

   for (uint i=0;i<aux.size();i++){
      index = -1;
      for (it = voisins[aux[i]].begin();it!=voisins[aux[i]].end() && index ==-1;it++)
         if (find(aux,*it)==-1) index = 0;
      if (index == 0) result.push_back(aux[i]);
   }
   return result;
}

vector<uint> lineExtraction(uint start, uint end, AimsSurface<3,Void> &gyrusSurf){
   vector<uint> result;
   AimsSurfaceTriangle gyrusMesh;
   gyrusMesh[0] = gyrusSurf;
   vector<set<uint> > voisins(SurfaceManip::surfaceNeighbours(gyrusMesh));
   set<uint>::iterator it;
   vector<uint> aux;
   int index;
   uint curr = start;

   result.push_back(curr);
   Point3df vecteur(gyrusSurf.vertex()[end][0]-gyrusSurf.vertex()[start][0],
      gyrusSurf.vertex()[end][1]-gyrusSurf.vertex()[start][1],
      gyrusSurf.vertex()[end][2]-gyrusSurf.vertex()[start][2]);

   for (uint i=0;i<gyrusSurf.vertex().size();i++){
      if (cross(*new Point3df(gyrusSurf.vertex()[i][0]-gyrusSurf.vertex()[start][0],
      gyrusSurf.vertex()[i][1]-gyrusSurf.vertex()[start][1],
      gyrusSurf.vertex()[i][2]-gyrusSurf.vertex()[start][2]),vecteur)[2]>0.0) aux.push_back(i);
   }

   for (uint i=0;i<aux.size();i++){
      index = -1;
      for (it = voisins[aux[i]].begin();it!=voisins[aux[i]].end() && index ==-1;it++)
         if (find(aux,*it)==-1) index = 0;
      if (index == 0) result.push_back(aux[i]);
   }
   return result;
}


uint getNearestPoint(const vector<uint> &vec, const Texture<double> &inTex, double value){

   // à partir d'une liste d'indices de points, d'une texture et d'une valeur donnée, renvoit le point dont la valeur sur la texture est la plus proche de la valeur cible
 
   uint best=0;
   for (uint i=0;i<vec.size();i++)
      if (abs(inTex.item(vec[i])-value)<abs(inTex.item(vec[best])-value)) best = i;
   return vec[best];
}


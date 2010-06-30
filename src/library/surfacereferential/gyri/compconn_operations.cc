#include <aims/io/io_g.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surface.h>
#include <cortical_surface/surfacereferential/gyri/vector_operations.h>
#include <cortical_surface/surfacereferential/gyri/vertices_operations.h>


using namespace std;
using namespace aims;

vector<vector<uint> > getComposantesConnexes(const set<uint> &v, const vector<set<uint> > &voisins){
   vector<set<uint> > result;
   vector<vector<uint> > result2;
   set<uint>::iterator it,it3,ite;
   for (it=v.begin();it!=v.end();it++){
      set<uint> v2;
      v2.insert(*it);
      result.push_back(v2);
   }
   
   for (uint i=0;i<result.size();i++){

      printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b     [%d]    ",getRealSize(result));
      fflush(NULL);
      for (ite=result[i].begin();ite!=result[i].end();ite++){
         for (it=voisins[*ite].begin(); it!=voisins[*ite].end();it++){
            if (v.find(*it)!=v.end()){
               for (uint k=0;k<result.size();k++){
                  if (k!=i && result[k].find(*it)!=result[k].end()) {
                     for (it3=result[k].begin();it3!=result[k].end();it3++){
                        result[i].insert(*it3);
                     }
                     result[k].clear();
                  }
               }
            }
         }         
      }
   }
   

   for (uint i=0;i<result.size();i++)
      if (result[i].size()!=0) result2.push_back(setToVector(result[i]));
       
   printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b     [%d]         \n",result2.size());
   fflush(NULL);
   return result2;   
}

vector<vector<uint> > getComposantesConnexes2(const vector<uint> &v, const vector<set<uint> > &voisins){

   vector<vector<uint> > result;   
   int test;
   set<uint>::iterator it;
   for (uint i=0;i<v.size();i++){
      vector<uint> v2;
      v2.push_back(v[i]);
      result.push_back(v2);
   }

   
   for (uint i=0;i<result.size();i++){
      if (test != -1) i=0;
      printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b     [%d]    ",result.size());
      fflush(NULL);
      test = -1;
      for (uint j=0;j<result[i].size() && test == -1;j++){
         for (it=voisins[result[i][j]].begin(); it!=voisins[result[i][j]].end() && test == -1;it++){
            test = -1;
            //printf("(%d,%d) %d %d | ", i,j,*it, voisins[result[i][j]].size());
            if (find(v,*it) != -1){
               //printf("\n%d\n", find(v,*it));
               for (uint k=0;k<result.size() && test == -1;k++){
                  if (k!=i && find(result[k], *it)!=-1) test = k;
               }
               if (test != -1){
                  for (uint m=0;m<result[test].size();m++){
                     result[i].push_back(result[test][m]);
                  }
                  result.erase(result.begin() + test);
               }
            }
         }         
      }
   }
   printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b     [%d]         \n",result.size());
   fflush(NULL);
   return result;   
}

vector<vector<uint> > getComposantesConnexes(short gyruslabel, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   set<uint> v;
   for (uint i=0;i<inTex.nItem();i++)
      if (inTex.item(i)==gyruslabel) v.insert(i);
   return getComposantesConnexes(v, voisins);
}


vector<uint> fusionComposantesConnexes(uint newpoint, const vector<vector<uint> > &compConn, const vector<set<uint> > &voisins){
   set<uint> result;
   set<uint>::iterator it;
   ASSERT(newpoint<voisins.size());
   for (it=voisins[newpoint].begin();it!=voisins[newpoint].end();it++)
      for (uint i=0;i<compConn.size();i++)
         if (find(compConn[i], *it)!=-1) result.insert(i);

   return setToVector(result);
}

vector<vector<uint> > fusionComposantesConnexes(uint comp1, uint comp2, const vector<vector<uint> > &compConn){
   vector<vector<uint> > result;
   ASSERT(comp1 < compConn.size() && comp2 < compConn.size());
   vector<uint> aux(compConn[comp1]);
   push_vector(aux, compConn[comp2]);
   for (uint i=0;i<compConn.size();i++)
      if (i==comp1) result.push_back(aux);
      else if (i==comp2);
      else result.push_back(compConn[i]);
      
   return result;
}

void raccomodage(vector<vector<uint> > &compConn, const vector<uint> &candidates, const vector<set<uint> > &voisins){
   uint test=0;
   for (uint i=0;i<candidates.size() && test==0;i++){
      vector<uint> fcc(fusionComposantesConnexes(candidates[i], compConn, voisins));
      if (fcc.size()==2){
         compConn[fcc[0]].push_back(candidates[i]);
         fusionComposantesConnexes(fcc[0],fcc[1],compConn);
         test=1;
      }
   }
   vector<uint> aux;
   for (uint i=0;i<compConn.size();i++)
      for (uint j=i+1;j<compConn.size();j++){
         aux.clear();
         push_vector(aux, compConn[i]);
         push_vector(aux, compConn[j]);
         if (getComposantesConnexes2(aux, voisins).size()==1) {
            compConn = fusionComposantesConnexes(i,j,compConn);
            j=i=compConn.size();
         }
      }
}

pair<vector<uint>, vector<uint> > sortRightLeft(AimsSurface<3,Void> &inMesh, const pair<vector<uint>, vector<uint> > &hautBas,
      const pair<vector<uint>, vector<uint> > &gaucheDroite, const vector<set<uint> > &voisins){
   double m=0.0,n=0.0;
   Point3df p1(0.0,0.0,0.0);
   Point3df p2(0.0,0.0,0.0);
   pair<vector<uint>, vector<uint> > result;
   /*
   Point3df o(0.0,0.0,0.0);
   
   Point3df v1(inMesh.vertex()[hautBas.first[0]][0]-inMesh.vertex()[hautBas.second[0]][0],
      inMesh.vertex()[hautBas.first[0]][1]-inMesh.vertex()[hautBas.second[0]][1],
      inMesh.vertex()[hautBas.first[0]][2]-inMesh.vertex()[hautBas.second[0]][2]);
   Point3df v2;
   if (hautBas.first.size()!=1){
      v2 = *new Point3df(inMesh.vertex()[hautBas.first[hautBas.first.size()-1]][0]-inMesh.vertex()[hautBas.second[0]][0],
         inMesh.vertex()[hautBas.first[hautBas.first.size()-1]][1]-inMesh.vertex()[hautBas.second[0]][1],
         inMesh.vertex()[hautBas.first[hautBas.first.size()-1]][2]-inMesh.vertex()[hautBas.second[0]][2]);
   }
   else {
      int test=-1;
      for (uint i=0;i<gaucheDroite.first.size() && test==-1;i++)
         if (voisins[hautBas.first[0]].find(gaucheDroite.first[i]) != voisins[hautBas.first[0]].end())
            test = i;
      printf("%d\n", test);      
      v2 = *new Point3df(inMesh.vertex()[gaucheDroite.first[test]][0]-inMesh.vertex()[hautBas.second[0]][0],
         inMesh.vertex()[gaucheDroite.first[test]][1]-inMesh.vertex()[hautBas.second[0]][1],
         inMesh.vertex()[gaucheDroite.first[test]][2]-inMesh.vertex()[hautBas.second[0]][2]);
   }
   Point3df v1v2(cross(v1,v2));
   printf("%.3f %.3f %.3f\n",v1v2[0],v1v2[1],v1v2[2]);
   
   for (uint i=0;i<gaucheDroite.first.size();i++){
      Point3df w(inMesh.vertex()[gaucheDroite.first[i]][0]-inMesh.vertex()[hautBas.second[0]][0],
            inMesh.vertex()[gaucheDroite.first[i]][1]-inMesh.vertex()[hautBas.second[0]][1],
            inMesh.vertex()[gaucheDroite.first[i]][2]-inMesh.vertex()[hautBas.second[0]][2]);
      Point3df vw(cross(v1,w));
      //m += vw.dot(v1v2);
      p1[0] += vw[0];
      p1[1] += vw[1];
      p1[2] += vw[2];
   }

   //m = m/gaucheDroite.first.size();
   p1[0] = p1[0] / gaucheDroite.first.size();
   p1[1] = p1[1] / gaucheDroite.first.size();
   p1[2] = p1[2] / gaucheDroite.first.size();

   printf("%.3f\n", m);
   
   for (uint i=0;i<gaucheDroite.second.size();i++){
      Point3df w(inMesh.vertex()[gaucheDroite.second[i]][0]-inMesh.vertex()[hautBas.second[0]][0],
            inMesh.vertex()[gaucheDroite.second[i]][1]-inMesh.vertex()[hautBas.second[0]][1],
            inMesh.vertex()[gaucheDroite.second[i]][2]-inMesh.vertex()[hautBas.second[0]][2]);
      Point3df vw(cross(v1,w));
      //n += vw.dot(v1v2);
   }

   n = n/gaucheDroite.second.size();
   printf("%.3f\n", n);
   */
   for (uint i=0;i<gaucheDroite.first.size();i++){
      p1[0] += inMesh.vertex()[gaucheDroite.first[i]][0];
      p1[1] += inMesh.vertex()[gaucheDroite.first[i]][1];
      p1[2] += inMesh.vertex()[gaucheDroite.first[i]][2];
   }
   p1[0] /= gaucheDroite.first.size();
   p1[1] /= gaucheDroite.first.size();
   p1[2] /= gaucheDroite.first.size();
   
   for (uint i=0;i<gaucheDroite.second.size();i++){
      p2[0] += inMesh.vertex()[gaucheDroite.second[i]][0];
      p2[1] += inMesh.vertex()[gaucheDroite.second[i]][1];
      p2[2] += inMesh.vertex()[gaucheDroite.second[i]][2];
   }
   p2[0] /= gaucheDroite.second.size();
   p2[1] /= gaucheDroite.second.size();
   p2[2] /= gaucheDroite.second.size();
   
   m = p1[1];
   n = p2[1];
   
   if (m>n) {
      result.first = *new vector<uint>(gaucheDroite.first);
      result.second = *new vector<uint>(gaucheDroite.second);
   }
   else if (n>m){
      result.first = *new vector<uint>(gaucheDroite.second);
      result.second = *new vector<uint>(gaucheDroite.first);
   }
   else {
      printf("Impossible de distinguer gauche et droite !\n");
      result.first = *new vector<uint>(gaucheDroite.first);
      result.second = *new vector<uint>(gaucheDroite.second);
   }   
   return result;
}

pair<vector<uint>, vector<uint> > getOppositeSides(pair<vector<uint>, vector<uint> > &hautBas,
               const vector<vector<uint> > &vertices, const vector<set<uint> > &voisins, const Texture<short> &inTex){

   vector<uint> hb,gd,suspects,temp1,temp2,temp,aux1;
   vector<vector<uint> > compConn;
   set<uint>::iterator it;
   push_vector(hb, hautBas.first);
   push_vector(hb, hautBas.second);
   
   gd = minusVector(vertices[1], hb);
   
   uint test=0;
   //compConn = getComposantesConnexes(hb, voisins);
   vector<set<uint> > voisins2,voisins3;
   for (uint i=0;i<voisins.size();i++){
      voisins2.push_back(*new set<uint>);
      if (find(gd,i)!=-1){
         for (it=voisins[i].begin();it!=voisins[i].end();it++){            
            if (find(gd,*it)!=-1){               
               voisins2[voisins2.size()-1].insert(*it);
            }
         }
      }                  
   }
   
   for (uint i=0;i<voisins2.size();i++){
      voisins3.push_back(*new set<uint>);
      for (it=voisins2[i].begin();it!=voisins2[i].end();it++){

         vector<uint> pts(getThirdPoints(i,*it,voisins));
         for (uint j=0;j<pts.size() && test==0;j++)
            if (inTex.item(pts[j])!=inTex.item(vertices[0][0])) test=1;
         if (test!=0){
            voisins3[voisins3.size()-1].insert(*it);
            test=0;                        
         }
      }
   }
   //   set<uint> gd;
   //   for (uint i=0;i<gd2.size();i++)
   //      gd.insert(gd2[i]);
   
   compConn = getComposantesConnexes2(gd, voisins3);

                       


   
   if (compConn.size()==1){
      
      test = 0;
      if (hautBas.first.size() != 0 && hautBas.second.size() != 0){
         
         if (hautBas.first.size() <= hautBas.second.size()) {
            temp1 = *new vector<uint>(hautBas.first);
            temp2 = *new vector<uint>(hautBas.second);
         }
         else {
            temp1 = *new vector<uint>(hautBas.second);
            temp2 = *new vector<uint>(hautBas.first);
         }
         
         while (compConn.size() != 2){
            test++;
            for (uint i=0;i<test;i++){
               for (uint j=0;j<temp1.size();j++){
                  for (it = voisins[temp1[j]].begin(); it != voisins[temp1[j]].end() ; it++){
                     if (find(vertices[1],*it)!=-1)
                        
                        aux1.push_back(*it);
                  }
               }
               push_vector(temp1,aux1);             
            }
            temp = *new vector<uint>(vertices[1]);
            temp = minusVector(temp, temp1);
            temp = minusVector(temp, temp2);
            
            compConn = getComposantesConnexes2(temp, voisins2);
         }
         if (hautBas.first.size() <= hautBas.second.size()) {
            hautBas.first = *new vector<uint>(temp1);
         }
         else {
            hautBas.second = *new vector<uint>(temp1);
         }
      }
   }

   

   if (compConn.size()==3) raccomodage(compConn, minusVector(vertices[1],gd), voisins);
   if (compConn.size()==3) {
      vector<uint> aux;
      push_vector(aux, compConn[0]);
      //push_vector(aux, compConn[1]);
      push_vector(aux, compConn[2]);
      printf("Attention : getOppositeSides renvoit 3 composantes connexes. (Verifier la parcellisation)\n");
      return *new pair<vector<uint>, vector<uint> >(aux, compConn[1]);
   }
   if (compConn.size()==1) {
      for (uint i=49;i<70;i++)
         printf("%d ", compConn[0][i]);
      printf("\n");      
      printf("Attention : getOppositeSides renvoit 1 composante connexe. (Verifier la parcellisation)\n");
      return *new pair<vector<uint>, vector<uint> >(compConn[0], *new vector<uint>);
   }
   if (compConn.size()==2) return *new pair<vector<uint>, vector<uint> >(compConn[0],compConn[1]);
   else {
      printf("getOppositeSides renvoit %d composantes connexes.\n", compConn.size());
      return *new pair<vector<uint>, vector<uint> >;
   }
}

vector<vector<uint> > nettoyerTaches(Texture<short> &inTex, const vector<set<uint> > &voisins){
   vector<vector<uint> > result(getComposantesConnexes(0,voisins,inTex));
   vector<uint> v;
   set<uint>::iterator it;
   vector<uint> gyriVoisins;
   uint maxi=0;
   for (uint n=0;n<result.size();n++)
      if (result[maxi].size()<result[n].size()) maxi = n;

   // Elimination des taches (on remplit avec le premier voisin qu'on trouve)
   for (uint i=0;i<result.size();i++){
      gyriVoisins.clear();
      //if (result[i].size()<25){
      if (i!=maxi){
         for (uint j=0;j<result[i].size();j++){
            for (it = voisins[result[i][j]].begin();it!=voisins[result[i][j]].end();it++){
               if (inTex.item(*it) != 0) gyriVoisins.push_back(inTex.item(*it));
            }
         }
         for (uint j=0;j<result[i].size();j++){
            inTex.item(result[i][j]) = gyriVoisins[0];
         }
      }
   }
   
   printf("Nettoyage effectu�..\n");
   return result;
}


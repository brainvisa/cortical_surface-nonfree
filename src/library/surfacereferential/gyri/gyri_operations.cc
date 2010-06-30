#include <aims/mesh/texture.h>
#include <cortical_surface/surfacereferential/gyri/vector_operations.h>
#include <cortical_surface/surfacereferential/gyri/vertices_operations.h>

using namespace std;
using namespace aims;



vector<short> getNeighbours(short gyrus1, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   // cette m�thode renvoit les couleurs dont des points ont des voisins dans gyrus1
   vector<short> result;   
   set<uint>::iterator it;
   set<short> aux;
   set<short>::iterator it2;
   for (uint i=0;i<voisins.size();i++){
      if (inTex.item(i) == gyrus1)
         for (it=voisins[i].begin();it!=voisins[i].end();it++)
            if (inTex.item(*it)!=gyrus1) aux.insert(inTex.item(*it));
   }
   for (it2=aux.begin();it2!=aux.end();it2++){
      //printf(",%d\n",*it2);
      result.push_back(*it2);
      }
   return result;

}

vector<short> getNeighbouringLabels(short gyrus1, short gyrus2, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   // cette m�thode renvoit les couleurs dont il existe un point ayant pour voisins des points de deux couleurs donn�es
   // LES GYRI DOIVENT OBLIGATOIREMENT SE TOUCHER
   vector<short> result;
   set<short> aux;
   set<uint>::iterator it;
   set<short>::iterator it2;
   ASSERT(gyrus1!=gyrus2 && getIntersection(gyrus1,gyrus2,voisins,inTex).size()!=0);
   vector<short> gyriVois(getNeighbours(gyrus1, voisins, inTex));
   for (uint i=0;i< gyriVois.size();i++)
      if (gyriVois[i]!=gyrus2 && getRealIntersection(gyrus1, gyrus2, gyriVois[i], voisins, inTex).size()!=0)
         aux.insert(gyriVois[i]);
   if (aux.size()!=2) printf("[%d,%d]\n",gyrus1,gyrus2);
   for (it2=aux.begin();it2!=aux.end();it2++){
      if (aux.size()!=2) printf("/%d\n",*it2); 
      result.push_back(*it2);
   }
   return result;
}


vector<short> getNeighbouringLabels(uint n, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   set<short> result;
   set<uint>::iterator it;
   for (it = voisins[n].begin();it!=voisins[n].end();it++)
      if (inTex.item(*it) != inTex.item(n)) result.insert(inTex.item(*it));

   return setToVector(result);
}

vector<short> getNeighbouringLabels(const vector<uint> &v, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   set<short> result;
   vector<short> couleursVoisines;
   for (uint i=0;i<v.size();i++){
      couleursVoisines = getNeighbouringLabels(v[i], voisins, inTex);
      for (uint j=0;j<couleursVoisines.size();j++)
         result.insert(couleursVoisines[j]);
   }

   return setToVector(result);
}

vector<short> getCommonNeighbours(short gyrus1, short gyrus2, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   // cette m�thode renvoit les couleurs dont des points ont des voisins dans gyrus1, et gyrus2 (pas forc�ment un meme point
   // voisin � la fois de gyrus1 et gyrus2, non non..)
   vector<short> result;
   ASSERT(gyrus1!=gyrus2 && getIntersection(gyrus1,gyrus2,voisins,inTex).size()==0);
   set<uint>::iterator it;
   set<short> aux, aux2;
   set<short>::iterator it2;
   for (uint i=0;i<voisins.size();i++){
      if (inTex.item(i) == gyrus1)
         for (it=voisins[i].begin();it!=voisins[i].end();it++)
            if (inTex.item(*it)!=gyrus1) aux.insert(inTex.item(*it));
   }
   for (uint i=0;i<voisins.size();i++){
      if (inTex.item(i) == gyrus2)
         for (it=voisins[i].begin();it!=voisins[i].end();it++)
            if (inTex.item(*it)!=gyrus2 && aux.find(inTex.item(*it))!=aux.end()) aux2.insert(inTex.item(*it));
   }
   for (it2=aux2.begin();it2!=aux2.end();it2++){
      //printf(",%d\n",*it2);
      result.push_back(*it2);
      }
   return result;
}

vector<pair<short,short> > getInBetweenLabels(short gyrus1, short gyrus2, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   // cette m�thode renvoit les couleurs dont des points ont des voisins dans gyrus1, et/ou gyrus2
   vector<pair<short,short> > result;
   set<pair<short,short> > aux;
   ASSERT(gyrus1!=gyrus2 && getIntersection(gyrus1,gyrus2,voisins,inTex).size()==0);
   set<pair<short,short> >::iterator it;
   vector<short> gyriVoisins(getCommonNeighbours(gyrus1,gyrus2,voisins,inTex));
   for (uint i=0;i<gyriVoisins.size();i++)
      for (uint j=0;j<gyriVoisins.size();j++)
         if (j!=i
            && getIntersection(gyrus1, gyriVoisins[i], gyriVoisins[j], voisins, inTex).size()!=0
            && getIntersection(gyrus2, gyriVoisins[i], gyriVoisins[j], voisins, inTex).size()!=0){
            aux.insert(*new pair<short,short>(min(gyriVoisins[i],gyriVoisins[j]),max(gyriVoisins[i],gyriVoisins[j])));
         }

   for (it=aux.begin();it!=aux.end();it++)
      result.push_back(*it);
   return result;
}


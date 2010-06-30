#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <cortical_surface/surfacereferential/gyri/gyri_operations.h>
#include <cortical_surface/surfacereferential/gyri/vertices_operations.h>
#include <cortical_surface/surfacereferential/gyri/vector_operations.h>

using namespace std;
using namespace aims;


vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, short stop, const vector<set<uint> > &voisins,
   const Texture<short> &inTex){
   vector<short> result;
   short curr=start;
   short prev=forbidden;
   short aux;
   vector<short> gyriVoisins;
   while (curr!=stop){
      gyriVoisins = getNeighbouringLabels(gyruslabel, curr, voisins, inTex);
      ASSERT(gyriVoisins.size()==2 && (gyriVoisins[0]==prev || gyriVoisins[1]==prev));
      result.push_back(curr);
      aux = curr;
      if (gyriVoisins[0]==prev) curr=gyriVoisins[1];
      else if (gyriVoisins[1]==prev) curr=gyriVoisins[0];
      prev = aux;
   }
   return result;
}

vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, const pair<short,short> &stop, const vector<set<uint> > &voisins,
   const Texture<short> &inTex){
   vector<short> result;
   short curr=start;
   short prev=forbidden;
   short aux;
   vector<short> gyriVoisins;
   while (curr!=stop.first && curr!=stop.second){
      gyriVoisins = getNeighbouringLabels(gyruslabel, curr, voisins, inTex);
      ASSERT(gyriVoisins.size()==2 && (gyriVoisins[0]==prev || gyriVoisins[1]==prev));
      result.push_back(curr);
      aux = curr;
      if (gyriVoisins[0]==prev) curr=gyriVoisins[1];
      else if (gyriVoisins[1]==prev) curr=gyriVoisins[0];
      prev = aux;
   }
   return result;
}

vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, const vector<pair<short,short> > &stop, const vector<set<uint> > &voisins,
   const Texture<short> &inTex){
   vector<short> result;
   short curr=start;
   short prev=forbidden;
   short aux;
   vector<short> gyriVoisins;
   bool stoptest;
   ASSERT(stop.size()<=2);
   
   if (stop.size()==1){
      stoptest = (!((curr==stop[0].first && prev==stop[0].second) || (curr==stop[0].second && prev==stop[0].first)));
   }
   else {
      if (stop[1].second==-1)
         stoptest = (!((curr==stop[0].first && prev==stop[0].second) || (curr==stop[0].second && prev==stop[0].first) ||
                     (curr==stop[1].first)));
      else
         stoptest = (!((curr==stop[0].first && prev==stop[0].second) || (curr==stop[0].second && prev==stop[0].first) ||
                     (curr==stop[1].first && prev==stop[1].second) || (curr==stop[1].second && prev==stop[1].first)));
   }
   
   while (stoptest){
      gyriVoisins = getNeighbouringLabels(gyruslabel, curr, voisins, inTex);
      ASSERT(gyriVoisins.size()==2 && (gyriVoisins[0]==prev || gyriVoisins[1]==prev));
      result.push_back(curr);
      aux = curr;
      if (gyriVoisins[0]==prev) curr=gyriVoisins[1];
      else if (gyriVoisins[1]==prev) curr=gyriVoisins[0];
      prev = aux;
      
      if (stop.size()==1){
         stoptest = (!((curr==stop[0].first && prev==stop[0].second) || (curr==stop[0].second && prev==stop[0].first)));
      }
      else {
         if (stop[1].second==-1)
            stoptest = (!((curr==stop[0].first && prev==stop[0].second) || (curr==stop[0].second && prev==stop[0].first) ||
                        (curr==stop[1].first)));
         else
            stoptest = (!((curr==stop[0].first && prev==stop[0].second) || (curr==stop[0].second && prev==stop[0].first) ||
                        (curr==stop[1].first && prev==stop[1].second) || (curr==stop[1].second && prev==stop[1].first)));
      }
   }

   return result;
}


uint lookUpIntersectionCase(short gyruslabel, short left, short central, short right, const vector<set<uint> > &voisins,
   const Texture<short> &inTex){
   if (getIntersection(gyruslabel, central, voisins, inTex).size()!=0){
      bool GCL = (getRealIntersection(gyruslabel, central, left, voisins, inTex).size()!=0);
      bool GCR = (getRealIntersection(gyruslabel, central, right, voisins, inTex).size()!=0);
      if (GCL && GCR){      
         return 0;
      }
      else if (GCL && !GCR){
         if (getIntersection(central, right, voisins, inTex).size()!=0){
            return 1;
         }
         else if (getIntersection(gyruslabel, right, voisins, inTex).size()!=0){
            return 2;
         }
         else {
            if (getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1){
               return 3;
            }
            else {
               for (uint i=0;i<getInBetweenLabels(gyruslabel, right, voisins, inTex).size();i++)
                  printf("!%d, %d\n",getInBetweenLabels(gyruslabel, right, voisins, inTex)[i].first,getInBetweenLabels(gyruslabel, right, voisins, inTex)[i].second );
               return 255;
            }
         }
      }
      else if (!GCL && GCR){
         if (getIntersection(central, left, voisins, inTex).size()!=0){
            return 4;
         }
         else if (getIntersection(gyruslabel, left, voisins, inTex).size()!=0){
            return 5;
         }
         else {
            if (getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1){
               return 6;
            }
            else {
               return 254;
            }
         }
      }
      else if (!GCL && !GCR){
         bool GL = (getIntersection(gyruslabel, left, voisins, inTex).size()!=0);
         bool GR = (getIntersection(gyruslabel, right, voisins, inTex).size()!=0);
         if (GL && GR) return 7;
         else if (GL && !GR) {
            if (getIntersection(central, right, voisins, inTex).size()!=0){
               if (getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1)
                  return 8;
               else
                  return 253;
            }
            else {
               if (getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1){
                  return 9;
               }
               else {
                  return 252;
               }
            }
         }
         else if (!GL && GR){
            if (getIntersection(central, left, voisins, inTex).size()!=0){
               if (getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1)
                  return 10;
               else
                  return 251;
            }
            else {
               if (getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1){
                  return 11;
               }
               else {
                  return 250;
               }
            }            
         }
         else if (!GL && !GR){
            bool LC = (getIntersection(central, left, voisins, inTex).size()!=0);
            bool RC = (getIntersection(central, right, voisins, inTex).size()!=0);
            if (LC && RC) return 12;
            else if (LC && !RC){
               if (getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1){
                  if (getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1)
                     return 13;
                  else
                     return 249;
               }
               else {
                  for (uint i=0;i<getInBetweenLabels(gyruslabel, right, voisins, inTex).size();i++)
                     printf("!%d, %d\n",getInBetweenLabels(gyruslabel, right, voisins, inTex)[i].first,getInBetweenLabels(gyruslabel, right, voisins, inTex)[i].second );
                  return 248;
               }               
            }
            else if (!LC && RC){
               if (getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1){
                  if (getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1)
                     return 14;
                  else
                     return 253;
               }
               else {
                  return 247;
               }
            }
            else if (!LC && !RC){
               if (getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1
                     && getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1){
                  return 15;
               }
               else {
                  return 246;
               }
            }
         }
      }
   }
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   else {
      vector<pair<short,short> > gyriVoisins(getInBetweenLabels(gyruslabel, central, voisins, inTex));
      if (gyriVoisins.size()>1){
         uint test=0;
         for (uint i=0;i<gyriVoisins.size() && test==0;i++)
            if (gyriVoisins[i].first == min(left,right) && gyriVoisins[i].second == max(left,right)) test=1;
         if (test==1) {
            gyriVoisins.clear();
            return 20;
         }
         else {
            return 195;
         }
      }
      else if (gyriVoisins.size()==1 && gyriVoisins[0].first==min(left,right) && gyriVoisins[0].second==max(left,right)){
         return 21;
      }
      else if (gyriVoisins.size()==1) {
         bool GL = (getIntersection(gyruslabel, left, voisins, inTex).size()!=0);
         bool GR = (getIntersection(gyruslabel, right, voisins, inTex).size()!=0);
         if (GL && GR) return 22;
         else if (GL && !GR){
            if (getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1)
               return 23;
            else
               return 200;
         }
         else if (!GL && GR){
            if (getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1)
               return 24;
            else
               return 205;
         }
         else if (!GL && !GR){
            if (getInBetweenLabels(gyruslabel, left, voisins, inTex).size()==1 &&
                           getInBetweenLabels(gyruslabel, right, voisins, inTex).size()==1 )
               return 25;
            else
               return 210;
         }
         else return 215;
      }
      else if (gyriVoisins.size()==0){
         bool GLR = (getRealIntersection(gyruslabel, left, right, voisins, inTex).size()!=0);
         bool CLR = (getRealIntersection(central, left, right, voisins, inTex).size()!=0);
         if (GLR && !CLR){
            bool CL = (getIntersection(central, left, voisins, inTex).size()!=0);
            bool CR = (getIntersection(central, right, voisins, inTex).size()!=0);
            if (CR && !CL){
               return 26;
            }
            else if (CL && !CR){
               return 27;
            }
            else if (getInBetweenLabels(central, right, voisins, inTex).size()==1 &&
                     !(getInBetweenLabels(central, left, voisins, inTex).size()==1)){
                  return 28;   
            }
            else if (getInBetweenLabels(central, left, voisins, inTex).size()==1 &&
                     !(getInBetweenLabels(central, right, voisins, inTex).size()==1)){
                  return 29;   
            }
            else return 220;
         }
         else if (CLR && !GLR){
            bool GL = (getIntersection(gyruslabel, left, voisins, inTex).size()!=0);
            bool GR = (getIntersection(gyruslabel, right, voisins, inTex).size()!=0);
            if (GR && !GL){
               gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
               if (gyriVoisins.size()==1){
                  return 30;
               }
               else return 222;
            }
            else if (GL && !GR){
               gyriVoisins = getInBetweenLabels(gyruslabel, right, voisins, inTex);
               if (gyriVoisins.size()==1){
                  return 31;
               }
               else return 224;
            }
            else if (getInBetweenLabels(central, right, voisins, inTex).size()==1
                     && !(getInBetweenLabels(central, left, voisins, inTex).size()==1)){
                  return 32;
            }
            else if (getInBetweenLabels(central, left, voisins, inTex).size()==1
                     && !(getInBetweenLabels(central, right, voisins, inTex).size()==1)){
                  return 33;
            }
            else return 230;
         }
         else return 235;
      }
   }
   return 256;
}

vector<uint> substractIntersections(const vector<uint> &intersection, short gyruslabel, const vector<short> &parcours, const vector<set<uint> > &voisins, const Texture<short> &inTex){
   // cette m�thode vise � retirer les intersections aux extr�mit�s du segment s�lectionn� : le vecteur intersection est un
   // ensemble de points s�lectionn�s � partir du vecteur de labels parcours. L'id�e est d'en retirer les intersections des
   // extr�mit�s (pour �viter les �ventuelles fusions entre haut et bas dans des cas d�favorables).
   vector<uint> result;
   short gyrus;
   vector<short> gyriVoisins(getNeighbouringLabels(gyruslabel, parcours[0], voisins, inTex));
   ASSERT(gyriVoisins.size()==2);
   if (parcours.size()>1){
      //printf(">1\n");
      if (gyriVoisins[0]==parcours[1]) gyrus = gyriVoisins[1]; else gyrus = gyriVoisins[0];
      vector<uint> minus;
      minus = getIntersection(gyruslabel, parcours[0], gyrus, voisins, inTex);
      result = minusVector(intersection, minus);
      gyriVoisins = getNeighbouringLabels(gyruslabel, parcours[parcours.size()-1], voisins, inTex);
      ASSERT(gyriVoisins.size()==2);
      if (gyriVoisins[0]==parcours[parcours.size()-2]) gyrus = gyriVoisins[1]; else gyrus = gyriVoisins[0];
      minus = getIntersection(gyruslabel, parcours[parcours.size()-1], gyrus, voisins, inTex);
      result = minusVector(result, minus);
   }
   else if (parcours.size()==1){
      vector<uint> minus;
      //printf("=1 %d %d %d\n", parcours[0], gyriVoisins[0], gyriVoisins[1]);
      minus = getIntersection(gyruslabel, gyriVoisins[0], voisins, inTex);
      result = minusVector(intersection, minus);
      minus = getIntersection(gyruslabel, gyriVoisins[1], voisins, inTex);
      result = minusVector(result, minus);
   }
   if (result.size()==0)
      return intersection;
   else 
      return result;
   
}

vector<uint> getIntersection(uint code, short gyruslabel, short left, short central, short right, const vector<set<uint> > &voisins,
                                    const Texture<short> &inTex){
   vector<uint> result, aux;
   vector<short> parcours, sideGyri;
   vector<pair<short,short> > gyriVoisins,gyriVoisinsL,gyriVoisinsR,stop;
   set<uint>::iterator it,it2;
   short forbidden, gyrus;
   switch(code){
      case 0:
         parcours.push_back(central);
         result = getIntersection(gyruslabel, central, voisins, inTex);
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;
      break;
      case 1:
         parcours.push_back(central);
         result = getIntersection(gyruslabel, central, voisins, inTex);
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;
      break;
      case 2:
         parcours = parcoursPerim(gyruslabel, central, left, right, voisins,inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;
      break;
      case 3:
         gyriVoisins = getInBetweenLabels(gyruslabel,right,voisins,inTex);
         ASSERT(gyriVoisins.size()==1);         
         parcours = parcoursPerim(gyruslabel, central, left, gyriVoisins[0], voisins, inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;
      break;
      case 4:
         return getIntersection(gyruslabel, central, voisins, inTex);
      break;
      case 5:
         parcours = parcoursPerim(gyruslabel, central, right, left, voisins,inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);   
         return result;
      break;
      case 6:
         gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         parcours = parcoursPerim(gyruslabel, central, right, gyriVoisins[0], voisins, inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;         
      break;
      case 7:
         sideGyri = getNeighbouringLabels(gyruslabel, central, voisins, inTex);
         ASSERT(sideGyri.size()==2);
         parcours.push_back(central);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[0], central, *new pair<short,short>(left,right), voisins,inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[1], central, *new pair<short,short>(left,right), voisins,inTex));
         
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         //printf("\n");
         //print_vector(parcours);
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         //print_vector(getNeighbouringLabels(result,voisins, inTex));
         return result;
      break;
      case 8:
         gyriVoisins = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);         
         if (gyriVoisins[0].first == central) forbidden = gyriVoisins[0].second;
         else forbidden = gyriVoisins[0].first;         
         parcours = parcoursPerim(gyruslabel, central, forbidden, left, voisins, inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);   
         return result;
      break;
      case 9:
         sideGyri = getNeighbouringLabels(gyruslabel, central, voisins, inTex);
         ASSERT(sideGyri.size()==2);
         gyriVoisins = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         stop.push_back(gyriVoisins[0]);
         stop.push_back(*new pair<short,short>(left,-1));
         parcours.push_back(central);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[0], central, stop, voisins, inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[1], central, stop, voisins, inTex));
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result; 
      break;
      case 10:
         gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         if (gyriVoisins[0].first == central) forbidden = gyriVoisins[0].second;
         else forbidden = gyriVoisins[0].first;         
         parcours = parcoursPerim(gyruslabel, central, forbidden, right, voisins, inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;      
      break;
      case 11:
         sideGyri = getNeighbouringLabels(gyruslabel, central, voisins, inTex);
         ASSERT(sideGyri.size()==2);
         gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         stop.push_back(gyriVoisins[0]);
         stop.push_back(*new pair<short,short>(right,-1));
         parcours.push_back(central);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[0], central, stop, voisins, inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[1], central, stop, voisins, inTex));
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result; 
      break;
      case 12:
         return getIntersection(gyruslabel, central, voisins, inTex);
      break;   
      case 13:
         gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         if (gyriVoisins[0].first == central) forbidden = gyriVoisins[0].second;
         else forbidden = gyriVoisins[0].first;
         gyriVoisins = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);     
         parcours = parcoursPerim(gyruslabel, central, forbidden, gyriVoisins[0], voisins, inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;
      break;
      case 14:
         gyriVoisins = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         if (gyriVoisins[0].first == central) forbidden = gyriVoisins[0].second;
         else forbidden = gyriVoisins[0].first;
         gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);     
         parcours = parcoursPerim(gyruslabel, central, forbidden, gyriVoisins[0], voisins, inTex);
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result;      
      break;
      case 15:
         sideGyri = getNeighbouringLabels(gyruslabel, central, voisins, inTex);
         ASSERT(sideGyri.size()==2);
         gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         stop.push_back(gyriVoisins[0]);
         gyriVoisins = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         stop.push_back(gyriVoisins[0]);
         parcours.push_back(central);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[0], central, stop, voisins, inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, sideGyri[1], central, stop, voisins, inTex));
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result; 
      break;
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
      case 20:
         result = getIntersection(gyruslabel, left, right, voisins, inTex);
         for (uint i=0;i<result.size();i++)
            for (it=voisins[result[i]].begin();it!=voisins[result[i]].end();it++)            
               if (inTex.item(*it)==gyruslabel)
                  for (it2=voisins[*it].begin();it2!=voisins[*it].end();it2++)
                     if (inTex.item(*it2)!=gyruslabel) aux.push_back(*it);
         push_vector(result,aux);
            
         return result;
      break;
      case 21:
         result = getRealIntersection(gyruslabel, left, right, voisins, inTex);
         /*for (uint i=0;i<result.size();i++)
            for (it=voisins[result[i]].begin();it!=voisins[result[i]].end();it++)            
               if (inTex.item(*it)==gyruslabel)
                  for (it2=voisins[*it].begin();it2!=voisins[*it].end();it2++)
                     if (inTex.item(*it2)!=gyruslabel) aux.push_back(*it);
         push_vector(result,aux);*/
            
         return result;
      break;
      case 22:
         gyriVoisins = getInBetweenLabels(gyruslabel, central, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);         
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].first, gyriVoisins[0].second,
               *new pair<short,short>(left,right), voisins, inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].second, gyriVoisins[0].first,
               *new pair<short,short>(left,right), voisins, inTex));
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result; 
      break;
      case 23:
         gyriVoisins = getInBetweenLabels(gyruslabel, central, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         gyriVoisinsR = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisinsR.size()==1);
         stop.push_back(gyriVoisinsR[0]);
         stop.push_back(*new pair<short,short>(left,-1));         
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].first, gyriVoisins[0].second, stop, voisins, inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].second, gyriVoisins[0].first, stop, voisins, inTex));
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result; 
      break;
      case 24:
         gyriVoisins = getInBetweenLabels(gyruslabel, central, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         gyriVoisinsL = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisinsL.size()==1);
         stop.push_back(gyriVoisinsL[0]);
         stop.push_back(*new pair<short,short>(right,-1));
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].first, gyriVoisins[0].second, stop, voisins, inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].second, gyriVoisins[0].first, stop, voisins, inTex));
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result; 
      break;
      case 25:
         gyriVoisins = getInBetweenLabels(gyruslabel, central, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         gyriVoisinsL = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisinsL.size()==1);
         stop.push_back(gyriVoisinsL[0]);
         gyriVoisinsR = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisinsR.size()==1);
         stop.push_back(gyriVoisinsR[0]);
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].first, gyriVoisins[0].second, stop, voisins, inTex));
         parcours = reverseVector(parcours);
         push_vector(parcours, parcoursPerim(gyruslabel, gyriVoisins[0].second, gyriVoisins[0].first, stop, voisins, inTex));
         for (uint i=0;i<parcours.size();i++)
            push_vector(result,getIntersection(gyruslabel,parcours[i],voisins,inTex));
         result = substractIntersections(result, gyruslabel, parcours, voisins, inTex);
         return result; 
      break;
      case 26:
         return getIntersection(gyruslabel, left, right, voisins, inTex);
      break;
      case 27:
         return getIntersection(gyruslabel, left, right, voisins, inTex);
      break;
      case 28:
         return getIntersection(gyruslabel, left, right, voisins, inTex);
      break;
      case 29:
         return getIntersection(gyruslabel, left, right, voisins, inTex);
      break;
      case 30:
         gyriVoisins = getInBetweenLabels(gyruslabel, left, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         gyrus=-1;
         if (gyriVoisins[0].first == right) gyrus = gyriVoisins[0].second;
         else gyrus = gyriVoisins[0].first;
         ASSERT(gyrus!=-1);
         return getIntersection(gyruslabel, right, gyrus, voisins, inTex);
      break;
      case 31:
         gyriVoisins = getInBetweenLabels(gyruslabel, right, voisins, inTex);
         ASSERT(gyriVoisins.size()==1);
         gyrus=-1;
         if (gyriVoisins[0].first == left) gyrus = gyriVoisins[0].second;
         else gyrus = gyriVoisins[0].first;
         ASSERT(gyrus!=-1);
         return getIntersection(gyruslabel, left, gyrus, voisins, inTex);
      break;
      case 32:
         gyriVoisins = getInBetweenLabels(central, right, voisins, inTex);
         return getIntersection(gyruslabel, gyriVoisins[0].first, gyriVoisins[0].second, voisins, inTex);
      break;
      case 33:
         gyriVoisins = getInBetweenLabels(central, left, voisins, inTex);
         return getIntersection(gyruslabel, gyriVoisins[0].first, gyriVoisins[0].second, voisins, inTex);
      break;
      default:
         return result;
   }                                    

   return result;
   
}


#include <stdio.h>
#include <aims/vector/vector.h>
#include <cortical_surface/surfacereferential/gyri/vector_operations.h>

using namespace std;


vector<string> readGyriToTextureFile(const char *path){
   FILE *f = fopen(path,"r");   
   vector<string> tab;
   char s[50];
   string str;
   while (fgets(s,50,f)!=NULL){
      str = s;
      str.erase(0,str.find('\t')+1);
      str.erase(str.end()-1);
      tab.push_back(str);
   }
   return tab;
}


pair<vector<vector<uint> >, vector<vector<uint> > > getDiffusionModel(const char *path){
   FILE *f = fopen(path,"r");   
   vector<vector<uint> > tab;
   vector<vector<uint> > tab2;
   char s[50];
   string str;
   uint i;
   int step = -1;
   while (fgets(s,50,f)!=NULL){
      str = s;
      if (step == -1 && (int) str.find("[Gyri]\n")!=-1) {
         step = 0;
         ASSERT(fgets(s,50,f)!=NULL);
         str = s;
      }
      else if (step == 0 && (int) str.find("[Constraints]\n")!=-1) {
         step = 1;
         ASSERT(fgets(s,50,f)!=NULL);
         str = s;
      }      
      if (step == 0){
         vector<uint> v;
         for (uint k=0;k<7;k++){
            valueOf(str,i);
            v.push_back(i);
            str.erase(0,str.find('\t')+1);
         }
         tab.push_back(v);
      }
      else if (step == 1){
         vector<uint> v;
         for (uint k=0;k<4;k++){
            valueOf(str,i);
            v.push_back(i);
            str.erase(0,str.find('\t')+1);
         }
         tab2.push_back(v);
      }      
   }
   return *new pair<vector<vector<uint> >,vector<vector<uint> > >(tab,tab2);
}

vector<uint> gyrToTexCorres(const vector<string> &tab){
   vector<uint> result;
   for (uint i=0;i<tab.size();i++)
      if (tab[i] == "background") result.push_back(0);
      else if (tab[i] == "gyrus14_right") result.push_back(2);
      else if (tab[i] == "Frontal-Superior_right") result.push_back(3);
      else if (tab[i] == "-Lingual_right") result.push_back(4);
      else if (tab[i] == "gyrus12_right") result.push_back(6);
      else if (tab[i] == "gyrus19_right") result.push_back(7);
      else if (tab[i] == "gyrus8_right") result.push_back(8);
      else if (tab[i] == "Frontal-Inferior_right") result.push_back(9);
      else if (tab[i] == "Temporal-Superior_right") result.push_back(10);
      else if (tab[i] == "gyrus20_right") result.push_back(11);
      else if (tab[i] == "gyrus10_right") result.push_back(12);
      else if (tab[i] == "Pre-Central_right") result.push_back(13);
      else if (tab[i] == "Post-Central_right") result.push_back(14);
      else if (tab[i] == "gyrus15_right") result.push_back(15);
      else if (tab[i] == "Frontal-Middle_right") result.push_back(16);
      else if (tab[i] == "Orbital_right") result.push_back(17);
      else if (tab[i] == "Temporal-Inferior_right") result.push_back(18);
      else if (tab[i] == "Temporal-Middle_right") result.push_back(19);
      else if (tab[i] == "Cuneus_right") result.push_back(1);
      else if (tab[i] == "gyrus14_left") result.push_back(2);
      else if (tab[i] == "Frontal-Superior_left") result.push_back(3);
      else if (tab[i] == "Lingual_left") result.push_back(4);
      else if (tab[i] == "gyrus12_left") result.push_back(6);
      else if (tab[i] == "gyrus19_left") result.push_back(7);
      else if (tab[i] == "gyrus8_left") result.push_back(8);
      else if (tab[i] == "Frontal-Inferior_left") result.push_back(9);
      else if (tab[i] == "Temporal-Superior_left") result.push_back(10);
      else if (tab[i] == "gyrus20_left") result.push_back(11);
      else if (tab[i] == "gyrus10_left") result.push_back(12);
      else if (tab[i] == "Pre-Central_left") result.push_back(13);
      else if (tab[i] == "Post-Central_left") result.push_back(14);
      else if (tab[i] == "gyrus15_left") result.push_back(15);
      else if (tab[i] == "Frontal-Middle_left") result.push_back(16);
      else if (tab[i] == "Orbital_left") result.push_back(17);
      else if (tab[i] == "Temporal-Inferior_left") result.push_back(18);
      else if (tab[i] == "Temporal-Middle_left") result.push_back(19);
      else if (tab[i] == "Cuneus_left") result.push_back(1);
      else if (tab[i] == "Machin") result.push_back(5);
      else if (tab[i] == "Nawak") result.push_back(20);
      else result.push_back(255);
   return result;
}

vector<uint> getGyriToTextureCorres(const char *path){
   return gyrToTexCorres(readGyriToTextureFile(path));
}

pair<Point3d, Point3d> getGyrusModel(uint gyrus, const vector<vector<uint> > &diffMod){
   for (uint i=0;i<diffMod.size();i++){
      ASSERT(diffMod[i].size()==7);
      if (diffMod[i][0] == gyrus)
         return *new pair<Point3d, Point3d>(*new Point3d(diffMod[i][1],diffMod[i][2],diffMod[i][3]),
               *new Point3d(diffMod[i][4],diffMod[i][5],diffMod[i][6]));
   }
   return *new pair<Point3d, Point3d>;    
}

vector<vector<short> > getGyrusConstraints(uint gyrus, const vector<vector<uint> > &diffMod){
   vector<vector<short> > result;
   for (uint i=0;i<diffMod.size();i++){
      ASSERT(diffMod[i].size()==4);
      if (diffMod[i][0] == gyrus){
         vector<short> v;
         for (uint j=1;j<4;j++)            
            v.push_back((short)diffMod[i][j]);
         result.push_back(v);
      }
         
   }
   return result;

}


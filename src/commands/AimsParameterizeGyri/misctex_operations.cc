#include <aims/io/io_g.h>
#include <aims/mesh/texture.h>

using namespace std;
using namespace aims;



void makeGenericTexture(Texture<short> &inTex, const vector<uint> &corres){
   for (uint i=0;i<inTex.nItem();i++){
      inTex.item(i) = corres[inTex.item(i)];
   }
}

void createAuxTexture(uint size,const vector<pair<pair<vector<uint>, vector<uint> >, pair<vector<uint>, vector<uint> > > > &points,
    char* path){
   Writer<TimeTexture<short> > wTex(path);
   TimeTexture<short> tex(points.size(),0);
   for (uint i=0;i<points.size();i++)
      for (uint j=0;j<size;j++)
         tex[i].push_back(0);

   for (uint i=0;i<points.size();i++){
      for (uint j=0;j<points[i].first.first.size();j++)
          tex[i].item(points[i].first.first[j]) = 2;
      for (uint j=0;j<points[i].first.second.size();j++)
          tex[i].item(points[i].first.second[j]) = 4;
      for (uint j=0;j<points[i].second.first.size();j++)
          tex[i].item(points[i].second.first[j]) = 6;
      for (uint j=0;j<points[i].second.second.size();j++)
          tex[i].item(points[i].second.second[j]) = 8;
   }
   wTex.write(tex);
}


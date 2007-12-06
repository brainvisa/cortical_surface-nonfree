#include <aims/io/io_g.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace std;

TimeTexture<short> getTimeTex(const Texture<short> &inTex){
   TimeTexture<short> tex(0,0);
   for (uint i=0;i<inTex.nItem();i++)
      tex[0].push_back(inTex.item(i));
   return tex;
}
TimeTexture<float> getTimeTex(const Texture<double> &inTex){
   TimeTexture<float> tex(0,0);
   for (uint i=0;i<inTex.nItem();i++)
      tex[0].push_back((float)inTex.item(i));
   return tex;
}
TimeTexture<float> getTimeTex(const Texture<float> &inTex){
   TimeTexture<float> tex(0,0);
   for (uint i=0;i<inTex.nItem();i++)
      tex[0].push_back((float)inTex.item(i));
   return tex;
}
void writeTexture(const Texture<short> &inTex, char *path){
   Writer<TimeTexture<short> > writer(path);
   writer.write(getTimeTex(inTex));
}
void writeTexture(const TimeTexture<short> &inTex, char *path){
   Writer<TimeTexture<short> > writer(path);
   writer.write(inTex);
}
void writeTexture(const Texture<double> &inTex, char *path){
   Writer<TimeTexture<float> > writer(path);
   writer.write(getTimeTex(inTex));
}
void writeTexture(const Texture<float> &inTex, char *path){
   Writer<TimeTexture<float> > writer(path);
   writer.write(getTimeTex(inTex));
}

void addToTexture(const Texture<double> &inTex, char *path/*, vector<uint> &vertices, vector<uint> &corr*/){
   Reader<TimeTexture<float> > reader(path);
   TimeTexture<float> aux;
   reader.read(aux);
   TimeTexture<float> aux2(aux.size()+1,aux[0].nItem());
   for (uint i=0;i<aux.size();i++)
      for (uint j=0;j<aux[0].nItem();j++)
         aux2[i].item(j) = aux[i].item(j);
   for (uint j=0;j<aux[0].nItem();j++)
      aux2[aux2.size()-1].item(j) = inTex.item(j);         
   /*for (uint j=0;j<aux[0].nItem();j++)
      aux2[aux2.size()-1].item(j) = 0;
   for (uint j=0;j<vertices.size();j++)
      aux2[aux2.size()-1].item(vertices[j]) = inTex.item(corr[vertices[j]]);*/
   Writer<TimeTexture<float> > writer(path);
   writer.write(aux2);
}


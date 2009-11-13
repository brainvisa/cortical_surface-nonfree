 
#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/representation.h>

using namespace aims;
using namespace carto;
using namespace std;
 
  
pair<Point2df, Point2df> getBoundingBox(set<int> &nodes_list, TimeTexture<float> &lat, TimeTexture<float> &lon){
  Point2df bbmin, bbmax;
  bbmin[0] = 181.0;
  bbmin[1] = 361.0;
  bbmax[0] = -1.0;
  bbmax[1] = -1.0;
  set<int>::iterator it;
  pair<Point2df, Point2df> bb;
  for (it = nodes_list.begin() ; it != nodes_list.end() ; it ++){
    if (lat[0].item(*it) < bbmin[0])
      bbmin[0]=lat[0].item(*it);
    if (lon[0].item(*it) < bbmin[1])
      bbmin[1]=lon[0].item(*it);
    
    if (lat[0].item(*it) > bbmax[0])
      bbmax[0]=lat[0].item(*it);
    if (lon[0].item(*it) > bbmax[1])
      bbmax[1]=lon[0].item(*it);
  }
  
  if (bbmax[1] > 300.0 && bbmin[1] < 60.0) {
    for (uint i=0;i<nodes_list.size();i++){
      if (lon[0].item(*it) >300.0 && lon[0].item(*it) < bbmax[1])
        bbmax[1]=lon[0].item(*it);
      if (lon[0].item(*it) < 60.0 && lon[0].item(*it) > bbmin[1])
        bbmin[1]=lon[0].item(*it);
    }
  }
  
  bb.first = bbmin;
  bb.second = bbmax;
  return bb;
}
 

AimsSurfaceTriangle getFlatMap(vector<set<int> > &nodes_lists, TimeTexture<float> &lat, TimeTexture<float> &lon, TimeTexture<float> &tex){
  AimsSurfaceTriangle objects;
  for (uint i=0;i<nodes_lists.size();i++){
    if (nodes_lists[i].size()!=0){
      pair<Point2df,Point2df> bb(getBoundingBox(nodes_lists[i],lat,lon));
      assert(bb.first[0]<=bb.second[0] || !(cout << bb.first[0] << " /\\" << bb.second[0] << endl));
      assert(bb.first[1]<=bb.second[1]|| !(cout << bb.first[1] << " /\\" << bb.second[1] << endl));
      float area = (bb.second[0]-bb.first[0])*(bb.second[1]-bb.first[1]);
      if(area<1000.0){
        tex[0].push_back(area);
        tex[0].push_back(area);
        tex[0].push_back(area);
        tex[0].push_back(area);
//         cerr << tex[0].nItem() << " " << flush;
        objects[0].vertex().push_back(Point3df(bb.first[0],bb.first[1],0.001));
        objects[0].vertex().push_back(Point3df(bb.first[0],bb.second[1],0.002));
        objects[0].vertex().push_back(Point3df(bb.second[0],bb.second[1],0.003));
        objects[0].vertex().push_back(Point3df(bb.second[0],bb.first[1],0.0005));
//         cerr << objects[0].vertex().size() << endl;
        objects[0].polygon().push_back(AimsVector<uint,3>(objects[0].vertex().size()-4,objects[0].vertex().size()-3,objects[0].vertex().size()-2));
        objects[0].polygon().push_back(AimsVector<uint,3>(objects[0].vertex().size()-2,objects[0].vertex().size()-1,objects[0].vertex().size()-4));
      }
    }
  }

  return objects;
}


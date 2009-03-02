#include <cstdlib>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/distancemap/meshdistance_d.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace aims::meshdistance;

set<uint> nearest_vertices(Point3df pf, AimsSurfaceTriangle &mesh, float rayon){
   set<uint> v;
   for (uint i=0;i<mesh.vertex().size();i++){
      Point3df t(mesh[0].vertex()[i]);
      if (abs(t[0] - pf[0]) < rayon && abs(t[1] - pf[1]) < rayon && abs(t[2] - pf[2]) < rayon)
         v.insert(i);
   }
   return v;
}

pair<int,float> plus_proche_point_normal(Point3df p, Point3df n, AimsSurfaceTriangle &mesh, set<uint> &vertices){
   float dist=0.0,min = 5000.0,minscal=5000.0,maxscal=0.0; // un cerveau de 5 m, ça va comme limite ?
   int index = -1;
   set<uint>::iterator it;
   for (it=vertices.begin();it!=vertices.end();it++){
      float dot = n.dot(mesh[0].normal()[*it]);
      if (dot>maxscal) maxscal=dot;
      if (dot<minscal) minscal=dot;
   }
   for (it=vertices.begin();it!=vertices.end();it++){
      Point3df t(mesh[0].vertex()[*it]);
      float dot = n.dot(mesh[0].normal()[*it]);
      dist = Point3df(p-t).norm();
      if (dist < min && dot > minscal + (maxscal-minscal)/100.0*30.0 ) {
         index = *it;
         min = dist;
      }
   }
   return pair<int,float>(index, min);
}

pair<int,float> plus_proche_point_normal(Point3df p, Point3df n, AimsSurfaceTriangle &mesh){
   
   set<uint> vertices(nearest_vertices(p,mesh,10.0));
   return plus_proche_point_normal(p,n,mesh,vertices);
   
}

pair<Point3df,bool> isInsideTriangle(const Point3df &pt, const Point3df &vertex0, const Point3df &vertex1, const Point3df &vertex2)
{
   Point3df res(0.0,0.0,0.0);
   Point3df u = vertex1 - vertex0;
   Point3df v = vertex2 - vertex0;
   Point3df w = pt - vertex0;

   float uu = u.dot(u);
   float uv = u.dot(v);
   float vv = v.dot(v);
   float wu = w.dot(u);
   float wv = w.dot(v);
   float d = uv * uv - uu * vv;

   float invD = 1 / d;
   float s = (uv * wv - vv * wu) * invD;
   float t = (uv * wu - uu * wv) * invD;

   if (s < 0 || s > 1)
      return pair<Point3df,bool>(res,false);

   if (t < 0 || (s + t) > 1)
      return pair<Point3df,bool>(res,false);
   res = Point3df(vertex0);
   u *= s;
   v *= t;
   res += u;
   res += v;
   return pair<Point3df,bool>(res,true);
}

pair<Point3df, int> plus_proche_point_sur_triangle(Point3df p, Point3df n, AimsSurfaceTriangle &mesh, set<uint> &vertices, vector<set<uint> > &voisins){
   // soit p le noeud dont on cherche le noeud le plus proche, n la normale à ce noeud, mesh le maillage où l'on va chercher le noeud le plus proche,
   // vertices un ensemble de noeuds pour cibler la recherche sur mesh, voisins les voisins de premier ordre de mesh
   float dist=0.0,min=5000.0; // un cerveau de 5 m, ça va comme limite ?
   int index1 = -1,index2 = -1;
   
   set<uint>::iterator it,it2;

   index1 = plus_proche_point_normal(p, n, mesh, vertices).first;

   vector<Point3df> projetes;
   Point3df p1(mesh[0].vertex()[index1]);
   
   for (it=voisins[index1].begin();it!=voisins[index1].end();it++){
      for (it2=voisins[index1].begin();it2!=voisins[index1].end();it2++){
         if (*it > *it2){
            Point3df p2(mesh[0].vertex()[*it]);
            Point3df p3(mesh[0].vertex()[*it2]);
            pair<Point3df,bool> res(isInsideTriangle(p, p1,p2,p3));
            if (res.second)
               projetes.push_back(res.first);
         }
      }
   }
   if (projetes.size() != 0){
      min = 5000.0;
      for (uint i=0;i<projetes.size();i++){
         Point3df t(projetes[i]);
         dist = Point3df(p-t).norm();
         if (dist < min) {
            index2 = i;
            min = dist;
         }
      }
      return pair<Point3df, int>(projetes[index2],index1);
   }
   else
      return pair<Point3df, int>(p1,index1);

}

pair<Point3df, int> plus_proche_point_sur_triangle(Point3df p, Point3df n, AimsSurfaceTriangle &mesh, vector<set<uint> > &voisins){
   set<uint> vertices(nearest_vertices(p,mesh,10.0));
   return plus_proche_point_sur_triangle(p,n,mesh,vertices,voisins);
}


vector<set<uint> > readVoisinsFromDisk(string path){
   ifstream instream(path.data(), ios::out);
   string s,str; stringstream sstr; istringstream is;
   instream >> s;
   uint t;
   int index;
   vector<set<uint> > extvoisins;
   uint i=0;
   cout << "Reading neighbours table file..." << endl;
   while (!instream.eof()){
      cout << "\b\b\b\b\b\b\b" << i << flush;
      instream >> s;
      extvoisins.push_back(*(new set<uint>));
      do {
         is.clear(); is.str(s);  is >> t;
         extvoisins[extvoisins.size()-1].insert(t);
         index = s.find(";");
         s.erase(0,index+1);
      } while (index!=-1);
      instream >> s;
      i++;
   }
   cout << endl << " done" << endl;
   return extvoisins;
}

void writeVoisinsToDisk(const string &path, vector<set<uint> > &extvoisins4){
   set<uint>::iterator ite;
   cout << "Writing neighbours table.. " << endl;
   ofstream outstream(path.data(), ios::out);
   for (uint i=0;i<extvoisins4.size();i++){
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/"<< extvoisins4.size()-1 << flush;
      outstream << "[" << i << "]" << endl;
      outstream << *(extvoisins4[i].begin());
      ite=extvoisins4[i].begin();
      for (ite++;ite!=extvoisins4[i].end();ite++)
         outstream << ";" << *ite;
      outstream << endl;
   }
   cout << endl << " done" << endl;
}

vector<set<uint> > compute_neighbours_order(AimsSurfaceTriangle extmesh, uint order){
   set<uint>::iterator ite,ite1;
   vector<set<uint> > extvoisins4,extvoisins;
   
   // on part d'aucun fichier existant : on commence par extraire les voisins d'ordre 1 puis ca continue en doublant à chaque itération

   extvoisins4
       = *(new vector<set<uint> >(SurfaceManip::surfaceNeighbours(extmesh)));
   extvoisins
       = *(new vector<set<uint> >(SurfaceManip::surfaceNeighbours(extmesh)));

   cout << "Creating neighbours table.." << endl;
   cout << "1st order done" << endl;
   for (uint j=1;j<order;j++){
      cout << " " << j+1 << "th order" << endl;
      vector<set<uint> > extvoisins5(extvoisins4);
      for (uint i=0;i<extvoisins4.size();i++){
         cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << extvoisins4.size()-1 << flush;
         for (ite=extvoisins4[i].begin();ite!=extvoisins4[i].end();ite++){
            for (ite1=extvoisins[*ite].begin();ite1!=extvoisins[*ite].end();ite1++){
               extvoisins5[i].insert(*ite1);
            }
         }
      }
      extvoisins4 = *(new vector<set<uint> >(extvoisins5));
      cout << endl << " done" << endl;
   }
   
   return extvoisins4;
   
}


pair<AimsSurfaceTriangle, TimeTexture<float> > build_median_surface(AimsSurfaceTriangle &intmesh, AimsSurfaceTriangle &extmesh, vector<set<uint> > &extvoisins4, int op){
   
   AimsSurfaceTriangle matching_mesh(intmesh), median(intmesh);
   uint intsize = intmesh[0].vertex().size();   uint extsize = extmesh[0].vertex().size();
   vector<uint> intcorr(intsize);   vector<uint> extcorr(extsize);
   // intcorr[i] désigne le noeud de la surface externe correspondant au noeud i du maillage interne

   vector<set<uint> > intvoisins(SurfaceManip::surfaceNeighbours(intmesh));
   vector<set<uint> > extvoisins(SurfaceManip::surfaceNeighbours(extmesh));
   
   set<uint>::iterator ite,ite1;
   
   vector<uint> current,processed,neighbours;
   vector<uint> classe(intsize);   for (uint i=0;i<intsize;i++)  classe[i]=0;
   
   uint curr=0;
   
   current.push_back(curr);
   classe[curr]=1;
   Point3df p(intmesh[0].vertex()[curr]);
   Point3df pn(intmesh[0].normal()[curr]);
   pair<Point3df,int> currRes(plus_proche_point_sur_triangle(p,pn,extmesh,extvoisins));
   intcorr[curr] = currRes.second;
   matching_mesh[0].vertex()[curr] = currRes.first;
   
   Texture<float> aux,tex1a;   for (uint i=0;i<intsize;i++)  aux.push_back(0.0);   aux.item(curr) = 100.0;   tex1a = MeshDistance(intmesh[0], aux, false);

   vector<uint>::iterator it, best;
   set<uint>::iterator sit;

   TimeTexture<float> width;
   for (uint i=0;i<intsize;i++){
      width[0].push_back(-1.0);
   }
   width[0].item(curr) = Point3df(p-currRes.first).norm();

   cout << "Propagating the front.." << endl;
   
   while (processed.size() <= intsize){
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << processed.size() << "/" << intsize << "("<< current.size()<< ")" <<  flush;
      
      float min=10000.0;
      for (it=current.begin();it!=current.end();it++){
         if (tex1a.item(*it) < min) {
            best = it;
            min = tex1a.item(*it);
         }
      }
      processed.push_back(*best);
      curr = *best;
      if (current.size()!=0)
         current.erase(best);

      classe[curr]=2;
      neighbours.clear();
      for (sit=intvoisins[curr].begin();sit!=intvoisins[curr].end();sit++){
         if (classe[*sit]==0){
            neighbours.push_back(*sit);
            classe[*sit]=3;
         }
      }
      set<uint> voisins(extvoisins4[intcorr[curr]]);
      for (it=neighbours.begin();it!=neighbours.end();it++){
         Point3df n(intmesh[0].vertex()[*it]);
         pair<Point3df,int> res(plus_proche_point_sur_triangle(n,intmesh[0].normal()[*it],extmesh,voisins,extvoisins));
         if (Point3df(n-res.first).norm() > 5.0 || ((intmesh[0].normal()[*it].dot(res.first-n) < 0.0 && op==0) || (intmesh[0].normal()[*it].dot(res.first-n) > 0.0 && op==1))) {
            res = plus_proche_point_sur_triangle(n,intmesh[0].normal()[*it],extmesh,extvoisins);
         }
         matching_mesh[0].vertex()[*it] = res.first;
         intcorr[*it] = res.second;
         current.push_back(*it);
         classe[*it]=1;
         width[0].item(*it) = Point3df(n-res.first).norm();
         if (width[0].item(*it) > 20.0) width[0].item(*it) = 7.20;

      }
      neighbours.clear();
   }

   cout << endl << " done" << endl;
   cout << "Creating the surface.." << flush;
   for (uint i=0;i<intmesh[0].vertex().size();i++){
      for (uint j=0;j<3;j++){
         median[0].vertex()[i][j] = (intmesh[0].vertex()[i][j] + matching_mesh[0].vertex()[i][j])/2.0;
      }
   }
   cout << " done" << endl;
   pair<AimsSurfaceTriangle, TimeTexture<float> > output(median, width);
   return output;
}
 

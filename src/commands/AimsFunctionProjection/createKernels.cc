#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/distancemap/meshdistance_d.h>
#include <aims/distancemap/meshvoronoi.h>


using namespace aims;
using namespace carto;
using namespace std;
using namespace aims::meshdistance;


inline float calcule_distance(const Point3df &p, const Point3df &t){
   Point3df aux(p-t);
   return aux.dnorm();
}

inline float calcule_distance(const Point3df &p, const Point3d &t){
   Point3df aux(p);
   aux[0] -= t[0]; aux[1] -= t[1];  aux[2] -= t[2];   
   return aux.dnorm();
}

inline float geod_weight_function(const float &d, const float &d0){
   float weight;
//    float d0=6.0;
   if (d>d0 || d<0.0) weight = 0.0;
   else weight = -1/d0 * d + 1.0;
   return weight;
}

inline float cortical_distance(const Point3d &nv3, const Point3df &vsize, const Point3df &v, const Point3df &n, const float &width, const float &seuil){
   float weight;
   Point3df d(nv3[0]*vsize[0]-v[0],nv3[1]*vsize[1]-v[1],nv3[2]*vsize[2]-v[2]);
   float distance = calcule_distance(Point3df(nv3[0]*vsize[0],nv3[1]*vsize[1],nv3[2]*vsize[2]),v);
   if (d.dot(n) > 0.0) {
      if (distance < width) weight = 1.0;
      else weight = max(-0.95/seuil * (distance-width) + 0.95, 0.0);
   }
   else weight = max(-0.95/seuil * distance + 0.95, 0.0);
   return weight;
}

string auxmaskpath2;
void getMaskPath(string path){
   auxmaskpath2 = path;
}

float cortical_distance_via_tex(Point3df v, short vertex, Point3d nv3, Point3df vsize, Point3df n){
   
   TimeTexture<float> tex;
   Reader<TimeTexture<float> > r(auxmaskpath2);
   r.read(tex);
   float width = tex[0].item(vertex);
   float seuil = 2.0;
         
   return cortical_distance(nv3, vsize, v, n, width, seuil);
}

inline vector<uint> nearest_vertices(Point3df pf, AimsSurfaceTriangle &mesh, float rayon){
   vector<uint> v;
   for (uint i=0;i<mesh.vertex().size();i++){
      Point3df t(mesh[0].vertex()[i]);
      if (abs(t[0] - pf[0]) < rayon && abs(t[1] - pf[1]) < rayon && abs(t[2] - pf[2]) < rayon)
         v.push_back(i);
   }
   return v;
}

pair<int,float> plus_proche_point(Point3df p, AimsSurfaceTriangle &mesh){
   
   float dist=0.0,min = 5000.0; // un cerveau de 5 m, � va comme limite ?
   int index = -1;
   vector<uint> vertices(nearest_vertices(p, mesh, 5.0));
   for (uint i=0;i<vertices.size();i++){
      // Hypoth�e : les points p et q sont �moins de 5 mm de la surface corticale
      Point3df t(mesh[0].vertex()[vertices[i]]);
      dist = calcule_distance(p,t);
      if (dist < min && !(mesh[0].normal()[i][0] == mesh[0].normal()[i][1] && mesh[0].normal()[i][0] == mesh[0].normal()[i][2] && mesh[0].normal()[i][0] == 0.0)) {
         index = vertices[i];
         min = dist;
      }
   }
   return pair<int,float>(index, dist);
}

struct ltstr{
   bool operator()(Point3d p, Point3d q) const   {
      bool res = false;
      if (p[0]<q[0]) res = true;
      else if (p[0]==q[0] && p[1]<q[1]) res = true;
      else if (p[0]==q[0] && p[1]==q[1] && p[2]<q[2]) res = true;
      return res;
   }
};

void compute_kernel(AimsData<float> &kernels, uint time, AimsSurfaceTriangle &mesh, uint node, const vector<map<uint,float> > &voisins2, AimsData<long> &vertex, AimsData<short> &classe, const float &geod_decay, const float &norm_decay, Point3df &vsize, int size){

   map<uint,float>::iterator hit,ite;
   set<Point3d,ltstr>::iterator it;
   set<Point3d, ltstr> neighbours;
   map<float, Point3d> current;
   map<float, Point3d>::iterator aux;
   float distance,min;
   float wg,wn,sum=0.0;
   int x,y,z;
   uint processed= 0;
   Texture<float> tex1a;

   Point3df &p = mesh[0].vertex()[node];     // le noeud dont on calcule le noyau de convolution
   Point3d nv((int)((p[0]+vsize[0]/2.0)/vsize[0]),(int)((p[1]+vsize[1]/2.0)/vsize[1]),(int)((p[2]+vsize[2]/2.0)/vsize[2])); // le voxel contenant le noeud p
  
   current[0.0] = nv;
   Point3d c(nv[0]-(int)(size/2), nv[1]-(int)(size/2), nv[2]-(int)(size/2));   // le voxel c est le voxel situ�"en haut �gauche" du noyau de convolution
   uint a1=nv[0] - c[0],a2=nv[1] - c[1],a3=nv[2] - c[2];                       // on s'en sert principalement pour se placer dans le r��entiel du noyau
   vertex(a1,a2,a3,0) = node;          // le noeud qui correspond au voxel du centre est �idemment le noeud dont on calcule le noyau de convolution
   classe(a1,a2,a3,0) = 1;             // on remplit la case correspondante dans la matrice classe avec "1:voxel en cours de traitement"
   // et on remplit kernel avec le poids normal, consid�ant que le poids g�d�ique est �al �1
   kernels(a1,a2,a3,time) = 1.0;

   Point3d cour=nv;
   uint totalsize = size*size*size;
   
   while (!((*current.begin()).first==1.0) && processed != (uint)totalsize){   // tant que le nombre des voxels trait� est diff�ent du nombre de voxels total
                                                                        // si le poids maximal dans le front courant est �al �z�o alors c'est fini !
//       if (processed.size()%1000 == 0) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << processed.size() << "/" << size*size*size << ";" << current.size() << " " << (*current.begin()).first << flush ; //<< endl;
      // 1�e �ape on voit quel voxel du front courant est le plus pr� de la surface
      aux=current.begin();
      processed++;
      Point3d curr = (*aux).second;
      current.erase(aux);
      cour = curr;
      classe(curr[0]-c[0],curr[1]-c[1],curr[2]-c[2],0) = 2;
      neighbours.clear();
      map<uint,float> vois1 = voisins2[vertex(curr[0] - c[0],curr[1] - c[1],curr[2] - c[2],0)];  // on r�up�e le set de voisins dans lequel on va chercher le point
      
      for (x=-1;x<=1;x++){
         for (y=-1;y<=1;y++){
            for (z=-1;z<=1;z++){
               Point3d n(curr[0]+x,curr[1]+y,curr[2]+z);       // n est un voxel voisin (26-voisinage) du voxel courant
               if (n[0]-c[0]<size && n[1]-c[1]<size && n[2]-c[2]<size && n[0]-c[0]>=0 && n[1]-c[1]>=0 && n[2]-c[2]>=0 && classe(n[0]-c[0],n[1]-c[1],n[2]-c[2],0) == 0 ){
                  neighbours.insert(n);
                  classe(n[0]-c[0],n[1]-c[1],n[2]-c[2],0) = 3;    // on r�up�e les voisins qui n'ont pas ��trait� et on les marque comme �ant �traiter
               }
            }
         }
      }

      //3�e �ape on prend chaque voisin et on cherche dans les alentours de son surfacepoint correspondant le point le plus proche
      // le vecteur voxelcenter.surfacepoint est le plus grand.
      // ensuite on met �jour le voisin en question et on le passe dans current
      for (it=neighbours.begin();it!=neighbours.end();it++){
         distance=1000.0, min=1000.0;
         a1=(*it)[0]-c[0]; a2=(*it)[1]-c[1]; a3=(*it)[2]-c[2];
         hit = vois1.begin();
         for (ite=vois1.begin();ite!=vois1.end();ite++){
           Point3df &pt = mesh[0].vertex()[(*ite).first];
            distance = calcule_distance(pt,Point3df((*it)[0]*vsize[0],(*it)[1]*vsize[1],(*it)[2]*vsize[2]) );
            if (distance < min){
               min = distance;
               hit = ite;
            }
         }

         vertex(a1,a2,a3,0) = (*hit).first;    // le noeud associ�est celui dont on a retenu la distance
 
         classe(a1,a2,a3,0) = 1;       // on marque ces voxels comme �ant de classe "1: voxel en cours de traitement"

         wn = cortical_distance(*it,vsize,mesh[0].vertex()[(*hit).first],mesh[0].normal()[(*hit).first], 3.0, norm_decay); 
         wg = geod_weight_function((*hit).second, geod_decay);

         kernels(a1,a2,a3,time) = wg*wn ;
         current[1.0-wg*wn] = (*it);       // ce voxel trait�est ajout��la liste des voxels "courants", i.e. au front de propagation
         
         sum += kernels(a1,a2,a3,time);

      }
   }

   // 5�e �ape : normalisation au sein d'un seul noyau
   for (x=0;x<size;x++)
      for (y=0;y<size;y++)
         for (z=0;z<size;z++){
            kernels(x,y,z,time) /= sum;
         }
}

int kernel_index;
void get_kernelindex(int index){
   kernel_index = index;
}

AimsData<float> fast_marching_kernels(string meshpath, int size, Point3df vsize, float geod_decay, float norm_decay){
   int operation = 2;
   if (kernel_index == -1) operation = 0;
   Reader<AimsSurfaceTriangle> r(meshpath);
   AimsSurfaceTriangle mesh;
   r.read(mesh);
   uint nb_nodes = mesh[0].vertex().size();
   if (operation == 2) nb_nodes = 1;
   
   AimsData<float> kernel(size,size,size,nb_nodes);
   kernel.setSizeXYZT(vsize[0],vsize[1],vsize[2],1.0);

   // On pr�are voisins2 qui contient la liste des voisins de deuxi�e ordre de chaque noeud du maillage
   vector<set<uint> > voisins(SurfaceManip::surfaceNeighbours(mesh));
   set<uint>::iterator it,it2,it3,it4,it_min;
   set<uint> vois2;
   vector<map<uint,float> > voisins2(voisins.size());
   cerr << "Computing the 2nd-order neighbors :" << endl; 
   for (uint i=0;i<voisins.size();i++){
     if (i%1000==0) cout << "\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << voisins.size() << flush;
     for (it=voisins[i].begin();it!=voisins[i].end();it++){ // on parcourt les voisins de i
       voisins2[i][*it] = calcule_distance(mesh[0].vertex()[i], mesh[0].vertex()[*it]);
       vois2.clear();

       for (it2=voisins[*it].begin();it2!=voisins[*it].end();it2++){ // on parcourt les voisins de *it
         if (voisins[i].find(*it2) == voisins[i].end() && *it2 != i){ // il ne faut pas que *it2 soit un voisin de premier ordre de i ni i
           vois2.insert(*it2);
         }
       }
       for (it2=vois2.begin();it2!=vois2.end();it2++){
         float distance=1000.0;
         for (it3=voisins[*it2].begin();it3!=voisins[*it2].end();it3++){
           it4 = voisins[i].find(*it3);
           if (it4 != voisins[i].end()){ // on considère les voisins de *it2 qui sont aussi voisins de i
             distance = min(distance, calcule_distance(mesh[0].vertex()[*it2], mesh[0].vertex()[*it4]) + calcule_distance(mesh[0].vertex()[*it4], mesh[0].vertex()[i]));
           }
         }
         voisins2[i][*it2] = distance;
       }
     }
   }
   cout << "\b\b\b\b\b\b\b\b\b\b\b\b" << voisins.size() << "/" << voisins.size() << endl;
   cerr << "Computing kernel for mesh \"" << meshpath << "\", size = " << size << " / voxel size = " << vsize << endl;

   // On pr�are les diff�entes matrices qui peuvent
   AimsData<long> vertex(size,size,size,1);
   vertex.setSizeXYZT(vsize[0],vsize[1],vsize[2],1.0);
   AimsData<short> classe(size,size,size,1);
   classe.setSizeXYZT(vsize[0],vsize[1],vsize[2],1.0);
   
   for (int x=0;x<size;x++)
      for (int y=0;y<size;y++)
         for (int z=0;z<size;z++){
            vertex(x,y,z,0) = 0;
         }

         for (uint i=0;i<nb_nodes;i++){
            if (i%100==0)
               cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << nb_nodes << flush ;
            for (int x=0;x<size;x++)
               for (int y=0;y<size;y++)
                  for (int z=0;z<size;z++){
                     classe(x,y,z,0)=0;
                     kernel(x,y,z,i)=0.0;
                  }
            if (operation == 2)               
               compute_kernel(kernel, i, mesh, kernel_index, voisins2, vertex, classe, geod_decay, norm_decay,vsize, kernel.dimX());
            else 
              compute_kernel(kernel, i, mesh, i, voisins2, vertex, classe, geod_decay, norm_decay, vsize, kernel.dimX());
         }
   
         cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << nb_nodes << "/" << nb_nodes << endl;
         return kernel;
}

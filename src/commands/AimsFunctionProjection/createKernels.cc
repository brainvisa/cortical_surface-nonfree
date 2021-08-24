#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <aims/data/data.h>



using namespace aims;
using namespace carto;
using namespace std;


inline float calcule_distance ( const Point3df &p, const Point3df &t ) {
    Point3df aux(p-t);
    return aux.dnorm();
}

inline float calcule_distance ( const Point3df &p, const Point3d &t ) {
    Point3df aux(p);
    aux[0] -= t[0];
    aux[1] -= t[1];
    aux[2] -= t[2];
    return aux.dnorm();
}

inline float geod_weight_function ( const float &d, const float &d0 ) {
    float weight;
    if ( d>d0 || d<0.0 )
        weight = 0.0;
    else weight = -1/d0 * d + 1.0;
    return weight;
}

inline float cortical_distance ( const Point3d &nv3,
                                 const Point3df &vsize,
                                 const Point3df &v,
                                 const Point3df &n,
                                 const float &width,
                                 const float &seuil ) {
    float weight;
    Point3df d ( nv3[0] * vsize[0] - v[0], nv3[1] * vsize[1] - v[1], nv3[2] * vsize[2] - v[2] );
    float distance = calcule_distance ( Point3df ( nv3[0] * vsize[0], nv3[1] * vsize[1], nv3[2] * vsize[2] ), v );
    if ( d.dot(n) > 0.0 ) {
        if ( distance < width ) weight = 1.0;
        else weight = std::max( -0.95 / seuil * (distance-width) + 0.95, 0.0);
    }
    else weight = std::max( -0.95 / seuil * distance + 0.95, 0.0 );
    return weight;
}

std::string auxmaskpath2;

void getMaskPath ( std::string path ) {
    auxmaskpath2 = path;
}

float cortical_distance_via_tex ( Point3df v, short vertex, Point3d nv3, Point3df vsize, Point3df n ) {

    TimeTexture<float> tex;
    Reader<TimeTexture<float> > r ( auxmaskpath2 );
    r.read(tex);
    float width = tex[0].item(vertex);
    float seuil = 2.0;

    return cortical_distance ( nv3, vsize, v, n, width, seuil );
}

inline std::vector<uint> nearest_vertices ( Point3df pf, AimsSurfaceTriangle &mesh, float rayon ) {
    std::vector<uint> v;
    for ( uint i = 0 ; i < mesh.vertex().size() ; i++ ) {
        Point3df t ( mesh[0].vertex()[i] );
        if ( abs(t[0] - pf[0]) < rayon
                && abs(t[1] - pf[1]) < rayon
                && abs(t[2] - pf[2]) < rayon )
            v.push_back(i);
    }
    return v;
}

std::pair<int, float> plus_proche_point ( Point3df p, AimsSurfaceTriangle &mesh ) {

    float dist = 0.0, min = 5000.0; // un cerveau de 5 m, � va comme limite ?
    int index = -1;
    std::vector<uint> vertices( nearest_vertices(p, mesh, 5.0) );
    for ( uint i = 0 ; i < vertices.size() ; i++ ) {
        // Hypoth�e : les points p et q sont �moins de 5 mm de la surface corticale
        Point3df t( mesh[0].vertex()[vertices[i]] );
        dist = calcule_distance ( p, t );
        if ( dist < min
                &&
             ! ( mesh[0].normal()[i][0] == mesh[0].normal()[i][1]
                  && mesh[0].normal()[i][0] == mesh[0].normal()[i][2]
                  && mesh[0].normal()[i][0] == 0.0 ) ) {
            index = vertices[i];
            min = dist;
        }
    }
    return std::pair<int,float>( index, dist );
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

void compute_kernel ( AimsData<float> &kernels,
                      uint time,
                      AimsSurfaceTriangle &mesh,
                      uint node,
                      const std::vector< std::map<uint,float> > &voisins2,
                      AimsData<long> &vertex,
                      AimsData<short> &classe,
                      const float &geod_decay,
                      const float &norm_decay,
                      Point3df &vsize,
                      int size)
{

    // cout << "compute_kernel " << node << endl;
    std::map< uint, float >::iterator hit, ite;
    std::set< Point3d, ltstr >::iterator it;
    std::set< Point3d, ltstr> neighbours;
    std::map< float, Point3d > current;
    std::map< float, Point3d >::iterator aux;
    float distance, min;
    float wg, wn, sum = 0.0;
    int x, y, z;
    uint processed = 0;
    Texture<float> tex1a;

    Point3df &p = mesh[0].vertex()[node];     // le noeud dont on calcule le noyau de convolution
    Point3d nv( (int) ( (p[0] + vsize[0]/2.0) / vsize[0]),
               (int) ( (p[1] + vsize[1]/2.0) / vsize[1]),
               (int) ( (p[2] + vsize[2]/2.0) / vsize[2])); // le voxel contenant le noeud p

    current[0.0] = nv;

    Point3d c ( nv[0]-(int)(size/2),
               nv[1]-(int)(size/2),
               nv[2]-(int)(size/2) );   // le voxel c est le voxel situ�"en haut �gauche" du noyau de convolution

    uint a1 = nv[0] - c[0],
        a2 = nv[1] - c[1],
        a3 = nv[2] - c[2];                       // on s'en sert principalement pour se placer dans le r��entiel du noyau

    vertex ( a1, a2, a3, 0 ) = node;          // le noeud qui correspond au voxel du centre est �idemment le noeud dont on calcule le noyau de convolution
    classe ( a1, a2, a3, 0 ) = 1;             // on remplit la case correspondante dans la matrice classe avec "1:voxel en cours de traitement"
    // et on remplit kernel avec le poids normal, consid�ant que le poids g�d�ique est �al �1
    kernels (a1, a2, a3, time) = 1.0;
    sum += 1.0;

    Point3d cour = nv;
    uint totalsize = size * size * size;

    while ( !current.empty() && current.begin()->first != 1.0 && processed != totalsize )
    {   // tant que le nombre des voxels trait� est diff�ent du nombre de voxels total
                                                                        // si le poids maximal dans le front courant est �al �z�o alors c'est fini !
        //       if (processed.size()%1000 == 0) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << processed.size() << "/" << size*size*size << ";" << current.size() << " " << (*current.begin()).first << flush ; //<< endl;
        // 1�e �ape on voit quel voxel du front courant est le plus pr� de la surface
        aux = current.begin();
        processed++;
//         std::cout << "processed: " << processed << " / " << totalsize << ", " << current.begin()->first << ", n: " << current.size() << std::endl;
        Point3d curr = (*aux).second;
        current.erase(aux);
        cour = curr;
        classe ( curr[0]-c[0], curr[1]-c[1], curr[2]-c[2], 0 ) = 2;
        neighbours.clear();
        std::map< uint, float > vois1 =
            voisins2 [ vertex(curr[0]-c[0], curr[1]-c[1], curr[2]-c[2], 0) ];  // on r�up�e le set de voisins dans lequel on va chercher le point

        for ( x = -1 ; x <= 1 ; x++ ) {
            for ( y = -1 ; y <= 1 ; y++ ) {
                for ( z = -1 ; z <= 1 ; z++ ) {
                    Point3d n ( curr[0]+x, curr[1]+y, curr[2]+z );       // n est un voxel voisin (26-voisinage) du voxel courant
                    if ( n[0]-c[0] < size
                           && n[1]-c[1] < size
                           && n[2]-c[2] < size
                           && n[0]-c[0] >= 0
                           && n[1]-c[1] >= 0
                           && n[2]-c[2] >= 0
                           && classe ( n[0]-c[0], n[1]-c[1], n[2]-c[2], 0 ) == 0 ) {
                                neighbours.insert(n);
                                    classe ( n[0]-c[0], n[1]-c[1], n[2]-c[2], 0 ) = 3;    // on r�up�e les voisins qui n'ont pas ��trait� et on les marque comme �ant �traiter
                    }
                }
            }
        }

        //3�e �ape on prend chaque voisin et on cherche dans les alentours de son surfacepoint correspondant le point le plus proche
        // le vecteur voxelcenter.surfacepoint est le plus grand.
        // ensuite on met �jour le voisin en question et on le passe dans current
//         cout << "neighbours: " << neighbours.size() << endl;
        for ( it = neighbours.begin() ; it != neighbours.end() ; it++ )
        {
            distance = 1000.0, min = 1000.0;
            a1 = (*it)[0]-c[0];
            a2 = (*it)[1]-c[1];
            a3 = (*it)[2]-c[2];
            hit = vois1.begin();
            for ( ite = vois1.begin() ; ite != vois1.end() ; ite++ ) {
                Point3df &pt = mesh[0].vertex()[(*ite).first];
                distance = calcule_distance ( pt,
                                              Point3df( (*it)[0] * vsize[0], (*it)[1] * vsize[1], (*it)[2] * vsize[2]) );
                if ( distance < min ) {
                    min = distance;
                    hit = ite;
                }
            }

            vertex ( a1, a2, a3, 0 ) = (*hit).first;    // le noeud associ�est celui dont on a retenu la distance

            classe ( a1, a2, a3, 0 ) = 1;       // on marque ces voxels comme �ant de classe "1: voxel en cours de traitement"

            wn = cortical_distance ( *it,
                                     vsize,
                                     mesh[0].vertex() [(*hit).first],
                                     mesh[0].normal() [(*hit).first],
                                     3.0,
                                     norm_decay );
            //std::cout << "wn:" << wn <<  " " << std::flush ;

            wg = geod_weight_function( (*hit).second, geod_decay );

            //std::cout << "dis:" << (*hit).second << " wg:" << wg <<  " " << std::flush ;

            kernels ( a1, a2, a3, time ) = wg * wn ;            

            //std::cout << "wg*wn:" << wg*wn <<  " " << std::flush ;


            current [ 1.0 - wg * wn ] = (*it);       // ce voxel trait�est ajout��la liste des voxels "courants", i.e. au front de propagation

            sum = sum + kernels ( a1, a2, a3, time );
            //std::cout << "sum:" << sum << " " << std::flush;

        }
    }
    std::cout << "SUM:" << sum << std::endl;
    // 5�e �ape : normalisation au sein d'un seul noyau
    if ( ! ( sum > 0.0) ) {
        std::cout << "sum: " << sum << std::endl;
        std::cout << "time: " << time << std::endl;        
    }
    assert( sum > 0.0 );
    float mini = 5.0, maxi = -5.0;
    for ( x = 0 ; x < size ; x++ )
        for ( y = 0 ; y < size ; y++ )
            for ( z = 0 ; z < size ; z++ ) {
                kernels ( x, y, z, time ) /= sum;
                if ( kernels ( x, y, z, time ) < mini )
                    mini = kernels ( x, y, z, time );
                if ( kernels ( x, y, z, time ) > maxi )
                    maxi = kernels ( x, y, z, time );
            }
}

int kernel_index;

void get_kernelindex ( int index ) {
    kernel_index = index;
}

std::map<uint, float> LocalMeshDistanceMap ( AimsSurface<3,Void> *mesh,
                                             const std::vector< std::set<unsigned> >    &neighbours,
                                             uint max_node,
                                             float ind ) {
    unsigned i;
    std::multimap<float,unsigned> front1, front2;
    std::multimap<float,unsigned> *cfront = &front1, *nfront = &front2, *tmpf;
    std::multimap<float,unsigned>::iterator    iv, fv;
    std::set<unsigned> neigh_local;
    std::set<unsigned>::const_iterator in, fn;
    float d, d2, l;
    Point3df pos;
    float dist = 0;
    std::map< uint, float > distances;
    distances[max_node] = 0.0;
    std::map<uint,float>::iterator it;
    for ( it = distances.begin() ; it != distances.end() ; it++ )
        if ( it->second == 0.0 )
            front1.insert( std::pair<float,unsigned>( 0, (*it).first ) );

    while( dist <= ind ) {
        nfront->clear();
        neigh_local.clear();

        for( iv = cfront->begin(), fv = cfront->end() ; iv != fv ; ++iv ) {

            i = (*iv).second;
            d = (*iv).first;
            for( in = neighbours[i].begin(), fn = neighbours[i].end() ; in != fn ; ++in ) {

                //d2 = tex.item( *in );
                if ( distances.find(*in) == distances.end() )
                    distances[*in] = FLT_MAX;
                d2 = distances[*in];

                pos = mesh->vertex()[i] - mesh->vertex()[*in];
                l = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );

                if( d2 > d + l ) {
                    //tex.item( *in ) = d + l;
                    distances[*in] = d+l;
                    neigh_local.insert( *in );

                    dist = distances [ *in ];
                }

            }

        }

        for ( in = neigh_local.begin(), fn = neigh_local.end() ; in != fn ; ++in )
            nfront->insert( std::pair<float,unsigned>( distances[*in], *in ) );

        tmpf = cfront;
        cfront = nfront;
        nfront = tmpf;
    }

    front1.clear();
    front2.clear();
    (*cfront).clear();
    (*nfront).clear();
    (*tmpf).clear();

    return ( distances );
}


AimsData<float> fast_marching_kernels ( std::string meshpath,
                                        int size,
                                        Point3df vsize,
                                        float geod_decay,
                                        float norm_decay )
{
    assert( geod_decay > 0.0 );
    assert( norm_decay > 0.0 );
    Reader<AimsSurfaceTriangle> r ( meshpath );
    AimsSurfaceTriangle mesh;
    r.read(mesh);

    uint nb_nodes = mesh[0].vertex().size();
    if ( kernel_index != -1 )
        nb_nodes = 1;

    AimsData<float> kernel ( size, size, size, nb_nodes );
    kernel.setSizeXYZT ( vsize[0], vsize[1], vsize[2], 1.0 );

    // On pr�are voisins2 qui contient la liste des voisins de deuxi�e ordre de chaque noeud du maillage
    std::vector< std::set<uint> > voisins ( SurfaceManip::surfaceNeighbours(mesh) );
    std::set<uint>::iterator it, it2, it3, it4, it_min;
    std::set<uint> vois2;
    std::vector<std::map<uint,float> > voisins2 ( voisins.size() );

//    std::cerr << "Computing the 2nd-order neighbors :" << std::endl;
//    for ( uint i = 0 ; i < voisins.size() ; i++ ) {
//        if ( i%1000 == 0 ) std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << voisins.size() << std::flush;
//        for ( it = voisins[i].begin() ; it != voisins[i].end() ; it++ ) { // on parcourt les voisins de i
//            voisins2[i][*it] = calcule_distance( mesh[0].vertex()[i], mesh[0].vertex()[*it] );
//            vois2.clear();
//
//            for ( it2 = voisins[*it].begin() ; it2 != voisins[*it].end() ; it2++ ) { // on parcourt les voisins de *it
//                if ( voisins[i].find(*it2) == voisins[i].end() && *it2 != i ) { // il ne faut pas que *it2 soit un voisin de premier ordre de i ni i
//                    vois2.insert(*it2);
//                }
//            }
//            for ( it2 = vois2.begin() ; it2 != vois2.end() ; it2++ ) {
//                float distance = 1000.0;
//                for ( it3 = voisins[*it2].begin() ; it3 != voisins[*it2].end() ; it3++ ) {
//                    it4 = voisins[i].find(*it3);
//                    if ( it4 != voisins[i].end() ) { // on considère les voisins de *it2 qui sont aussi voisins de i
//                        distance = std::min ( distance,
//                             calcule_distance ( mesh[0].vertex()[*it2],
//                                     mesh[0].vertex()[*it4]) + calcule_distance ( mesh[0].vertex()[*it4], mesh[0].vertex()[i] ) );
//                    }
//                }
//                voisins2[i][*it2] = distance;
//            }
//        }
//    }
//    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << voisins.size() << "/" << voisins.size() << std::endl;
    std::cerr << "Computing kernel for mesh " << meshpath << ", size = " << size << " / voxel size = " << vsize << std::endl;

    voisins2.clear();
    for ( uint i = 0 ; i < voisins.size() ; i++ ) {
        if ( i%1000 == 0 ) std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << voisins.size() << std::flush;
        voisins2.push_back(LocalMeshDistanceMap( &(mesh[0]), voisins, i, 2.5) );
        assert(voisins2[voisins2.size()-1].size() > 0);
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << voisins.size() << "/" << voisins.size() << std::endl;


    // On pr�are les diff�entes matrices qui peuvent
    AimsData<long> vertex ( size, size, size, 1 );
    vertex.setSizeXYZT ( vsize[0], vsize[1], vsize[2], 1.0 );
    AimsData<short> classe ( size, size, size, 1 );
    classe.setSizeXYZT ( vsize[0], vsize[1], vsize[2], 1.0 );

    for ( int x = 0 ; x < size ; x++ )
        for ( int y = 0 ; y < size ; y++ )
            for ( int z = 0 ; z < size ; z++ ) {
                vertex ( x, y, z, 0 ) = 0;
            }

    for ( uint i = 0 ; i < nb_nodes ; i++ )
    {
        if ( i%100 == 0 )
            std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << nb_nodes << std::flush ;
        for ( int x = 0 ; x < size ; x++ )
            for ( int y = 0 ; y < size ; y++ )
                for ( int z = 0 ; z < size ; z++ ) {
                    classe ( x, y, z, 0 ) = 0;
                    kernel ( x, y, z, i ) = 0.0;
                }
        if ( kernel_index != -1 )
        {
          if( i == kernel_index )
            compute_kernel( kernel, i, mesh, kernel_index, voisins2, vertex, classe, geod_decay, norm_decay, vsize, kernel.dimX() );
        }
        else
            compute_kernel( kernel, i, mesh, i, voisins2, vertex, classe, geod_decay, norm_decay, vsize, kernel.dimX() );
    }

    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << nb_nodes << "/" << nb_nodes << std::endl;

    /*
    float mini = 5.0, maxi = -5.0;
    for ( uint x = 0 ; x < size ; x++ )
        for ( uint y = 0 ; y < size ; y++ )
            for ( uint z = 0 ; z < size ; z++ ) {
                if ( kernel ( x, y, z, 3114 ) < mini )
                    mini = kernel ( x, y, z, 3114 );
                if ( kernel ( x, y, z, 3114 ) > maxi )
                    maxi = kernel ( x, y, z, 3114 );
            }
    */

    if ( kernel_index != -1 )
    {
        std::cout << " node coordinates : (index = " << kernel_index << ") " << mesh[0].vertex()[kernel_index] << std::endl;
        Point3df &p = mesh[0].vertex()[kernel_index];     // le noeud dont on calcule le noyau de convolution
        Point3d nv( (int) ( (p[0] + vsize[0]/2.0) / vsize[0]),
                   (int) ( (p[1] + vsize[1]/2.0) / vsize[1]),
                   (int) ( (p[2] + vsize[2]/2.0) / vsize[2])); // le voxel contenant le noeud p

        Point3d c ( (nv[0]-(int)(size/2)) * vsize[0],
                   (nv[1]-(int)(size/2)) * vsize[1],
                   (nv[2]-(int)(size/2)) * vsize[2]);
        std::cout << " transfo to get in the kernel ref : " << c << std::endl;

    }
    return kernel;
}

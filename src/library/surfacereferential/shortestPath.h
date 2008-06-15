#ifndef AIMS_SHORTESTPATH_H
#define AIMS_SHORTESTPATH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

namespace aims
{
#define NOEUD struct _noeud
#define ARC struct _arc

NOEUD
{
     int nom;
     int marque;
     ARC *arcs;                    /* liste d'arcs */
     NOEUD *suivant;               /* dans le graphe */
};

ARC
{
     NOEUD *noeud;
     int cout;
     ARC *suivant;
};

typedef struct
{
     NOEUD *noeuds;
     int n_noeuds;
} GRAPHE;

#define CHEMIN struct _chemin

struct _chemin
{
     NOEUD *noeud;
     int cout;
     CHEMIN *suivant;
};

template<typename Val> class GraphPath
{

     public:
     
     GraphPath() {plus_court=NULL; longueur =0;}

     NOEUD * nouveau_noeud(const int nom);
     NOEUD * insere_noeud(GRAPHE *graphe, NOEUD *noeud);
     NOEUD * trouve_noeud(const GRAPHE *graphe, int nom);
     NOEUD * ajoute_noeud(GRAPHE *graphe, const int nom);
     void ajoute_arc(NOEUD *depart, const NOEUD *arrivee, const int cout);
     GRAPHE * lecture_graphe(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value);
     void copie_chemin(const CHEMIN *chemin);
     void traite_chemin(const CHEMIN *chemin);
     void tous_chemins(const GRAPHE *graphe, NOEUD *depart,const NOEUD *arrivee, CHEMIN *chemin, const int cout);
     TimeTexture<Val> imprime_plus_court(int size, Val value);
     TimeTexture<Val> process(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr);
     float getLongueur(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr);
     TimeTexture<Val> cleanPath(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh);
     
     protected:
     
     CHEMIN *plus_court;
     int longueur;
};

/* Fabrique un nouveau noeud */
template<typename Val>
  struct aims::_noeud * GraphPath<Val>::nouveau_noeud(const int nom)
{
     NOEUD *noeud= (NOEUD *) malloc(sizeof(NOEUD));
     if (noeud != NULL)
     {
          noeud->nom= nom;
          noeud->arcs= NULL;
          noeud->suivant= NULL;
          noeud->marque=0;
     }

     return noeud;
}

/* Insere un noeud dans la liste de noeuds du graphe */
template<typename Val>
  struct aims::_noeud * GraphPath<Val>::insere_noeud(GRAPHE *graphe, NOEUD *noeud)
{
     /* Insertion */
     noeud->suivant= graphe->noeuds;
     graphe->noeuds= noeud;
     return noeud;
}

template<typename Val>
  struct aims::_noeud * GraphPath<Val>::trouve_noeud(const GRAPHE *graphe, const int nom)
{
     NOEUD *noeud;

     for (noeud= graphe->noeuds; noeud != NULL; noeud= noeud->suivant)
          if (nom == noeud->nom)
               break;

     return noeud;
}

/* Cherche si un noeud d'apres son nom, en creee un nouveau s'il
n'existe pas */
template<typename Val>
  struct aims::_noeud * GraphPath<Val>::ajoute_noeud(GRAPHE *graphe, const int nom)
{
     NOEUD *noeud= trouve_noeud(graphe, nom);

     if (noeud == NULL)
     {
          noeud= nouveau_noeud(nom);
          if (noeud != NULL)
          {
               graphe->n_noeuds++;
               return insere_noeud(graphe, noeud);
          }
     }

     return noeud;
}

template<typename Val>
void GraphPath<Val>::ajoute_arc(NOEUD *depart, const NOEUD *arrivee, const int cout)
{
     ARC *arc= (ARC *) malloc(sizeof(ARC));
     if (arc != NULL)
     {
          arc->noeud= (NOEUD *) arrivee;
          arc->cout= cout;
          /* Insertion de l'arc dans la liste des arcs */
          arc->suivant= depart->arcs;
          depart->arcs= arc;
     }
}

template<typename Val>
GRAPHE * GraphPath<Val>::lecture_graphe(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value)
{
     GRAPHE *graphe= (GRAPHE *) malloc(sizeof(GRAPHE));

     std::vector<std::set<uint> > neigh;
     std::set<uint>::const_iterator itneigh;
     neigh = SurfaceManip::surfaceNeighbours( initmesh );
     int numero_ligne =0;
     //std::cout<<"YEAH1"<<std::endl;
     if (graphe != NULL)
     {
          graphe->noeuds= NULL;
          graphe->n_noeuds= 0;
          for(unsigned i=0; i<tex[0].nItem(); i++)
          {
               itneigh=neigh[i].begin();
               for(;itneigh!=neigh[i].end();++itneigh)
               {
                    if(tex[0].item(*itneigh)==value)
                    {
                         int nom1;
                         int nom2;
                         int cout_ch;
                         nom1=i;
                         nom2=(*itneigh);
                         cout_ch=1;

                         numero_ligne++;

                         NOEUD *n1, *n2;
                         int cout=0;

                         n1= ajoute_noeud(graphe, nom1);
                         n2= ajoute_noeud(graphe, nom2);
                         ajoute_arc(n1, n2, cout);
                    }
               }
          }
     }
     return graphe;
}


// GLOBAL VARIABLE (!!!!!!) WAS THERE

template<typename Val>
void GraphPath<Val>::copie_chemin(const CHEMIN *chemin)
{
     longueur= 0;
     while (chemin)
     {
          /* Copie des etapes */
          plus_court[longueur]= *chemin;
          chemin= chemin->suivant;
          longueur++;
     }
}

template<typename Val>
void GraphPath<Val>::traite_chemin(const CHEMIN *chemin)
{
     if (longueur == 0 || chemin->cout < plus_court->cout)
          copie_chemin(chemin);
}

template<typename Val>
void GraphPath<Val>::tous_chemins(const GRAPHE *graphe, NOEUD *depart,
     const NOEUD *arrivee, CHEMIN *chemin, const int cout)
{
     CHEMIN etape;
     etape.noeud= depart;
     etape.suivant= chemin;
     etape.cout= cout;

     if (depart == arrivee)
     {
          traite_chemin(&etape);
     }
     else
     {
          ARC *arc;

          depart->marque= 1;

          for (arc= depart->arcs; arc != NULL; arc= arc->suivant)
          {
               if (arc->noeud->marque == 0)
               {
                    if (longueur == 0 || etape.cout < plus_court->cout)
                         tous_chemins(graphe, arc->noeud, arrivee, &etape, cout + arc->cout);
               }
          }

          depart->marque= 0;
/*          std::cout << "+" << std::flush;*/
     }
}

template<typename Val>
TimeTexture<Val> GraphPath<Val>::imprime_plus_court(int size, Val value)
{
     int index;
     //std::cout<<"ImprimePlusCourt et size="<<size<<std::endl;
     TimeTexture<Val> result(1,size);
     for(int i=0; i<size; i++)
          result[0].item(i)=0;

     //std::cout<<"ImprimePlusCourt 2 et longueur="<<longueur<<std::endl;
     //std::cout<<plus_court[longueur-1].cout<<std::endl;

     for (index= longueur - 1; index >= 0; index--)
     {
          //std::cout<<plus_court[index].noeud->nom<<std::endl;
          //std::cout<<" - "<<plus_court[index].cout<<std::endl;

          result[0].item( plus_court[index].noeud->nom )=value;
     }
     //std::cout<<"voil�"<<std::endl;
     return result;
}

template<typename Val>
TimeTexture<Val> GraphPath<Val>::process(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr)
{
     TimeTexture<Val> texDirty, texFinal;
     GRAPHE *graphe= lecture_graphe(tex, initmesh, value);
     NOEUD *depart, *arrivee;
     plus_court= (CHEMIN *) calloc(graphe->n_noeuds, sizeof(CHEMIN));

     longueur= 0;

     // special case that generates a bug : the two points are neighbors

     std::vector<std::set<uint> > neigh;
     std::set<uint> depNeigh;
     neigh = SurfaceManip::surfaceNeighbours( initmesh );
     depNeigh=neigh[dep];
     if (depNeigh.find(arr) != depNeigh.end())
     {
/*          std::cout << "\t\t\t Shortest::process -> neighbors !!!!" << std::endl;*/
          uint ns=tex[0].nItem();
          texFinal=TimeTexture<Val>(1, ns);
          for (uint i=0; i<ns; i++)
          {
               if ((i==arr) || (i==dep))
                    texFinal[0].item(i)=(Val) 1;
               else
                    texFinal[0].item(i)=(Val) 0;
          }
     }
     else
     {
     depart= trouve_noeud(graphe, dep);
     arrivee= trouve_noeud(graphe, arr);
/*     std::cout << "\t\t\t Shortest::process -> Tous chemins" << std::endl;*/
     tous_chemins(graphe, depart, arrivee, NULL, 0);
    
//      std::cout << "\t\t\t Shortest::process -> Imprime plus court" << std::endl;

     texDirty=imprime_plus_court(tex[0].nItem(), value);

//      std::cout << "\t\t\t Shortest::process -> Clean" << std::endl;

     delete plus_court;
     texFinal=cleanPath(texDirty, initmesh);  // hack to solve a bug (Olivier)
//      // Shortest path sometimes include triangles. I have not programmed this and I cannot
//      // find the problem so I decided to work around it with a postprocessing of the texture
//      // obviously this is a dirty hack;
//

/*     std::cout << "\t\t\t Shortest::process -> OK" << std::endl;*/
     }

     return texFinal;
}

template<typename Val> float GraphPath<Val>::getLongueur(TimeTexture<Val> & tex, AimsSurfaceTriangle & initmesh, Val value, int dep, int arr)
{
     process(tex, initmesh, value, dep, arr);
     return(longueur);
}


template<typename Val> TimeTexture<Val> GraphPath<Val>::cleanPath(TimeTexture<Val> & path, AimsSurfaceTriangle & initmesh)
{

     uint i, j; 
     uint ns=path[0].nItem();
     TimeTexture<Val> newPath(1, ns);
     std::vector<std::set<uint> > neigh;
     std::set<uint> setPath;
     std::set<uint> remove;
     std::map<uint, int> nbIn;
     neigh = SurfaceManip::surfaceNeighbours( initmesh );
     std::set<uint> listN, listN2;
     std::set<uint>::iterator itN, itN2; 
     
     
//      std::cout << "DEBUG : Start cleaning" << std::endl;

     for (i=0; i<ns; i++)
          newPath[0].item(i)=path[0].item(i);
          
//      std::cout << "DEBUG : step1" << std::endl;

     for (i=0; i<ns; i++)
     {
          if (path[0].item(i)!=0)
          {
               setPath.insert(i);
               listN=neigh[i];
               itN=listN.begin();
               int countN=0;
               for ( ; itN!=listN.end(); ++itN)
               {
                    if (path[0].item(*itN)!=0)
                         countN++;
               }
               nbIn[i]=countN;
          }
     }
     
//      std::cout << "DEBUG : step2" << std::endl;

     std::set<uint>::iterator setN=setPath.begin();
     for ( ; setN!=setPath.end(); ++setN)
     {
          i=*setN;
          if (nbIn[i]==3)
          {
               listN=neigh[i];
               itN=listN.begin();
               for ( ; itN!=listN.end(); ++itN)
               {
                    j=*itN;
                    if ((nbIn.find(j)!=nbIn.end()) && (nbIn[j]==3))
                    {
                         listN2=neigh[j];
                         itN2=listN2.begin();
                         for ( ; itN2!=listN2.end(); ++itN2)
                         {
                              if ((nbIn.find(*itN2)!=nbIn.end()) && (nbIn[*itN2]==2) && (listN.find(*itN2)!=listN.end()))
                                   remove.insert(*itN2);
                         }
                    }
               }
          }
     }
     
//      std::cout << "DEBUG : step3" << std::endl;

     std::set<uint>::iterator removeIt;
     for (removeIt=remove.begin(); removeIt!=remove.end(); ++removeIt)
     {
          newPath[0].item(*removeIt)=(Val) 0;
     }
     
//      std::cout << "DEBUG : writing cleaning results" << std::endl;
//      Writer<TimeTexture<Val> >  tex1W( "/home/olivier/beforeClean.tex" );
//      tex1W.write( path );
//      Writer<TimeTexture<Val> >  tex2W( "/home/olivier/afterClean.tex" );
//      tex2W.write( newPath );     
     return newPath;
}


}
#endif

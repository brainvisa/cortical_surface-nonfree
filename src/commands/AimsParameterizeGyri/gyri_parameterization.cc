#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/curv.h>
#include <aims/io/reader.h>
#include <stdio.h>
#include "verif_operations.h"
#include "compconn_operations.h"
#include "intersec_operations.h"
#include "misctex_operations.h"
#include "model_operations.h"
#include "vertices_operations.h"
#include "mesh_operations.h"
#include "constraints_operations.h"
#include "vector_operations.h"

using namespace std;
using namespace aims;


TimeTexture<float> GyriParamTexture(AimsSurface<3,Void> &inMesh, Texture<short> &inTex, const Texture<float> &spmTex,
      const vector<uint> &corres, const pair<vector<vector<uint> >,vector<vector<uint> > > &diffMod, uint option, float criter, float dt){

   bool step = true;
   cerr << "... verifying mesh and texture are matching ..." << flush;
   if (verifMeshTexMatch(inMesh, inTex))
      cout << " done" << endl;
   else {
      step = false;
      cerr << endl << "Attention : mesh and texture do not match." << endl;
   }

   cout << "... verifying table and texture are matching ..." << flush;
   if (verifTableTexMatch(inTex, corres))
      cout << " done" << endl;
   else{
      step = false;
      cerr << endl << "Attention : table and texture do not match." << endl;
   }

   AimsSurfaceTriangle mesh;
   mesh[0] = inMesh;
   vector<std::set<uint> > voisins(SurfaceManip::surfaceNeighbours(mesh));
   nettoyerTaches(inTex,voisins);


   cerr << "... verifying gyral parcellation matches intersections model ..." << flush;
   if (verifIntersectionsModel(voisins, inTex, corres))
      cout << " done" << endl;
   else {
      //step = false;
      cerr << endl << "Attention : the gyral parcellation does not match the intersections model." << endl;
   }

   TimeTexture<float> final(3,0);

   if (step){

         cerr << "DEBUG : ENTERING THE LOOP" << endl;

         makeGenericTexture(inTex,corres);

         vector<pair<pair<vector<uint>, vector<uint> >, pair<vector<uint>, vector<uint> > > > points;

         AimsSurface<3, Void> mesh_base = inMesh;
         map<unsigned, set<pair<unsigned,float> > > poids = AimsMeshWeightFiniteElementLaplacian (mesh_base, 0.98);


         for (uint j=0;j<inTex.nItem();j++){
            final[0].push_back((float)inTex.item(j));
            final[1].push_back(0.0);
            final[2].push_back(0.0);
         }
         vector<vector<uint> > intersMod(diffMod.first);
         vector<vector<uint> > constrMod(diffMod.second);

         cerr << "DEBUG : SECOND STEP" << endl;

         for (uint j=0; j<intersMod.size();j++){
            uint code;
            short g = intersMod[j][0];
            printf("Gyrus %d :\n", g);
            pair<vector<uint>, vector<uint> > hautBas, gaucheDroite;
            Point3d gHaut(getGyrusModel(g, intersMod).first);
            Point3d gBas(getGyrusModel(g, intersMod).second);

            cerr << "... verifying existence of gyri refered by the model ..." << flush;
            if (verifGyriExistence(g, gHaut, gBas, corres)){
               cerr << " done" << endl;

               vector<vector<uint> > vertices(sortVertices(g, voisins, inTex));
               vector<uint> corr;
               AimsSurfaceTriangle gyrusMesh = getGyrusMesh(inMesh, vertices[0], corr);
               map<unsigned, set<pair<unsigned,float> > > poidsGyrus = getGyrusWeight(poids,vertices[0],corr);

               code = lookUpIntersectionCase(g, gHaut[0], gHaut[1], gHaut[2], voisins,inTex);
               printf("Code intersections : (%d) et ",code);
               hautBas.first = getIntersection(code,g, gHaut[0], gHaut[1], gHaut[2], voisins,inTex);

               code = lookUpIntersectionCase(g, gBas[0], gBas[1], gBas[2], voisins,inTex);
               printf("(%d)\n",code);

               hautBas.second = getIntersection(code,g, gBas[0], gBas[1], gBas[2], voisins, inTex);
               if (!(hautBas.first.empty() || hautBas.second.empty())){
                  gaucheDroite = getOppositeSides(hautBas, vertices , voisins, inTex);
                  gaucheDroite = sortRightLeft(inMesh, hautBas, gaucheDroite, voisins);
                  pair<pair<vector<uint>, vector<uint> >, pair<vector<uint>, vector<uint> > > intersection(hautBas, gaucheDroite);
                  points.push_back(intersection);
               }
               else {
                  if (hautBas.first.empty()) printf("Attention : détection intersection \"haut\" foireuse.\n");
                  if (hautBas.second.empty()) printf("Attention : détection intersection \"bas\" foireuse.\n");
               }
               if (!(gaucheDroite.first.empty() || gaucheDroite.second.empty())){

                  vector<uint> haut(getCorresVector(hautBas.first,corr));
                  vector<uint> bas(getCorresVector(hautBas.second,corr));
                  vector<uint> gauche(getCorresVector(gaucheDroite.first, corr));
                  vector<uint> droite(getCorresVector(gaucheDroite.second, corr));

                  Texture<double> verticDiff, horizDiff;
                  AimsSurface<3,Void> flatMesh;

                  if (option != 5){

                     verticDiff= diffusion(poidsGyrus, gyrusMesh[0], haut, bas, *new vector<pair<vector<uint>,short> >, -1, vertices[0], corr, criter, dt);
                     horizDiff = diffusion(poidsGyrus, gyrusMesh[0], gauche,droite,*new vector<pair<vector<uint>,short> > , -1, vertices[0], corr, criter, dt);

                     for (uint i2=0;i2<vertices[0].size();i2++)
                        final[1].item(vertices[0][i2]) = (float) verticDiff.item(corr[vertices[0][i2]]);
                     for (uint i2=0;i2<vertices[0].size();i2++)
                        final[2].item(vertices[0][i2]) = (float) horizDiff.item(corr[vertices[0][i2]]);

                     flatMesh = getFlatMesh(gyrusMesh[0], vertices[0], corr, final);
                  }

                  pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > constraints;

                  switch(option){
                     case 1:
                     constraints = getConstraints(g,constrMod,hautBas,gaucheDroite,corr,voisins,inTex,gyrusMesh,verticDiff,horizDiff);
                     break;
                     case 2:
                     constraints = getConstraints(g,constrMod,hautBas,gaucheDroite,corr,voisins,inTex,flatMesh,verticDiff,horizDiff);
                     break;
                     case 3:
                     constraints = getConstraints(g,constrMod,hautBas,gaucheDroite,corr,voisins,inTex,flatMesh,verticDiff,horizDiff,spmTex);
                     break;
                     case 4:
                     constraints =
                     getConstraints(g,constrMod,hautBas,gaucheDroite,corr,voisins,inTex,flatMesh,gyrusMesh,verticDiff,horizDiff,spmTex);
                     break;
                     case 5:
                     constraints = getConstraints(g, constrMod, hautBas, gaucheDroite, corr, voisins, inTex);

                     break;
                  }

                  if (constraints.first.size()!=0 || option == 5)
                     verticDiff = diffusion(poidsGyrus, gyrusMesh[0], haut, bas, constraints.first, 50, vertices[0], corr, criter, dt);
                  else
                  printf("Pas de contrainte !\n");

                  if (constraints.second.size()!=0 || option == 5){
                     horizDiff = diffusion(poidsGyrus, gyrusMesh[0], gauche,droite, constraints.second, 50, vertices[0], corr, criter, dt);
                  }

          /*
                  vector<uint> line(lineExtraction(constr[0],getNearestPoint(droite,verticDiff,40.0),flatMesh));

                  for (uint r=0;r<verticDiff.nItem();r++)
                     verticDiff.item(r)=0;
                  for (uint r=0;r<line.size();r++)
                     verticDiff.item(line[r])=r;
                  writeTexture(verticDiff, "/opt/goperto/testTex.tex");
                  for (uint i2=0;i2<vertices[0].size();i2++)
                     final[1].item(vertices[0][i2]) = (float) verticDiff.item(corr[vertices[0][i2]]);
                  writeTexture(final[1], "/opt/goperto/auxTex.tex");

                  constraints.clear();
                  constraints.push_back(* new pair<vector<uint>,short>(line, 40));

                  verticDiff = diffusion(poidsGyrus, gyrusMesh[0], haut, bas, constraints, vertices[0], corr);
                  //horizDiff = diffusion(poidsGyrus, gyrusMesh[0], gauche, droite, *new vector<pair<vector<uint>,short> >,vertices[0],  corr);*/

                  for (uint i2=0;i2<vertices[0].size();i2++)
                     final[1].item(vertices[0][i2]) = (float) verticDiff.item(corr[vertices[0][i2]]);
                  for (uint i2=0;i2<vertices[0].size();i2++)
                     final[2].item(vertices[0][i2]) = (float) horizDiff.item(corr[vertices[0][i2]]);
               }
               else {
                  printf("Attention : détection intersections \"gaucheDroite\" foireuse.\n");
               }
            }
            else {
               printf("Attention : VerifGyriExistence failed\n");
            }

            printf("\n");
         }
   }



   return final;
}

TimeTexture<float> GyriParamTexture(string meshfile, string intexfile, string gyrifile, string diffmodfile, string spmtexfile, uint constraint_method, float criter, float dt){



  //
  // read triangulation
  //
  cout << "reading triangulation   : " << flush;
  AimsSurfaceTriangle inMesh;
  Reader<AimsSurfaceTriangle> triR( meshfile );
  triR >> inMesh;
  cout << "done" << endl;

  //
  // read input texture
  //
  cout << "reading texture   : " << flush;
  TimeTexture<short> inTex;
  Reader<TimeTexture<short> > texR( intexfile );
  texR >> inTex;
  cout << "done" << endl;

  //
  // read gyri to texture table
  //

  vector<uint> corres(getGyriToTextureCorres(gyrifile.data()));

  //
  // read diffusion model
  //

  pair<vector<vector<uint> >,vector<vector<uint> > > diffMod(getDiffusionModel(diffmodfile.data()));

  //
  // read spmtexfile
  //

 /*
  cout << "reading texture   : " << flush;
  Reader<TimeTexture<float> > spmRdr(spmtexfile);
*/
  TimeTexture<float> spmTex;
/*  spmRdr.read(spmTex);
  cout << "done" << endl;
*/


  cout << "mesh vertices : " << inMesh[0].vertex().size() << endl;
  cout << "mesh polygons : " << inMesh[0].polygon().size() << endl;
  cout << "texture dim   : " << inTex[0].nItem() << endl;

   return GyriParamTexture(inMesh[0], inTex[2], spmTex[0], corres, diffMod, constraint_method, criter, dt);
}




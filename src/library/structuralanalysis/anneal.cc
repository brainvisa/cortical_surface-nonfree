#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include <cortical_surface/structuralanalysis/anneal.h>

using namespace aims;
using namespace carto;
using namespace std;



Anneal::Anneal(Graph &primal){
  MinimizationSetup(primal);

}

void Anneal::Step(vector<int> &random, long double temp, uint &mod){
  
    long double somme=0.0;
    int old;
    mod=0;
    set<uint>::iterator it;

    cout << "STEP " << flush;

    // Iterating On The Vector Of Sites (In Randomized Order) (long loop)
    for (uint i = 0 ; i < random.size() ; i++){
        somme=0.0;

        old = sites[random[i]]->label;

        vector<int> zoneLab;
//         cout << "[" << i << "] " << flush;

        // Selecting The Few Labels That The Current Site Can Take (short loop)
        for ( it = listeZones[ random[i] ].begin() ; it != listeZones[ random[i] ].end() ; it++)
            zoneLab.push_back( *it ); 
        
        vector<long double> globalenergieslabels( zoneLab.size() ),
                            expenergies( zoneLab.size() ),
                            total( zoneLab.size() );

        // Iterating On These Labels
        for ( uint k = 0 ; k < zoneLab.size() ; k++ ) {

//             cout << zoneLab[k] << " " << flush;
            sites[ random[i] ]->label = zoneLab[ k ];
            globalenergieslabels[ k ] = energy;

            int nclsim1 = 0, nclsim2 = 0, nbips1 = 0, nbips2 = 0;

            // Iterating On The Cliques The Current Site Is Involved In
            for ( uint n = 0 ; n < cliquesDuSite[ random[i] ].size() ; n++ ) {
                uint aux = cliquesDuSite[ random[i] ][ n ];

                // Updating The System Depending On The Clique Type
                if ( cliques[ aux ].type == DATADRIVEN || cliques[ aux ].type == INTRAPRIMALSKETCH ){

                    globalenergieslabels[ k ] += cliques[ aux ].updateEnergy( random[ i ], old, false, nbsujets );
                    
                }
                else if ( cliques[aux].type == SIMILARITY ){
                    
                    globalenergieslabels[ k ] += cliques[aux].updateEnergy( random[i], old, false, nbsujets );
                    uint index = 0;
                    if ( cliques[aux].blobs[0]->index == (uint) random[i] )
                        index = 1;
                    else if ( cliques[aux].blobs[1]->index == (uint) random[i] )
                        index = 0;
                    else {
                        ASSERT(false);
                    }
                    if ( cliques[aux].blobs[index]->label == zoneLab[k] && zoneLab[k] != 0 ) nclsim1++;
                    if ( cliques[aux].blobs[index]->label == old && old != 0 ) nclsim2++;
                    
                }
            }
            // Iterating On The Intra Primal Sketch Cliques
            for ( uint n = 0 ; n < ipscliques.size() ; n++ ){
                
                uint aux = ipscliques[n];
                if ( cliques[aux].blobs[0]->subject != sites[ random[i] ]->subject ){
                    if ( zoneLab[k] != 0 )
                        nbips1 += cliques[aux].labelscount[ zoneLab[k] ];
                    if (old != 0)
                        nbips2 += cliques[aux].labelscount[ old ];
                }
                
            }



            total[k] = (nbips1-nclsim1 - (nbips2-nclsim2));


        //       globalenergieslabels[k] += Clique::intrapsweight * total[k];

            somme += exp(-globalenergieslabels[k]/temp);
            if (isnan(somme)) cout << "#####################################" << endl;
        }

        long double somme2=0.0;
        uint acc;
        if (somme > exp(700.0)) {
        acc=0;
        for (uint k=0;k<zoneLab.size();k++)
            if (globalenergieslabels[k]<globalenergieslabels[acc]) acc = k;

        }
        else {
        if (somme > exp(700.0)) cout << "dist:[";
        for (uint k=0;k<zoneLab.size();k++){
            somme2 += exp(-globalenergieslabels[k]/temp)/somme;
            expenergies[k] = somme2;
            if (somme > exp(700.0)){
            cout << somme2 << " " ;
            }
        }
        if (somme > exp(700.0)) cout << "]";

        long double tirage = ((long double)UniformRandom() * 1.0);

        for (acc=0;acc<expenergies.size() && expenergies[acc]<tirage;acc++){} //cout << globalenergieslabels[acc] << " " ;
        }
        if (old != zoneLab[acc]) {
        sites[random[i]]->label=zoneLab[acc];
        for (uint m=0;m<cliquesDuSite[random[i]].size();m++){

            energy +=cliques[cliquesDuSite[random[i]][m]].updateEnergy(random[i],old,true,nbsujets);
        }


    //       energy += Clique::intrapsweight*total[acc];
        mod++;
        }
        else {
        sites[random[i]]->label = old;
        }

    }


}

void Anneal::Run(int verbose){


//   for (uint k=0;k<cliques.size();k++){
//     cliques[k].updateLabelsCount();
//     cliques[k].computeEnergy(true,nbsujets);
//   }
// 
//   for (uint i=0;i<cliques.size();i++)
//     if (cliques[i].type == INTRAPRIMALSKETCH)
//       ipscliques.push_back(i);
// 
//   cout << ipscliques.size() << " cliques intraps" << endl;
// 
//   energy = getTotalEnergy();
// 
//   cout << "energie initiale : " << energy << endl;
//   SummaryLabels();


  vector<int> indices_start;
  for(uint i=0;i<sites.size();i++)
    indices_start.push_back(i);

  long double temp=300.0;

  uint mod=1, ite=0, nb_under_threshold=0,test=1;

  cout.precision(2);

  FILE * f1, *f;
  if (recuitpath!="")
    f1 = fopen (recuitpath.data(),"w");
  if (energypath!=""){
      f = fopen (energypath.data(),"a");
      fprintf(f, "== DEBUT NOUVEAU RECUIT ==\n");
  }
  //   cin >> test;
  //   test=0;

  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }

  if (run==1){
      while (nb_under_threshold<5 || mod!=0){ 
      //    while (temp>200.0){

        if (mod!=0) nb_under_threshold=0;
        else nb_under_threshold++;
        cout << " T=" << temp << " it="<< ite++ << " " << flush ;

        if (recuitpath!=""){
          for (uint i0=0;i0<sites.size();i0++){
            fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
          }
//         double sim0 = getTypeEnergy(SIMILARITY), intraps0=getTypeEnergy(INTRAPRIMALSKETCH), lower0=getTypeEnergy(BESTLOWERSCALE),dd0=getTypeEnergy(DATADRIVEN);
        
//         ASSERT(sim0+dd0+intraps0+lower0==energy);
//         fprintf(f, "%3lf %3lf %3lf %3lf %3lf\n", (float)dd0, (float)sim0, (float) intraps0, (float) lower0, (float)energy);

          fprintf(f1, "\n");
        }
        vector<int> indices(indices_start);
        vector<int> random;
        cout << "test" << flush;
        for (uint i=0;i<sites.size();i++){
          int index = (int)(UniformRandom() * indices.size());
          random.push_back(indices[index]);
          indices.erase(indices.begin()+index);
        }
        cout << "__ok__" << flush;
        ASSERT(random.size() == sites.size());
        Step(random, temp, mod);

        cout << " chg:" << mod << " " << flush;


        if (verbose == 1) ShortSummaryLabels();
//         double everif= getTotalEnergy();
        cout << " E=" << energy << endl; //" Everif=" << everif << endl;
//         if (pow(everif-energy,2)<0.01) cout << "ok"<<endl; else printf("no %.3lf\n",(double)(everif-energy));

        temp = temp*0.99;

    }
  }

  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }
  energy = getTotalEnergy();

  cout << "energie finale : " << energy << endl;
  ShortSummaryLabels();
  if (energypath!=""){
    fprintf(f, "%3lf\nFIN RECUIT\n", (float)energy);
    fclose(f);
  }
  if (recuitpath!="")
    fclose(f1);
  //   for (uint i=0;i<sites.size();i++)
  //     sites[i]->label = (uint)sites[i]->t;
}

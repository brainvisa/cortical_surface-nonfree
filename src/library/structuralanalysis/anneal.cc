#include <aims/math/random.h>
#include <cortical_surface/structuralanalysis/anneal.h>

using namespace aims;
using namespace std;

void Anneal::Step ( std::vector<int> &random, double temp, uint &mod ) {

    double somme = 0.0;
    int old;
    mod = 0;
    std::set<uint>::iterator it;
    double sum_update = 0.0;
    
    // Iterating On The Vector Of Sites (In Randomized Order) (long loop)
    for ( uint i = 0 ; i < random.size() ; i++ ) {
        somme = 0.0;

        old = sites[random[i]]->label;

        std::vector<int> zoneLab;
//         cout << "[" << i << "] " << flush;

        // Selecting The Few Labels That The Current Site Can Take (short loop)
        for ( it = listeZones[ random[i] ].begin() ; it != listeZones[ random[i] ].end() ; it++)
            zoneLab.push_back( *it );

        std::vector<double> globalenergieslabels( zoneLab.size() ),
                            expenergies( zoneLab.size() ),
                            total( zoneLab.size() );

        // Iterating On These Labels
        for ( uint k = 0 ; k < zoneLab.size() ; k++ ) {

            // cout << zoneLab[k] << " " << flush;
            sites[ random[i] ]->label = zoneLab[ k ];
            globalenergieslabels[ k ] = energy;

            int nclsim1 = 0, nclsim2 = 0;

            // Iterating On The Cliques The Current Site Is Involved In
            for ( uint n = 0 ; n < cliquesDuSite[ random[i] ].size() ; n++ ) {
                uint aux = cliquesDuSite[ random[i] ][ n ];

                // Updating The System Depending On The Clique Type
                if ( cliques[ aux ].type == DATADRIVEN ||
                        cliques[ aux ].type == INTRAPRIMALSKETCH ||
                        cliques[ aux ].type == BESTLOWERSCALE ) {

                    globalenergieslabels[ k ] += cliques[ aux ].updateEnergy( random[ i ], old, false, nbsujets );

                }
                else if ( cliques[aux].type == SIMILARITY ) {

                    globalenergieslabels[ k ] += cliques[aux].updateEnergy( random[i], old, false, nbsujets );
                    uint index = 0;
                    if ( cliques[aux].blobs[0]->index == (uint) random[i] )
                        index = 1;
                    else if ( cliques[aux].blobs[1]->index == (uint) random[i] )
                        index = 0;
                    else 
                        ASSERT(false);
                    if ( cliques[aux].blobs[index]->label == zoneLab[k] && zoneLab[k] != 0 ) nclsim1++;
                    if ( cliques[aux].blobs[index]->label == old && old != 0 ) nclsim2++;

                }
            }
            // Iterating On The Intra Primal Sketch Cliques
            uint nbips1 = 0, nbips2 = 0;
            for ( uint n = 0 ; n < ipscliques.size() ; n++ ) {

                uint aux = ipscliques[n];
                if ( cliques[aux].blobs[0]->subject != sites[ random[i] ]->subject ) {
                    if ( zoneLab[k] != 0 )
                        nbips1 += cliques[aux].labelscount[ zoneLab[k] ];
                    if (old != 0)
                        nbips2 += cliques[aux].labelscount[ old ];
                }
            }

            total[k] = ( nbips1 - nclsim1 - (nbips2 - nclsim2) );

            somme += exp ( - globalenergieslabels[k] / temp );
            assert ( !isnan(somme) || ( std::cout << globalenergieslabels[k] << "/" << temp << std::endl && false ) );
        }

        double somme2 = 0.0;
        uint acc;
        if ( somme > exp(700.0) ) {
            acc = 0;
            for ( uint k = 0 ; k < zoneLab.size() ; k++ )
				if ( globalenergieslabels[k] < globalenergieslabels[acc] ) acc = k;
		}
        else {
            for ( uint k = 0 ; k < zoneLab.size() ; k++ ) {
                somme2 += exp ( - globalenergieslabels[k] / temp ) / somme;
                expenergies[k] = somme2;
            }

			double tirage = (double) UniformRandom() * 1.0 ;

			for ( acc = 0 ;
				  acc < expenergies.size() && expenergies[acc] < tirage ;
				  acc++ )    { }  //cout << globalenergieslabels[acc] << " " ;
        }
        if ( old != zoneLab[acc] ) {

            sites[random[i]]->label = zoneLab[acc];
            for ( uint m = 0 ; m < cliquesDuSite[random[i]].size() ; m++ ) {
                if ( cliques [cliquesDuSite[random[i]][m]].type == DATADRIVEN ||
                    cliques [cliquesDuSite[random[i]][m]].type == SIMILARITY ||
                    cliques [cliquesDuSite[random[i]][m]].type == BESTLOWERSCALE ||
                    cliques [cliquesDuSite[random[i]][m]].type == INTRAPRIMALSKETCH ) {
                        
                        double update = cliques [cliquesDuSite[random[i]][m]].updateEnergy ( random[i], old, true, nbsujets );
                        sum_update += update;
                        energy += update;
                }
            }
            mod++;
        }
        else {
            sites[random[i]]->label = old;
        }
    }
}

void Anneal::Run ( int verbose ){

    std::vector< int >  indices_start;
    for( uint i = 0 ; i < sites.size() ; i++ )
        indices_start.push_back(i);

    double temp = 300.0;

    uint mod = 1, ite = 0, nb_under_threshold = 0, test = 1;

    std::cout.precision(2);

    FILE * f1, *f;
    if ( labelsPath != "" )
        f1 = fopen ( labelsPath.data(), "w" );
    if ( energyPath != "" ) {
        f = fopen ( energyPath.data(), "a" );
        fprintf(f, "== DEBUT NOUVEAU RECUIT ==\n");
    }

    for ( uint k = 0 ; k < cliques.size() ; k++ ){
        cliques[k].updateLabelsCount();
        cliques[k].computeEnergy(true, nbsujets);
    }

    while ( nb_under_threshold < 5 || mod != 0 ) {

        if ( mod != 0 )
            nb_under_threshold = 0;
        else
            nb_under_threshold++;

        std::cout << " T=" << temp << " it="<< ite++ << " " << std::flush ;

        if ( this->labelsPath != "" ) {
            for ( uint i0 = 0 ; i0 < sites.size() ; i0++ ) {
                fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
            }
            fprintf(f1, "\n");
        }
        std::vector< int > indices( indices_start );
        std::vector< int > random;

        for ( uint i = 0 ; i < sites.size() ; i++ ) {
            int index = (int)(UniformRandom() * indices.size());
            random.push_back(indices[index]);
            indices.erase(indices.begin()+index);
        }
        ASSERT ( random.size() == sites.size() );
        Step(random, temp, mod);

        std::cout << " chg:" << mod << " " << std::flush;

        if (verbose == 1) ShortSummaryLabels();
        double everif = getTotalEnergy();
        std::cout << " E=" << energy << std::endl << " Everif=" << everif << endl;

        temp = temp * 0.99;

    }

    for ( uint k = 0 ; k < cliques.size() ; k++ ) {
        cliques[k].updateLabelsCount();
        cliques[k].computeEnergy(true, nbsujets);
    }
    energy = getTotalEnergy();

    std::cout << "final energy : " << energy << std::endl;
    ShortSummaryLabels();
    if ( this->energyPath != "" ) {
        fprintf(f, "%3lf\nFIN RECUIT\n", (float)energy);
        fclose(f);
    }
    if ( this->labelsPath != "" )
        fclose(f1);

}

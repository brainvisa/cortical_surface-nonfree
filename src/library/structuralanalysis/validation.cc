#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include <cortical_surface/structuralanalysis/validation.h>

using namespace aims;
using namespace carto;
using namespace std;


void StructuralAnalysis_Validation::printFile ( std::vector<double> &samples, FILE *f ) {
    for ( uint j = 0 ; j < samples.size() ; j++ ) {
        // cout << (float)(j)*(float)step+(float)mini << " " << histo[j] << endl;
        // if (samples[j] > 10000.0) cout << "attention outlier " ;
        fprintf(f, "%lf\n", samples[j]);
    }
    std::cout << std::endl;
}

std::vector<int> StructuralAnalysis_Validation::getCompConn ( SurfaceBased_StructuralAnalysis &ssb, 
                                                              std::vector<uint> &indicesCliques, 
                                                              std::set<uint> &listeSites ) {
    std::vector<int> comp ( ssb.sites.size() );
    uint blob0, blob1;
    Site *s0, *s1;
    int lcomp, nbcomp, aux;
    int label0,label1;
    std::set<uint>::iterator it;
    
    for ( uint i = 0 ; i < ssb.sites.size() ; i++ )
        comp[i] = -1;
    
    lcomp = 0;
    nbcomp = 0;

    for ( it = listeSites.begin() ; it != listeSites.end() ; it++ )
        comp[*it] = 0;

    for ( uint i = 0 ; i < indicesCliques.size() ; i++ )
        if ( ssb.cliques[indicesCliques[i]].type == SIMILARITY ) {
            s0 = ssb.cliques[indicesCliques[i]].blobs[0];
            s1 = ssb.cliques[indicesCliques[i]].blobs[1];
            blob0 = s0->index;
            blob1 = s1->index;
            comp[blob0] = 0;
            comp[blob1] = 0;
        }
    
    for ( uint i = 0 ; i < indicesCliques.size() ; i++ ) {
        if ( ssb.cliques[indicesCliques[i]].type == SIMILARITY ) {

            s0 = ssb.cliques[indicesCliques[i]].blobs[0];
            s1 = ssb.cliques[indicesCliques[i]].blobs[1];
            blob0 = s0->index;
            blob1 = s1->index;

            if ( comp[blob0] == 0 ) {
                if ( comp[blob1] == 0 ) {
                    lcomp++;
                    nbcomp++;
                    comp[blob0] = lcomp;
                    comp[blob1] = lcomp;
  //             cout << "a"<< nbcomp << "-" ;
                }
                else if ( comp[blob1] > 0 ) {
                    comp[blob0] = comp[blob1];
                }
            }
            else if ( comp[blob1] == 0 ) {
                comp[blob1] = comp[blob0];
            }
            else if ( comp[blob0] == comp[blob1] ) {
                // clique cyclique ne pas s'inquiéter mais répéter 100 fois rapidement
            }
            else if ( comp[blob0] > 0 && comp[blob1] > 0 ) { // fusion
                label0 = comp[blob0];
                label1 = comp[blob1];
                if ( label0 > label1 ) {
                    aux = label0;
                    label0 = label1;
                    label1 = aux;
                }
                //           cout << nbcomp << ";" << lcomp << "(" << label1 << "-" << label0 << ")-";
                for ( uint i0 = 0 ; i0 < comp.size() ; i0++ ) {
                    if ( comp[i0] == label1 )
                        comp[i0] = label0;
                }
                for ( uint i0 = 0 ; i0 < comp.size() ; i0++ ){
                    if ( comp[i0] == lcomp ) {
                        comp[i0] = label1;
                        //               cout << "o" ;
                    }
                }
                lcomp--;
                nbcomp--;
            }
        }
    }

    std::vector<uint> zeros;
    for ( uint i = 0 ; i < comp.size() ; i++ )
        if ( comp[i] == 0 ) 
            zeros.push_back(i);
    for ( uint i = 0 ; i < zeros.size() ; i++ ) {
        lcomp++;
        nbcomp++;
        comp[zeros[i]] = lcomp;
    }
    //   cout << "nbcomp:" << nbcomp << endl;
    return comp;
}

std::vector<std::set<uint> > StructuralAnalysis_Validation::getCompConnVector ( std::vector<int> &comp){
    std::vector< std::set<uint> > cc;
    uint cpt = 0, cpt2, nbsites = 0;
    for ( uint i = 0 ; i < comp.size() ; i++ )
        if ( comp[i] >= 0 )  
            nbsites++;
    //   cout << nbsites << endl;
    for ( uint i0 = 0 ; cpt < nbsites ; i0++ ) {
        cpt2 = 0;
//     cout << i0 << " " << flush;
        cc.push_back( std::set<uint>() );
        for ( uint i = 0 ; i < comp.size() ; i++ )
            if ( comp[i] == (int) i0 ) {
                cpt2++; 
                cpt++;
                cc[i0].insert(i);
            }
    }
    //   cout << endl;
    return cc;
}

uint StructuralAnalysis_Validation::nbcombinaisons( SurfaceBased_StructuralAnalysis &ssb,
                                                    std::set<uint> &graphe, uint card){
    uint nbfinal = 0;
    std::set<uint>::iterator it2;
    std::set<unsigned short int>::iterator it;
    std::set<std::set<unsigned short int> >::iterator itt;
    std::set<std::set<unsigned short int> > composantes;
    std::set<std::set<unsigned short int> > composantes_np1;
    uint taille = 2;

    for ( it2 = graphe.begin() ; it2 != graphe.end() ; it2++ ) {
        // cout << "\b\b\b\b\b\b\b\b\b\b\b\b\bit="<< *it << flush;
        std::set<unsigned short int> composante;
        composante.insert(*it2);
        composantes_np1.insert(composante);
    }
    std::cout << "nb taille 1 :" << composantes_np1.size() << std::endl;

    while (taille <=card && composantes_np1.size() < 1000000) {
        composantes.clear();
        composantes = set<set<unsigned short int> >(composantes_np1);
        composantes_np1.clear();
        //     cout << "composantes:"<<composantes.size() << endl;
        uint aux=0;
        for (itt = composantes.begin();itt!=composantes.end() && composantes_np1.size() < 1000000;itt++){
            if (aux++%1000==0) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << aux << "/" << composantes.size()<<"("<< composantes_np1.size()<<")" << flush;
            set<unsigned short int> voisins;
            //       set<unsigned short int> currcomp(*itt);
            for (it = (*itt).begin();it!= (*itt).end(); it++){
                for (uint i=0;i<ssb.cliquesDuSite[*it].size() ;i++){
                    uint currclique = ssb.cliquesDuSite[*it][i];
                    if (ssb.cliques[currclique].type==SIMILARITY){
                        //         cout << "last="<<last << " i="<<i << " "<< flush; ;
                        if ((*itt).find(ssb.cliques[currclique].blobs[0]->index) == (*itt).end() && graphe.find(ssb.cliques[currclique].blobs[0]->index)!=graphe.end())
                            voisins.insert(ssb.cliques[currclique].blobs[0]->index);
                        if ((*itt).find(ssb.cliques[currclique].blobs[1]->index) == (*itt).end() && graphe.find(ssb.cliques[currclique].blobs[1]->index)!=graphe.end())
                            voisins.insert(ssb.cliques[currclique].blobs[1]->index);
                    }
                }
            }
            //       cout << "voisins.size() = " << voisins.size() << endl;
            for (it = voisins.begin();it !=voisins.end();it++){
                set<unsigned short int> comp_np1((*itt));
                comp_np1.insert(*it);
                composantes_np1.insert(comp_np1);
            }
        }
        std::cout << "nb taille="<<taille << " :" << composantes_np1.size() << std::endl;
        if ( card < 5 ) {
            for ( itt =composantes_np1.begin() ; itt!=composantes_np1.end() ; itt++ ) {
                for ( it = (*itt).begin() ; it!= (*itt).end() ; it++ )
                    cout << *it << "-" << flush;
                cout << " " << flush;
            }
            cout << endl;
        }
        taille++;
    }
    return composantes_np1.size();
}

double StructuralAnalysis_Validation::WalshTest( std::vector<double> &samplesdist, int r){

  std::sort(samplesdist.begin(), samplesdist.end());
  double c =(double) ceil(sqrt(2.F*samplesdist.size())); uint k=r+c; double b2 = 1.0/0.05;
  double a = (1.0 + sqrt(b2) * sqrt((c-b2)/(c-1)))/(c-b2-1.0);
//   Point2df res(samplesdist[0]-(1+a)*samplesdist[1]+a*samplesdist[k-1],samplesdist[samplesdist.size()-1]-(1+a)*samplesdist[samplesdist.size()-2]+a*samplesdist[samplesdist.size()-k]);
// Xr - (1+a)Xr+1 + aXk < 0
  return samplesdist[r-1] - (1+a)*samplesdist[r] + a*samplesdist[k-1];
  return samplesdist[samplesdist.size()-r]-(1+a)*samplesdist[samplesdist.size()-r-1]+a*samplesdist[samplesdist.size()-k];
//   return res;
}

std::vector<double> StructuralAnalysis_Validation::getCaracSample ( SurfaceBased_StructuralAnalysis &ssb,
                                                                    std::vector<uint> &composante ) {
    double tmoy = 0.0, 
           rec = 0.0, 
           sum = 0.0, 
           compac = 0.0, 
           Ttest = 0.0, 
           compaccent;
    uint nbblobsrec = 0;
    Point3df barycentre ( 0.0, 0.0, 0.0 );
    for ( uint k = 0 ; k < composante.size() ; k++ ) {
        tmoy += ssb.sites[composante[k]]->tValue;
        compac += ssb.sites[composante[k]]->t;
    }
    compac /= composante.size();
    for ( uint k = 0 ; k < composante.size() ; k++ )
        sum += pow ( ssb.sites[composante[k]]->tValue - compac, 2 );
    Ttest = sqrt ( (float) ( composante.size() * ( composante.size() - 1 ) ) ) * compac / sqrt(sum);
    
    for ( uint k = 0 ; k < composante.size() - 1 ; k++ ) {
        for ( uint m = k + 1 ; m < composante.size() ; m++ ) {
            Point3df bbmax1 = ssb.sites[composante[k]]->boundingbox_max, 
                     bbmax2 = ssb.sites[composante[m]]->boundingbox_max;
            Point3df bbmin1 = ssb.sites[composante[k]]->boundingbox_min, 
                     bbmin2 = ssb.sites[composante[m]]->boundingbox_min;
            uint no_overlap = 1;
            float reco = getOverlap ( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
            if ( no_overlap == 0 )
                rec += reco;
            nbblobsrec++;
        }
        barycentre += ssb.sites[composante[k]]->gravitycenter;
    }
    barycentre /= composante.size();
    
    double distance_moy = 0.0;
    
    for ( uint k = 0 ; k < composante.size() - 1 ; k++ ) { 
        Point3df aux = ssb.sites[composante[k]]->gravitycenter - barycentre;
        distance_moy += aux.norm();
    }
    distance_moy /= composante.size();
    
    std::vector<double> sample;
    sample.push_back ( compac );
    sample.push_back ( Ttest );
    sample.push_back ( rec / (double) nbblobsrec );
    sample.push_back ( ssb.getClusterEnergy(composante) );
    sample.push_back ( distance_moy );
    return sample;
}

std::vector<double> StructuralAnalysis_Validation::getBackup ( SurfaceBased_StructuralAnalysis &ssb, 
                                                               std::vector<uint> &composante ) {
    double tmoy = 0.0, 
           rec = 0.0, 
           sum = 0.0, 
           compac = 0.0, 
           Ttest = 0.0;
    uint nbblobsrec = 0;
    std::vector<double> sample;

    for ( uint k = 0 ; k < composante.size() ; k++ ) {
        tmoy += ssb.sites[composante[k]]->t;
        sample.push_back ( ssb.sites[composante[k]]->tValue );
        compac += ssb.sites[composante[k]]->tValue;
    }
    compac /= composante.size();
    for ( uint k = 0 ; k < composante.size() ; k++ )
        sum += pow ( ssb.sites[composante[k]]->tValue - compac, 2 );
    Ttest = sqrt ( (float) ( composante.size() * ( composante.size() - 1 ) ) ) * compac / sqrt(sum);
    for ( uint k = 0 ; k < composante.size() - 1 ; k++ )
        for ( uint m = k + 1 ; m < composante.size() ; m++ ) {
            Point3df bbmax1 = ssb.sites[composante[k]]->boundingbox_max, 
                     bbmax2 = ssb.sites[composante[m]]->boundingbox_max;
            Point3df bbmin1 = ssb.sites[composante[k]]->boundingbox_min, 
                     bbmin2 = ssb.sites[composante[m]]->boundingbox_min;

            uint no_overlap = 1;
            float reco = getOverlap ( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
            if ( no_overlap == 0 ) {
                rec += reco;
                sample.push_back ( reco );
                nbblobsrec++;
            }
        }
    sample.push_back(compac);
    sample.push_back(sqrt(sum));
    return sample;
}

double StructuralAnalysis_Validation::getRankPercentile ( std::vector< double > &samplesT,
                           std::set<uint>  &activblobs, 
                           SurfaceBased_StructuralAnalysis &ssb ) {
    std::vector<uint> composante;
    std::vector<double> samplesTri ( samplesT );
    std::set <uint>::iterator it;

    for ( it = activblobs.begin() ; it != activblobs.end() ; it++ )
        composante.push_back( *it );
    std::vector<double> sample ( getCaracSample( ssb, composante) );
    samplesTri.push_back( sample[0] );
    std::sort( samplesTri.begin(), samplesTri.end() );
    double r;
    for ( r = 0 ;
          r < (int) samplesTri.size() && pow( samplesTri[r] - sample[0], 2) > 0.001 ;
          r++ ) { }
    assert ( r != samplesTri.size() );
    std::cout << "perc. " << ":" << r / (double) samplesTri.size() * 100.0 << " " << sample[0] << std::endl;
    return r / (double) samplesTri.size() * 100.0;
}

std::vector<uint> StructuralAnalysis_Validation::getComponent ( SurfaceBased_StructuralAnalysis &ssb,
                                 std::set<uint> &activblobs, 
                                 std::set< uint > &forbidden,
                                 std::vector< int > &tirage ) {
    // cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << j << flush;
    // création du cluster
    uint startsite;
    std::vector< uint > composante;
    std::set <uint>::iterator it;

    std::set<uint> dejapris = std::set< uint > ( forbidden );
    do {
        startsite = (float) UniformRandom() * ssb.sites.size();
        }
    while ( tirage[startsite] < 1 ||
            dejapris.find( startsite ) != dejapris.end() );

    composante.push_back( startsite );
    
    dejapris.insert( startsite );
    
    std::set<uint> voisins;
    while ( composante.size() < activblobs.size() ) {
    
        uint test = composante[ composante.size() - 1 ];
        for ( uint m = 0 ; m < ssb.cliquesDuSite[test].size() ; m++ ) {
            uint currclique = ssb.cliquesDuSite[test][m];
            if ( ssb.cliques[currclique].type == SIMILARITY ) {
                if ( dejapris.find(ssb.cliques[currclique].blobs[0]->index) == dejapris.end() )
                    voisins.insert( ssb.cliques[currclique].blobs[0]->index );
                if ( dejapris.find(ssb.cliques[currclique].blobs[1]->index) == dejapris.end() )
                    voisins.insert( ssb.cliques[currclique].blobs[1]->index );
            }
        }
        assert( voisins.size() != 0 );
    
        uint random = (float) UniformRandom() * voisins.size();
        it = voisins.begin();
        for ( uint n = 0 ; n < random ; n++, it++ ) { }
        if ( dejapris.find(*it) == dejapris.end() ) {
          composante.push_back(*it);
          dejapris.insert(*it);
          voisins.erase(it);
        }
    }
    return composante;
}

void StructuralAnalysis_Validation::saveSignificanceInfo ( uint label, 
                                                           SurfaceBased_StructuralAnalysis &swc, 
                                                           std::vector<std::vector<uint> > & composantes, 
                                                           FILE *f ) {
    

    for ( uint i = 0 ; i < composantes.size() ; i ++) {
        std::vector<double> sample( getCaracSample(swc, composantes[i]) );
        //samplesCarac.push_back( sample );
        fprintf(f, "%d, %.3f, %.3f, %.3f, %.3f\n", label, sample[3], sample[0], sample[2], sample[3] );
        printf("%d, %.3f, %.3f, %.3f, %.3f\n", label, sample[3], sample[0], sample[2], sample[3] );
        //samplesT.push_back( sample[0] );
        //samplesRec.push_back( sample[2] );
        //if ( sample[3] < 0.0 )
        //    samplesNeg.push_back( sample[0] );
        //samplesDispersion.push_back( sample[4] );
    }
}

std::vector<std::vector<uint> > StructuralAnalysis_Validation::ValidAround ( uint label, 
                                                                             SurfaceBased_StructuralAnalysis &ssb,
                                                                             uint number_of_samples ) {
    std::vector<std::vector<uint> > composantes;
    std::set <uint>::iterator it;
  
    std::set<uint> activblobs;  // activblobs contains the different sites indices carrying a positive label 
    for ( uint j = 0 ; j < ssb.sites.size() ; j++ )
        if ( ssb.sites[j]->label == (uint) label )
            activblobs.insert( j );

    std::vector <uint> cliquesV;
    std::set <uint> sitesV;

    for ( uint j = 0 ; j < ssb.cliques.size() ; j++ )
        cliquesV.push_back( j );
    for ( uint j = 0 ; j < ssb.sites.size() ; j++ )
        sitesV.insert( j );      

    std::vector< int > ccc = getCompConn( ssb, cliquesV, sitesV );
    std::vector< std::set<uint> > cc = getCompConnVector( ccc );
    uint startsite;

    std::set< uint > activblobsglobal;
    for ( it = activblobs.begin() ; it != activblobs.end() ; it++ )
         activblobsglobal.insert( *it );      


    std::cout << "ABG=" << activblobsglobal.size() << std::endl;
    
    std::set< uint > forbidden, autorized ( activblobs );

    std::cout << std::endl << activblobs.size() << endl;
    for ( it = activblobs.begin() ; it != activblobs.end() ; it++ )
        std::cout << *it << " ";
    std::cout << endl;

    // on va tirer au sort des clusters de taille activblobs.size()
    std::vector< int > select_cc, tirage( ssb.sites.size() );
    for ( uint k = 0 ; k < cc.size() ; k++ )
        if ( cc[k].size() >= activblobs.size() )
            select_cc.push_back( k );
    for ( uint k = 0 ; k < ssb.sites.size() ; k++ ) {
        tirage[k] = -1;
        if ( cc [ccc[k]].size() >= activblobs.size() )
            tirage[k] = ccc[k];
    }
    
 
    for ( it = activblobs.begin() ; it != activblobs.end() ; it++) {
        uint test = *it;
        for ( uint m = 0 ; m < ssb.cliquesDuSite[test].size() ; m++ ) {
            uint currclique = ssb.cliquesDuSite[test][m];
            if ( ssb.cliques[currclique].type == SIMILARITY ) {
                autorized.insert ( ssb.cliques[currclique].blobs[0]->index );
                autorized.insert ( ssb.cliques[currclique].blobs[1]->index );
            }
        }
    }
    std::cout << "restent :"<< ssb.sites.size() - forbidden.size() << std::endl;
    //            std::vector< std::vector<double> > samplesCarac;
    for ( uint j = 0 ; j < number_of_samples ; j++ )
        composantes.push_back( getComponent ( ssb, activblobs, forbidden, tirage ) );
    
    return composantes;

}

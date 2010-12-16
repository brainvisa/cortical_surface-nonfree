#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include <cortical_surface/structuralanalysis/minimization.h>
#include <cortical_surface/structuralanalysis/texturetoblobs.h>
#include <float.h>


using namespace aims;
using namespace carto;
using namespace std;

float max ( float a, float b ) {
  if ( a > b )
      return a;
  else
      return b;
  return a;
}
float min ( float a, float b ) {
  if ( a < b )
      return a;
  else
      return b;
  return a;
}

struct ltstr_vec
{
    bool operator ( ) ( const std::vector<uint> p1,
                      const std::vector<uint> p2) {
        assert( p1.size() == p2.size() );
        for ( uint i = 0 ; i < p1.size() ; i++ ) {
            if ( p1[p1.size()-1-i] < p2[p1.size()-1-i] )
                return true;
            else if ( p1[p1.size()-1-i]>p2[p1.size()-1-i] )
                return false;
        }
    }
};

//###############################################################################################


void SurfaceBased_StructuralAnalysis::noLabelsZones ( int number_of_labels ) {/* vector<pair<Point2df,Point2df> > &labelsZones,
                                                           vector<set<uint> > &zonesListesBlobs,
                                                           vector<set<uint> > &listeZones ){*/
    listeZones = std::vector<std::set<uint> > (sites.size());
//     zonesListesBlobs = vector<set<uint> >(labelsZones.size()+1);
//     labelsZones = vector<pair<Point2df,Point2df> > ();
//     for (uint it=0;it<20;it++){
//         pair<Point2df, Point2df> zone;
//                 zone.first = Point2df(0.0, 0.0);
//                 zone.second = Point2df(180.0,360.0);
//         labelsZones.push_back(zone);
//     }

    for ( uint i = 0 ; i < number_of_labels ; i++ )
        labels.push_back(i);

    for ( uint j = 0 ; j < sites.size() ; j++ ) {
//         listeZones[j].insert(0);
//         zonesListesBlobs[0].insert(j);
        for ( uint k = 0 ; k < number_of_labels ; k++ ) {
            listeZones[j].insert(k);
//             zonesListesBlobs[k+1].insert(j);
        }
    }
}

//###############################################################################################


void SurfaceBased_StructuralAnalysis::regionLabelsZones () {/* vector<pair<Point2df,Point2df> > &labelsZones,
                                                           vector<set<uint> > &zonesListesBlobs,
                                                           vector<set<uint> > &listeZones ){*/
    uint i = 0;
    labelsZones = std::vector<std::pair<Point2df,Point2df> > ();

    for ( uint it = 0 ; it < 3 ; it++ )
        for ( float zonelat = -20 ; zonelat < 180.0 ; zonelat += 36.0 )
            for ( float zonelon = -20 ; zonelon < 360.0 ; zonelon += 72.0 ) {

                std::pair<Point2df, Point2df> zone;
                zone.first = Point2df( max(zonelat, 0.0), max(zonelon, 0.0) );
                zone.second = Point2df( min(zonelat + 76.0, 180.0), min(zonelon + 112.0,360.0) );
                // std::cout << i << " " << zone.first[0] << ";" << zone.first[1] << " " << zone.second[0] << ";" << zone.second[1] << std::endl;
                labelsZones.push_back(zone);
                i++;
            }

    std::cout << labelsZones.size() << " zones" << std::endl;
    for ( i = 0 ; i < labelsZones.size() + 1 ; i++ )
        labels.push_back(i);
    std::vector<int> zonescount, labelscount;

    for ( i = 0 ; i < labelsZones.size() + 1 ; i++ ) {
        zonescount.push_back(0);
        labelscount.push_back(0);
    }

    zonesListesBlobs = std::vector< std::set<uint> >( labelsZones.size() + 1 );
    listeZones = std::vector<std::set<uint> > (sites.size());

    for ( uint j = 0 ; j < sites.size() ; j++ ) {
        uint count = 0;
        listeZones[j].insert(0);
        zonesListesBlobs[0].insert(j);
        for ( uint k = 0 ; k < labelsZones.size() ; k++ ) {
            Point3df bbmin1 = sites[j]->boundingbox_min,
                     bbmax1 = sites[j]->boundingbox_max;
            uint no_overlap = 0;
            getOverlap ( bbmin1,
                         bbmax1,
                         Point3df ( labelsZones[k].first[0], labelsZones[k].first[1], 0.0 ),
                         Point3df ( labelsZones[k].second[0],labelsZones[k].second[1], 0.0 ),
                         &no_overlap ) ;
            if ( no_overlap == 0 ) {
                zonescount[k+1]++;
                zonesListesBlobs[k+1].insert(j);
                listeZones[j].insert(k+1);
                count++;
            }
        }
        labelscount[count]++;
        assert ( listeZones[j].size() != 0 );
    }

    for ( i = 0 ; i < labelsZones.size() + 1 ; i++ )
        std::cout << zonescount[i] << " ";
    std::cout << std::endl;
    for ( i = 0 ; i < labelsZones.size() + 1 ; i++ )
        std::cout << labelscount[i] << " ";
    std::cout << std::endl;

}

//###############################################################################################


// Initialization Gets All Labels Set to 0 (if initLabel Set to True), Counts Labels
//   Occurences, Initial Energy,

void SurfaceBased_StructuralAnalysis::Initialization ( bool initLabel ) {

    std::cout << "Initialization..." << std::endl;
    if ( initLabel ) {
      std::cout << " All labels set to " << labels[0] << std::endl;

      for ( uint i = 0 ; i < sites.size() ; i++ )
        sites[i]->label = labels[0];
    }
    else {
      std::cout << "No label initialization" << std::endl;
      for ( uint i = 0 ; i < sites.size() ; i++ )
          std::cout << sites[i]->label << " " << std::flush;
    }

    for ( uint k = 0 ; k < cliques.size() ; k++ ){
        cliques[k].updateLabelsCount();
        cliques[k].updateSubjectsCount();
        cliques[k].computeEnergy( true, nbsujets );
    }
    energy = getTotalEnergy();
    std::cout << "initial energy : " << energy << std::endl;
    SummaryLabels();
}

double SurfaceBased_StructuralAnalysis::getClusterEnergy ( std::vector<uint> &composante ) {
    std::vector<int> old;
    for ( uint i = 0 ; i < sites.size() ; i++ ) {
        old.push_back ( sites[i]->label );
        sites[i]->label = 0;
    }
    for ( uint i = 0 ; i < composante.size() ; i++ )
        sites[composante[i]]->label = 1;

    double energy = getTotalEnergy();
    for ( uint i =  0; i < sites.size() ; i++ )
        sites[i]->label = old[i];
    return energy;
}

double SurfaceBased_StructuralAnalysis::getLabelEnergy(int label, int type){
    bool test = true;
    double energydd = 0.0, energyls = 0.0, energysim = 0.0, energy = 0.0;
    uint nclsim = 0, nbips = 0;
    std::vector<int> bysub( nbsujets );

    for ( uint i = 0 ; i < cliques.size() ; i++ ) {
        cliques[i].updateLabelsCount();
        cliques[i].updateSubjectsCount();
        cliques[i].computeEnergy ( true, nbsujets );
        if ( type == UNKNOWN || cliques[i].type == type ) {
            test = true;
            for ( uint j = 0 ; j < cliques[i].blobs.size() && test == true ; j++ )
                if ( cliques[i].blobs[j]->label != label )
                    test = false;
            if (test) {
                if ( cliques[i].type == DATADRIVEN ) {
                    energydd += cliques[i].energie;
                }
                else if ( cliques[i].type == BESTLOWERSCALE ) {
                    energyls += cliques[i].energie;
                }
                else if ( cliques[i].type == SIMILARITY )  {
                    nclsim++;
                    energysim += cliques[i].energie;
                }

    //         if (false && cliques[i].type==SIMILARITY) std::cout << "(("<<cliques[i].rec << "(" << cliques[i].blobs[0]->label << "(" <<  cliques[i].blobs[0]->node << ")" << "-" << cliques[i].blobs[1]->label <<"(" <<  cliques[i].blobs[1]->node << ")" << ")=>" << cliques[i].energie << ")) ";
            }
    //       if (cliques[i].type == INTRAPRIMALSKETCH && cliques[i].labelscount[label] > 1 && label != 0) {energy += (cliques[i].labelscount[label]-1)*Clique::getIntraPSWeight();}
        }
    }
    uint i, j;
    for ( uint n = 0 ; n < ipscliques.size() ; n ++ ) {
          bysub[n] = cliques[ipscliques[n]].labelscount[label];
    //       std::cout << bysub[n] << " " ;
    }
    //     std::cout << std::endl;
    uint nb = 0;
    for ( i = 0 ; i < nbsujets - 1 ; i++ ) {
        for ( j = i + 1 ; j < nbsujets ; j++ ) {
            nb += bysub[i]*bysub[j];
        }
    }
    uint nb2 = 0;
    for ( i = 0 ; i < nbsujets ; i++ )
        if ( bysub[i] > 1 )
            nb2 += bysub[i] - 1;

    //     std::cout << nclsim << "|"<< energydd << " " << energysim << "|";
    nbips += nb;
    // std::cout << Clique::intrapsweight*(nbips-nclsim) << " ";
    ASSERT( nbips>=nclsim || (std::cout << nbips << ">=" << nclsim << std::endl && false) );
    //   energy += Clique::intrapsweight*(nbips-nclsim);
    energy += Clique::intrapsweight * nb2 * nbsujets;
    energy += energydd + energyls;
    energy += energysim;

    return energy;
}

double SurfaceBased_StructuralAnalysis::getTypeEnergy(int type){ // RETOURNE L'ENERGIE PAR TYPE DE CLIQUE
    double energy = 0.0;
    for ( uint i = 0 ; i < cliques.size() ; i++ )
        if ( cliques[i].type == type )
            energy += cliques[i].energie;
    return energy;
}

double SurfaceBased_StructuralAnalysis::getTotalEnergy ( ) {

    double energy = 0.0;
    int nclsim = 0, nbips = 0, nb2 = 0;
    std::vector<int> bysub( nbsujets );

    for ( uint i = 0 ; i < cliques.size() ; i++ ) {
        cliques[i].updateLabelsCount();
        cliques[i].updateSubjectsCount();
        if ( cliques[i].type == DATADRIVEN || cliques[i].type == BESTLOWERSCALE ) {
            energy += cliques[i].computeEnergy( true, nbsujets );
        }
        else if ( cliques[i].type == SIMILARITY ) {
            energy += cliques[i].computeEnergy( true, nbsujets );
            if ( cliques[i].blobs[0]->label == cliques[i].blobs[1]->label && cliques[i].blobs[0]->label != 0 )
                nclsim++;
        }
    }

    for ( uint k = 1 ; k < labels.size() ; k++ ) {

        for ( uint n = 0 ; n < ipscliques.size() ; n++ )
            bysub[n] = cliques[ipscliques[n]].labelscount[labels[k]];
            // std::cout << bysub[n] << " " ;
        // std::cout << std::endl;
        uint nb = 0;
        for ( uint i = 0 ; i < nbsujets - 1 ; i++ )
            for ( uint j = i + 1 ; j < nbsujets ; j++ )
                nb += bysub[i] * bysub[j];
        //  std::cout << "nb"<< nb << " ";
        nbips += nb;
        for ( uint i = 0 ; i < nbsujets ; i++ )
            if ( bysub[i] > 1 )
                nb2 += bysub[i] - 1;
    }

    ASSERT( nbips >= nclsim || (std::cout << nbips << ">=" << nclsim << std::endl && false ));
    //std::cout << Clique::intrapsweight << " " << nb2 << " " << nbsujets << std::endl;
    energy += Clique::intrapsweight * nb2 * nbsujets;

    return energy;
}



void SurfaceBased_StructuralAnalysis::SummaryLabels ( ) {

    double Edd, Els, Esim, energy, Eintra;
    std::cout << labels[0] << ":";
    uint nblabel = 0;
    for ( uint il = 0 ; il < cliques.size() ; il++ ) {
        nblabel += cliques[il].labelscount[0];
    }
    std::cout << nblabel << " - ";
    std::vector<uint> nblab;

    double Etotal = 0.0;
    std::cout << std::endl << std::endl;
//   FILE * f;
//   if (energyPath!=""){
//    f = fopen (energyPath.data(),"a");
//    fprintf(f, "== SUMMARYLABELS ==\n");
//   }
    for ( uint lab = 1 ; lab < labels.size() ; lab++ ) {
        Edd = 0.0; Els = 0.0; Esim = 0.0; Eintra = 0.0;

        energy = 0.0;
        int nclsim = 0, nbips = 0;
        std::vector<int> bysub ( nbsujets );

        for ( uint i = 0 ; i < cliques.size() ; i++ ) {
            cliques[i].updateLabelsCount();
            cliques[i].updateSubjectsCount();

            if ( cliques[i].type == DATADRIVEN && cliques[i].blobs[0]->label == labels[lab] ) {
                energy += cliques[i].computeEnergy ( true, nbsujets );
                Edd += cliques[i].computeEnergy ( true, nbsujets );
            }
            else if ( cliques[i].type == BESTLOWERSCALE && cliques[i].blobs[0]->label == labels[lab] ) {
                energy += cliques[i].computeEnergy ( true, nbsujets );
                Els += cliques[i].computeEnergy ( true, nbsujets );
            }
            else if ( cliques[i].type == SIMILARITY
                    && cliques[i].blobs[0]->label == cliques[i].blobs[1]->label
                    && cliques[i].blobs[0]->label== labels[lab] ) {
                energy += cliques[i].computeEnergy ( true, nbsujets );
                Esim += cliques[i].computeEnergy ( true, nbsujets );
                if ( cliques[i].blobs[0]->label == cliques[i].blobs[1]->label
                        && cliques[i].blobs[0]->label != 0 )
                    nclsim++;
            }
        }
        for ( uint n = 0 ; n < ipscliques.size() ; n++ ) {
            bysub[n] = cliques[ipscliques[n]].labelscount[labels[lab]];
        }

        uint nb = 0;
        for ( uint k = 0 ; k < nbsujets - 1 ; k++ ) {
            for ( uint j = k + 1 ; j < nbsujets ; j++ ) {
                nb += bysub[k] * bysub[j];
            }
        }
        nbips += nb;

        uint nb2 = 0.0;
        for ( uint i = 0 ; i < nbsujets ; i++ )
            if ( bysub[i] > 1 )
                nb2 += bysub[i] - 1;

        ASSERT( nbips >= nclsim || (std::cout << nbips << ">=" << nclsim << std::endl && false ));
        energy += Clique::intrapsweight * nb2 * nbsujets;
        Eintra += Clique::intrapsweight * nb2 * nbsujets;

        Etotal += energy;

        double diff = Edd + Els + Esim + Eintra - energy;
        ASSERT ( fabs(diff) < 0.000000001);

        nblabel = 0;
        for ( uint il = 0 ; il < ipscliques.size() ; il++ ) {
            nblabel += cliques[ipscliques[il]].labelscount[labels[lab]];
        }
        if ( nblabel != 0 ) {
            std::cout << "label " << labels[lab] << " : " << cliques[globalclique].subjectscount[ labels[lab] ].size() << \
            " subjects (" << Edd << ";" << Esim << ";" << Els << ";" <<Eintra<<") " << energy << " " << std::flush;

            std::cout << nblabel << " - " << std::endl;
            // fprintf(f,"label %d : (%3lf;%3lf) %i;\n", labels[lab],(double)Edd,(double)Esim,nblabel);
        }
        nblab.push_back(nblabel);

    }
    std::cout << " ";
    for ( uint il = 0 ; il < nblab.size() ; il++ )
        if ( nblab[il] != 0 )
            std::cout << "<<" << nblab[il] << ">>-";
        else
            std::cout << nblab[il] << "-" ;
    for ( uint i = 0 ; i < sites.size() ; i++ ) {
        for ( uint j = 0 ; j < nblab.size() ; j++ ) {
            if ( sites[i]->label != 0 )
                sites[i]->label_occur_number = nblab[sites[i]->label-1];
            else
                sites[i]->label_occur_number = 0;
        }
    }
    std::cout <<"\b ";
    std::cout << "Etot=" << Etotal << std::endl;
//   if (energyPath!="")
//     fclose(f);

}

void SurfaceBased_StructuralAnalysis::ShortSummaryLabels(){
    // Various energy values : "energy" is first computed for each label, "Etotal"
    //   being the sum of them.
    double Edd, Els, Esim, energy, Eintra, Etotal = 0.0;

    // First counting the zero label occurrences
    std::cout << labels[0] << ":";
    uint nblabel = 0;
    for ( uint il = 0 ; il < cliques.size() ; il++ )
        nblabel += cliques[il].labelscount[0];
    std::cout << nblabel << " - ";

    // Now dealing with the other labels
    std::vector<uint> nblab;
    for ( uint lab = 1 ; lab < labels.size() ; lab++ ) {
        Els = 0.0; Edd = 0.0; Esim = 0.0; Eintra = 0.0;
        energy = 0.0;
        int nclsim = 0;

        for ( uint i = 0 ; i < cliques.size() ; i++ ) {
            cliques[i].updateLabelsCount();
            cliques[i].updateSubjectsCount();
            if ( ( cliques[i].type == DATADRIVEN || cliques[i].type == BESTLOWERSCALE )
                    && cliques[i].blobs[0]->label == labels[lab] ) {
                energy += cliques[i].computeEnergy( true, nbsujets );
                if ( cliques[i].type == DATADRIVEN )
                    Edd += cliques[i].computeEnergy( true, nbsujets );
                else if ( cliques[i].type == BESTLOWERSCALE )
                    Els += cliques[i].computeEnergy( true, nbsujets );
            }
            else if ( cliques[i].type == SIMILARITY &&
                    cliques[i].blobs[0]->label == cliques[i].blobs[1]->label &&
                    cliques[i].blobs[0]->label == labels[lab] ) {
                energy += cliques[i].computeEnergy ( true, nbsujets );
                Esim += cliques[i].computeEnergy ( true, nbsujets );
                //std::cout << cliques[i].blobs[0]->index << "-"<< cliques[i].blobs[1]->index << "-" << cliques[i].blobs[0]->subject \
                    << "-" << cliques[i].blobs[1]->subject << std::endl;
                if (cliques[i].blobs[0]->label == cliques[i].blobs[1]->label && cliques[i].blobs[0]->label != 0)
                    nclsim++;
            }
        }
        assert( ipscliques.size() == nbsujets ) ;
        // Now counting the number of subjects in which the current label appears
        std::vector<int> bysub ( nbsujets );
        for ( uint i = 0 ; i < nbsujets ; i ++ ) {
            bysub[i] = cliques[ipscliques[i]].labelscount[ labels[lab] ];
            //std::cout << bysub[i] << "-" << std::flush;
        }
        //std::cout << std::endl;
        // nbips is how many intersubject pairs can be made with the sites
        //   carrying the current label
        uint nbips = 0;
        for ( uint i = 0 ; i < nbsujets - 1 ; i++ )
            for ( uint j = i + 1 ; j < nbsujets ; j++ )
                  nbips += bysub[i] * bysub[j];
        //std::cout << nbips << std::endl;

        uint nb = 0;
        for ( uint i = 0 ; i < nbsujets ; i++ )
            if ( bysub[i] > 1 ) nb += bysub[i] - 1;

        energy += Clique::intrapsweight * nb * nbsujets;
        Eintra += Clique::intrapsweight * nb * nbsujets;
        Etotal += energy;

        ASSERT( pow ( Edd + Esim + Eintra  + Els - energy, 2 ) < DBL_EPSILON );

        nblabel = 0;
        for ( uint il = 0 ; il < ipscliques.size() ; il++ ) {
            nblabel += cliques[ ipscliques[il] ].labelscount[ labels[lab] ];
        }
        if ( nblabel != 0 ) {
            std::cout << "L" << labels[lab] << ":" << cliques[globalclique].subjectscount[labels[lab]].size() << \
            " subjects (" << energy << "=" << Edd << "+" << Esim << "+" << Eintra << "+" << Els <<"):" << std::flush;
            std::cout << nblabel << " - ";
        }
        nblab.push_back ( nblabel );
        ASSERT( nbips >= nclsim || (std::cout << nbips << ">=" << nclsim << std::endl && false ) );

    }
    std::cout << " ";
    for ( uint il = 0 ; il < nblab.size() ; il++ )
        std::cout << nblab[il] << "-" ;
    std::cout << "\b ";
    std::cout << "Etot=" << Etotal << " "  ;

}



void SurfaceBased_StructuralAnalysis::setModelParameters ( float _ddweight,
                                                           float _intrapsweight,
                                                           float _simweight,
                                                           float _lsweight,
                                                           float _ddx1,
                                                           float _ddx2,
                                                           float _simx1,
                                                           float _simx2,
                                                           float _lsx1,
                                                           float _lsx2,
                                                           float _ddh ) {
    Clique::setParameters ( _ddweight, _intrapsweight, _simweight, _lsweight, _ddx1, _ddx2, _simx1, _simx2, _lsx1, _lsx2, _ddh );
}

void SurfaceBased_StructuralAnalysis::StoreToGraph ( Graph &primal ) {
    std::set<Vertex *>::iterator iv, jv;

    int index;
    std::string subject;
    primal.setProperty ( "global_energy", (float) energy );

    for ( iv = primal.vertices().begin() ; iv != primal.vertices().end() ; ++iv ) {
        std::string test;
        (*iv)->getProperty( "index", index );
        (*iv)->getProperty( "subject", subject );

        for ( uint i = 0 ; i < sites.size() ; i++ ) {

            if ( sites[i]->graph_index == (uint) index && sites[i]->subject == subject ) {
                std::ostringstream s;
                s << sites[i]->label ;
                (*iv)->setProperty( "label", s.str());
                (*iv)->setProperty( "name", s.str());
            }
        }
    }
}

void SurfaceBased_StructuralAnalysis::StoreSignificanceToGraph ( Graph &primal ) {
    std::set<Vertex *>::iterator iv, jv;
    int index, label_occur_number;
    float value;
    std::string subject;

    for ( iv = primal.vertices().begin() ; iv != primal.vertices().end() ; ++iv ) {
        std::string test;
        (*iv)->getProperty( "index", index );
        (*iv)->getProperty( "subject", subject );

        for ( uint i = 0 ; i < sites.size() ; i++ ) {

            if ( sites[i]->graph_index == (uint) index && sites[i]->subject == subject ) {
                // node = sites[i]->node;
                label_occur_number = sites[i]->label_occur_number;
                (*iv)->setProperty( "label_occur_number", label_occur_number);
                value = sites[i]->significance;
                (*iv)->setProperty( "significance", value);
                value = sites[i]->t_rankperc;
                (*iv)->setProperty( "t_rankperc", value);
                value = sites[i]->sim_rankperc;
                (*iv)->setProperty( "sim_rankperc", value);
            }
        }
    }
}

//void SurfaceBased_StructuralAnalysis::BuildSites ( Graph &primal, std::vector<Site *> sites ) {
//
//    std::set<Vertex *>::iterator iv;
//    std::vector<float> bc;
//    float tmin = -1.0,
//        tmax = -1.0,
//        trep = -1.0,
//        tvalue = -1.0,
//        t = -1.0;
//        //x = -1.0,
//        //y = -1.0;
//    int index = -1,
//        rank = -1;
//    std::string subject = "";
//    int newindex = 0;
//
//    std::vector< int > listepts;
//    std::vector< float > bbmax, bbmin;
//
//    for ( iv = primal.vertices().begin() ; iv != primal.vertices().end() ; ++iv ) {
//        std::string test;
//        (*iv)->getProperty( "index", index);
//        (*iv)->getProperty( "subject", subject);
//        (*iv)->getProperty( "label", test);
//        (*iv)->getProperty( "gravity_center", bc);
//        (*iv)->getProperty( "tmin", tmin);
//        (*iv)->getProperty( "tmax", tmax);
//        (*iv)->getProperty( "trep", trep);
//        (*iv)->getProperty( "t", t);
//        (*iv)->getProperty( "tValue", tvalue);
//        (*iv)->getProperty( "rank", rank);
//        (*iv)->getProperty( "nodes_list", listepts);
//        (*iv)->getProperty( "boundingbox_max", bbmax);
//        (*iv)->getProperty( "boundingbox_min", bbmin);
//
//        sites.push_back( new Site() );
//        Site *s = sites[sites.size()-1];
//
//        s->label = atoi(test.data());
//        (*iv)->setProperty("label", test);
//        (*iv)->setProperty( "name", test);
//        (*iv)->getProperty( "label", test);
//        (*iv)->getProperty( "name", test);
//        s->index = newindex++;
//        s->graph_index = index;
//        s->subject = subject;
//
//        if ( bc.size() == 3 ) {
//            s->gravitycenter = Point3df();
//            s->gravitycenter[0] = bc[0];
//            s->gravitycenter[1] = bc[1];
//            s->gravitycenter[2] = bc[2];
//        }
//        s->tmin = tmin;
//        s->tmax = tmax;
//        s->trep = trep;
//        s->rank = rank;
//        if ( bbmax.size() == 3 ) {
//            s->boundingbox_max[0] = bbmax[0];
//            s->boundingbox_max[1] = bbmax[1];
//            s->boundingbox_max[2] = bbmax[2];
//        }
//        if ( bbmin.size() == 3 ) {
//            s->boundingbox_min[0] = bbmin[0];
//            s->boundingbox_min[1] = bbmin[1];
//            s->boundingbox_min[2] = bbmin[2];
//        }
//
//        s->t = t;
//        s->tValue = tvalue;
//        for ( uint i = 0 ; i < listepts.size() ; i++ )
//            s->nodes_list.insert( listepts[i] );
//
//        (*iv)->setProperty( "sites_index", sites.size() - 1 );
//    }
//    //return sites;
//}

void SurfaceBased_StructuralAnalysis::ConvertSSBlobsToSites ( std::vector<surf::ScaleSpaceBlob *> &ssblobs, std::vector<Site *> &sites ) {

    for ( uint i = 0 ; i < ssblobs.size() ; i++ ) {

        sites.push_back(new Site());
        Site *s = sites[sites.size() - 1];        
        s->index = sites.size() - 1;
        s->graph_index = ssblobs[i]->index;
        s->subject = ssblobs[i]->subject;
        s->label = ssblobs[i]->label;
        s->tmin = ssblobs[i]->tmin;
        s->tmax = ssblobs[i]->tmax;
        s->t = ssblobs[i]->t;
        s->nodes_list = ssblobs[i]->nodes;

        if ( ssblobs[i]->nodes.size() == 1 ) {
            // if the ssblobs are provided with one point in vector raw_coordinates,
            //  it means that the graph was provided with a gravity_center for every ssb node
            int node = *(ssblobs[i]->nodes.begin());
            s->gravitycenter = Point3df ( ssblobs[i]->raw_coordinates[node][0],
                                          ssblobs[i]->raw_coordinates[node][1],
                                          ssblobs[i]->raw_coordinates[node][2] ) ;
        }

        std::pair<Point2df, Point2df> bb = ssblobs[i]->get2DBoundingBox();
        s->boundingbox_min = Point3df(bb.first[0], bb.first[1], 0);
        s->boundingbox_max = Point3df(bb.second[0], bb.second[1], 0);

    }

}

void SurfaceBased_StructuralAnalysis::GetSimilarityCliquesFromSSBCliques ( std::vector<surf::Clique> &ssbcliques,
                                std::vector<Site *> &sites,
                                std::vector<Clique> &cliques,
                                std::vector<std::vector<int> > &cliquesDuSite){

    cliquesDuSite = std::vector<std::vector<int> >( sites.size() );

    for ( uint i = 0 ; i < ssbcliques.size() ; i++ ) {
        Clique simc;
        simc.type = SIMILARITY;

        simc.similarity = ssbcliques[i].similarity;
        simc.distance = ssbcliques[i].distance;
        assert (simc.similarity == -1.0 || simc.distance == -1.0);

        surf::ScaleSpaceBlob *ssb1, *ssb2;
        int iSSB1, iSSB2;
        ssb1 = ssbcliques[i].ssb1;
        ssb2 = ssbcliques[i].ssb2;
        
        for ( iSSB1 = 0 ; iSSB1 < sites.size() && 
            !(sites[iSSB1]->graph_index == ssb1->index && sites[iSSB1]->subject == ssb1->subject ) ;
             iSSB1++ ) { }
        for ( iSSB2 = 0 ; iSSB2 < sites.size() && 
            !(sites[iSSB2]->graph_index == ssb2->index && sites[iSSB2]->subject == ssb2->subject ) ;
             iSSB2++ ) { }
        
        //iSSB1 = ssb1->index;
        //iSSB2 = ssb2->index;
        cliquesDuSite[ sites[iSSB1]->index ].push_back(i);
        cliquesDuSite[ sites[iSSB2]->index ].push_back(i);

        simc.blobs.push_back( sites[iSSB1] );
        simc.blobs.push_back( sites[iSSB2] );
        cliques.push_back(simc);
    }

    for ( uint i = 0 ; i < sites.size() ; i++ ) {
        for ( uint n = 0 ; n < cliquesDuSite[i].size() ; n++ ) {
            uint aux = cliquesDuSite[ i ][ n ];
            if ( cliques[aux].type == SIMILARITY ) {
                if ( cliques[aux].blobs[0]->index == (uint) i ) { }
                else if (cliques[ aux ].blobs[1]->index == (uint) i ) { }
                else {
                    std::cout << i << " " << aux << " " << cliques[aux].type << " " << cliques[aux].blobs.size() << " " << cliques[aux].blobs[0]->index << " " << cliques[aux].blobs[1]->index << std::endl;
                    ASSERT(false);
                }
                assert ( cliques[aux].blobs[0]->subject != cliques[aux].blobs[1]->subject );
            }
        }
    }
    std::cout << cliques.size() << "cliques recovered from ssbcliques" << std::endl;

}

//##############################################################################

// This function takes the "ssblobs" vector and figures out which pairs of blobs
//  overlap. The resulting vector "cliques" associates to every relevant pair of
//  scale-space blobs (noted by their indices) its calculated spatial overlap.

std::vector<surf::Clique> SurfaceBased_StructuralAnalysis::BuildSimilarityCliques ( std::vector<surf::ScaleSpaceBlob *>   &ssblobs,
                                                 std::vector<std::vector<surf::GreyLevelBlob *> > &matchingblobs ) {

    std::vector<surf::Clique > cliques;
    matchingblobs = std::vector<std::vector<surf::GreyLevelBlob *> > (ssblobs.size());

    std::set<surf::GreyLevelBlob *>::iterator itB1, itB2;
    surf::GreyLevelBlob *b1max, *b2max;

    // Start of cliques construction
    for ( uint i = 0 ; i < ssblobs.size() - 1 ; i++ ) {
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << ssblobs.size() << "(" << cliques.size() << ")" << std::flush;
        for ( uint j = i + 1 ; j < ssblobs.size() ; j++ ) {

            // For every single pair of scale-space blobs, computes a maximal overlap
            //   between every possible pair of grey-level blobs.

            // We consider only pairs of scale-space blobs from different subjects.
            if ( ssblobs[i]->subject != ssblobs[j]->subject ) {

                float overmax = -1.0;
                for ( itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() ; itB1++ ) {
                    for ( itB2 = ssblobs[j]->blobs.begin() ; itB2 != ssblobs[j]->blobs.end() ; itB2++ ) {

                        // For every possible pair of grey-level blobs between these two scale-
                        //   space blobs, we figure out their possible spatial overlap.

                        // vector<int> listNodesB1(set2vector((*itB1)->nodes_list)),
                        //             listNodesB2(set2vector((*itB2)->nodes_list));
                        std::pair<Point2df, Point2df> bbi = (*itB1)->get2DBoundingBox(),
                            //getBoundingBox((*itB1)->nodes, data[ssblobs[i]->subject].lat, data[ssblobs[i]->subject].lon),
                                                bbj = (*itB2)->get2DBoundingBox();
                            //(*itB2)->nodes, data[ssblobs[j]->subject].lat, data[ssblobs[j]->subject].lon);
                        Point2df bbmin1 (bbi.first[0], bbi.first[1]),
                                bbmax1 (bbi.second[0], bbi.second[1]),
                                bbmin2 (bbj.first[0], bbj.first[1]),
                                bbmax2 (bbj.second[0], bbj.second[1]) ;
                        uint no_overlap = 2;
                        double overlap = TextureToBlobs::getOverlapMeasure( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
                        //cout << "bbi("<< (*itB1)->nodes.size() << "):" << bbi.first[0] << "-" << bbi.first[1] << " " <<
                        //bbi.second[0] << "-" << bbi.second[1] << " " <<
                        //"bbj("<< (*itB2)->nodes.size() << "):" << bbj.first[0] << " " << bbj.first[1] << " " <<
                        //bbj.second[0] << " " << bbj.second[1]  << endl;

                        if ( no_overlap == 0 ){

                            // If the current pair's overlap is maximal, then the glb indices are stored.

                            //cout << "bbi("<< (*itB1)->nodes.size() << "):" << bbi.first[0] << "-" << bbi.first[1] << " " <<
                            //   bbi.second[0] << "-" << bbi.second[1] << " " <<
                            //   "bbj("<< (*itB2)->nodes.size() << "):" << bbj.first[0] << " " << bbj.first[1] << " " <<
                            //   bbj.second[0] << " " << bbj.second[1] << " over:" << overlap << endl;
                            //cout << (*itB1)->scale << " " << (*itB2)->scale << endl;

                            if ( overlap > overmax ) {
                                overmax = overlap;
                                b1max = *itB1;
                                b2max = *itB2;
                            }
                        }

                    }
                }

                // Here all the possible glb pairs have been processed for the two current ssb
                if ( overmax > 0.0 &&
                        !((ssblobs[j]->tmin > ssblobs[i]->tmax) || (ssblobs[i]->tmin > ssblobs[j]->tmax)) ) {
                    // If the two scale-space blobs have at least one pair of grey-level
                    //   overlapping (bounding-boxes) (+ scales overlapping), then a clique
                    // is created between these two ssb and the max-overlapping pair of glb
                    // is stored in "matchingblobs".

                    cliques.push_back(surf::Clique(ssblobs[i], ssblobs[j], -1.0, overmax));
                    matchingblobs[i].push_back(b1max);
                    matchingblobs[j].push_back(b2max);
                    //cout << "max (" << ssblobs[i]->index <<","<< ssblobs[j]->index << ") between:" << b1max->index << " "
                    //<< b2max->index << " overmax:" << overmax << endl;
                    //cout << "scales: " << b1max->scale << " " << b1max->scale << endl;
                }
            }
            // The next pair of scale-space blobs will now be processed.
        }
    }
    std::cout << ssblobs.size() << "/" << ssblobs.size() << "(" << cliques.size() << ")" << std::endl;
    // Construction of a representation blob for each scale-space blob
    for ( uint i = 0 ; i < ssblobs.size() ; i++ ) {

        // For every scale-space blob, we create a representation blob
        //   from the set of grey-level blobs found to be max-matching
        //   with some others (from other scale-space blobs)
        std::set<uint>::iterator it;

        if ( matchingblobs[i].size() != 0 )
            std::cout << i << ":";

        // for (it = matchingblobs[i].begin() ; it != matchingblobs[i].end() ; it++){
        for ( uint j = 0 ; j < matchingblobs[i].size() ; j++ ) {
            std::set<int> blobNodes( matchingblobs[i][j]->nodes );
            ssblobs[i]->nodes.insert( blobNodes.begin(), blobNodes.end() );
            std::cout << ssblobs[i]->nodes.size() << " " << std::flush;
        }

        if (matchingblobs[i].size()!=0)
            std::cout << std::endl ;
    }
    return cliques;
}

////##############################################################################


std::vector<surf::Clique> SurfaceBased_StructuralAnalysis::BuildSimilarityCliques3D ( std::vector<surf::ScaleSpaceBlob *>   &ssblobs,
                                                   GroupData &data,
                                                   float threshold,
                                                   int type_distance ) {

    std::vector<surf::Clique > cliques;

    std::set<surf::GreyLevelBlob *>::iterator itB1, itB2;
    surf::GreyLevelBlob *b1max, *b2max;

    // Start of cliques construction
    if ( type_distance == DISTANCE_LATITUDES )
        std::cout << "DISTANCE_LATITUDES" << std::endl;
    else if ( type_distance == DISTANCE_LONGITUDES )
        std::cout << "DISTANCE_LONGITUDES" << std::endl;
    else if ( type_distance == DISTANCE_LATLON )
            std::cout << "DISTANCE_LATLON" << std::endl;
    else if ( type_distance == DISTANCE_3DEUCLIDIAN )
        std::cout << "DISTANCE_3DEUCLIDIAN" << std::endl;

    for ( uint i = 0 ; i < ssblobs.size() - 1 ; i++ ) {

        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << ssblobs.size() << "(" << cliques.size() << ")" << std::flush;
        for ( uint j = i + 1 ; j < ssblobs.size() ; j++ ) {

            // For every single pair of scale-space blobs, computes a maximal overlap
            //   between every possible pair of grey-level blobs.

            // We consider only pairs of scale-space blobs from different subjects.
            if ( ssblobs[i]->subject != ssblobs[j]->subject ) {

                float distance=-1.0;
                assert( ssblobs[i]->blobs.size() == 1 );
                itB1 = ssblobs[i]->blobs.begin();
                b1max = *itB1;
                itB2 = ssblobs[j]->blobs.begin();
                b2max = *itB2;
                int max1 = b1max->getMaximumNode(*(data[ssblobs[i]->subject]->tex));
                int max2 = b2max->getMaximumNode(*(data[ssblobs[j]->subject]->tex));

                if ( type_distance == DISTANCE_LATITUDES ) {
                    Point3df p1 (b1max->coordinates[max1][0], b1max->coordinates[max1][1], 0.0);
                    Point3df p2 (b2max->coordinates[max2][0], b2max->coordinates[max2][1], 0.0);
                    Point3df p = p1-p2;
                    distance = sqrt(p[0]*p[0]);
                }
                else if ( type_distance == DISTANCE_LONGITUDES ) {
                    Point3df p1 (b1max->coordinates[max1][0], b1max->coordinates[max1][1], 0.0);
                    Point3df p2 (b2max->coordinates[max2][0], b2max->coordinates[max2][1], 0.0);
                    Point3df p = p1-p2;
                    distance = sqrt(p[1]*p[1]);
                }
                else if ( type_distance == DISTANCE_LATLON ) {
                    Point3df p1 (b1max->coordinates[max1][0], b1max->coordinates[max1][1], 0.0);
                    Point3df p2 (b2max->coordinates[max2][0], b2max->coordinates[max2][1], 0.0);
                    Point3df p = p1-p2;
                    distance = sqrt(p[0]*p[0] + p[1]*p[1]);
                }
                else if ( type_distance == DISTANCE_3DEUCLIDIAN ) {
                    Point3df p1 (b1max->raw_coordinates[max1][0], b1max->raw_coordinates[max1][1], b1max->raw_coordinates[max1][2]);
                    Point3df p2 (b2max->raw_coordinates[max2][0], b2max->raw_coordinates[max2][1], b2max->raw_coordinates[max2][2]);
                    Point3df p = p1-p2;
                    distance = p.norm();
                }

                // std::cout << distance << " " << std::flush;

                if ( distance < threshold ) {
                    cliques.push_back( surf::Clique(ssblobs[i], ssblobs[j], distance, -1.0 ) );
                }
            }
            // The next pair of scale-space blobs will now be processed.
        }
    }
    std::cout << std::endl;
    std::cout << ssblobs.size() << "/" << ssblobs.size() << "(" << cliques.size() << ")" << std::endl;
    return cliques;
}

//##############################################################################

















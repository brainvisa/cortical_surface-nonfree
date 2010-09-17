#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/minimization.h>


using namespace aims;
using namespace carto;
using namespace std;

float Clique::ddweight, Clique::intrapsweight, Clique::simweight, Clique::lsweight, Clique::ddx2, Clique::ddx1, Clique::ddh, Clique::globalweight;
void Clique::setParameters( float _ddweight, float _intrapsweight, float _simweight, float _lsweight, float _ddx2, float _ddx1, float _ddh, float _globalweight ){
  ddweight=_ddweight; intrapsweight = _intrapsweight; simweight=_simweight; lsweight=_lsweight; ddx2 =_ddx2;  ddx1 = _ddx1; ddh=_ddh;
  globalweight=_globalweight;
}

void Clique::updateLabelsCount ( ) {
    if ( type == INTRAPRIMALSKETCH ) {
        labelscount = std::map<int, uint>();
        for ( uint i = 0 ; i < blobs.size() ; i++ ) {
            if ( labelscount.find(blobs[i]->label) == labelscount.end() )
                labelscount[blobs[i]->label] = 1;
            else
                labelscount[blobs[i]->label]++;
        }
        std::map<int, uint>::iterator it;
        uint chksum = 0;
    //     cout << blobs[0]->subject << endl;
        for ( it = labelscount.begin() ; it != labelscount.end() ; it++ ) {
    //       cout << it->first << " " << it->second << "-";
            chksum += (*it).second;
        }
    //     cout  << endl;
        ASSERT( chksum == blobs.size() );
    }
}

void Clique::updateSubjectsCount ( ) {
    if ( type == GLOBAL ) {
        subjectscount = std::map<int, std::set<std::string> >();

        for ( uint i = 0 ; i < blobs.size() ; i++ ) {
            if ( subjectscount.find( blobs[i]->label) == subjectscount.end() ) {
                subjectscount[blobs[i]->label] = std::set<std::string>();
            }
            subjectscount[blobs[i]->label].insert(blobs[i]->subject);
        }

        std::map<int, std::set<std::string> >::iterator itset;
        for ( itset = subjectscount.begin() ; itset != subjectscount.end() ; itset++ ) {
            std::cout << (*itset).second.size() << " " << std::flush;
        }
        std::cout << std::endl;

    }
}


double getOverlap(Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap){

  float overlap_x,overlap_y,aux;
  double rec=0.0;

  if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) < 0.0001) {bbmax1[0] += 0.5; /*cout << "bbmax10+ ";*/}
  if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) < 0.0001) {bbmax1[1] += 0.5; /*cout << "bbmax11+ ";*/}
  if (sqrt(pow(bbmin2[0]-bbmax2[0],2)) < 0.0001) {bbmax2[0] += 0.5; /*cout << "bbmax20+ ";*/}
  if (sqrt(pow(bbmin2[1]-bbmax2[1],2)) < 0.0001) {bbmax2[1] += 0.5; /*cout << "bbmax21+ ";*/}
  // if (bbmin1[1]>bbmax1[1] && bbmin2[1] < bbmax2[1] ) {//alors i a bouclé autour de 360/0
  if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) >300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) <300){
    //  cout << "i boucle lon " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
    //  ASSERT(bbmin1[1]>bbmax2[1]);
      if (360-bbmax2[1]<bbmin2[1]){
        aux = bbmax1[1];
        bbmax1[1] = bbmin1[1] + 360.0;
        bbmin1[1] = aux;
      }
      else {
        aux = bbmin1[1];
        bbmin1[1] = bbmax1[1] - 360.0;
        bbmax1[1] = aux;
      }
  }
  //  else if (bbmin1[1]<bbmax1[1] && bbmin2[1] > bbmax2[1] ) {//alors j a bouclé autour de 360/0
  else if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) <300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) >300){
  //  cout << "j boucle lon " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
  //  ASSERT(bbmin2[1]>bbmax1[1]);
      if (360-bbmax1[1]<bbmin1[1]){
        aux = bbmax2[1];
        bbmax2[1] = bbmin2[1] + 360.0;
        bbmin2[1] = aux;
      }
      else {
        aux = bbmin2[1];
        bbmin2[1] = bbmax2[1] - 360.0;
        bbmax2[1] = aux;
      }
  }
  // else if (bbmin1[1]>bbmax1[1] && bbmin2[1]>bbmax2[1] ) {//alors i&j ont bouclé
  else if (sqrt(pow(bbmin1[1]-bbmax1[1],2)) >300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) >300){
  // cout << "i et j bouclent lon " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
        aux = bbmin1[1];
        bbmin1[1] = bbmax1[1] - 360.0;
        bbmax1[1] = aux;
        aux = bbmin2[1];
        bbmin2[1] = bbmax2[1] - 360.0;
        bbmax2[1] = aux;
  }
// on s'occupe de la latitude
// if (bbmin1[0]>bbmax1[0] && bbmin2[0] < bbmax2[0] ) {//alors i a bouclé autour de 360/0
if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) >150 && sqrt(pow(bbmin2[0]-bbmax2[0],2)) <150){
//  cout << "i boucle lat" << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//  ASSERT(bbmin1[0]>bbmax2[0]);
      if (180-bbmax2[0]<bbmin2[0]){
        aux = bbmax1[0];
        bbmax1[0] = bbmin1[0] + 180.0;
        bbmin1[0] = aux;
      }
      else {
        aux = bbmin1[0];
        bbmin1[0] = bbmax1[0] - 180.0;
        bbmax1[0] = aux;
      }
  }
//  else if (bbmin1[0]<bbmax1[0] && bbmin2[0] > bbmax2[0] ) {//alors j a bouclé autour de 360/0
else if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) <150 && sqrt(pow(bbmin2[0]-bbmax2[0],2)) >150){
//  cout << "j boucle lat" << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//  ASSERT(bbmin2[0]>bbmax1[0]);
      if (180-bbmax1[0]<bbmin1[0]){
        aux = bbmax2[0];
        bbmax2[0] = bbmin2[0] + 180.0;
        bbmin2[0] = aux;
      }
      else {
        aux = bbmin2[0];
        bbmin2[0] = bbmax2[0] - 180.0;
        bbmax2[0] = aux;
      }

  }
// else if (bbmin1[0]>bbmax1[0] && bbmin2[0]>bbmax2[0] ) {//alors i&j ont bouclé
else if (sqrt(pow(bbmin1[0]-bbmax1[0],2)) >150 && sqrt(pow(bbmin2[0]-bbmax2[0],2)) >150){
// cout << "i et j bouclent lat" << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
        aux = bbmin1[0];
        bbmin1[0] = bbmax1[0] - 360.0;
        bbmax1[0] = aux;
        aux = bbmin2[0];
        bbmin2[0] = bbmax2[0] - 360.0;
        bbmax2[0] = aux;
  }
  // prétraitements effectués on calcule le recouvrement
  // if (*no_overlap==0) cout << "rec: " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
  float margin=2.0;
  bbmin1[0]-=margin;
  bbmin2[0]-=margin;
  bbmin1[1]-=margin/2.0;
  bbmin2[1]-=margin/2.0;
  bbmax1[0]+=margin;
  bbmax2[0]+=margin;
  bbmax1[1]+=margin/2.0;
  bbmax2[1]+=margin/2.0;
  *no_overlap=0;
  if (bbmin1[0]<=bbmin2[0])
    if (bbmax1[0]<bbmin2[0]) *no_overlap=1;
  else overlap_x= (bbmax2[0] < bbmax1[0] ? bbmax2[0] : bbmax1[0]) - bbmin2[0] ;
  else
    if (bbmax2[0]<bbmin1[0]) *no_overlap=1;
  else overlap_x= (bbmax1[0] < bbmax2[0] ? bbmax1[0] : bbmax2[0]) - bbmin1[0];

  if (*no_overlap==0) {

    if (bbmin1[1]<=bbmin2[1])
      if (bbmax1[1]<bbmin2[1])
          *no_overlap=1;
      else
          overlap_y= (bbmax2[1] < bbmax1[1] ? bbmax2[1] : bbmax1[1]) - bbmin2[1];
    else
      if (bbmax2[1]<bbmin1[1])
          *no_overlap=1;
      else
          overlap_y= (bbmax1[1] < bbmax2[1] ? bbmax1[1] : bbmax2[1]) - bbmin1[1];

    if (*no_overlap==0) {

      rec=overlap_x*overlap_y;
      double div=( ((bbmax1[0]-bbmin1[0])*(bbmax1[1]-bbmin1[1]))
            + ((bbmax2[0]-bbmin2[0])*(bbmax2[1]-bbmin2[1])));

  // if (*no_overlap==0 && rec > div/2.0) {cout << "rec: " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl; cout << rec << " " << div << " " << 2*rec/div << endl;}
  // assert((rec > div/2.0) || !(cout << bbmax1[0] << " " << bbmax1[1] << " " << bbmax1[2] << " "
      rec=2 * rec / div;
  // cout << rec << " " ;
  // assert(rec<1.0 || !(cout << "!" << rec << "!"<< flush));

    }

  }

  return rec;
}


void BuildMaximalOrderCliques ( vector<Site *> &sites,
                               vector<vector<int> > &cliquesDuSite,
                               vector<Clique> &cliques ) {

    vector< Clique > intraps;
    vector< string > subjects;

    for ( uint i = 0 ; i < sites.size() ; i++ ) {
        uint j = 0;

        for ( ; j < subjects.size() && subjects[j] != sites[i]->subject ; j++ ) { }

        if ( j == subjects.size() ) {
            subjects.push_back( sites[i]->subject );
            intraps.push_back( Clique() );
            intraps[j].type = INTRAPRIMALSKETCH;

        }

        intraps[j].blobs.push_back( sites[i] );

    }

    for ( uint i = 0 ; i < intraps.size() ; i++ ) {
        cliques.push_back( intraps[i] );
        for ( uint j = 0 ; j < intraps[i].blobs.size() ; j++ ) {
            cliquesDuSite[ intraps[i].blobs[j]->index ].push_back( cliques.size() - 1 );
        }
    }

}


void BuildDataDrivenCliques( vector<Site *> &sites,
                                  vector<vector<int> > &cliquesDuSite,
                                  vector<Clique> &cliques ) {

  for ( uint i = 0 ; i < sites.size() ; i++ ){
      Clique c;
      c.type = DATADRIVEN; /*ls.type = BESTLOWERSCALE;*/
      cliquesDuSite[ sites[i]->index ].push_back( cliques.size() );
      c.blobs.push_back( sites[i] );
      cliques.push_back( c );
        // cliquesDuSite[sites[i]->index].push_back(cliques.size());
        // ls.blobs.push_back(sites[i]);
        // cliques.push_back(ls);
  }
}

void BuildGlobalClique( vector<Site *> &sites,
                                  vector<vector<int> > &cliquesDuSite,
                                  vector<Clique> &cliques ) {
    Clique c;
    c.type = GLOBAL;

    for ( uint i = 0 ; i < sites.size() ; i++ ) {

        cliquesDuSite[ sites[i]->index ].push_back( cliques.size() );
        c.blobs.push_back( sites[i] );
    }
    cliques.push_back( c );


}

//##############################################################################

void SurfaceBased_StructuralAnalysis::GetSimilarityCliquesFromSSBCliques ( std::vector<surf::SSBClique> &ssbcliques,
                                std::vector<Site *> &sites,
                                std::vector<Clique> &cliques,
                                std::vector<std::vector<int> > &cliquesDuSite){

    cliquesDuSite = std::vector<std::vector<int> >( sites.size() );

    for ( uint i = 0 ; i < ssbcliques.size() ; i++ ) {
       Clique simc;
       simc.type = SIMILARITY;
       simc.rec = ssbcliques[i].similarity;
       surf::ScaleSpaceBlob *ssb1, *ssb2;
       int iSSB1, iSSB2;
       ssb1 = ssbcliques[i].ssb1;
       ssb2 = ssbcliques[i].ssb2;
       iSSB1 = ssb1->index;
       iSSB2 = ssb2->index;

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
                    cout << i << " " << aux << " " << cliques[aux].type << " " << cliques[aux].blobs.size() << " " << cliques[aux].blobs[0]->index << " " << cliques[aux].blobs[1]->index << endl;
                    ASSERT(false);
                }
            }
        }
    }
    std::cout << cliques.size() << "cliques recovered from ssbcliques" << std::endl;

}

//##############################################################################

// This function takes the "ssblobs" vector and figures out which pairs of blobs
//  overlap. The resulting vector "cliques" associates to every relevant pair of
//  scale-space blobs (noted by their indices) its calculated spatial overlap.

std::vector<surf::SSBClique> SurfaceBased_StructuralAnalysis::BuildSimilarityCliques ( std::vector<surf::ScaleSpaceBlob *>   &ssblobs,
                                                 std::vector<std::vector<surf::GreyLevelBlob *> > &matchingblobs ) {

    std::vector<surf::SSBClique > cliques;
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

                        pair<Point2df,Point2df> bbi = (*itB1)->get2DBoundingBox(),
                            //getBoundingBox((*itB1)->nodes, data[ssblobs[i]->subject].lat, data[ssblobs[i]->subject].lon),
                                                bbj = (*itB2)->get2DBoundingBox();
                            //(*itB2)->nodes, data[ssblobs[j]->subject].lat, data[ssblobs[j]->subject].lon);

                        Point2df bbmin1 (bbi.first[0], bbi.first[1]),
                                bbmax1 (bbi.second[0], bbi.second[1]),
                                bbmin2 (bbj.first[0], bbj.first[1]),
                                bbmax2 (bbj.second[0], bbj.second[1]) ;

                        uint no_overlap = 2;
                        double overlap = getOverlapMeasure( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
                // cout << "bbi("<< (*itB1)->nodes.size() << "):" << bbi.first[0] << "-" << bbi.first[1] << " " <<
                // bbi.second[0] << "-" << bbi.second[1] << " " <<
                // "bbj("<< (*itB2)->nodes.size() << "):" << bbj.first[0] << " " << bbj.first[1] << " " <<
                // bbj.second[0] << " " << bbj.second[1]  << endl;

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

                if ( overmax > 0.10 &&
                        !((ssblobs[j]->tmin > ssblobs[i]->tmax) || (ssblobs[i]->tmin > ssblobs[j]->tmax)) ) {

                    // If the two scale-space blobs have at least one pair of grey-level
                    //   overlapping (bounding-boxes) (+ scales overlapping), then a clique
                    // is created between these two ssb and the max-overlapping pair of glb
                    // is stored in "matchingblobs".



                    cliques.push_back(surf::SSBClique(ssblobs[i], ssblobs[j], overmax));
                    matchingblobs[i].push_back(b1max);
                    matchingblobs[j].push_back(b2max);
                    //cout << "max (" << ssblobs[i]->index <<","<< ssblobs[j]->index << ") between:" << b1max->index << " "
                    //    << b2max->index << " overmax:" << overmax << endl;
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
            cout << endl ;

    }

    return cliques;
}

////##############################################################################


std::vector<surf::SSBClique> SurfaceBased_StructuralAnalysis::BuildSimilarityCliques3D ( std::vector<surf::ScaleSpaceBlob *>   &ssblobs,
                                                   GroupData &data,
                                                   float threshold,
                                                   int type_distance ) {

    std::vector<surf::SSBClique > cliques;

    std::set<surf::GreyLevelBlob *>::iterator itB1, itB2;
    surf::GreyLevelBlob *b1max, *b2max;

    // Start of cliques construction
    if ( type_distance == DISTANCE_LATITUDES ) 
        std::cout << "DISTANCE_LATITUDES" << std::endl;
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
                    distance = sqrt(10*p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
                }
                else if ( type_distance == DISTANCE_3DEUCLIDIAN ) {
                    Point3df p1 (b1max->raw_coordinates[max1][0], b1max->raw_coordinates[max1][1], b1max->raw_coordinates[max1][2]);
                    Point3df p2 (b2max->raw_coordinates[max2][0], b2max->raw_coordinates[max2][1], b2max->raw_coordinates[max2][2]);
                    Point3df p = p1-p2;
                    distance = p.norm();
                }

                // std::cout << distance << " " << std::flush;

                if ( distance < threshold ) {
                    cliques.push_back( surf::SSBClique(ssblobs[i], ssblobs[j], distance ) );
                }
            }
            // The next pair of scale-space blobs will now be processed.
        }
    }
    std::cout << endl;
    std::cout << ssblobs.size() << "/" << ssblobs.size() << "(" << cliques.size() << ")" << std::endl;
    return cliques;
}

//##############################################################################


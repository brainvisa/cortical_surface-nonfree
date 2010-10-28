#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/cliques.h>


using namespace aims;
using namespace carto;
using namespace std;

float Clique::ddweight, 
      Clique::intrapsweight, 
      Clique::simweight, 
      Clique::lsweight, 
      Clique::ddx1, 
      Clique::ddx2, 
      Clique::simx1, 
      Clique::simx2, 
      Clique::ddh, 
      Clique::globalweight;

void Clique::setParameters ( float _ddweight, 
                             float _intrapsweight, 
                             float _simweight, 
                             float _lsweight, 
                             float _ddx1, 
                             float _ddx2,
                             float _simx1,
                             float _simx2,
                             float _ddh, 
                             float _globalweight ){
    ddweight = _ddweight; 
    intrapsweight = _intrapsweight; 
    simweight = _simweight; 
    lsweight = _lsweight;
    ddx1 = _ddx1; 
    ddx2 = _ddx2; 
    simx1 = _simx1;
    simx2 = _simx2;
    ddh = _ddh;
    globalweight = _globalweight;
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

//        std::map<int, std::set<std::string> >::iterator itset;
//        for ( itset = subjectscount.begin() ; itset != subjectscount.end() ; itset++ ) {
//            std::cout << (*itset).second.size() << " " << std::flush;
//        }
//        std::cout << std::endl;

    }
}


double getOverlap(Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap){

  float overlap_x,overlap_y,aux;
  double rec = 0.0;

  if ( sqrt(pow(bbmin1[0]-bbmax1[0],2)) < 0.0001 ) { bbmax1[0] += 0.5; /*cout << "bbmax10+ ";*/ }
  if ( sqrt(pow(bbmin1[1]-bbmax1[1],2)) < 0.0001 ) { bbmax1[1] += 0.5; /*cout << "bbmax11+ ";*/ }
  if ( sqrt(pow(bbmin2[0]-bbmax2[0],2)) < 0.0001 ) { bbmax2[0] += 0.5; /*cout << "bbmax20+ ";*/ }
  if ( sqrt(pow(bbmin2[1]-bbmax2[1],2)) < 0.0001 ) { bbmax2[1] += 0.5; /*cout << "bbmax21+ ";*/ }
  // if (bbmin1[1]>bbmax1[1] && bbmin2[1] < bbmax2[1] ) {//alors i a bouclé autour de 360/0
  if ( sqrt(pow(bbmin1[1]-bbmax1[1],2)) >300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) < 300 ) {
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
  else if ( sqrt(pow(bbmin1[1]-bbmax1[1],2)) <300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) > 300 ) {
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
  else if ( sqrt(pow(bbmin1[1]-bbmax1[1],2)) >300 && sqrt(pow(bbmin2[1]-bbmax2[1],2)) > 300 ) {
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
if ( sqrt(pow(bbmin1[0]-bbmax1[0],2)) >150 && sqrt(pow(bbmin2[0]-bbmax2[0],2)) < 150 ) {
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


void BuildDataDrivenCliques ( std::vector <Site *> &sites,
                              std::vector < std::vector<int> > &cliquesDuSite,
                              std::vector < Clique > &cliques ) {

    for ( uint i = 0 ; i < sites.size() ; i++ ){
        Clique c;
        c.type = DATADRIVEN;
        cliquesDuSite[ sites[i]->index ].push_back( cliques.size() );
        c.blobs.push_back( sites[i] );
        cliques.push_back( c );
    }
}

void BuildLowerScaleCliques ( std::vector < Site *> &sites,
                              std::vector < std::vector<int> > &cliquesDuSite,
                              std::vector < Clique > &cliques ) {

    for ( uint i = 0 ; i < sites.size() ; i++ ){
        Clique ls;
        ls.type = BESTLOWERSCALE;
        cliquesDuSite[ sites[i]->index ].push_back( cliques.size() );
        ls.blobs.push_back( sites[i] );
        cliques.push_back( ls);
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

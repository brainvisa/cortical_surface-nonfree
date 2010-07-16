#include <aims/getopt/getopt2.h>
#include <cortical_surface/structuralanalysis/cliques.h>

using namespace aims;
using namespace carto;
using namespace std;

float Clique::ddweight, Clique::intrapsweight, Clique::simweight, Clique::lsweight, Clique::ddx2, Clique::ddx1, Clique::ddh;
void Clique::setParameters(float _ddweight, float _intrapsweight, float _simweight, float _lsweight, float _ddx2, float _ddx1, float _ddh){
  ddweight=_ddweight; intrapsweight = _intrapsweight; simweight=_simweight; lsweight=_lsweight; ddx2 =_ddx2;  ddx1 = _ddx1; ddh=_ddh;
}

void Clique::updateLabelsCount(){
  if (type == INTRAPRIMALSKETCH){
    labelscount = map<int, uint>();
    for (uint i=0;i<blobs.size();i++){
      if (labelscount.find(blobs[i]->label) == labelscount.end())
        labelscount[blobs[i]->label] = 1;
      else
        labelscount[blobs[i]->label]++;
    }
    map<int,uint>::iterator it;
    uint chksum=0;
//     cout << blobs[0]->subject << endl;
    for (it=labelscount.begin();it!=labelscount.end();it++){
//       cout << it->first << " " << it->second << "-";
      chksum += (*it).second;
    }
//     cout  << endl;
    ASSERT(chksum == blobs.size());
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


void ConstruireCliquesIntraPS( vector<Site *> &sites,
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


void ConstruireCliquesDataDriven( vector<Site *> &sites,
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


void ConstruireCliquesSimilarity(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, vector<Clique> &cliques){
  uint temp=0,temp2=0,temp3=0,temp4=0;
  double rec;
  set<uint> sitesaux, sitesindexes;
  float cliques_thresh=40.0;
  cout << "Similarity cliques distance limit : << " << cliques_thresh << " >>" << endl;
  for (uint i=0;i<sites.size();i++){
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << sites.size() << "(" << sites[i]->rank << "-" << sites[i]->subject << "-" << sites[i]->nodes_list.size() << ") " << flush;
    for (uint j=i;j<sites.size();j++) {
      if (sites[i]->subject != sites[j]->subject) {

//         uint cpt0=0;
//         for (it=sitesnodes.begin();*it!=sites[i]->node;it++,cpt0++){}
//         if (distmaps[sites[i]->node].find(sites[j]->node) == distmaps[sites[i]->node].end())
//           rec= 600.0;
//         else 
//           rec = distmaps[sites[i]->node][sites[j]->node];
//         rec = distmaps[sites[i]->node][j];
        
        // DISTANCE AVEC SURFACES
//           (*jv)->getProperty( "boundingbox_max", bbmax_2);
//           (*jv)->getProperty( "boundingbox_min", bbmin_2);
          Point3df bbmax1=sites[i]->boundingbox_max, bbmax2=sites[j]->boundingbox_max;
          Point3df bbmin1=sites[i]->boundingbox_min, bbmin2=sites[j]->boundingbox_min;
         // on s'occupe de la longitude
uint no_overlap=1;
        rec = getOverlap(bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap);
//         if (rec >0.1) cout << "REC: " << bbmin1[0] << " " << bbmin1[1] << " " << bbmax1[0] << " " << bbmax1[1] << " " << bbmin2[0] << " " << bbmin2[1] << " " << bbmax2[0] << " " << bbmax2[1] << " " << endl;
//         rec = sqrt(pow(meshes[sites[0]->subject].vertex()[sites[i]->node][0]-meshes[sites[0]->subject].vertex()[sites[j]->node][0],2)+pow(meshes[sites[0]->subject].vertex()[sites[i]->node][1]-meshes[sites[0]->subject].vertex()[sites[j]->node][1],2)+pow(meshes[sites[0]->subject].vertex()[sites[i]->node][2]-meshes[sites[0]->subject].vertex()[sites[j]->node][2],2));  // USING EUCLIDEAN DISTANCE
//         cout << i<< ";" << j << ";" << sites[i]->node << "-" << sites[j]->node << "rec:" << rec << endl;
        if (rec < 0.000000001) {sitesaux.insert(sites[i]->node); sitesaux.insert(sites[j]->node); sitesindexes.insert(i); sitesindexes.insert(j);}
        
        if ((no_overlap==0 && rec < cliques_thresh) && !((sites[j]->tmin > sites[i]->tmax) || (sites[i]->tmin > sites[j]->tmax)) /*&& (sites[i]->tValue * sites[j]->tValue > 2.0)*/) {
//           if (rec < 0.2) cout << "BLIP " << endl;
          Clique simc;
          simc.type = SIMILARITY;
          simc.rec = rec;
//           cout << simc.rec << endl ;
          cliquesDuSite[sites[i]->index].push_back(cliques.size());
          cliquesDuSite[sites[j]->index].push_back(cliques.size());

          simc.blobs.push_back(sites[i]);
          simc.blobs.push_back(sites[j]);

//           cout << "(("<<simc.rec << " (lab:" << simc.blobs[0]->label << " t:"<< simc.blobs[0]->t << "(nod:" <<  simc.blobs[0]->node << " suj:"<<simc.blobs[0]->subject<<")" << "-lab:" << simc.blobs[1]->label << " t:"<< simc.blobs[1]->t << "(nod:" <<  simc.blobs[1]->node << " suj:"<<simc.blobs[1]->subject<<")" << ")=>" << simc.energie << ")) " << endl;
          cliques.push_back(simc);
  //           }
        }
      }
    }
  }
  cout << sitesaux.size() << " " << sitesindexes.size() << endl;
  cout << "TEMP :" << temp<< " " << temp2 << " " << temp3 << " " << temp4 << " "<<temp+temp4<< " ";
    for (uint i=0;i<cliques.size();i++){
      if (cliques[i].type == SIMILARITY){
          assert(cliques[i].blobs[0]->index <sites.size());
          assert(cliques[i].blobs[1]->index <sites.size());
      }
    } 


}



vector<Clique> ConstruireCliques(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite)  {
  
  vector<Clique> cliques;

  cliquesDuSite = vector<vector<int> >(sites.size());

  ConstruireCliquesIntraPS(sites, cliquesDuSite, cliques);
  
  ConstruireCliquesDataDriven(sites, cliquesDuSite, cliques);

  ConstruireCliquesSimilarity(sites, cliquesDuSite, cliques);
    
  return cliques;
}


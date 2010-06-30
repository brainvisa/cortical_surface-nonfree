
#include <aims/data/data_g.h>

using namespace std;

int find(const vector<uint> &v, uint item){
   int c=-1;
   for (uint i=0; i<v.size() && c==-1 ; i++){
      if (item == v[i]) c=(int)i;
   }
   return c;
}

int find(const vector<short> &v, short item){
   int c=-1;
   for (uint i=0; i<v.size() && c==-1 ; i++){
      if (item == v[i]) c=(int)i;
   }
   return c;
}

void push_vector(vector<uint> &v, const vector<uint> &w){
   for (uint i=0;i<w.size();i++)
      v.push_back(w[i]);   
}

void push_vector(vector<short> &v, const vector<short> &w){
   for (uint i=0;i<w.size();i++)
      v.push_back(w[i]);   
}

vector<short> reverseVector(const vector<short> &v){
   vector<short> result;
   for (uint i=0;i<v.size();i++)
      result.push_back(v[v.size()-1-i]);
   return result;
}
   
vector<uint> setToVector(const set<uint> &v){
   set<uint>::iterator it;
   vector<uint> result;
   for (it = v.begin();it!=v.end();it++)
      result.push_back(*it);

   return result;
}

vector<uint> setToVector(const set<int> &v){
    set<int>::iterator it;
    vector<uint> result;
    for (it = v.begin();it!=v.end();it++)
        result.push_back(*it);
    
    return result;
}

vector<short> setToVector(const set<short> &v){
   set<short>::iterator it;
   vector<short> result;
   for (it = v.begin();it!=v.end();it++)
      result.push_back(*it);

   return result;
}

vector<uint> getCorresVector(const vector<uint> &inVector, const vector<uint> &corres){
   vector<uint> outVector(inVector.size());
   for (uint i=0;i<inVector.size();i++)
      outVector[i]=corres[inVector[i]];
   return outVector;
}

vector<uint> minusVector(const vector<uint> &v, const vector<uint> &w){
   vector<uint> result;

   for (uint i=0;i<v.size();i++)
      if (find(w,v[i])==-1) result.push_back(v[i]);

   return result;
}

vector<short> minusVector(const vector<short> &v, const vector<short> &w){
   vector<short> result;

   for (uint i=0;i<v.size();i++)
      if (find(w,v[i])==-1) result.push_back(v[i]);

   return result;
}

vector<uint> intersVector(const vector<uint> &v, const vector<uint> &w){
   vector<uint> result;

   for (uint i=0;i<v.size();i++)
      if (find(w,v[i])!=-1) result.push_back(v[i]);

   return result;
}

int getRealSize(const vector<set<uint> > &v){
   int result=0;
   for (uint i=0;i<v.size();i++)
      if (v[i].size() != 0) result++;
   return result;
}

bool valueOf(const std::string &s, uint &obj){
   std::istringstream is(s);
   return is >> obj;
}

string f2str(float f){
   stringstream sstr;  string s;   sstr << f;   sstr >> s;    return s;
}

string i2str(uint f){
   stringstream sstr;  string s;   sstr << f;   sstr >> s;    return s;
}

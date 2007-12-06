#ifndef AIMS_PARAMETERIZEGYRI_VECTOR_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_VECTOR_OPERATIONS_H


using namespace std;

int find(const vector<uint> &v, uint item);

int find(const vector<short> &v, short item);

void push_vector(vector<uint> &v, const vector<uint> &w);

void push_vector(vector<short> &v, const vector<short> &w);

vector<short> reverseVector(const vector<short> &v);
   
vector<uint> setToVector(const set<uint> &v);

vector<short> setToVector(const set<short> &v);

vector<uint> getCorresVector(const vector<uint> &inVector, const vector<uint> &corres);

vector<uint> minusVector(const vector<uint> &v, const vector<uint> &w);

vector<short> minusVector(const vector<short> &v, const vector<short> &w);

vector<uint> intersVector(const vector<uint> &v, const vector<uint> &w);

int getRealSize(const vector<set<uint> > &v);

bool valueOf(const std::string &s, uint &obj);

string f2str(float f);

string i2str(uint f);

#endif


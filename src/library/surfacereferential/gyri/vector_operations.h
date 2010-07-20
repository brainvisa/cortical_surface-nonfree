#ifndef AIMS_PARAMETERIZEGYRI_VECTOR_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_VECTOR_OPERATIONS_H


int find(const std::vector<uint> &v, uint item);

int find(const std::vector<short> &v, short item);

void push_vector(std::vector<uint> &v, const std::vector<uint> &w);

void push_vector(std::vector<short> &v, const std::vector<short> &w);

std::vector<short> reverseVector(const std::vector<short> &v);
   
std::vector<uint> setToVector(const std::set<uint> &v);

std::vector<uint> setToVector(const std::set<int> &v);

std::vector<short> setToVector(const std::set<short> &v);

std::vector<uint> getCorresVector(const std::vector<uint> &inVector, const std::vector<uint> &corres);

std::vector<uint> minusVector(const std::vector<uint> &v, const std::vector<uint> &w);

std::vector<short> minusVector(const std::vector<short> &v, const std::vector<short> &w);

std::vector<uint> intersVector(const std::vector<uint> &v, const std::vector<uint> &w);

int getRealSize(const std::vector<std::set<uint> > &v);

bool valueOf(const std::string &s, uint &obj);

std::string f2str(float f);

std::string i2str(uint f);

#endif


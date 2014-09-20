#ifndef GIBBS_UTIL_H
#define GIBBS_UTIL_H

#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

template <typename T> 
vector<T>* string_to_vector(string s, T (*convert)(const char*) );

template <typename T>
string vector_to_string(vector<T> *v);

#endif

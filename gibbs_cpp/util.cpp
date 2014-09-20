#include <vector>
#include <iostream>
#include <sstream>
#include <cstring>

using namespace std;

template <typename T> 
vector<T>* string_to_vector(string s, T (*convert)(const char*) ) {
    vector<T> *v = new vector<T>();

    const char *cstr = s.c_str();
    char *pch = strtok((char*)cstr, "[], ");
    while (pch != NULL) {
        //cout << "element: " << pch << endl;
        v->push_back( (*convert)(pch) );

        pch = strtok(NULL, "[], ");
    }

    return v;
}

template <typename T>
string vector_to_string(vector<T> *v) {
    ostringstream oss;

    oss << "[";
    for (unsigned int i = 0; i < v->size(); i++) {
        if (i > 0) oss << ",";
        oss << v->at(i);
    }
    oss << "]";

    return oss.str();
}

template 
string vector_to_string<long double>(vector<long double> *v);


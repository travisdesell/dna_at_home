#include <iostream>
#include <vector>
#include <sstream>

#ifndef GIBBS_ARGUMENTS_H
#define GIBBS_ARGUMENTS_H

using namespace std;

template <typename T>
void get_argument_vector(vector<string> arguments, string argument, bool required, vector<T> &results) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            i++;
            while (arguments.at(i).substr(0,2).compare("--") != 0) {
                T result;
                if ( !(stringstream(arguments.at(i++)) >> result) ) {
                    cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                    exit(0);
                }
                results.push_back(result);
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(0);
    }

    cerr << "parsed argument '" << argument << "' successfully:";
    for (unsigned int i = 0; i < results.size(); i++) {
        cerr << " " << results.at(i);
    }
    cerr << endl;
}

template <>
void get_argument_vector<string>(vector<string> arguments, string argument, bool required, vector<string> &results) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            i++;
            while (arguments.at(i).substr(0,2).compare("--") != 0) {
                results.push_back(arguments.at(i++));
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(0);
    }

    cerr << "parsed argument '" << argument << "' successfully:";
    for (unsigned int i = 0; i < results.size(); i++) {
        cerr << " " << results.at(i);
    }
    cerr << endl;
}


template <typename T>
void get_argument(vector<string> arguments, string argument, bool required, T &result) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            if ( !(stringstream(arguments.at(++i)) >> result) ) {
                cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                exit(0);
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(0);
    }

    cerr << "parsed argument '" << argument << "' successfully: " << result << endl;
}

template <>
void get_argument<string>(vector<string> arguments, string argument, bool required, string &result) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            result = arguments.at(++i);
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(0);
    }

    cerr << "parsed argument '" << argument << "' successfully: " << result << endl;
}


template <typename T1, typename T2>
void get_arguments(vector<string> arguments, string argument, bool required, T1 &result1, T2 &result2) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            if ( !(stringstream(arguments.at(++i)) >> result1) ) {
                cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                exit(0);
            }
            if ( !(stringstream(arguments.at(++i)) >> result2) ) {
                cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                exit(0);
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(0);
    }
    cerr << "parsed argument '" << argument << "' successfully: " << result1 << " " << result2 << endl;
}

bool argument_exists(vector<string> arguments, string argument) {
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            cerr << "parsed argument '" << argument << "' successfully." << endl;
            return true;
        }
    }
    return false;
}


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

#endif

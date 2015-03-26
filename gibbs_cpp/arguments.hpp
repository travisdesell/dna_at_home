#ifndef GIBBS_ARGUMENTS_H
#define GIBBS_ARGUMENTS_H

#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

template <typename T>
bool get_argument_vector(vector<string> arguments, string argument, bool required, vector<T> &results) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            i++;
            while (i < arguments.size() && arguments.at(i).substr(0,2).compare("--") != 0) {
                T result;
                if ( !(stringstream(arguments.at(i++)) >> result) ) {
                    cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                    exit(1);
                }
                results.push_back(result);
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(1);
    }

    if (found) {
        cerr << "parsed argument '" << argument << "' successfully:";
        for (unsigned int i = 0; i < results.size(); i++) {
            cerr << " " << results.at(i);
        }
        cerr << endl;
    }
    return found;
}

template <>
bool get_argument_vector<string>(vector<string> arguments, string argument, bool required, vector<string> &results) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            i++;
            while (i < arguments.size() && arguments.at(i).substr(0,2).compare("--") != 0) {
                results.push_back(arguments.at(i++));
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(1);
    }

    if (found) {
        cerr << "parsed argument '" << argument << "' successfully:";
        for (unsigned int i = 0; i < results.size(); i++) {
            cerr << " " << results.at(i);
        }
        cerr << endl;
    }
    return found;
}


template <typename T>
bool get_argument(vector<string> arguments, string argument, bool required, T &result) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            if ( !(stringstream(arguments.at(++i)) >> result) ) {
                cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                exit(1);
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(1);
    }

    if (found) {
        cerr << "parsed argument '" << argument << "' successfully: " << result << endl;
    }
    return found;
}

template <>
bool get_argument<string>(vector<string> arguments, string argument, bool required, string &result) {
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
        exit(1);
    }

    if (found) {
        cerr << "parsed argument '" << argument << "' successfully: " << result << endl;
    }
    return found;
}


template <typename T1, typename T2>
bool get_arguments(vector<string> arguments, string argument, bool required, T1 &result1, T2 &result2) {
    bool found = false;
    for (unsigned int i = 0; i < arguments.size(); i++) {
        if (argument.compare(arguments.at(i)) == 0) {
            if ( !(stringstream(arguments.at(++i)) >> result1) ) {
                cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                exit(1);
            }
            if ( !(stringstream(arguments.at(++i)) >> result2) ) {
                cerr << "ERROR: invalid argument '" << argument << "': " << arguments.at(i) << endl;
                exit(1);
            }
            found = true;
            break;
        }
    }

    if (required && !found) {
        cerr << "ERROR: argument '" << argument << "' required and not found." << endl;
        exit(1);
    }

    if (found) {
        cerr << "parsed argument '" << argument << "' successfully: " << result1 << " " << result2 << endl;
    }
    return found;
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

#endif

#ifndef STRUCTS_H
#define STRUCTS_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

#define ALPHABET_LENGTH 4
#define ALPHABET "ACGT"

#define A_POSITION 0
#define C_POSITION 1
#define G_POSITION 2
#define T_POSITION 3

#define A_POSITION_REVERSE 3
#define C_POSITION_REVERSE 2
#define G_POSITION_REVERSE 1
#define T_POSITION_REVERSE 0

#ifndef _DEBUG_
#define at(x) operator[](x)
#endif

#endif

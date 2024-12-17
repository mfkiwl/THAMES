/**
@file JSONParser.h
@brief Method declaration for JSONParser portable JSON parser.

*/
#ifndef JSONPARSERH
#define JSONPARSERH

#include "global.h"
#include "utils.h"
#include "valid.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <typeinfo>
#include <variant>
#include <vector>

using namespace std;

/**
@struct PoreSizeVolume
@brief Volume fraction of a sub-voxel pore of a giveen effective diameter
*/

struct PoreSizeVolume {
  double diam;
  double volume;
  double volfrac;
};

/**
@class JSONNode
@brief Handles the tracking of phases and communications between GEM and THAMES.
**/
class JSONNode {
  enum class Type { OBJECT, LIST, STRING, NUMBER, BOOLEAN, NULL_TYPE };
};
#endif

#ifndef MAPPING_H
#define MAPPING_H

#include <util/containers/DArray.h> 

using namespace Util;

class clusterInfo;

void mapping (clusterInfo* step0, clusterInfo* step1, double cutoff_C);

DArray<int> reinitializeArray (int size);

#endif

#ifndef MAP_CLUSTER_H
#define MAP_CLUSTER_H

#include <util/containers/DArray.h>              // member template
#include <util/containers/GArray.h>              // member template
#include <iostream>

using namespace Util;

void mapping 
      // Work array to calculate contribution of each cluster
      // This array would be of size = number of clusters at
      // previous timestep
      DArray<int > workArray;



#endif

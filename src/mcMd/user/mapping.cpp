#ifndef MAPPING_CPP
#define MAPPING_CPP

#include "mcMd/user/mapping.h"
#include "mcMd/user/clusterInfo.h"
#include <util/containers/DArray.h>              // member template

using namespace Util;

void mapping (clusterInfo* step0, clusterInfo* step1, int cutoff_U, double cutoff_C)
{
   // The following lines declare workArrays to find 
   // the contributions. It is important to remember
   // that the indices of this array would correspond
   // to the array clusterIds of the class clusterInfo.
   // So to find out the respective clusterId you will 
   // have to use that member of the clusterInfo class.
   
   // This is a workArray to analyze the contribution
   // of clusters at step0 to step1
   DArray < int> contribution0;
   // This is a workArray to analyze the contribution
   // of clusters at step0 to step1
   DArray < int> contribution1;
   
   contribution0.allocate (step0->nClusters);
   contribution1.allocate (step1->nClusters);




   // Remember to set nCluster here
}

#endif

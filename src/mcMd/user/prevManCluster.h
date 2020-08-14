#ifndef MAN_CLUSTER_H
#define MAN_CLUSTER_H

#include <util/containers/DArray.h>              // member template
#include <util/containers/GArray.h>              // member template


using namespace Util;

class manCluster{

   public:

      void readStep();

      void mapClusters();

      void writeStep();

      // Consists of all the information on clusters and comprising molecules
      // at the respective timesteps
      GArray< GArray<int > > prevStep;
      GArray< GArray<int > > nextStep;

      // Cluster IDs presently in use
      GArray< int > clusterIds;

      // Max cluster ID in use
      //int maxClusterID;

      // number of samples examined
      int iSample;    

      // Work array to calculate contribution of each cluster
      DArray<int> workArray;

};

#endif

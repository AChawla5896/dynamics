#ifndef CLUSTER_INFO_H
#define CLUSTER_INFO_H

#include <util/containers/DArray.h>              // member template
#include <util/containers/GArray.h>              // member template

#include <vector>
#include <iterator>
#include <fstream>
#include <iostream>

using namespace Util;

class clusterInfo{

   public:

      // readStep0 function will also calculate the number of molecules 
      // in the simulation and return it. Thereafter, the DArray, 
      // whichCluster will be allocated. This needs to be called only
      // once and in the beginning.
      // After this only call readStep.
      // cutoff_U is required as a parameter so as to call oragnize.
      // See organize for its details. 
      void readStep0 (std::istream& in, int cutoff_U);
      void readStep (std::istream& in, int cutoff_U);

      // This function will organize the read in clusters so as to have
      // the cluster [0] consist of clusters having aggregation number 
      // less than equal to cutoff_U. Rest of the clusters will remain 
      // untouched. nClusters, clusterId and whichCluster should be changed 
      // accordingly.
      void organize_clus (std::vector< std::vector<int > > clustersRead, int cutoff_U);

      void writeStep (std::ostream& out);

      // Will find the number of lines in the output file
      // and also update nClusters
      void countLines (std::istream& in);

      // Will be used to update the clusterID from prevID to 
      // newID after mapping part of the cluster
      // Things to update: clusterIds, whichCluster
      void updateClusterId (int prevId, int newId);

      // Need to clear all the std::vectors after writing new cluster outputs
      // Do not clear the allocation of DArray (but can think of assingning 
      // a random value). Call before you read next time step (to be safe).
      void clear();

      // Consists of all the information on clusters and comprising molecules
      // Outer vector size: nClusters
      // Inner vector size: aggregation number of each of the clusters
      // clusters[0] has the list of all the molecules part of clusters which 
      // have an aggregation number less than equal to cutoff_U.
      // See: organize
      std::vector< std::vector<int > > clusters;

      // Cluster IDs presently in use. Need to update this when you have 
      // mapped the clusters
      // Size: nClusters
      std::vector<int > clusterIds;

      // This size of this array needs to be equal to the number of molecules
      // Will be used to get the clusterId using the molecule
      DArray<int > whichCluster;

      // total number of clusters at that timestep, update when you find 
      // the total number of lines
      int nClusters;

      // Total number of molecules in the simulation
      int nMolecules;

      // Max cluster ID in use
      // remember to change maxClusterId when mapping
      static int maxClusterId;

};

#endif

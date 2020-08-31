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
      // This function also resets isProcessed
      // cutoff_U is required as a parameter so as to call oragnize.
      // See organize for its details. 
      void readStep0 (std::istream& in, int cutoff_U);
      void readStep (std::istream& inClusters, std::istream& inCOMs, std::istream& inMoments, int cutoff_U);

      // This function will organize the read in clusters so as to have
      // the cluster [0] consist of clusters having aggregation number 
      // less than equal to cutoff_U. Rest of the clusters will remain 
      // untouched. nClusters, clusterId and whichCluster should be changed 
      // accordingly.
      void organize (std::vector< std::vector<int > > clustersRead, std::vector< std::vector<double > > COMsRead, 
                        std::vector< std::vector<double > > momentsRead, int cutoff_U);

      // This function will be used to find the index of the given clusterId
      // in the vector cluster Id
      // During the mapping process (i.e. when the whole process is not 
      // complete), it is possible that two clusters have the same cluster Ids.
      // So an optional vector notAccept is also taken as an input. This vector
      // would consist of the indices which are already processed and hence not 
      // accepted as the correct output of the cluster Index function
      int clusterIndex (int Id);

      void writeStep (std::ostream& outClusters, std::ostream& outCOMs, std::ostream& outMoments);

      // Will find the number of lines in the output file
      // and also update nClusters
      void countLines (std::istream& in);

      // Will be used to update the clusterID from prevID to 
      // newID after mapping part of the cluster
      // Things to update: clusterIds, whichCluster
      void updateClusterId (int prevId, int newId);

      // Need to clear all the std::vectors after writing new cluster outputs
      // Do not clear the allocation of DArray (but can think of assingning 
      // a random value).
      // Call before you read next time step (to be safe).
      void clear();

      // Consists of all the information on clusters and comprising molecules
      // Outer vector size: nClusters
      // Inner vector size: aggregation number of each of the clusters
      std::vector< std::vector<int > > clusters;

      // Consists of all the information on COMs
      // Outer vector size: nClusters
      // Inner vector size: 3
      std::vector< std::vector<double > > COMs;

      // Consists of all the information on momentsTensors
      // Outer vector size: nClusters
      // Inner vector size: 9
      std::vector< std::vector<double > > moments;

      // Consists of the list of all the molecules part of clusters which 
      // have an aggregation number less than equal to cutoff_U.
      // See: organize
      std::vector<int > melt;

      // Cluster IDs presently in use. Need to update this when you have 
      // mapped the clusters
      // melt clusterId is not present in this vector. It is always denoted by 
      // clusterId = 0
      // Size: nClusters
      // Would initially range from 1 - (nClusters)
      std::vector<int > clusterIds;

      // This size of this array needs to be equal to the number of molecules
      // Will be used to get the clusterId using the molecule
      // Size: nMolecules
      DArray<int > whichClusterId;

      // This vector will keep track of the molecules which have had their 
      // cluster Id updated.
      // At the end of mapping all the molecules would have had their Id 
      // updated.
      // Size: nMolecules
      DArray<bool > isProcessed;

      // total number of clusters at that timestep (does not include melt),
      // Therefore, total number of clusters including melt will be 
      // nClusters + 1
      // update when you find 
      // the total number of lines
      int nClusters;

      // Total number of molecules in the simulation
      int nMolecules;

      // Max cluster ID in use
      // remember to change maxClusterId when mapping
      static int maxClusterId;

};

#endif

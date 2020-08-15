#ifndef CLUSTER_INFO_CPP
#define CLUSTER_INFO_CPP	

#include "mcMd/user/clusterInfo.h"

#include <util/containers/DArray.h>              // member template
#include <util/containers/GArray.h>              // member template
#include <util/containers/ArrayIterator.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace Util;

void clusterInfo::readStep0 (std::istream& in, int cutoff_U) 
{
   // The clusterId will be read as string and then 
   // converted to integer. These are the two corresponding 
   // variables for it. Will be used to set clusterIds, i.e.
   // the cluster IDs presently in use.  
   std::string clusterIdSt;
   int clusterIdInt;

   // The aggregation number will be read as string and then 
   // converted to integer. These are the two corresponding 
   // variables for it.
   std::string clusterAggSt;
   int clusterAggInt;

   // Work Array to store in read in clusters
   std::vector< std::vector<int > > clustersRead;

   // Work Array to store in read in clusterIds
   // You might not need this.
   std::vector<int > clusterIdsRead;

   // Will be used to keep track of total number of molecules
   // present in the simulation 
   // The numbering of the molecules starts from 0 in the 
   // output
   nMolecules = -1;

   // Calling countLines to count the number of clusters
   // in the file. nClusters is updated.
   countLines (in); 

   // A while loop can go to the end of file. This loop 
   // will fill in the corresponding values.

   int iCluster = 0;
   std::string line;

   while (getline( in, line )) {

      std::istringstream linestream( line );

      // Reads in the clusterId from the file and sets it in the
      // vector clusterIds.
      std::getline(linestream, clusterIdSt, '(' );
      clusterIdInt = stoi(clusterIdSt);
      clusterIdsRead.push_back (clusterIdInt);

      std::getline(linestream, clusterAggSt, ')' );
      clusterAggInt = stoi(clusterAggSt);
         
      clustersRead.push_back ( std::vector<int>( std::istream_iterator<int>(linestream), 
                                  std::istream_iterator<int>()) );

      iCluster++;

   }

   // Resetting the ifstream to the beginning of the file
   in.clear();
   in.seekg (0, std::ios::beg);    

   // Assert iCluster is equal to nClusters
   if (iCluster != nClusters) {
      std::cout << "iCluster is not equal to nCluster" << std::endl;
   }

   // orgainze the clusters so as to have clusters[0] have a list of molecules
   // in clusters having aggregation number less than equal to a cutoff
  //  organize (clustersRead, cutoff_U);

   // Finding the total number of molecules in the simulation
   // also allocating the DArray, whichCluster 
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      for (int iMol = 0; iMol < clusters [iCluster].size(); iMol++) {
         if ( clusters [iCluster] [iMol] > nMolecules ) {
            nMolecules = clusters [iCluster] [iMol];
         }
      }
   }
   // Incrementing because the numbering of molecules starts from 0
   nMolecules++;
   whichCluster.allocate(nMolecules);
 
   // Setting the DArray which cluster.
   // Will be used in mapping
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      for (int iMol = 0; iMol < clusters [iCluster].size(); iMol++) {
         whichCluster [ clusters [iCluster] [iMol] ] = clusterIds [iCluster];
      }
   }
}

void clusterInfo::readStep (std::istream& in, int cutoff_U)
{
   // The clusterId will be read as string and then 
   // converted to integer. These are the two corresponding 
   // variables for it. Will be used to set clusterIds, i.e.
   // the cluster IDs presently in use.  
   std::string clusterIdSt;
   int clusterIdInt;

   // The aggregation number will be read as string and then 
   // converted to integer. These are the two corresponding 
   // variables for it.
   std::string clusterAggSt;
   int clusterAggInt;

   // Work Array to store in read in clusters
   std::vector< std::vector<int > > clustersRead;

   // Work Array to store in read in clusterIds
   // You might not need this.
   std::vector<int > clusterIdsRead;

   // Calling countLines to count the number of clusters
   // in the file. nClusters is updated.
   countLines (in);  

   // A while loop can go to the end of file. This loop 
   // will fill in the corresponding values.

   int iCluster = 0;
   std::string line;

   while (getline( in, line )) {

      std::istringstream linestream( line );

      // reads in the clusterId from the file and sets it in the
      // vector clusterIds.
      std::getline(linestream, clusterIdSt, '(' );
      clusterIdInt = stoi(clusterIdSt);
      clusterIdsRead.push_back (clusterIdInt);

      std::getline(linestream, clusterAggSt, ')' );
      clusterAggInt = stoi(clusterAggSt);
      
      clustersRead.push_back ( std::vector<int>( std::istream_iterator<int>(linestream), 
                                  std::istream_iterator<int>()) );

      iCluster++;

   }

   // Resetting the ifstream to the beginning of the file
   in.clear();
   in.seekg (0, std::ios::beg);

   // Assert iCluster is equal to nClusters
   if (iCluster != nClusters) {
      std::cout << "iCluster is not equal to nCluster" << std::endl;
   }

   // orgainze the clusters so as to have clusters[0] have a list of molecules
   // in clusters having aggregation number less than equal to a cutoff
   //organize (clustersRead, cutoff_U);

   // Setting the DArray which cluster.
   // Will be used in mapping
   // Assumes it has already been allocated
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      for (int iMol = 0; iMol < clusters [iCluster].size(); iMol++) {
         whichCluster [ clusters [iCluster] [iMol] ] = clusterIds [iCluster]; 
      }
   }
}

void organize (std::vector< std::vector<int > > clustersRead, int cutoff_U)
{
//   // Initialize the first element as -1 and then delete it after reorganizing
//   // Change nClusters
//   // Have clusterIds go from 0 - (nClusters-1)
//   
//   // Will be deleted after reorganizing
     clusters [0][0] = -1;
// 
//   for (int iCluster = 0; iCluster < clustersRead.size(); iCluster++) {
//      if (clustersRead [iCluster].size() <= cutoff_U) {
//         //clusters[0].insert (clusters[0].end(), clustersRead[iCluster].begin(), 
//         //                       clustersRead[iCluster].end());
//      }
//      else {
//         clusters.push_back(clustersRead[iCluster]);
//      }
//   }
//
//   //clusters[0].erase(0);
//
//   // Resetting nClusters after organizing
//   nClusters = clusters.size();
//
//   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
//      clusterIds.push_back (iCluster);     
//   }
}

void clusterInfo::writeStep (std::ostream& out)
{
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      out << clusterIds [iCluster] << "\t";
      out << "(" << clusters [iCluster].size() << ")" << "\t";
      for (int iMol = 0; iMol < clusters [iCluster].size(); iMol++) {
         out << clusters [iCluster].at(iMol) << "\t";
      }
      out << "\n";
   }   
}

void clusterInfo::clear()
{  
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      clusters [iCluster].clear();
   }
   clusters.clear();
   clusterIds.clear();

   nClusters = -1;
}

void clusterInfo::updateClusterId (int prevId, int newId)
{
   // Relies on the fact that before updating the clusterIds, the clusterIds
   // will range from 0 - (nClusters-1). So the prevID will be a number in this
   // range directly equal to the vector index.
   //
   // Updating: clusterIds, whichCluster 
   clusterIds [prevId] = newId;
   for (int iMol = 0; iMol < nMolecules; iMol++) {
      if (whichCluster [iMol] == prevId) {
         whichCluster [iMol] = newId;
      }
   }    

}

void clusterInfo::countLines (std::istream& in)
{
   int count = 0;
   std::string line;

   while (getline (in,line)) {
      count++;
   }

   nClusters = count;
   // Resetting the ifstream to the beginning of the file
   in.clear();
   in.seekg (0, std::ios::beg); 
}


int clusterInfo::maxClusterId = 0;


#endif

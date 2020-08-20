#ifndef CLUSTER_INFO_CPP
#define CLUSTER_INFO_CPP	

#include "mcMd/user/clusterInfo.h"

#include <util/containers/DArray.h>              // member template
#include <util/containers/GArray.h>              // member template
#include <util/containers/ArrayIterator.h>
#include <iostream>
#include <algorithm>
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
      std::cout << "Read file error: iCluster is not equal to nCluster" << std::endl;
   }

   // Orgainze the clusters so as to have clusters[0] have a list of molecules
   // in clusters having aggregation number less than equal to a cutoff
   organize (clustersRead, cutoff_U);

   // Finding the total number of molecules in the simulation
   // also allocating the DArray, whichClusterId 
   for (int iMol = 0; iMol < melt.size(); iMol++) {
      if ( melt [iMol] > nMolecules ) {
         nMolecules = melt [iMol];
      }
   }
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      for (int iMol = 0; iMol < clusters [iCluster].size(); iMol++) {
         if ( clusters [iCluster] [iMol] > nMolecules ) {
            nMolecules = clusters [iCluster] [iMol];
         }
      }
   }


   // Incrementing because the numbering of molecules starts from 0
   nMolecules++;
   whichClusterId.allocate(nMolecules);
   isProcessed.allocate(nMolecules);
 
   // Setting the DArray which cluster.
   // Will be used in mapping
   for (int iMol = 0; iMol < melt.size(); iMol++) {
      whichClusterId [ melt [iMol] ] = 0;
      isProcessed [ melt [iMol] ] = 0;
   }
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      for (int iMol = 0; iMol < clusters [iCluster].size(); iMol++) {
         whichClusterId [ clusters [iCluster] [iMol] ] = clusterIds [iCluster];
         isProcessed [ clusters [iCluster] [iMol] ] = 0;
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
      std::cout << "Read file error: iCluster is not equal to nCluster" << std::endl;
   }

   // Orgainze the clusters so as to have clusters[0] have a list of molecules
   // in clusters having aggregation number less than equal to a cutoff
   organize (clustersRead, cutoff_U);

   // Setting the DArray which cluster.
   // Will be used in mapping
   // Assumes it has already been allocated
   for (int iMol = 0; iMol < melt.size(); iMol++) {
      whichClusterId [ melt [iMol] ] = 0;
      isProcessed [ melt [iMol] ] = 0;
   }
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      for (int iMol = 0; iMol < clusters [iCluster].size(); iMol++) {
         whichClusterId [ clusters [iCluster] [iMol] ] = clusterIds [iCluster];
         isProcessed [ clusters [iCluster] [iMol] ] = 0; 
      }
   }

}

void clusterInfo::organize (std::vector< std::vector<int > > clustersRead, int cutoff_U)
{
   // clusterIds go from 1 - (nClusters)
   // melt has clusterId=0
  
   for (int iCluster = 0; iCluster < clustersRead.size(); iCluster++) {
      if (clustersRead [iCluster].size() <= cutoff_U) {
         melt.insert (melt.end(), clustersRead[iCluster].begin(), 
                                clustersRead[iCluster].end());
      }
      else {
         clusters.push_back(clustersRead[iCluster]);
      }
   }

   // Resetting nClusters after organizing
   nClusters = clusters.size();

   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      clusterIds.push_back (iCluster + 1); 
          
   }
}

int clusterInfo::clusterIndex (int Id) 
{
   int index;

   // Corresponding to melt
   if (Id == 0) {
      index = -1;
   }
   // Corresponding to being in a cluster
   else {
      auto it = std::find (clusterIds.begin(), clusterIds.end(), Id);
      if (it != clusterIds.end()) {
         index = std::distance (clusterIds.begin(), it);
      }
      else {
         std::cout<<"Cluster Id not found in the vector clusterIds"<<std::endl;
         index = -5; 
      }  
   }

   return index;
}

void clusterInfo::writeStep (std::ostream& out)
{
   // Printing the molecules making up the melt
   out << "0" << "  ";
   out << "(" << melt.size() << ")" << "\t";
   for (int iMol = 0; iMol < melt.size(); iMol++) {
      out << melt.at(iMol) << "  ";
   }
   out << "\n";
 
   int max = -1;
   int index = -1;
   int value = -1;
  
   max = * (std::max_element(clusterIds.begin(), clusterIds.end()));
   value = max;

   std::cout<<"***********Writing************"<<std::endl;
//   std::cout<<"max = "<<max<<std::endl; 

   // Will check if the cluster has already been printed or not
   std::vector<bool > print;
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      print.push_back(0);
   }

   // Printing whole clusters in ascending order with respect to 
   // cluster Ids
   for (int iPrint = 0; iPrint < print.size(); iPrint++) {
      for (int iCluster = 0; iCluster < nClusters; iCluster++) {
         if (print [iCluster] == 0) {
  //          std::cout<<"clusterIds ["<<iCluster<<"] = "<<clusterIds [iCluster]<<std::endl;
            if (clusterIds [iCluster] <= value) {
               value = clusterIds [iCluster];
            }
  //          std::cout<<"value = "<<value<<std::endl;
         }
      } 
      index = clusterIndex (value);   
      out << clusterIds [index] << "  ";
      out << "(" << clusters [index].size() << ")" << "\t";
      for (int iMol = 0; iMol < clusters [index].size(); iMol++) {
         out << clusters [index].at(iMol) << "  ";
      }   
      out << "\n";
      print [index] = 1;
      value = max;
   }
}

void clusterInfo::clear()
{  
   for (int iCluster = 0; iCluster < nClusters; iCluster++) {
      clusters [iCluster].clear();
   }
   clusters.clear();
   melt.clear();
   clusterIds.clear();

   nClusters = -1;
}

void clusterInfo::updateClusterId (int prevId, int newId)
{
   // Relies on the fact that before updating the clusterIds, the clusterIds
   // will range from 1 - (nClusters). 
   // So the prevId will be present at index (prevId-1)
   // Updating: clusterIds, whichClusterId
   
   std::cout<<"---------updating ID----------"<<std::endl;
   std::cout<<"prevId In: "<<prevId<<std::endl;
   std::cout<<"newId In: "<<newId<<std::endl;

 
   std::cout<<"Corresponding clusterId: "<<clusterIds [prevId - 1]<<std::endl;
   clusterIds [prevId - 1] = newId;
   for (int iMol = 0; iMol < nMolecules; iMol++) {
      if (whichClusterId [iMol] == prevId && isProcessed [iMol] == 0) {
//         std::cout<<"whichClusterId ["<<iMol<<"] = "<<whichClusterId [iMol]<<std::endl;
         whichClusterId [iMol] = newId;
//         std::cout<<"whichClusterId ["<<iMol<<"] = "<<whichClusterId [iMol]<<std::endl;
         isProcessed [iMol] = 1;
      }
   }    
   std::cout<<"------------------------------"<<std::endl;
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

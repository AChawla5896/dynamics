#ifndef MAPPING_CPP
#define MAPPING_CPP

#include "mcMd/user/mapping.h"
#include "mcMd/user/clusterInfo.h"
#include <util/containers/DArray.h>              // member template

#include <vector>

using namespace Util;

void maxId (clusterInfo* step0, clusterInfo* step1)
{
   if (step0->maxClusterId < step1->nClusters) {
      step0->maxClusterId = step1->nClusters;
   }
}


void reinitializeArray (DArray <double> * in, int size)
{
   for (int i = 0; i < size; i++) {
      (*in) [i]= 0.0;
   }
}

void mapping (clusterInfo* step0, clusterInfo* step1, double cutoff_P, double cutoff_F, 
              std::vector<double> * tally, std::ostream& summary)
{
   /* We will start by mapping the step0 clusters to 
    * step1. In this mapping we will only be able to 
    * identify the fusion events. In this case, a 
    * micelle at step1 will show a contribution <cutoff_P. 
    * This would mean the identity of the micelle was not 
    * "preserved". Having identified such a micelle at 
    * step1 we will find micelles at step0 having a 
    * contribution > cutoff_F. This process helps us 
    * identify the micelles at step0 undergoing fusion 
    * and also in identifying the micelle they form at 
    * step1. In the process of mapping micelles from 
    * step0 to step1 we can also identify the stepwise
    * association process. For this, a micelle at step1
    * should have a contribution > cutoff_P from 
    * cluster 0 (melt) at step0. 
    *
    * Next we will start mapping clusters from step1 to 
    * step0. This will help us identify both fission and
    * stepwise dissociation. A similar procedure/reasoning
    * will follow.
    *
    * The third and the final step pertains to identifying 
    * the chain insertion/expulsion events. (Hopefully) As
    * we have now covered all the processes which can result
    * in birth/death of a micelle, the number of micelles
    * not analyzed till now will be the same for both the 
    * steps. So now we will only consider those micelles 
    * which have not been "processed". (Hopefully), we 
    * should find only one micelle at step1 which maps to
    * step0. The condition for mapping that will be used 
    * is that the contribution needs to be > cutoff_P.
    * This would imply that the identity of micelle is 
    * preserved. Checking this condition from 
    * step0 -> step1 is necessary. However, this condition
    * will also be verfied from step1 -> step0 just to be
    * sure. Then by using the whichClusterId Array we will 
    * be able to identify the chain insertion/expulsion
    * events.
    */

   // The following lines declare workArrays to find 
   // the contributions. It is important to remember
   // that the size of these arrays is equal to 
   // nClusters + 1. This is because nClusters only 
   // account for proper aggregates and not the melt.
   // Index 0 of these arrays would represent the melt.
   // Rest of the indices would represent proper 
   // aggregates. To get their clusterId you would have
   // to call clusterId ("DArray index" - 1). If "DArray 
   // index" is 0, that means the melt
   
   // This is a workArray to analyze the contribution
   // of clusters at step0 to step1
   DArray < double> contribution0;
   // This is a workArray to analyze the contribution
   // of clusters at step1 to step0
   DArray < double> contribution1;

   // These vectors would keep track of the clusters which 
   // have been mapped at the respective time step. The 
   // sizes of these vectors will be stepi->nClusters.
   std::vector<bool > processed0;
   std::vector<bool > processed1;
 
   // These vectors will keep track of the clusters which are 
   // preserved and which are not
   std::vector<bool > preserved0;
   std::vector<bool > preserved1;

   // The maxContributioni vector would store the maximum 
   // contribution (ratio) a cluster has from the other step.
   // The index vector would store the DArray index for 
   // which this maxContribution occured.
   // Size: stepi->nClusters
   std::vector<double > maxContribution0;
   std::vector<double > maxContribution1;
   std::vector<int > index0;
   std::vector<int > index1;

   // Size of the respective DArrays and allocating them.
   int nClusters0 = (step0->nClusters);
   int nClusters1 = (step1->nClusters);
   
   // contribution1 will be used when we map from 
   // step0 -> step1. As it represents the clusters at
   // step0 its size will be determined from nClusters 
   // at step0.
   contribution0.allocate (nClusters1 + 1);
   contribution1.allocate (nClusters0 + 1);

   reinitializeArray (&contribution0, (nClusters1 + 1));
   reinitializeArray (&contribution1, (nClusters0 + 1));

   // Initializing the vectors processed maxContribution 
   // and index.
   for (int i = 0; i < nClusters0; i++) {
      processed0.push_back (0); 
      preserved0.push_back (1);
      maxContribution0.push_back (-1.0);
      index0.push_back (-1);  
   }
   for (int i = 0; i < nClusters1; i++) {
      processed1.push_back (0); 
      preserved1.push_back(1);
      maxContribution1.push_back (-1.0);
      index1.push_back (-1);  
   } 

   // These variables would keep track of maximum and the 
   // sum of the respective contribution arrays.
   // Required only in the loops.
   int max0 = -1;
   int max1 = -1;
   int sum0 = 0;
   int sum1 = 0;

   // The selected molecule from the selected cluster
   int selectMol = -1;
   // Cluster Id of the selected molecule at the other step
   int ClusId = -5;
   // Position in the DArray
   int pos = -1;
   
   // This variable will be used to count the number of 
   // clusters involved in the fusion/fission process
   int nDynamic = 0;

   // FIRST STEP: Identify fusion and stepwise association
   //
   for (int iCluster = 0; iCluster < nClusters1; iCluster++) {
   //   std::cout<<"iCluster :"<<iCluster<<std::endl;
      for (int iMol = 0; iMol < (step1->clusters) [iCluster].size(); iMol++) {
         selectMol = (step1->clusters) [iCluster] [iMol];
         ClusId = step0->whichClusterId[selectMol];
         // Added one because position 0 represents melt over here
         // and clusterIndex would return -1 for melt
         pos = step0->clusterIndex(ClusId) + 1;
         contribution1 [pos]++;
         sum1++;
         if (contribution1 [pos] > max1) {
            max1 = contribution1 [pos];
            index1 [iCluster] = pos;
         } 
      }
      maxContribution1 [iCluster] = (double) max1/ (double) sum1;

      // Identifying fusion events
      if (maxContribution1 [iCluster] < cutoff_P) {
         summary<<"Fusion : ";
         for (int i = 0; i < nClusters0 + 1; i++) {
            contribution1 [i] = contribution1[i]/(double) sum1;
            if (contribution1 [i] >= cutoff_F) {
     //          std::cout<<"A cluster satisfied fusion conditions"<<std::endl;
               if (i != 0) {
                  // These are the micelles which are fusing to form
                  // the micelle at step1.
                  processed0 [i - 1] = 1;
                  preserved0 [i - 1] = 0;
                  nDynamic++;
                  summary<<step0->clusterIds[i-1]<<" ";
     //             std::cout<<"contribution1 :"<<contribution1 [i]<<std::endl;
     //             std::cout<<"Cluster id: "<<step0->clusterIds[i-1]<<std::endl;
     //             std::cout<<"nDynamic: "<<nDynamic<<std::endl;
               }
               else {
                  std::cout << "Algorithmic error: Micelle fusing with melt"
                                << std::endl;
               }
            }
     //       std::cout<<"contribution1 ["<<i<<"]"<<contribution1 [i]<<std::endl;
         }
         if (nDynamic >= 2) {
            step0->maxClusterId++;
      //      std::cout<<"Max Cluster ID :"<<step0->maxClusterId<<std::endl;
      //      std::cout<<"Final nDynamic: "<<nDynamic<<std::endl;
      //      std::cout<<"Prev Id: "<<step1->clusterIds [iCluster]<<std::endl;
      //      std::cout<<"New Id: "<<step0->maxClusterId<<std::endl;
            step1->updateClusterId(step1->clusterIds [iCluster], step0->maxClusterId);
            (* tally) [5]++; 
            processed1 [iCluster] = 1;
            preserved1 [iCluster] = 0;
         }
         else {
            std::cout << "Algorithmic error: Less than two micelles fusing"
                          << std::endl;
         }
         summary<<"= "<<step0->maxClusterId<<std::endl<<std::endl;
      }
      else {
         // Micelle forming from melt
         if (index1 [iCluster] == 0) {
     //       std::cout<<"A cluster forming from the melt"<<std::endl;
            summary<<"Stepwise association : ";
            step0->maxClusterId++;
     //       std::cout<<"Max Cluster ID :"<<step0->maxClusterId<<std::endl;
     //       std::cout<<"Prev Id: "<<step1->clusterIds [iCluster]<<std::endl;
    //        std::cout<<"New Id: "<<step0->maxClusterId<<std::endl;
   //         std::cout<<"iCluster : "<<iCluster<<std::endl;
    //        std::cout<<"Max contribution :"<<maxContribution1 [iCluster]<<std::endl;
            step1->updateClusterId(step1->clusterIds [iCluster], step0->maxClusterId);
            summary<<step0->maxClusterId<<std::endl<<std::endl;
            (* tally) [2]++; 
            processed1 [iCluster] = 1; 
            preserved1 [iCluster] = 0;
         }
      }

      reinitializeArray (&contribution1, (nClusters0 + 1));
      max1 = -1;
      sum1 = 0;
      nDynamic = 0;
   }   

   // SECOND STEP: Identify fission and stepwise dissociation
   //
   for (int iCluster = 0; iCluster < nClusters0; iCluster++) {
      if (processed0 [iCluster] == 0) {
      //   std::cout<<"iCluster :"<<iCluster<<std::endl;
         for (int iMol = 0; iMol < (step0->clusters) [iCluster].size(); iMol++) {
            selectMol = (step0->clusters) [iCluster] [iMol];
            ClusId = step1->whichClusterId[selectMol];
            // Added one because position 0 represents melt over here
            // and clusterIndex would return -1 for melt
            pos = step1->clusterIndex(ClusId) + 1;
            contribution0 [pos]++;
            //std::cout<<"Here...iMol = "<<iMol<<std::endl;
            //std::cout<<"step1->whichClusterId["<<selectMol<<"] ="<< step1->whichClusterId[selectMol] <<std::endl;
            sum0++;
            if (contribution0 [pos] > max0) {
               max0 = contribution0 [pos];
               index0 [iCluster] = pos;
        //       std::cout<<"pos = "<<pos<<std::endl;
            }
         }
         maxContribution0 [iCluster] = (double) max0/ (double) sum0;
         // Identifying fission events
         if (maxContribution0 [iCluster] < cutoff_P) {
            summary<<"Fission : "<< step0->clusterIds [iCluster] <<" = ";
            for (int i = 0; i < nClusters1 + 1; i++) {
               contribution0 [i] = (double) contribution0[i]/(double) sum0;
               if (contribution0 [i] >= cutoff_F) {
          //        std::cout<<"A cluster satisfied fission conditions"<<std::endl;
                  // These are the micelles which are splitting to form
                  // the micelles at step 1
                  if (i != 0 ) {
                     processed1 [i - 1] = 1;
                     preserved1 [i - 1] = 0;
                     step0->maxClusterId++;
          //           std::cout<<"contribution0 :"<<contribution0 [i]<<std::endl;
          //           std::cout<<"Max Cluster ID :"<<step0->maxClusterId<<std::endl;
          //           std::cout<<"Cluster number: "<<i<<std::endl;
          //           std::cout<<"Prev Id: "<<step1->clusterIds [i-1]<<std::endl;
          //           std::cout<<"New Id: "<<step0->maxClusterId<<std::endl;
                     step1->updateClusterId(step1->clusterIds [i-1], step0->maxClusterId);
                     summary<<step0->maxClusterId<<" ";
                     nDynamic++;
          //           std::cout<<"nDynamic: "<<nDynamic<<std::endl;
                  }
                  else {
                     std::cout <<"Algorithmic error: Micelle splitting to melt"
                                   << std::endl;
                  }
               }   
            }
            summary<<std::endl;   
            if (nDynamic >= 2) {
               processed0 [iCluster] = 1;
               preserved0 [iCluster] = 0;
               (* tally) [4]++;
            }   
            else {
               std::cout << "Algorithmic error: Splitting to less than two micelles"
                             << std::endl;
            }
            summary<<std::endl;
         }
         else {
            // Micelle dissociating into melt
            if (index0 [iCluster] == 0) {
               (* tally) [3]++;
               processed0 [iCluster] = 1;
               preserved0 [iCluster] = 0;
               summary<<"Stepwise dissociation : "<<step0->clusterIds[iCluster]
                         <<std::endl<<std::endl;
       //        std::cout<<"Micelle dissociating into the melt"<<std::endl;
       //        std::cout<<"clusterId :"<<step0->clusterIds[iCluster]<<std::endl;
       //        std::cout<<"Max Contribution :"<<maxContribution0 [iCluster]<<std::endl ;
            }
         }

         reinitializeArray (&contribution0, (nClusters1 + 1));
         max0 = -1;
         sum0 = 0;
         nDynamic = 0;
      }
   }

   // THIRD STEP: Identifying the chain insertion/expulsion events
   //
   int countUnprocessed0 = 0;
   int countUnprocessed1 = 0;
 
   // These variables will be used when updating the clusterId at
   // step1
   int oldId = -1;
   int newId = -1;
   int thisId = -1;
 
   for (int iCluster = 0; iCluster < nClusters0; iCluster++) {
      if (processed0 [iCluster] == 0) {
         countUnprocessed0++;
      }
   }
   for (int iCluster = 0; iCluster < nClusters1; iCluster++) {
      if (processed1 [iCluster] == 0) {
         countUnprocessed1++;
      }   
   }

   if (countUnprocessed0 != countUnprocessed1) {
      std::cout<<"Algorithmic error: unequal unprocessed clusters"
                  <<" after processing birth and death"<<std::endl;
      std::cout<<"Step0 has "<<countUnprocessed0<<std::endl;
      std::cout<<"Step1 has "<<countUnprocessed1<<std::endl;
   }
   else {
      //std::cout<<"Step0 has "<<countUnprocessed0<<std::endl;
      //std::cout<<"Step1 has "<<countUnprocessed1<<std::endl;

      for (int iCluster = 0; iCluster < nClusters1; iCluster++) {
         //std::cout << "iCluster = "<< iCluster<<std::endl;
         //std::cout << "processed1 ["<<iCluster<<"] = "<<processed1 [iCluster]<<std::endl;
         //std::cout << "maxContribution1 ["<<iCluster<<"] = "<<maxContribution1 [iCluster]<<std::endl;
         //std::cout << "index1 = " << (index1 [iCluster]) <<std::endl;
         if (processed1 [iCluster] == 0) {
            if (maxContribution1 [iCluster] >= cutoff_P) {
               oldId = step1->clusterIds[iCluster];
               newId = step0->clusterIds[(index1 [iCluster] - 1)];
         //      temp = step0->clusterIndex(newId);
         //      std::cout<<"oldId_map =" <<oldId<<std::endl;
         //      std::cout<<"newId_map =" <<newId<<std::endl;
         //      std::cout<<"maxContribution0 ["<<temp<<"]"<<maxContribution0[temp]<<std::endl;
         //      std::cout<<"index0 =" << index0[temp]<<std::endl;
               if (index0 [step0->clusterIndex(newId)] == (iCluster+1)) {
                  step1->updateClusterId(oldId, newId);
                  processed1 [iCluster] = 1;
                  processed0 [index1 [iCluster] - 1] = 1;
               }
               else {
                  std::cout<<"Algorithmic error: Id is not reversibly preserved"
                               <<std::endl;
               }
            }
            else {
               std::cout<<"Algorithmic error: Third step shows unpreserved cluster"
                            <<std::endl;
            }        
         }
      }

 
      for (int iCluster = 0; iCluster < nClusters1; iCluster++) {
         if (preserved1 [iCluster] == 1) {
            thisId = step1->clusterIds[iCluster];
            if (maxContribution1 [iCluster] != 1) {
               summary<<"Chain Insertion : ";
               for (int iMol = 0; iMol < step1->clusters[iCluster].size(); iMol++) {
                  selectMol = step1->clusters[iCluster][iMol];
                  ClusId = step0->whichClusterId[selectMol];
                  if(ClusId != thisId) {
                     summary<<selectMol<<"\t"<<ClusId<<" -> "<<thisId<<std::endl;
                     summary<<"                  ";
                     (* tally) [0]++;
                  }   
               }
               summary<<std::endl;   
            }

            if (maxContribution0 [(index1 [iCluster] - 1)] != 1) {
               summary<<"Chain Expulsion : ";
               for (int iMol = 0; iMol < step0->clusters[(index1 [iCluster] - 1)].size(); iMol++) {
                  selectMol = step0->clusters[(index1 [iCluster] - 1)][iMol];
                  ClusId = step1->whichClusterId[selectMol];
                  if (ClusId != thisId) {
                     summary<<selectMol<<"\t"<<thisId<<" -> "<<ClusId<<std::endl;
                     summary<<"                  ";
                     (* tally) [1]++;
                  }
               }
               summary<<std::endl;
            }
         }
      }      
 
   }


   // Ensure that count unprocessed is zero for both

   contribution0.deallocate();
   contribution1.deallocate();

}

#endif

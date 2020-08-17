#ifndef MAPPING_CPP
#define MAPPING_CPP

#include "mcMd/user/mapping.h"
#include "mcMd/user/clusterInfo.h"
#include <util/containers/DArray.h>              // member template

#include <vector>

using namespace Util;

void reinitializeArray (DArray <double> * in, int size)
{
   for (int i = 0; i < size; i++) {
      (*in) [i]= 0.0;
   }
}

void mapping (clusterInfo* step0, clusterInfo* step1, double cutoff_P, double cutoff_F, 
              std::vector<int> * tally, std::ostream& out)
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
    * sure. Then by using the whichCluster Array we will 
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
      maxContribution0.push_back (-1.0);
      index0.push_back (-1);  
   }
   for (int i = 0; i < nClusters1; i++) {
      processed1.push_back (0); 
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
   
   // This variable will be used to count the number of 
   // clusters involved in the fusion/fission process
   int nDynamic = 0;

   // FIRST STEP: Identify fusion and stepwise association
   //
   for (int iCluster = 0; iCluster < nClusters1; iCluster++) {
      for (int iMol = 0; iMol < (step1->clusters) [iCluster].size(); iMol++) {
         selectMol = (step1->clusters) [iCluster] [iMol];
         contribution1 [step0->whichCluster[selectMol]]++;
         sum1++;
         if (contribution1 [step0->whichCluster[selectMol]] > max1) {
            max1 = contribution1 [step0->whichCluster[selectMol]];
            index1 [iCluster] = step0->whichCluster[selectMol];
         } 
      }
      maxContribution1 [iCluster] = (double) max1/ (double) sum1;

      // Identifying fusion events
      if (maxContribution1 [iCluster] < cutoff_P) {
         for (int i = 0; i < nClusters0 + 1; i++) {
            contribution1 [i] = contribution1[i]/(double) sum1;
            if (contribution1 [i] >= cutoff_F) {
               if (i != 0) {
                  // These are the micelles which are fusing to form
                  // the micelle at step1.
                  processed0 [i - 1] = 1;
                  nDynamic++;
               }
               else {
                  std::cout << "Algorithmic error: Micelle fusing with melt"
                                << std::endl;
               }
            }
         }
         if (nDynamic >= 2) {
            step0->maxClusterId++;
            step1->updateClusterId(step1->clusterIds [iCluster], step0->maxClusterId);
            (* tally) [5]++; 
            processed1 [iCluster] = 1;
         }
         else {
            std::cout << "Algorithmic error: Less than two micelles fusing"
                          << std::endl;
         }
      }
      else {
         // Micelle forming from melt
         if (index1 [iCluster] == 0) {
            step0->maxClusterId++;
            step1->updateClusterId(step1->clusterIds [iCluster], step0->maxClusterId);
            (* tally) [2]++; 
            processed1 [iCluster] = 1; 
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
         for (int iMol = 0; iMol < (step0->clusters) [iCluster].size(); iMol++) {
            selectMol = (step0->clusters) [iCluster] [iMol];
            contribution0 [step1->whichCluster[selectMol]]++;
            sum0++;
            if (contribution0 [step1->whichCluster[selectMol]] > max0) {
               max0 = contribution0 [step1->whichCluster[selectMol]];
               index0 [iCluster] = step1->whichCluster[selectMol];
            }
         }
         maxContribution0 [iCluster] = (double) max0/ (double) sum0;

         // Identifying fission events
         if (maxContribution0 [iCluster] < cutoff_P) {
            for (int i = 0; i < nClusters1 + 1; i++) {
               contribution0 [i] = contribution0[i]/(double) sum0;
               if (contribution0 [i] >= cutoff_F) {
                  // These are the micelles which are splitting to form
                  // the micelles at step 1
                  if (i != 0 ) {
                     processed1 [i - 1] = 1;
                     step0->maxClusterId++;
                     step1->updateClusterId(step1->clusterIds [i-1], step0->maxClusterId);
                     (* tally) [4]++;
                     nDynamic++;
                  }
                  else {
                     std::cout <<"Algorithmic error: Micelle splitting to melt"
                                   << std::endl;
                  }
               }   
            }   
            if (nDynamic >= 2) {
               processed0 [iCluster] = 1;
            }   
            else {
               std::cout << "Algorithmic error: Splitting to less than two micelles"
                             << std::endl;
            }
         }
         else {
            // Micelle dissociating into melt
            if (index0 [iCluster] == 0) {
               (* tally) [3]++;
               processed0 [iCluster] = 1;
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







   // Chain insertion events at step1
   
   // Reset max and sum at the end of this loop
   for (int iCluster = 0; iCluster < (step1->nClusters); iCluster++) {
      for (int iMol = 0; iMol < (step1->clusters) [iCluster].size(); iMol++) {
         selectMol = (step1->clusters) [iCluster] [iMol];
         contribution0[step0->whichCluster[selectMol]]++;
         sum0++;
         if (contribution0[step0->whichCluster[selectMol]] > max0) {
            max0 = contribution0[step0->whichCluster[selectMol]];
            index0 = step0->whichCluster[selectMol];
         }
      }
      if (((double)max0/(double)sum0) >= cutoff_P) {

         /* Need to have a check over here to see if it is 
          * not a fission event. In a fission event, micelle 
          * at step0 will split and form two or more micelles 
          * at step1. Both of these micelles at step1 will 
          * show the same step0 micelle having the dominant 
          * contribution (and hence clearing the cutoff). 
          * However, such a problem would not arise in a fusion
          * event. In a fusion event, two or more micelles at 
          * step0 would make up a single micelle at step1. This
          * would result in that micelle not satisfying the 
          * above condition on contribution. The micelles at 
          * step0 would have approximately equal contributions
          * to the micelle at step1 (e.g. a contribution of 
          * about 50% is expected, from each of the micelles at
          * step0, if two micelles at step0 form a micelle at 
          * step1).
          */ 


         if (index0 != 0) {
            // index0 - 1 because DArray index 0 represents the melt
            // Therefore the clusters start from DArray index 1
            step1->updateClusterId(step1->clusterIds [iCluster], step0->clusterIds[(index0-1)]);
            out << "Chain Insertion :" << "\t" ;
            for (int i = 0; i < size0; i++) {
               if (i!=index0) {
                  out << "(" << contribution0[i] << ")" << "   " << i << "->" <<index0;
                  out << std::endl;
                  out << "                 " << "\t"; 
               }
            }
            out << std::endl;
         }
         else {
            // Will have to max Cluster ID
            // Identify the stepwise association process
            step0->maxClusterId++;
            step1->updateClusterId(step1->clusterIds [iCluster], step0->maxClusterId);
            out << "Unimers associated to form micelle " << step0->maxClusterId << std::endl;
         }
      }
      else {
         // Put code here to identify fusion
         // Idnetify fission using a different technique
         // Will have to update maxClusterId over here

      }
      reinitializeArray (&contribution0, size0);
      max0 = -1;
      sum0 = 0;
   }


}

#endif

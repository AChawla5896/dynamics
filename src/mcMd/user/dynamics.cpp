#include "mcMd/user/clusterInfo.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main (){

   std::ifstream commands ("commands");
   /* Format of commands file:
   *  Input Prefix: IPREFIX
   *  Output Prefix: OPREFIX
   *  Start: int (t0)
   *  Increment: int (del_T)
   *  Final: int (tf)
   *  cutoff unimers: int (cutoff to allocate the cluster to the set of unimers)
   *  cutoff combine: double (cutoff to identify fision and fusion processes)
   */
 
   std::string l1 ("Input Prefix"); 
   std::string IPre;
   std::string l2 ("Output Prefix");
   std::string OPre;
   std::string l3 ("Start");
   int T0;
   std::string l4 ("Increment");
   int delT;
   std::string l5 ("Final");
   int Tf;
   std::string l6 ("Cutoff unimers");
   int cutoff_U;
   std::string l7 ("Cutoff combine");
   double cutoff_C;

   // Reading commands file
   std::string line;
   std::string head;
   std::string ParaInput;
   
   while ( std::getline( commands, line ) ) {

      std::stringstream linestream(line);
      std::getline (linestream, head, ':');
      if (head == l1){
         std::getline(linestream, IPre);
         IPre.erase(0,1);
         std::cout<<l1<<": "<<IPre<<std::endl;
      }
      else if (head == l2) {
         std::getline(linestream, OPre);
         OPre.erase(0,1);
         std::cout<<l2<<": "<<OPre<<std::endl;
      } 
      else if (head == l3) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         T0 = stoi (ParaInput);
         std::cout<<l3<<": "<<T0<<std::endl;
      }
      else if (head == l4) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         delT = stoi (ParaInput);
         std::cout<<l4<<": "<<delT<<std::endl;
      }
      else if (head == l5) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         Tf = stoi (ParaInput);
         std::cout<<l5<<": "<<Tf<<std::endl;
      }  
      else if (head == l6) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         cutoff_U = stoi (ParaInput);
         std::cout<<l6<<": "<<cutoff_U<<std::endl;
      }  
      else if (head == l7) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         cutoff_C = stod (ParaInput);
         std::cout<<l7<<": "<<cutoff_C<<std::endl;
      }  
      else {
         std::cout << "Wrong input format"<<std::endl;
      }

   }

   // Setting the file streams and cluster variables
   std::ifstream inFileName0 (IPre+to_string(T0));
   std::ifstream inFileName1 (IPre+to_string(T0+delT));
   std::ofstream outFileName0 (OPre+to_string(T0));
   std::ofstream outFileName1 (OPre+to_string(T0+delT));
   clusterInfo step0, step1, temp;
 
   /* Read the first two time steps using readStep0 so as to allocate 
   * the DArrays
   */
   step0.readStep0(inFileName0);
   step1.readStep0(inFileName1);

   // Update maxClusterId. This is a static variable.
   // After this it needs to be updated in the mapping function
   step0.maxClusterId = step0.nClusters;

   // Add a MAPPING command over here
   
   // Writing these timesteps to new output files
   step0.writeStep(outFileName0);
   step1.writeStep(outFileName1);
   
   // Clearing the input and output streams for step0. Reinitialize 
   // using the open function of these streams.
   inFileName0.close(); 
   outFileName0.close();

   // iFile is the string that needs to be appended to the filename
   // iSample is the number of samples that have been analyzed
   // 2 files have already been mapped using readStep0
   int iSample = 2;

   for (int iFile = (T0+(2*delT)); iFile <= Tf; iFile+=delT) {

      // Clearing the IO variables for step 1
      inFileName1.close();
      outFileName1.close();

      // Copying step1 to step0 after clearing step0.
      step0.clear();
      step0 = step1;
      step1.clear();

      // Reinitializing the IO variable
      inFileName1.open(IPre+to_string(iFile));
      outFileName1.open(OPre+to_string(iFile));
 
      // Reading the next cluster file
      step1.readStep(inFileName1);

      // MAPPING from step0 to step1

      // Writing the step1 file 
      step1.writeStep(outFileName1);

      iSample++;

   }

   std::cout<<"Total samples analyzed: " << iSample << std::endl;
   
   std::cout<<"Is this it....?"<< std::endl;


   // Do the initial allocation of maxClusterId using the nClusters after
   // calling readStep0 for step0
   // Start incrementing it in mapping thereafter

   // Take the initial timestep, final timestep and the increment as an input

   // Then have a for loop which reads all the files one by one, maps the 
   // clusters and then writes it to a seperate folder. You might also need 
   // to have the path of this folder as an input. 

   // need to have a function in this file to map the clusters
 
   // might also need to write a function to switch step0 and step1 as 
   // when you simply copy, the inequality between the sizes of the two 
   // GArrays/ DArrays might turn up as an error.

   return 0;
}

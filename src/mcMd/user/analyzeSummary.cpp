#include "mcMd/user/readSample.h"
#include <util/containers/RingBuffer.h> //member

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>


using namespace Util;
using namespace std;


void addNext (int Nh, RingBuffer <readSample> *histories, readSample temp) 
{
   (*histories) [Nh-1].clear();
   (*histories).append (temp);
}

void check (int Nh, RingBuffer <readSample> *histories)
{ 
   // Finding reverted fusion in the history
   for (int analyze = 2; analyze <= Nh; analyze++) {
      for (int m = 0; m < (*histories) [Nh-1].fusion.size(); m++) {
         if ((*histories) [Nh-1].fusion [m] [3] == 0) {
            for (int n = 0; n < (*histories) [Nh-analyze].fission.size(); n++) {
               // Checking if the formed micelle in fusion is being used in a fission 
               // at some timestep
               //std::cout<<"histories ["<<Nh-1<<"].fusion ["<<m<<"][2] ="<<histories [Nh-1].fusion [m][2]<<std::endl;
               if ((*histories) [Nh-analyze].fission[n][0] == (*histories) [Nh-1].fusion [m][2]) {
                  if ((*histories) [Nh-analyze].fission[n][3] == 0) {
                     (*histories) [Nh-analyze].fission[n][3] = 1;
                     (*histories) [Nh-1].fusion [m][3] = 1;
                  }   
                  else {
                     std::cout <<"Algorithmic error: The selected history has already been reverted"
                                << std::endl;
                  }
               }
            }
         }
      }
   }

   // Finding reverted fission in the history
   for (int analyze = 2; analyze <= Nh; analyze++) {
      for (int m = 0; m < (*histories) [Nh-1].fission.size(); m++) {
         if ((*histories) [Nh-1].fission [m] [3] == 0) {
            for (int n = 0; n < (*histories) [Nh-analyze].fusion.size(); n++) {
               // Checking if the formed micelles in fission, fuse together in a fusion 
               // at some timestep
               if ((*histories) [Nh-analyze].fusion[n][0] == (*histories) [Nh-1].fission [m][1]
                         || (*histories) [Nh-analyze].fusion[n][1] == (*histories) [Nh-1].fission [m][1]) {
                  if ((*histories) [Nh-analyze].fusion[n][0] == (*histories) [Nh-1].fission [m][2]
                            || (*histories) [Nh-analyze].fusion[n][1] == (*histories) [Nh-1].fission [m][2]) {
                     if ((*histories) [Nh-analyze].fusion[n][3] == 0) {
                        (*histories) [Nh-analyze].fusion[n][3] = 1;
                        (*histories) [Nh-1].fission [m][3] = 1;
                     }
                     else {
                        std::cout <<"Algorithmic error: The selected history has already been reverted"
                                    << std::endl;
                     }
                  }
               }
            }
         }
      }
   }
}

void print (RingBuffer <readSample> *histories, int Nh, std::ostream& out, std::vector<double> * tally)
{
   out << "SAMPLE "<< ": " << (*histories) [Nh-1].iSample << std::endl;
   for (int m = 0; m < (*histories) [Nh-1].fission.size(); m++) {
      if ((*histories) [Nh-1].fission [m] [3] == 0) {
         out << "Fission "<< ": ";
         out << (*histories) [Nh-1].fission [m] [0] << " = ";
         out << (*histories) [Nh-1].fission [m] [1] <<" "<<(*histories) [Nh-1].fission [m] [2] 
                << std::endl;
         (*tally) [4]++;
      }   
   }   
    
   for (int n = 0; n < (*histories) [Nh-1].fusion.size(); n++) {
      if ((*histories) [Nh-1].fusion [n] [3] == 0) {
         out << "Fusion "<< ": ";
         out << (*histories) [Nh-1].fusion [n] [0] <<" "<<(*histories) [Nh-1].fusion [n] [1]; 
         out << " = " << (*histories) [Nh-1].fusion [n] [2] 
                 << std::endl;
         (*tally) [5]++;
      }   
   }   
   out << std::endl << std::endl;
}

int main (){

   std::ifstream commands ("commands");
   /* Format of commands file:
   *  Input Prefix: IPREFIX
   *  Output Prefix: OPREFIX
   *  Start: int (t0)
   *  Increment: int (del_T)
   *  Final: int (tf)
   *  Cutoff unimers: int (cutoff to allocate the cluster to the set of unimers)
   *  Cutoff for preserved micelle: double (cutoff to identify if the micelle 
   *                                        preserved its identity or not)
   *  Cutoff for fission/fusion: double (cutoff to identify which of the micelles 
   *                                     fused/split)                              
   *  Number of histories: int (number of histories to use while analyzing if 
   *                            there are repeated fission or fusion of the 
   *                            same cluster)  
   */

   // summary file outputs the all the dynamic processes taking place from the 
   // previous time step -> next time step
 
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
   std::string l7 ("Cutoff for preserved micelle");
   double cutoff_P;
   std::string l8 ("Cutoff for fission/fusion");
   double cutoff_F;
   std::string l9 ("Number of histories");
   int Nh;

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
         cutoff_P = stod (ParaInput);
         std::cout<<l7<<": "<<cutoff_P<<std::endl;
      }
      else if (head == l8) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         cutoff_F = stod (ParaInput);
         std::cout<<l8<<": "<<cutoff_F<<std::endl;
      }
      else if (head == l9) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         Nh = stoi (ParaInput);
         std::cout<<l9<<": "<<Nh<<std::endl;
      }  
      else {
         std::cout << "Wrong input format"<<std::endl;
      }
   }

   std::cout<<"======================================"<<std::endl;

   // RingBuffer to store all the histories
   RingBuffer <readSample > histories;
   histories.allocate(Nh);
   // Temporary variable to store data. Will be appended to 
   // histories.
   readSample temp;

   std::ifstream summary ("summary");
   std::ofstream out ("analyzedSummary");

   /* This vector will keep track of all the dynamic processes
    * taking place in the micelles
    * Index 0: Total chain insertion events
    * Index 1: Total chain expulsion events
    * Index 2: Total stepwise association events
    * Index 3: Total stepwise dissociation events
    * Index 4: Total fission events
    * Index 5: Total fusion events
    */
   std::vector<double > tally;
   // Initializing this vector to 0
   for (int i = 0; i < 6; i++) {
      tally.push_back (0);
   }  

   // header strings in summary file
   std::string newSample ("SAMPLE ");
   std::string fission ("Fission ");
   std::string fusion ("Fusion ");
   std::string chainInsert ("Chain Insertion ");
   std::string chainExpul ("Chain Expulsion ");
   std::string associate ("Stepwise association ");
   std::string dissociate ("Stepwise dissociation ");
   std::string emptySpace ("                  ");

   std::vector <int > readVec;
   std::string clusInput;
   int countEmpty = 0;

   while ( getline (summary, line ) ) {
      std::istringstream summaryStream(line);
      std::getline (summaryStream, head, ':'); 
      if (head == newSample){
         std::getline(summaryStream, ParaInput);
         ParaInput.erase(0,1);
         temp.iSample = stoi (ParaInput);
         std::cout<<newSample<<": "<<temp.iSample<<std::endl;
         countEmpty = 0;
         while ( getline (summary, line ) ) {
            std::istringstream summaryStream(line);
            std::getline (summaryStream, head, ':');

            if (head == fusion){
               std::getline (summaryStream, clusInput, '=');
               std::istringstream vecStream (clusInput);
               readVec.insert (readVec.begin(), std::istream_iterator<int>(vecStream), 
                         std::istream_iterator<int>() );
               std::getline (summaryStream, clusInput);
               clusInput.erase(0,1);  
               readVec.insert (readVec.end(), stoi(clusInput));
               readVec.insert (readVec.end(), 0);
               temp.fusion.push_back (readVec);
               readVec.clear();
               countEmpty = 0;
            }
            else if (head == fission) {
               std::getline (summaryStream, clusInput, '=');
               clusInput.erase(0,1);  
               readVec.insert (readVec.begin(), stoi(clusInput));
               std::getline (summaryStream, clusInput);
               std::istringstream vecStream (clusInput);
               readVec.insert (readVec.end(), std::istream_iterator<int>(vecStream), 
                         std::istream_iterator<int>() );
               readVec.insert (readVec.end(), 0); 
               temp.fission.push_back (readVec);
               readVec.clear();
               countEmpty = 0;
            }
            else if (head == chainInsert) {
               countEmpty = 0;
               break;
            }
            else if (head == chainExpul) {
               countEmpty = 0;
               break;
            }
            else if (head == associate) {
               countEmpty = 0;
               break;
            } 
            else if (head == dissociate) {
               countEmpty = 0;
               break;
            } 
            else if (head == emptySpace || head.empty()) {
               countEmpty++;
            }
            else if (countEmpty == 3) {
               break;
            }   
         }

         //for (int i = 0; i < temp.fusion.size(); i++) { 
         //   for (int j = 0; j < temp.fusion [i].size(); j++) {
         //      std::cout<<temp.fusion [i] [j] <<"  ";
         //   }  
         //   std::cout<<endl;
         //}

         //for (int i = 0; i < temp.fission.size(); i++) {
         //   for (int j = 0; j < temp.fission [i].size(); j++) {
         //      std::cout<<temp.fission [i] [j] <<"  ";
         //   }
         //   std::cout<<endl;
         //}

         if(histories.isFull()) {
            check (Nh, & histories);
            print (& histories, Nh, out, & tally);       
            addNext (Nh, & histories, temp);
         }
         else {
            histories.append(temp);
         }

         temp.clear();

      }

   } 

   check (Nh, & histories);
   print (& histories, Nh, out, & tally);

   out<<"============================================="<<std::endl<<std::endl;
   out<<"Total samples analyzed: " << histories [Nh-1].iSample << std::endl;

   for (int i = 0; i < tally.size(); i++) {
      tally [i] = tally[i]/((double)histories [Nh-1].iSample);
   }
   out<<"Chain insertion events per step : "<< tally [0] << setprecision(8) << std::endl;
   out<<"Chain expulsion events per step : "<< tally [1] << setprecision(8) << std::endl;
   out<<"Stepwise association events per step : "<< tally [2] << setprecision(8) << std::endl;
   out<<"Stepwise dissociation events per step : "<< tally [3] << setprecision(8) << std::endl;
   out<<"Fission events per step : "<< tally [4] << setprecision(8) << std::endl;
   out<<"Fusion events per step : "<< tally [5] << setprecision(8) << std::endl;


}

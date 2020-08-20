#ifndef MAPPING_H
#define MAPPING_H

#include <util/containers/DArray.h> 
#include <vector>

namespace Util {template <typename Data> class DArray;}

class clusterInfo;

void maxId (clusterInfo* step0, clusterInfo* step1);

void reinitializeArray (Util::DArray <double> * in, int size);

void mapping (clusterInfo* step0, clusterInfo* step1, double cutoff_P, double cutoff_F, 
              std::vector<double> * tally, std::ostream& summary);

#endif

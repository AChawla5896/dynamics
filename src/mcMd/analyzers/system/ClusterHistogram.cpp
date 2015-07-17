#ifndef MCMD_CLUSTER_HISTOGRAM_CPP
#define MCMD_CLUSTER_HISTOGRAM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ClusterHistogram.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>

#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>
#include <sstream>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   ClusterHistogram::ClusterHistogram(System& system) 
    : SystemAnalyzer<System>(system),
      identifier_(system),
      hist_(),
      outputFile_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_(),
      histMin_(),
      histMax_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("ClusterHistogram"); }

   /// Read parameters from file, and allocate arrays.
   void ClusterHistogram::readParameters(std::istream& in) 
   {  
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      read<int>(in, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }
      
      read<double>(in, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }
      identifier_.setup(speciesId_, atomTypeId_, cutoff_);

      read<int>(in,"histMin", histMin_);
      read<int>(in,"histMax", histMax_);
      hist_.setParam(histMin_, histMax_);
      hist_.clear();
    
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void ClusterHistogram::loadParameters(Serializable::IArchive& ar)
   {
      // Load (everything but accumulators_)
      Analyzer::loadParameters(ar);

      loadParameter<int>(ar,"speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      loadParameter<int>(ar, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }

      loadParameter<double>(ar, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }
      identifier_.setup(speciesId_, atomTypeId_, cutoff_);

      loadParameter<int>(ar,"histMin", histMin_);
      loadParameter<int>(ar,"histMax", histMax_);
      hist_.setParam(histMin_, histMax_);
      hist_.clear();

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void ClusterHistogram::save(Serializable::OArchive& ar)
   {
      Analyzer::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & histMin_;
      ar & histMax_;
      ar & hist_;
   }

   /*
   * Clear accumulators.
   */
   void ClusterHistogram::setup() 
   {  
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
      hist_.clear();
      nSample_ = 0;
   }

   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void ClusterHistogram::sample(long iStep) 
   { 
      if (isAtInterval(iStep)) {
         identifier_.identifyClusters();
         for (int i = 0; i < identifier_.nCluster(); i++) {
             hist_.sample(identifier_.cluster(i).size());
         }
         ++nSample_;
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ClusterHistogram::output() 
   {
      // Write parameter file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Write histogram output
      fileMaster().openOutputFile(outputFileName(".hist"), outputFile_);
      hist_.output(outputFile_);
      outputFile_.close();
   }

}
#endif 

#ifndef DDMD_TRAJECTORY_WRITER_H
#define DDMD_TRAJECTORY_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Based class for classes that write trajectories to a single file.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class TrajectoryWriter : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      TrajectoryWriter(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~TrajectoryWriter()
      {} 
   
      /**
      * Read parameters and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Clear nSample counter.
      */
      virtual void clear();
  
      /**
      * Dump configuration to file
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);

      /**
      * Close ouput file.
      */
      virtual void output();

   protected:

      /**
      * Write data that should appear once, at beginning of the file. 
      *
      * Called by sample when iStep == 0.
      */
      virtual void writeHeader(std::ostream& out, long iStep) = 0;

      /**
      * Write data that should appear in every frame.
      * 
      * Called by sample on every step.
      */
      virtual void writeFrame(std::ostream& out, long iStep) = 0;

   private:
 
      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;
   
      /// Has readParam been called?
      long isInitialized_;
   
   };

}
#endif 

#ifndef EXCHANGER_TEST_H
#define EXCHANGER_TEST_H

#include <ddMd/communicate/Domain.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/BondStorage.h>
#include <ddMd/configIo/ConfigIo.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/communicate/Exchanger.h>
#include <util/random/Random.h>
#include <util/mpi/MpiLogger.h>

#ifdef  UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

using namespace Util;
using namespace DdMd;

class ExchangerTest: public ParamFileTest<Exchanger>
{
private:

   Boundary boundary;
   Domain domain;
   Buffer buffer;
   AtomStorage atomStorage;
   BondStorage bondStorage;
   ConfigIo configIo;
   Random random;
   int atomCount;

public:

   void setUp()
   {
      std::ifstream configFile;

      // Set connections between atomDistributors
      domain.setBoundary(boundary);
      configIo.associate(domain, boundary, atomStorage, bondStorage, buffer);
      object().associate(domain, boundary, atomStorage, bondStorage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      configIo.setParamCommunicator(communicator());
      random.setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      openFile("in/Exchanger");

      domain.readParam(file());
      atomStorage.readParam(file());
      bondStorage.readParam(file());
      buffer.readParam(file());
      configIo.readParam(file());
      random.readParam(file());

      // Finish reading parameter file
      closeFile();

      configIo.readConfig("in/config");

      object().allocate();
      object().setPairCutoff(0.5);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.

      // Check that all atoms are accounted for after distribution.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (domain.gridRank() == 0) {
         //std::cout << std::endl;
         // std::cout << "Total atom count (post-distribute) = " 
         //          << nAtomAll << std::endl;
         atomCount = nAtomAll;
      }

   }

   void displaceAtoms(double max)
   {
      double min = -max;
      AtomIterator atomIter;
      for(int i = 0; i < 3; i++) {
         atomStorage.begin(atomIter);
         for ( ; !atomIter.atEnd(); ++atomIter) {
            atomIter->position()[i] += random.uniform(min, max);
         }
      }
   }

   virtual void testDistribute()
   { 
      printMethod(TEST_FUNC); 
   }

   void testAtomExchange()
   {
      printMethod(TEST_FUNC);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      // Check that all atoms are within the processor domain.
      AtomIterator  atomIter;
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      MpiLogger logger;

      #if 0
      logger.begin();
      std::cout << "Processor: " << myRank << ", Post-distribute nAtom = "
                << atomStorage.nAtom() << std::endl;
      logger.end();
      #endif

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  false));

      double range = 0.4;
      displaceAtoms(range);
      object().exchangeAtoms();

      // Check that all atoms are accounted for after exchange.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         //std::cout << "Total atom count (post atom exchange) = " 
         //<< nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  false));

      #if 0
      //Print the number of atoms with each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank << " : Post-exchange Atoms count = "
                << atomStorage.nAtom() << std::endl;
      logger.end();
      #endif

   }

   void testGhostExchange()
   {
      printMethod(TEST_FUNC);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator  atomIter;
      GhostIterator ghostIter;

      double range = 0.4;
      displaceAtoms(range);
      
      object().exchange();

      // Check that all atoms are accounted for after ghost exchange.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all local atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  false));

      #if 0
      MpiLogger logger;

      //Print number of atoms on each processor after the ghost exchange.
      logger.begin();
      std::cout << "Processor " << myRank 
                << " : Post-ghost exchange Atom  count = "
                << atomStorage.nAtom() << std::endl;
      logger.end();

      // Print number of ghosts on each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank 
                << " : Post-ghost exchange Ghost count = "
                << atomStorage.nGhost() << std::endl;
      logger.end();
      #endif

   }


   void testGhostUpdate()
   {
      printMethod(TEST_FUNC);

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator   atomIter;
      GhostIterator  ghostIter;
      DArray<Vector> ghostPositions;

      double range = 0.4;
      displaceAtoms(range);
      object().exchange();

      // Record number of atoms and ghosts after exchange
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Update ghost positions
      object().update();

      // Check number of atoms and ghosts on each processor is unchanged.
      TEST_ASSERT(nAtom == atomStorage.nAtom());
      TEST_ASSERT(nGhost == atomStorage.nGhost());

      // Check that all atoms are accounted for after atom and ghost exchanges.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  false));

   }

   void testGhostUpdateCycle()
   {
      printMethod(TEST_FUNC);

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator   atomIter;
      GhostIterator  ghostIter;

      double range = 0.4;
      displaceAtoms(range);

      object().exchange();
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      object().update();
      TEST_ASSERT(nAtom == atomStorage.nAtom());
      TEST_ASSERT(nGhost == atomStorage.nGhost());

      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  false));

      range = 0.1;
      for (int i=0; i < 3; ++i) {
         displaceAtoms(range);

         for (int j=0; j < 3; ++j) {
            object().update();
            TEST_ASSERT(nGhost == atomStorage.nGhost());
            TEST_ASSERT(nAtom == atomStorage.nAtom());
            displaceAtoms(range);
         }

         object().exchange();
         nAtom  = atomStorage.nAtom();
         nGhost = atomStorage.nGhost();

         // Check that all atoms are within the processor domain.
         atomStorage.begin(atomIter);
         for ( ; !atomIter.atEnd(); ++atomIter) {
            TEST_ASSERT(domain.isInDomain(atomIter->position()));
         }

         // Check that all ghosts are outside the processor domain.
         atomStorage.begin(ghostIter);
         for ( ; !ghostIter.atEnd(); ++ghostIter) {
            TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
         }

         TEST_ASSERT(atomStorage.isValid());
         TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(),
                                         false)); 
      }

   }

};

TEST_BEGIN(ExchangerTest)
TEST_ADD(ExchangerTest, testDistribute)
TEST_ADD(ExchangerTest, testAtomExchange)
TEST_ADD(ExchangerTest, testGhostExchange)
TEST_ADD(ExchangerTest, testGhostUpdate)
TEST_ADD(ExchangerTest, testGhostUpdateCycle)
TEST_END(ExchangerTest)

#endif /* EXCHANGER_TEST_H */

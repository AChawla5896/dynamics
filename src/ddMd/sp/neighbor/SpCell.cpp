#ifndef DDMD_SP_CELL_CPP
#define DDMD_SP_CELL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpCell.h"

namespace DdMd
{

   using namespace Util;

   SpCell::SpCell()
    : begin_(0),
      offsetsPtr_(0),
      nextCellPtr_(0),
      nAtom_(0),
      atomCapacity_(0)
   {}

   SpCell::~SpCell()
   {
      delete offsetsPtr_;
   }

   void SpCell::setOffsetArray(SpCell::OffsetArray& offsets)
   {  offsetsPtr_ = &offsets; }

   void SpCell::setNextCell(SpCell& nextCell)
   {  nextCellPtr_ = &nextCell; }

   void SpCell::setLastCell()
   {  nextCellPtr_ = 0; }

   /*
   * Fill an array with pointers to SpCellAtom objects in this cell and neighbors.
   *
   * Upon return, the NeighborArray neighbors contains pointers to all of the
   * atoms this cell and neighboring cells.  The first nAtom() elements are the
   * the atoms in this cell.
   */
   void SpCell::getNeighbors(NeighborArray &neighbors) const
   {
      // Preconditions
      assert(offsetsPtr_);

      SpCellAtom* atom;
      int offset;

      neighbors.clear();

      // add neighbors from this cell
      atom = begin_;
      while (atom < begin_ + nAtom_) {
         neighbors.append(atom);
         atom++;
      }

      for (int i = -1; i <= 1; i++)
         for (int j = -1; j <= 1; j++)
            for (int k = -1; k <= 1; k++) {
               offset = 0;
               // we already counted this cell
               if (!(i == 0 && j == 0 && k == 0)) {
                  offset += ((i == 0) ? 0 : ((i > 0) ? (*offsetsPtr_)[0].first : (*offsetsPtr_)[0].second));
                  offset += ((j == 0) ? 0 : ((j > 0) ? (*offsetsPtr_)[1].first : (*offsetsPtr_)[1].second));
                  offset += ((k == 0) ? 0 : ((k > 0) ? (*offsetsPtr_)[2].first : (*offsetsPtr_)[2].second));

                  if ((this +offset)->id_ >= id_) {
                     atom = (this + offset)->begin_;
                     while (atom < (this + offset)->begin_ + (this + offset)->nAtom_) {
                        neighbors.append(atom);
                        atom++;
                     }
                  }
               }
            }

   }

}
#endif

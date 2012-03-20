#ifndef ORTHORHOMBIC_BOUNDARY_CPP
#define ORTHORHOMBIC_BOUNDARY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthorhombicBoundary.h"
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/random/Random.h>
#include <util/math/Constants.h>
#include <util/math/feq.h>
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Util
{

   /* 
   * Default constructor.
   */
   OrthorhombicBoundary::OrthorhombicBoundary() 
    : OrthoRegion(),
      lattice_(Orthorhombic)
   {
      for (int i = 0; i < Dimension; ++i) {

         bravaisBasisVectors_.append(Vector::Zero);
         bravaisBasisVectors_[i][i] = lengths_[i];

         reciprocalBasisVectors_.append(Vector::Zero);
         reciprocalBasisVectors_[i][i] = 2.0*Constants::Pi/lengths_[i];

      }
   }

   /* 
   * Set box lengths and then call reset.
   */
   void OrthorhombicBoundary::setLengths(const Vector &lengths) 
   {  
      maxima_  = lengths; 
      lattice_ = Orthorhombic;
      reset(); 
   }

   /* 
   * Set box lengths and then call reset.
   */
   void OrthorhombicBoundary::setTetragonalLengths(double ab, double c) 
   {  
      maxima_[0] = ab;
      maxima_[1] = ab;
      maxima_[2] = c;
      lattice_ = Tetragonal;
      reset(); 
   }

   /* 
   * Set box lengths and call reset.
   */
   void OrthorhombicBoundary::setCubicLengths(double a) 
   {  
      maxima_[0] = a;
      maxima_[1] = a;
      maxima_[2] = a;
      lattice_ = Cubic;
      reset(); 
   }

   /* 
   * Reset all quantities that depend on unit cell lengths.
   */
   void OrthorhombicBoundary::reset()
   {
      resetRegion();
      for (int i = 0; i < Dimension; ++i) {
         bravaisBasisVectors_[i][i] = lengths_[i];
         reciprocalBasisVectors_[i][i] = 2.0*Constants::Pi/lengths_[i];
      }
   }

   /* 
   * Generate a random position within the box
   *
   * \param random random number generator object
   * \param r      Vector of random coordinates
   */
   void OrthorhombicBoundary::randomPosition(Random &random, Vector &r) const 
   {
     for (int i=0; i < Dimension; ++i) {
        r[i] = random.uniform(minima_[i], maxima_[i]);
     }
   }

   /* 
   * Check consistency of data.
   */
   bool OrthorhombicBoundary::isValid() 
   {  
      OrthoRegion::isValid(); 
      for (int i = 0; i < Dimension; ++i) {
         if (!feq(minima_[i], 0.0))
            UTIL_THROW("minima_[i] != 0");
      }
      return true;
   }

   /* 
   * Input a OrthorhombicBoundary from an istream, without line breaks.
   */
   std::istream& operator>>(std::istream& in, OrthorhombicBoundary &boundary)
   {
      LatticeSystem lattice;
      in >> lattice;
      if (lattice == Orthorhombic) {
         in >> boundary.maxima_;
      } else
      if (lattice == Tetragonal) {
         double ab, c;
         in >> ab >> c;
         boundary.maxima_[0] = ab;
         boundary.maxima_[1] = ab;
         boundary.maxima_[2] = c;
      } else 
      if (lattice == Cubic) {
         double a; 
         in >> a; 
         boundary.maxima_[0] = a;
         boundary.maxima_[1] = a;
         boundary.maxima_[2] = a;
      } else {
         UTIL_THROW("Lattice must be orthorhombic, tetragonal or cubic");
      }
      boundary.lattice_ = lattice;
      boundary.reset();
      return in;
   }

   /* 
   * Output an OrthorhombicBoundary to an ostream, without line breaks.
   */
   std::ostream& 
   operator<<(std::ostream& out, const OrthorhombicBoundary &boundary) 
   {
      out << boundary.lattice_ << "   ";
      if (boundary.lattice_ == Orthorhombic) {
         out << boundary.lengths_;
      } else
      if (boundary.lattice_ == Tetragonal) {
         out << Dbl(boundary.lengths_[0]);
         out << Dbl(boundary.lengths_[2]);
      } else 
      if (boundary.lattice_ == Cubic) {
         out << Dbl(boundary.lengths_[0]);
      }
      return out;
   }
 
   #ifdef UTIL_MPI
   template <>
   void send<Util::OrthorhombicBoundary>(MPI::Comm& comm, 
             Util::OrthorhombicBoundary& data, int dest, int tag)
   {
      Vector lengths = data.lengths();
      send<Vector>(comm, lengths, dest, tag);
   }

   template <>
   void recv<Util::OrthorhombicBoundary>(MPI::Comm& comm, 
             Util::OrthorhombicBoundary& data, int source, int tag)
   {
      Vector lengths;
      recv<Vector>(comm, lengths, source, tag);
      data.setLengths(lengths);
   }

   template <>
   void bcast<Util::OrthorhombicBoundary>(MPI::Intracomm& comm, 
              Util::OrthorhombicBoundary& data, int root)
   {
      Vector lengths; 
      int    rank = comm.Get_rank();
      if (rank == root) 
         lengths = data.lengths();
      bcast<Vector>(comm, lengths, root);
      if (rank != root) 
         data.setLengths(lengths);
   }

   /*
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<Util::OrthorhombicBoundary>::type = MPI::BYTE;
   bool MpiTraits<Util::OrthorhombicBoundary>::hasType = false;
   #endif

}

#endif
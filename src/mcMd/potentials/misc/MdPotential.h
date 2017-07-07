#ifndef MCMD_MD_POTENTIAL_H
#define MCMD_MD_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>               // base class
#include <mcMd/potentials/misc/EnergyCalculator.h>   // base class
#include <mcMd/potentials/misc/StressCalculator.h>   // base class

namespace McMd
{

   using namespace Util;

   /**
   * Potential for an MD simulation.
   *
   * \ingroup McMd_Potential_Module
   */
   class MdPotential : public virtual ParamComposite, 
                       public virtual EnergyCalculator, 
                       public virtual StressCalculator
   {

   public:

      /**
      * Destructor (does nothing)
      */
      virtual ~MdPotential();

      /**
      * Add forces from this potential to all atomic forces.
      */
      virtual void addForces() = 0;

   protected:

      /**
      * Constructor.
      *
      * Derived class constructor must set hasStress true or false.
      */
      MdPotential(bool createsStress = true);

   };

} 
#endif

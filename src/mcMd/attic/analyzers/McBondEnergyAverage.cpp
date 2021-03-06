/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McBondEnergyAverage.h"                  
#include <mcMd/analyzers/base/AverageAnalyzer.tpp>
#include <mcMd/potentials/bond/BondPotential.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McBondEnergyAverage::McBondEnergyAverage(McSystem& system)
    : AverageAnalyzer<McSystem>(system)
   {  setClassName("McBondEnergyAverage"); }

   /*
   * Evaluate total bond energy.
   */
   void McBondEnergyAverage::compute() 
   {  value_ = system().bondPotential().energy(); }

}

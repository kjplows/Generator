//____________________________________________________________________________
/*!

\class    genie::PhaseSpaceIteratorModel

\brief    A "bogus" model that scans over a (pN, Eb) space uniformly and returns one nucleon
          from each point in that space.

\author   John Plows

\created  January 5th, 2026

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _PHASESPACEITERATORMODEL_H_
#define _PHASESPACEITERATORMODEL_H_

#include <TH2D.h>
#include <vector>
#include <algorithm>
#include "Physics/NuclearState/NuclearModelI.h"

// This class should:
/*
 *  --> Accept a TH2D with bin edges on (pN, Eb)
 *  --> Construct a vector with (pN, Eb) pairs
 *  --> Implement a Next() function that is called in GenerateNucleon() to pick the next nucleon
 *      Next() should return a suitable end type when end is reached so while loops work
 *      as it will be called upon in the integrator.
 */

namespace genie {

class PhaseSpaceIteratorModel : public NuclearModelI {

public:
  PhaseSpaceIteratorModel();
  PhaseSpaceIteratorModel(string config);
  ~PhaseSpaceIteratorModel() = default;

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- Override the default functions for generating a nucleon at a specific (pN, Eb)
  bool GenerateNucleon(const double pN, const double Eb, const Target & tgt) const;
  double Prob(const double p, const double w, const Target & tgt) const;
  
  //-- implement the NuclearModelI interface
  bool GenerateNucleon(const Target &) const { return GenerateNucleon(0.0, 0.0); }
  NuclearModel_t ModelType (const Target &) const { return kNucmPhaseSpaceIterator; }
  
  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);
  
  //-- Sets the xsec with error
  void SetXSec(const double xs, const double err=0.) const {
    fPhaseSpace->SetBinContent( fX, fY, xs );
    fPhaseSpace->SetBinError( fX, fY, err );
  }

  //-- Getters to interface with callers
  int X() const { return fX; }
  int Y() const { return fY; }
  const TH2D & PhaseSpace() const { return *fPhaseSpace; }

  // -- Iterate over histogram bins
  bool Next() const;
  
 protected:
  //-- Note that LoadConfig() will be responsible for initialising the histogram.
  void   LoadConfig (void);
  
 private:

  // -- Data members

  // -- iterators to bin number
  mutable int fX; ///< Points to bin number along X dimension of phase space TH2D
  mutable int fY; ///< Points to bin number along Y dimension of phase space TH2D

  // -- Phase space histogram. This Model owns the pointer and only it can make changes to the histogram.
  mutable std::unique_ptr<TH2D> fPhaseSpace; ///< Phase space in (pN, Eb) to iterate over. Content and error is xsec

}; // class PhaseSpaceIteratorModel

} // namespace genie

#endif // #ifndef _PHASESPACEITERATORMODEL_H_

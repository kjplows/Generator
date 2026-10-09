#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*!

\class    genie::MarleyGenerator

\brief    Simulate events using an interface to the external MARLEY
          generator for low-energy neutrino interactions.

          MARLEY generates the *complete* event: the primary reaction and the
          full de-excitation cascade of the residual nucleus. The entire
          MARLEY history (excited residues, ejected nucleons and gammas, final
          residue) is copied into the GHEP record, with mother/daughter links
          and GENIE particle statuses, and with momenta converted from MeV
          to GeV. As for GENIE's own nuclear remnant, the *final* residual
          nucleus is recorded as a hadronic blob (kPdgHadronicBlob with status
          kIStFinalStateNuclearRemnant); the earlier de-exciting residues keep
          their nuclear PDG codes.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  July 18, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MARLEY_GENERATOR_H
#define _MARLEY_GENERATOR_H

#include "TLorentzVector.h"

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Physics/MARLEY/MarleyInterface.h"

namespace genie {

class Interaction;

class MarleyGenerator : public EventRecordVisitorI {

public :

  MarleyGenerator();
  MarleyGenerator(string config);
  ~MarleyGenerator();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord* event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry& config);
  void Configure(string config);

  /// Appends one MARLEY particle to the GHEP record (converting MeV -> GeV),
  /// registers it as a daughter of its mother(s), and returns its index
  int AddMarleyParticle( GHepRecord* event,
    const HepMC3::GenParticle& part, int mom1, int mom2,
    GHepStatus_t status, const TLorentzVector& v4 ) const;

private:

  void LoadConfig(void);

  const genie::MarleyInterface* fMARLEY;

};

}      // genie namespace
#endif // _MARLEY_GENERATOR_H
#endif // __GENIE_MARLEY_ENABLED__

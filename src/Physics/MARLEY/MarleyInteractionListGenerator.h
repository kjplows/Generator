#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*!

\class    genie::MarleyInteractionListGenerator

\brief    Creates a list of the MARLEY interactions (ProcessInfo: kScMarley)
          that can be simulated for a given initial state.

          An interaction is listed only if the MARLEY generator held by the
          algorithm's MarleyInterface (i.e., the one configured by the
          MARLEY job file for this channel) has at least one reaction
          matching the probe, the nuclear target, and the interaction type
          (CC or NC) requested for this list generator. This way the
          MARLEY channel only appears for the initial states that MARLEY
          can actually handle (e.g., nu_e + 40Ar for CC).

          Configurable parameters:
            is-CC   bool  (default false)  create MARLEY-CC interactions
            is-NC   bool  (default false)  create MARLEY-NC interactions
            MarleyAlg  alg                  genie::MarleyInterface instance

\created  October 2026
*/
//____________________________________________________________________________

#ifndef _MARLEY_INTERACTION_LIST_GENERATOR_H_
#define _MARLEY_INTERACTION_LIST_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class MarleyInterface;

class MarleyInteractionListGenerator : public InteractionListGeneratorI {

public :

  MarleyInteractionListGenerator();
  MarleyInteractionListGenerator(string config);
  ~MarleyInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig(void);

  bool fIsCC;
  bool fIsNC;
  const MarleyInterface * fMARLEY;
};

}      // genie namespace

#endif // _MARLEY_INTERACTION_LIST_GENERATOR_H_
#endif // __GENIE_MARLEY_ENABLED__

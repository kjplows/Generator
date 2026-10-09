#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

// Standard library includes
#include <cassert>

// GENIE includes
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InteractionType.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Interaction/ScatteringType.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/MARLEY/MarleyInteractionListGenerator.h"
#include "Physics/MARLEY/MarleyInterface.h"

// MARLEY includes
#include "marley/Generator.hh"
#include "marley/Reaction.hh"

using namespace genie;

//___________________________________________________________________________
MarleyInteractionListGenerator::MarleyInteractionListGenerator() :
InteractionListGeneratorI("genie::MarleyInteractionListGenerator"),
fIsCC(false), fIsNC(false), fMARLEY(nullptr)
{

}
//___________________________________________________________________________
MarleyInteractionListGenerator::MarleyInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::MarleyInteractionListGenerator", config),
fIsCC(false), fIsNC(false), fMARLEY(nullptr)
{

}
//___________________________________________________________________________
MarleyInteractionListGenerator::~MarleyInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * MarleyInteractionListGenerator::CreateInteractionList(
  const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  if ( !fIsCC && !fIsNC ) {
    LOG("IntLst", pWARN) << "Neither is-CC nor is-NC is set. Returning a NULL"
      << " InteractionList for init-state: " << init_state.AsString();
    return 0;
  }

  const int probe_pdg = init_state.ProbePdg();
  const int tgt_pdg   = init_state.TgtPdg();

  // Check whether the MARLEY generator for this channel knows about any
  // nuclear reaction of the requested type for this probe and target
  using PT = marley::Reaction::ProcessType;
  const marley::Generator * mgen = fMARLEY->GetMarleyGenerator();

  bool found = false;
  for ( const auto& react : mgen->get_reactions() ) {

    if ( react->pdg_a() != probe_pdg ) continue;
    if ( react->atomic_target().pdg() != tgt_pdg ) continue;

    const PT pt = react->process_type();
    const bool is_cc = ( pt == PT::NeutrinoCC_Discrete
      || pt == PT::NeutrinoCC_Continuum
      || pt == PT::AntiNeutrinoCC_Discrete
      || pt == PT::AntiNeutrinoCC_Continuum );
    const bool is_nc = ( pt == PT::NC_Discrete || pt == PT::NC_Continuum );

    if ( (fIsCC && is_cc) || (fIsNC && is_nc) ) {
      found = true;
      break;
    }
  }

  if ( !found ) {
    LOG("IntLst", pINFO) << "MARLEY has no matching reaction. Returning a"
      << " NULL InteractionList for init-state: " << init_state.AsString();
    return 0;
  }

  InteractionList * intlist = new InteractionList;

  InteractionType_t itype = fIsCC ? kIntWeakCC : kIntWeakNC;
  ProcessInfo proc_info( kScMarley, itype );
  Interaction * interaction = new Interaction( init_state, proc_info );

  // MARLEY does not use a hit nucleon: it treats the nucleus as a whole
  intlist->push_back( interaction );
  return intlist;
}
//___________________________________________________________________________
void MarleyInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyInteractionListGenerator::LoadConfig(void)
{
  GetParamDef( "is-CC", fIsCC, false );
  GetParamDef( "is-NC", fIsNC, false );

  fMARLEY = dynamic_cast<const MarleyInterface*>( this->SubAlg("MarleyAlg") );
  assert( fMARLEY );
}
//___________________________________________________________________________

#endif // __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*!

\class    genie::DummyLLPInteractionListGenerator

\brief    

\author   John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford

	  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	  University of Liverpool & STFC Rutherford Appleton Laboratory

\created  February 9th, 2022

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DUMMY_LLP_INTERACTION_LIST_GENERATOR_H_
#define _DUMMY_LLP_INTERACTION_LIST_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class DummyLLPInteractionListGenerator : public InteractionListGeneratorI {

public :
  DummyLLPInteractionListGenerator();
  DummyLLPInteractionListGenerator(string config);
 ~DummyLLPInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;
};

}      // genie namespace
#endif // _DUMMY_LLP_INTERACTION_LIST_GENERATOR_H_

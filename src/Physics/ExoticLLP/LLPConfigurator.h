//----------------------------------------------------------------------------
/*!
  
  An Algorithm that constructs an ExoticLLP instance with knowledge of its
  own mass and visible states it couples with.

\class    genie::bsm::LLPConfigurator

\brief    Constructs concrete instance of ExoticLLP

\author   John Plows <kjplows \at liverpool.ac.uk>
          University of Liverpool

\created  May 10th, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _LLP_CONFIGURATOR_H_
#define _LLP_CONFIGURATOR_H_

// -- C++ includes
#include <iterator>
#include <sstream>
#include <utility>

// -- GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

#include "Physics/ExoticLLP/LLPChannelCalculatorI.h"
#include "Physics/ExoticLLP/ExoticLLP.h"

namespace genie {
  
  class PDGLibrary;

  namespace llp {

    class ExoticLLP;
    class LLPConfigurator : public ChannelCalculatorI {

    public:

      LLPConfigurator();
      LLPConfigurator(string name);
      LLPConfigurator(string name, string config);
      ~LLPConfigurator();

      // -- implement the ChannelCalculatorI interface

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

      // inline function to allow access to the underlying LLP
      inline ExoticLLP RetrieveLLP() const { return fLLP; }

    private:

      void LoadConfig(void);

      mutable bool fIsConfigLoaded;

      mutable ExoticLLP fLLP; //! The concrete LLP instance.
      mutable double fMass; //! LLP mass in MeV

    }; // class LLPConfigurator

  } // namespace llp

} // namespace genie

#endif // #ifndef _LLP_CONFIGURATOR_H_

//----------------------------------------------------------------------------
/*!

\class    genie::llp::Decayer

\brief    Responsible for generating decays of particles, with or without angular correlations.
          Singleton class.

\author   John Plows <kjplows \at liverpool.ac.uk>
          University of Liverpool

\created  August 7th, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _LLP_DECAYER_H_
#define _LLP_DECAYER_H_

#include <cmath>
#include <cassert>

#include <TGenPhaseSpace.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

#include "Physics/ExoticLLP/ExoticLLP.h"
#include "Physics/ExoticLLP/LLPConfigurator.h"

typedef std::vector< genie::GHepParticle > LorentzMap; //! To store the 4-momentum of each product

namespace genie {

  class PDGCodeList;
  class GHepParticle;

  namespace llp {

    class ExoticLLP;
    struct ModeObject;

    class Decayer {

    public:

      //! Access instance
      static Decayer * Instance();

      //! Pass the Decayer a PDGCodeList to use
      void SetProducts( genie::PDGCodeList pdgv ) const;
      //! And pass it a boost vector to use
      void SetBoost( TVector3 boost_vec ) const;

      //! Perform an unpolarised decay
      bool UnpolarisedDecay( bool fudge = false ) const;

      //! Get the results back
      LorentzMap GetResults() const;
      LorentzMap GetRestFrameResults() const;

      double GetMasslessEnergy() const;

      void ClearEvent() const;

    private:

      Decayer();
      Decayer( const Decayer & dec );
      virtual ~Decayer();

      mutable double fLLPMass; //! This is very important so we will keep it.
      mutable genie::PDGCodeList fPDGCodeList;
      mutable LorentzMap fParticles; //! Lab frame
      mutable LorentzMap fParticles_rest;

      mutable double fMasslessEnergy; //! needed for acceptance calculations. See LLPFluxCreator

      mutable TVector3 fBoostVec; //! Boost vector of the lab frame

      static Decayer * fInstance;
      bool fInitialized; //! Done initialising singleton?

      struct Cleaner {
	void DummyMethodAndSilentCompiler() {} // see RandomGen.h
	~Cleaner() {
	  if( Decayer::fInstance != 0 ) {
	    delete Decayer::fInstance;
	    Decayer::fInstance = 0;
	  }
	}
      }; // struct Cleaner

      friend struct Cleaner;

      ClassDef(Decayer, 1)

    }; // class Decayer
    

  } // namespace llp

} // namespace genie

#endif // #ifndef _LLP_DECAYER_H_

//----------------------------------------------------------------------------
/*!

  This is a module for GENIE to read in hadron flux ntuples and construct LLP fluxes
  on the fly. 
  Assumes the hard geometrical calculations have already been run, and saved into a FluxContainer
  object which allows the module to constrain kinematics.
  
  Core loop: + Get ancestry information
	     + Assume decay to LLP
	     + Calculate LLP production mode based on parameter space read from config
	     + Calculate kinematics of LLP.
	     + Update the FluxContainer with constrained kinematics.

\class      genie::llp::FluxCreator

\brief      Calculates LLP production kinematics.
            Is a concrete implementation of the FluxRecordVisitorI interface

\author     John Plows <kplows \at liverpool.ac.uk>
            University of Liverpool

\created    June 24th, 2024

\cpright    Copyright (c) 2003-2024, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

 */
//----------------------------------------------------------------------------

#ifndef _LLP_FLUXCREATOR_H_
#define _LLP_FLUXCREATOR_H_

// -- C++ includes
#include <array>
#include <cassert>
#include <iomanip> // for momentum balance stream
#include <map>
#include <list>
#include <sstream>
#include <unordered_map>

// -- ROOT includes
#include "TChain.h"
#include "TDecayChannel.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

// -- GENIE includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

#include "Physics/ExoticLLP/LLPFluxRecordVisitorI.h"
#include "Physics/ExoticLLP/LLPFluxContainer.h"
#include "Physics/ExoticLLP/ExoticLLP.h"
#include "Physics/ExoticLLP/LLPConfigurator.h"
#include "Physics/ExoticLLP/AliasedBranch.h"
#include "Physics/ExoticLLP/VolumeSeeker.h"
#include "Physics/ExoticLLP/Decayer.h"

namespace genie{

  class PDGCodeList;
  class GHepParticle;

  namespace llp{
    
    struct ModeObject;
    class ExoticLLP;
    class FluxContainer;
    
    class FluxCreator : public FluxRecordVisitorI {

    public:

      FluxCreator();
      FluxCreator(string name);
      FluxCreator(string name, string config);
      ~FluxCreator();

      //-- implement the FluxRecordVisitorI interface
      void ProcessEventRecord(GHepRecord * event_rec) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);
      
      // set the input path for a flux
      void SetInputFluxPath( std::string finpath ) const;
      // get N(flux input entries)
      int GetNFluxEntries() const;
      // set first entry for read-in from chain
      void SetFirstFluxEntry( int iFirst ) const;

      // update flux info
      void UpdateFluxInfo( FluxContainer info ) const;
      // get flux info
      FluxContainer RetrieveFluxInfo() const;

    private:

      void LoadConfig(void);

      void SetCurrentEntry( int iCurr ) const;

      // init
      void OpenFluxInput( std::string finpath ) const;
      std::list<TString> RecurseOverDir( std::string finpath ) const;
      void InitialiseTree() const;
      void InitialiseMeta() const;

      // returns LLP 4-momentum from random decay in same frame as p4_parent
      TLorentzVector LLPEnergy() const;

      // There is no longer a choice at FluxCreator to point to random place in BBox.
      // This is handled by VolumeSeeker and passed to the FluxContainer

      // No reading in BRs either. That's taken up by the ExoticLLP instance.

      /*
       * There are three components for the acceptance calculation. 
       * 1) Lab-frame angular size, handled by VolumeSeeker.
       * 2) Acceptance correction, == 1 identically for massless.
       * 3) Boost correction factor, gives rest-frame angular size as S |-> B^2 S
       
       * The reason the acceptance correction needs to be done first is because first you choose the 
       * rest-frame emission angle (forward or backward?), and then you get the energy based on that.
       */
      
      // The acceptance correction. Returns the LLP energy lab frame (first) and the correction (second).
      std::pair< double, double > CalculateAcceptanceCorrection( TLorentzVector p4_LLP, 
								 double zeta ) const;
      // helper methods for the above method.
      // Need to reduce that...
      double AccCorr_Sqrt( double thetalab, double mass, double EPar, double MPar, double ENu ) const;
      double AccCorr_Denom( double thetalab, double mass, double EPar, double MPar, double ENu ) const;
      double AccCorr_SolnArgs( double thetalab, double mass, double EPar, double MPar, double ENu,
			       bool isPos ) const;
      double AccCorr_Solution( double thetalab, double mass, double EPar, double MPar, double ENu,
			       bool isPos ) const;
      double Forwards_Fcn( double Theta, TLorentzVector p4par, TLorentzVector p4LLP ) const; // deg, deg
      double Inverted_Fcn( double theta, TLorentzVector p4par, TLorentzVector p4LLP, bool backwards ) const; // deg, deg
      double Derivative( double theta, TLorentzVector p4par, TLorentzVector p4LLP ) const;
      // helper method
      static double labangle( double * x, double * par ); // function formula for correction. x is IN DEGREES. par is [E_parent, p(x,y,z)_parent, p_llp_rest, E_llp_rest (GeV)]
      
      // And a function that just returns the energy in lab frame from rest frame emission angle
      double CalculateLabFrameEnergy( double Theta, TLorentzVector p4par, TLorentzVector p4LLP ) const;

      // A function that calculates the rest-frame kinematics of the entire particle stack, given
      // the geometrical constraint of flux.
      std::vector<genie::GHepParticle> 
	FindStackKinematics( std::vector<genie::GHepParticle> random_system,
			     TVector3 desired_direction ) const;

      // some utilities for C++ strings, grabbed from https://www.stackoverflow.com/questions/874134/
      bool ends_with(const std::string & str, const std::string & suffix) const {
	return str.size() >= suffix.size() && str.compare( str.size() - suffix.size(), suffix.size(), suffix ) == 0;
      }
    
      bool begins_with(const std::string & str, const std::string & prefix) const {
	return str.size() >= prefix.size() && str.compare( 0, prefix.size(), prefix ) == 0;
      }

      // current path to keep track of what is loaded
      mutable std::string fCurrPath = ""; mutable bool fPathLoaded = false;
      // and which entry we're on
      mutable int iCurrEntry = 0;
      // which one was first?
      mutable int fFirstEntry = 0;
      // out of how many?
      mutable int fNEntries = 0;

      mutable TChain * ctree = 0, * cmeta = 0;
      mutable string ctree_name = "", cmeta_name = "";

      mutable double fMass; // LLP mass, GeV

      // tree variables. Always a pair of branch *value* and branch *alias*.
      // This is done so that the user can name their branches at runtime.

      // ctree
      mutable AliasedBranch<int>            m_ct_parent_pdg;    ///< PDG code of parent
      mutable std::string                   m_ct_parent_pdg_alias;
      mutable AliasedBranch<TLorentzVector> m_ct_parent_p4;     ///< Parent 4-momentum, NEAR [GeV]
      mutable AliasedBranch<TLorentzVector> m_ct_production_v4; ///< Production vertex, NEAR [m, ns]
      // cmeta
      mutable AliasedBranch<int>            m_cm_run_no;        ///< Simulation run number

      mutable genie::llp::FluxContainer fFluxInfo;
      
      mutable bool fIsConfigLoaded = false;

      // The underlying ExoticLLP instance... It's a particle type, basically.
      mutable ExoticLLP fExoticLLP;

      // A convenient unordered map to keep track of which production modes are from which parent
      mutable std::unordered_map< int, std::vector< ModeObject > > fGroupedModes;

    }; // class FluxCreator
      
  } // namespace llp
} // namespace genie


#endif // #ifndef _LLP_FLUXCREATOR_H_

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// Standard library includes
#include <algorithm>
#include <array>
#include <cstdlib>
#include <map>
#include <vector>

// GENIE includes
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/HepMC3Converter.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/MARLEY/MarleyGenerator.h"
#include "Physics/MARLEY/MarleyInterface.h"

#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/hepmc3_utils.hh"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
MarleyGenerator::MarleyGenerator() :
  EventRecordVisitorI( "genie::MarleyGenerator" )
{

}
//___________________________________________________________________________
MarleyGenerator::MarleyGenerator(string config) :
  EventRecordVisitorI( "genie::MarleyGenerator", config )
{

}
//___________________________________________________________________________
MarleyGenerator::~MarleyGenerator()
{

}
//___________________________________________________________________________
void MarleyGenerator::ProcessEventRecord(GHepRecord* event) const
{
  // Only used for the conversion of the NuHepMC particle status codes
  // (MARLEY and GENIE use the same ones for remnants, see HepMC3Converter)
  static genie::HepMC3Converter hepmc3_conv;

  genie::Interaction* inter = event->Summary();

  const InitialState& init_state = inter->InitState();
  int probe_pdg = init_state.ProbePdg();
  int tgt_pdg = init_state.TgtPdg();

  TLorentzVector* temp_probe_p4 = init_state.GetProbeP4( kRfLab );
  TLorentzVector probe_p4 = *temp_probe_p4;
  delete temp_probe_p4;

  // Probe kinetic energy (MeV)
  double probe_KE = ( probe_p4.E() - probe_p4.M() ) / genie::units::MeV;

  // Build a unit vector in the probe's direction of motion
  TVector3 dir = probe_p4.Vect().Unit();
  std::array< double, 3 > probe_dir = { dir.X(), dir.Y(), dir.Z() };

  // Get the 4-position of the interaction vertex (set by GENIE)
  TLorentzVector v4( *event->Probe()->X4() );

  // Create a new MARLEY event for the given initial state. MARLEY picks the
  // reaction and the kinematics, and (if enabled in its job file) simulates
  // the de-excitation of the residual nucleus.
  // The MARLEY generator is seeded from the GENIE seed (see MarleyInterface).
  marley::Generator* marley_gen = fMARLEY->GetMarleyGenerator();
  std::shared_ptr< HepMC3::GenEvent > marley_evt;
  try {
    marley_evt = marley_gen->create_event( probe_pdg, probe_KE, tgt_pdg,
      probe_dir );
  }
  // If MARLEY runs into a problem, convert its exception into a GENIE
  // EVGThreadException and back up the event generation thread
  catch ( const marley::Error& err ) {
    event->EventFlags()->SetBitNumber( genie::kGenericErr, true );
    genie::exceptions::EVGThreadException exception;
    exception.SetReason( err.what() );
    exception.SwitchOnStepBack();
    exception.SetReturnStep( 0 );
    throw exception;
  }

  const int probe_idx = event->ProbePosition();
  const int tgt_idx = event->TargetNucleusPosition();

  // Maps the (1-based) HepMC3 particle ID in the MARLEY event to the index of
  // the corresponding particle in the GHEP record. The MARLEY projectile and
  // target are identified with the probe and the target nucleus that GENIE
  // already put in the record, so they are not added again.
  std::map< int, int > marley_to_ghep;

  // Ions that MARLEY adds to the record (residues and, possibly, ejected light
  // ions), used below to identify the final residual nucleus
  struct IonInfo { int ghep_idx; int A; bool final_state; };
  std::vector< IonInfo > ions;

  // The particles are listed in the order that they were created, so a
  // particle always comes after all of its mothers
  for ( const auto& part : marley_evt->particles() ) {

    const int mstatus = part->status();

    if ( mstatus == marley_hepmc3::NUHEPMC_PROJECTILE_STATUS ) {
      marley_to_ghep[ part->id() ] = probe_idx;
      continue;
    }
    if ( mstatus == marley_hepmc3::NUHEPMC_TARGET_STATUS ) {
      marley_to_ghep[ part->id() ] = tgt_idx;
      continue;
    }

    // Find the mother(s) of this particle using its production vertex
    std::vector< int > moms;
    const auto& prod_vtx = part->production_vertex();
    if ( prod_vtx ) {
      for ( const auto& in : prod_vtx->particles_in() ) {
        auto it = marley_to_ghep.find( in->id() );
        if ( it != marley_to_ghep.end() ) moms.push_back( it->second );
      }
    }

    int mom1 = -1, mom2 = -1;
    const int abs_pdg = std::abs( part->pid() );
    if ( abs_pdg >= 11 && abs_pdg <= 16 ) {
      // The final-state primary lepton (charged lepton for CC, outgoing
      // neutrino for NC) is a daughter of the probe, as in other GENIE
      // channels
      mom1 = probe_idx;
    }
    else {
      // The probe is not a mother of anything but the primary lepton. A
      // particle that has no other mother descends from the target nucleus.
      moms.erase( std::remove( moms.begin(), moms.end(), probe_idx ),
        moms.end() );
      if ( moms.empty() ) mom1 = tgt_idx;
      else {
        mom1 = moms.front();
        if ( moms.size() > 1u ) mom2 = moms.back();
      }
    }

    GHepStatus_t status = hepmc3_conv.GetGHepParticleStatus( mstatus );
    int idx = this->AddMarleyParticle( event, *part, mom1, mom2, status, v4 );
    marley_to_ghep[ part->id() ] = idx;

    if ( genie::pdg::IsIon(part->pid()) ) {
      ions.push_back( { idx, genie::pdg::IonPdgCodeToA(part->pid()),
        mstatus == marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS } );
    }
  }

  // Put the final residual nucleus into a "hadronic blob", as GENIE does for
  // its own nuclear remnant (pseudo-particle with PDG code kPdgHadronicBlob
  // and status kIStFinalStateNuclearRemnant, carrying the remnant 4-momentum).
  // The earlier, de-exciting residues keep their nuclear PDG codes and
  // intermediate statuses, so the history of the cascade is preserved.
  //
  // The final residue is the heaviest ion in the final state. (Light ions
  // like alpha particles can be ejected during the de-excitation.) If MARLEY
  // did not mark any ion as final-state (e.g., a ground-state transition with
  // no de-excitation steps), use the last ion that it listed.
  int blob_idx = -1;
  int best_A = -1;
  for ( const auto& ion : ions ) {
    if ( ion.final_state && ion.A >= best_A ) {
      best_A = ion.A;
      blob_idx = ion.ghep_idx;
    }
  }
  if ( blob_idx < 0 && !ions.empty() ) blob_idx = ions.back().ghep_idx;

  if ( blob_idx >= 0 ) {
    GHepParticle* blob = event->Particle( blob_idx );
    blob->SetPdgCode( genie::kPdgHadronicBlob );
    blob->SetStatus( genie::kIStFinalStateNuclearRemnant );
  }
  else {
    LOG("Marley", pWARN) << "No residual nucleus found in the MARLEY event."
      << " No hadronic blob was added to the event record.";
  }

  // Set the final-state lepton's polarization
  genie::utils::SetPrimaryLeptonPolarization( event );
}
//___________________________________________________________________________
void MarleyGenerator::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyGenerator::Configure(string config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyGenerator::LoadConfig(void)
{
  fMARLEY = dynamic_cast<const MarleyInterface*>( this->SubAlg("MarleyAlg") );
}
//___________________________________________________________________________
int MarleyGenerator::AddMarleyParticle( GHepRecord* event,
  const HepMC3::GenParticle& part, int mom1, int mom2,
  GHepStatus_t status, const TLorentzVector& v4 ) const
{
  // Get the particle's 4-momentum, and convert to using GeV instead of MeV
  const HepMC3::FourVector& mom4 = part.momentum();
  TLorentzVector p4( mom4.px() * genie::units::MeV,
    mom4.py() * genie::units::MeV, mom4.pz() * genie::units::MeV,
    mom4.e() * genie::units::MeV );

  event->AddParticle( part.pid(), status, mom1, mom2, -1, -1, p4, v4 );
  const int idx = event->GetEntries() - 1;

  // Register this particle as a daughter of its mother(s). GHEP expects the
  // daughters of a particle to be contiguous, which is the case here because
  // MARLEY lists the products of each step together.
  for ( int m : { mom1, mom2 } ) {
    if ( m < 0 ) continue;
    GHepParticle* mother = event->Particle( m );
    if ( mother->FirstDaughter() < 0 ) mother->SetFirstDaughter( idx );
    mother->SetLastDaughter( idx );
  }

  return idx;
}
//___________________________________________________________________________

#endif // __GENIE_MARLEY_ENABLED__

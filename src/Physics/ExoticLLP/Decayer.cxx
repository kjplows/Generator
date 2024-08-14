//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kjplows \at liverpool.ac.uk>
 University of Liverpool

*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/Decayer.h"

using namespace genie;
using namespace genie::llp;

//____________________________________________________________________________
Decayer * Decayer::fInstance = 0;
double fLLPMass = 0;
double fMasslessEnergy = 0;
PDGCodeList fPDGCodeList;
LorentzMap fParticles;
LorentzMap fParticles_rest;
//____________________________________________________________________________
Decayer::Decayer()
{
  LOG("ExoticLLP", pDEBUG) << "Decayer late initialization";

  fInitialized = false;
  fInstance = 0;
  fLLPMass = 0.0;
  fMasslessEnergy = 0.0;

  if(fParticles.size() > 0) fParticles.clear(); 
  if(fParticles_rest.size() > 0) fParticles_rest.clear();
  if(fPDGCodeList.size() > 0) fPDGCodeList.clear();
}
//____________________________________________________________________________
Decayer::~Decayer()
{
  fInstance = 0;
  fInitialized = false;
  fLLPMass = 0.0;
  fMasslessEnergy = 0.0;

  if(fParticles.size() > 0) fParticles.clear(); 
  if(fParticles_rest.size() > 0) fParticles_rest.clear();
  if(fPDGCodeList.size() > 0) fPDGCodeList.clear();
}
//____________________________________________________________________________
void Decayer::ClearEvent() const
{
  if(fParticles.size() > 0) fParticles.clear(); 
  if(fParticles_rest.size() > 0) fParticles_rest.clear();
  if(fPDGCodeList.size() > 0) fPDGCodeList.clear();
}
//____________________________________________________________________________
Decayer * Decayer::Instance()
{
  if( fInstance == 0 ) {
    static Decayer::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new Decayer;
  }
  return fInstance;
}
//____________________________________________________________________________
bool Decayer::UnpolarisedDecay( bool fudge ) const
{
  // first, get the decay product masses
  std::vector<int>::const_iterator pdg_iter;
  int idec = 0;
  assert(fPDGCodeList.size() > 1 && "At least one daughter of 0th particle");
  double * mass = new double[fPDGCodeList.size()-1];
  double   sum  = 0.0;
  
  for( pdg_iter = fPDGCodeList.begin()+1; pdg_iter != fPDGCodeList.end(); ++pdg_iter ) {
    double m = PDGLibrary::Instance()->Find( *(pdg_iter) )->Mass();
    if( fudge && 
	std::abs(PDGLibrary::Instance()->Find( *(pdg_iter) )->PdgCode()) == kPdgLLP ) m = 0.0;
    mass[idec++] = m;
    sum += m;
  }

  LOG("ExoticLLP", pINFO)
    << "Decaying N = " << fPDGCodeList.size()-1 << " particles / total mass = " << sum;

  double parent_mass = PDGLibrary::Instance()->Find( *(fPDGCodeList.begin()) )->Mass();
  // First construct the rest-frame parent
  TLorentzVector parent_p4_rest( 0.0, 0.0, 0.0, parent_mass );

  // Set the decay
  TGenPhaseSpace PSGen;
  bool permitted = PSGen.SetDecay( parent_p4_rest, fPDGCodeList.size()-1, mass );
  if( !permitted ) {
    LOG("ExoticLLP", pERROR)
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << sum << "\n"
       << " Parent mass = " << parent_mass;
    // clean up and throw exception
    delete [] mass;
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Decay not permitted kinematically");
    exception.SwitchOnFastForward();
    throw exception;
  } // decay not allowed

  double wmax = -1.0;
  for( int i = 0; i < 200; i++ ) {
    double w = PSGen.Generate();
    wmax = std::max( wmax, w );
  }
  assert(wmax>0 && "Phase-space decayer works");
  wmax *= 2;

  LOG("ExoticLLP", pNOTICE)
     << "Max phase space gen. weight @ current system: " << wmax;

  // Do the actual decay;
  RandomGen * rnd = RandomGen::Instance();
  bool accept_decay = false;
  unsigned int itry = 0;

  while( ! accept_decay ) {
    itry++;
    if( itry > controls::kMaxUnweightDecayIterations ) {
      LOG("ExoticLLP", pWARN)
	<< "Couldn't generate an unweighted phase space decay after "
	<< itry << " attempts";
      permitted = true;
      return permitted;
    } // too many attempts
    
    double w = PSGen.Generate();
    if( w > wmax )
      LOG( "ExoticLLP", pWARN ) << "Decay weight = " << w << " > max decay weight = " << wmax;
    
    double gw = wmax * rnd->RndGen().Rndm();
    accept_decay = ( gw <= w );
  }

  // RETHERE: boost the particles into the lab frame if you've been given the boost factor!
  // Insert final state products into fParticles_rest;
  // Note the parent does NOT make it to this stage

  for( int idp = 0; idp < fPDGCodeList.size()-1; idp++ ) {
    TLorentzVector p4 = *(PSGen.GetDecay(idp));
    TLorentzVector v4(0.0, 0.0, 0.0, 0.0); // we don't really care about the v4 info at this stage
    GHepStatus_t ist = kIStStableFinalState;
    GHepParticle particle_in_stack( fPDGCodeList.at(idp+1), ist, 0, -1, -1, -1, p4, v4 );
    if( fudge && std::abs(particle_in_stack.Pdg()) == kPdgLLP ) {
      particle_in_stack.SetPdgCode( 0 ); // this is a special particle only meant for acceptance calcs
      fMasslessEnergy = particle_in_stack.E();
    }
    if( ! fudge )
      fParticles_rest.emplace_back( particle_in_stack );
  }

  LOG( "ExoticLLP", pDEBUG ) << "After a successful decay there are " << fParticles_rest.size()
			     << " particles in the stack";

  permitted = true;
  return permitted;
}
//____________________________________________________________________________
LorentzMap Decayer::GetResults() const
{
  return fParticles;
}
//____________________________________________________________________________
LorentzMap Decayer::GetRestFrameResults() const
{
  return fParticles_rest;
}
//____________________________________________________________________________
double Decayer::GetMasslessEnergy() const
{
  return fMasslessEnergy;
}
//____________________________________________________________________________
bool Decayer::PrepareDecay() const
{
  LOG("ExoticLLP", pFATAL) << "RETHERE.";
  return true;
}
//____________________________________________________________________________
void Decayer::SetProducts( PDGCodeList pdgv ) const
{ 
  LOG( "ExoticLLP", pDEBUG ) << "Input PDGCodeList has " << pdgv.size() << " elements";
  for( PDGCodeList::iterator itp = pdgv.begin(); itp != pdgv.end(); ++itp )
    {
      LOG("ExoticLLP", pDEBUG) << "Adding element " << *itp;
      fPDGCodeList.push_back( *itp );
    }
  LOG( "ExoticLLP", pDEBUG ) << "Decayer PDGCodeList now has " << fPDGCodeList.size() << " elements";
}

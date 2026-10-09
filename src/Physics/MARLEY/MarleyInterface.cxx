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
#include <cstdlib>

// GENIE includes
#include "Physics/MARLEY/MarleyInterface.h"

// MARLEY includes
#include "marley/JSONConfig.hh"

// GENIE includes
#include "Framework/Numerical/RandomGen.h"

// Standard library includes
#include <cstdint>
#include <functional>

using namespace genie;

//____________________________________________________________________________
MarleyInterface::MarleyInterface() : Algorithm( "genie::MarleyInterface" )
{

}
//____________________________________________________________________________
MarleyInterface::MarleyInterface(string config) :
Algorithm( "genie::MarleyInterface", config )
{

}
//____________________________________________________________________________
MarleyInterface::~MarleyInterface()
{

}
//____________________________________________________________________________
void MarleyInterface::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void MarleyInterface::Configure( string config )
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void MarleyInterface::LoadConfig(void)
{
  // Get configuration file name for initializing MARLEY
  std::string config_file_name;
  GetParam( "ConfigFileName", config_file_name ) ;

  // Initialize a new marley::Generator object
  std::string full_path = std::getenv( "GENIE" );
  full_path += "/data/evgen/marley/" + config_file_name;
  marley::JSONConfig jc( full_path );
  fMarleyGenerator = jc.create_generator();

  // Seed MARLEY from the GENIE random number seed (so that jobs with
  // different GENIE seeds do not produce identical MARLEY events). The name
  // of this configuration is mixed in so that several MARLEY generators
  // (e.g., MARLEY-CC and MARLEY-NC) in the same job are not correlated.
  std::uint64_t seed = static_cast< std::uint64_t >(
    genie::RandomGen::Instance()->GetSeed() );
  seed ^= static_cast< std::uint64_t >(
    std::hash< std::string >()( this->Id().Key() ) ) + 0x9e3779b97f4a7c15ULL
    + ( seed << 6 ) + ( seed >> 2 );
  fMarleyGenerator.reseed( seed );
}
//____________________________________________________________________________
marley::Generator* MarleyInterface::GetMarleyGenerator() const {
  return &fMarleyGenerator;
}
//____________________________________________________________________________

#endif // __GENIE_MARLEY_ENABLED__

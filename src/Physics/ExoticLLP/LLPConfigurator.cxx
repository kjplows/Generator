//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kjplows \at liverpool.ac.uk>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/LLPConfigurator.h"

using namespace genie;
using namespace genie::llp;

//____________________________________________________________________________
LLPConfigurator::LLPConfigurator() : 
  ChannelCalculatorI("genie::llp::LLPConfigurator")
{
  fIsConfigLoaded = false;
}
//____________________________________________________________________________
LLPConfigurator::LLPConfigurator(string name) : ChannelCalculatorI(name)
{
  fIsConfigLoaded = false;
}
//____________________________________________________________________________
LLPConfigurator::LLPConfigurator(string name, string config) : 
  ChannelCalculatorI(name, config)
{
  fIsConfigLoaded = false;
}
//____________________________________________________________________________
LLPConfigurator::~LLPConfigurator()
{

}
//____________________________________________________________________________
void LLPConfigurator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LLPConfigurator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LLPConfigurator::LoadConfig(void)
{
  if( fIsConfigLoaded ) return;

  this->GetParam( "Mass", fMass );

  // RETHERE Majorana vs Dirac?
  // RETHERE spin?
  
  // RETHERE this should read in the path to an XML and get the appropriate modes.
  // For now, just make a single production mode K+ --> X + e+ and decay mode X --> e- e+
  std::vector<int> prodDummy = { kPdgKP, kPdgLLP, kPdgPositron };
  std::vector<int> decDummy  = { kPdgLLP, kPdgElectron, kPdgPositron };

  std::pair< double, std::vector<int> > dprod_pair = std::make_pair( 1.0, prodDummy );
  std::pair< double, std::vector<int> > ddec_pair  = std::make_pair( 1.0, decDummy  );

  std::vector< std::pair< double, std::vector<int> > > dmodes = { dprod_pair, ddec_pair };

  ExoticLLP * ptLLP = new ExoticLLP( fMass, dmodes );
  
  fLLP = *ptLLP;
}
//____________________________________________________________________________

//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kjplows \at liverpool.ac.uk>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/LLPChannelCalculatorI.h"

using namespace genie;
using namespace genie::llp;

//____________________________________________________________________________
ChannelCalculatorI::ChannelCalculatorI() : Algorithm()
{

}
//____________________________________________________________________________
ChannelCalculatorI::ChannelCalculatorI(string name) : Algorithm(name)
{

}
//____________________________________________________________________________
ChannelCalculatorI::ChannelCalculatorI(string name, string config) : 
  Algorithm(name, config)
{

}
//____________________________________________________________________________

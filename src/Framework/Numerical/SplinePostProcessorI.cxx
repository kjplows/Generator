//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kplows \at liverpool.uk>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/Numerical/SplinePostProcessorI.h"

using namespace genie;

//____________________________________________________________________________
SplinePostProcessorI::SplinePostProcessorI() :
  Algorithm()
{

}
//____________________________________________________________________________
SplinePostProcessorI::SplinePostProcessorI(string name) :
  Algorithm(name)
{

}
//____________________________________________________________________________
SplinePostProcessorI::SplinePostProcessorI(string name, string config) :
  Algorithm(name, config)
{

}
//____________________________________________________________________________
SplinePostProcessorI::~SplinePostProcessorI()
{

}
//____________________________________________________________________________

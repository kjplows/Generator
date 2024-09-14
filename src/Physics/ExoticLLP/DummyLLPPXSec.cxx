//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/DummyLLPPXSec.h"

using namespace genie;

//____________________________________________________________________________
DummyLLPPXSec::DummyLLPPXSec() :
XSecAlgorithmI("genie::DummyLLPPXSec")
{

}
//____________________________________________________________________________
DummyLLPPXSec::DummyLLPPXSec(string config) :
XSecAlgorithmI("genie::DummyLLPPXSec", config)
{

}
//____________________________________________________________________________
DummyLLPPXSec::~DummyLLPPXSec()
{

}
//____________________________________________________________________________
double DummyLLPPXSec::XSec(const Interaction * , KinePhaseSpace_t ) const
{
  return 0;
}
//____________________________________________________________________________
double DummyLLPPXSec::Integral(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
bool DummyLLPPXSec::ValidProcess(const Interaction * ) const
{
  return true;
}
//____________________________________________________________________________

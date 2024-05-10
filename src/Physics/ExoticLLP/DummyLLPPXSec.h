//____________________________________________________________________________
/*!

\class    genie::DummyLLPPXSec

\brief

\author   

\created  May 05, 2009

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _DUMMYLLP_PXSEC_H_
#define _DUMMYLLP_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class DummyLLPPXSec : public XSecAlgorithmI {

public:
  DummyLLPPXSec();
  DummyLLPPXSec(string config);
 ~DummyLLPPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
};

}       // genie namespace
#endif  //

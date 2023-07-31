//____________________________________________________________________________
/*!

\class    genie::hnl::HNLProductionMode

\brief    Enumeration of HNL production modes.

\author   John Plows <komninos-john.plows \at physics.ox.ac.uk>

\created  May 06, 2022

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_PRODUCTION_MODE_H_
#define _HNL_PRODUCTION_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {
  namespace hnl {
    
    typedef enum t_HNLProd {

      // Using Table 4 from Coloma et al, EPJ C 81 (2021) 78
      // Only considering LFC modes
      
      kHNLProdNull          = -1,
      kHNLProdPion2Muon     = 0, // pi --> N + mu
      kHNLProdPion2Electron = 1, // pi --> N + e
      kHNLProdKaon2Muon     = 2, // K  --> N + mu
      kHNLProdKaon2Electron = 3, // K  --> N + e
      kHNLProdKaon3Muon     = 4, // K  --> N + mu   + pi0
      kHNLProdKaon3Electron = 5, // K  --> N + e    + pi0
      kHNLProdNeuk3Muon     = 6, // K0 --> N + mu   + pi
      kHNLProdNeuk3Electron = 7, // K0 --> N + e    + pi
      kHNLProdMuon3Numu     = 8, // mu --> N + numu + e
      kHNLProdMuon3Nue      = 9, // mu --> N + nue  + e
      //kHNLProdMuon3Nutau    = 10, // mu --> N + nutau + e (LFV!)
      kHNLProdTau2Rho       = 10, // \tau --> N + rho
      kHNLProdTau2Pion      = 11, // \tau --> N + pi
      kHNLProdTau3Numu      = 12, // \tau --> N + numu + mu
      kHNLProdTau3NutauMu   = 13, // \tau --> N + nutau + mu
      kHNLProdTau3Nue       = 14, // \tau --> N + nue + e
      kHNLProdTau3NutauE    = 15, // \tau --> N + nutau + e
      kHNLProdD2Tau         = 16, // D --> N + tau
      kHNLProdD3Muon        = 17, // D --> N + mu + K0
      kHNLProdD3Electron    = 18, // D --> N + e + K0
      kHNLProdD2Muon        = 19, // D --> N + mu
      kHNLProdD2Electron    = 20, // D --> N + e
      kHNLProdDs2Tau        = 21, // D_s --> N + tau
      kHNLProdDs2Muon       = 22, // D_s --> N + mu
      kHNLProdDs2Electron   = 23, // D_s --> N + e
      
    } HNLProd_t;

  } // namespace hnl
} // namespace genie

#endif // #ifndef _HNL_PRODUCTION_MODE_H_

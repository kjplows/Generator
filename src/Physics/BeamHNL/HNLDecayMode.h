//____________________________________________________________________________
/*!

\class    genie::hnl::HNLDecayMode

\brief    Enumeration of HNL decay modes.

\author   John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford

	  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	  University of Liverpool & STFC Rutherford Appleton Laboratory

\created  November 10, 2011

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_DECAY_MODE_H_
#define _HNL_DECAY_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {
  namespace hnl {

 typedef enum EHNLDecayMode {

   // Using list from Bondarenko et al, JHEP 11 (2018) 032, Table 5
   // "only channels with the branching ratio above 1% covering the
   //  HNL mass range up to 5 GeV"

   kHNLDcyNull      = -1,       // dummy
   kHNLDcyNuNuNu    = 0,	// N --> 3 nus. Summed over all flavours. [0 MeV]
   kHNLDcyNuEE      = 1,	// N --> \nu_{a}     e^{\mp}   e^{\pm}. W and Z interfere [1.02 MeV]
   kHNLDcyNuMuE     = 2,	// N --> \nu_{\mu/e} e^{\mp} \mu^{\pm}. Only W. Summed over nue and numu [106 MeV]
   kHNLDcyPi0Nu     = 3,	// N --> \pi^{0}   \nu (any kind) [135 MeV]
   kHNLDcyPiE       = 4,        // N --> \pi^{\pm}   e^{\mp} [140 MeV]
   kHNLDcyNuMuMu    = 5,	// N --> \nu_{a}   \mu^{\mp} \mu^{\pm}. W and Z interfere. [211 MeV]
   kHNLDcyPiMu      = 6,        // N --> \pi^{\pm} \mu^{\mp} [245 MeV]
   kHNLDcyPi0Pi0Nu  = 7,	// N --> \pi^{0}   \pi^{0} \nu (any kind) [268 MeV]
   kHNLDcyPiPi0E    = 8,	// N --> \pi^{\pm} \pi^{0}   e^{\mp} [275 MeV]
   kHNLDcyPiPiNu    = 9,        // N --> \pi^{\pm} + \pi^{\mp} \nu (any kind) [279 MeV]
   kHNLDcyPiPi0Mu   = 10,	// N --> \pi^{\pm} \pi^{0} \mu^{\mp} [380 MeV]
   kHNLDcyNu3Pi0_03 = 11,       // N --> \nu + (3\pi)^{0}, 3 pi^{0} s [405 MeV]
   kHNLDcyNu3Pi0_21 = 12,       // N --> \nu + (3\pi)^{0}, 1 pi^{0} and 1 each of \pi^{\pm} [405 MeV]
   kHNLDcyE3Pi_12   = 13,       // N --> e^{\mp} + (3\pi)^{\pm}, 1 \pi^{\pm} and 2 \pi^{0} [410 MeV]
   kHNLDcyE3Pi_30   = 14,       // N --> e^{\mp} + (3\pi)^{\pm}, 2 \pi^{\pm} and 1 \pi^{\mp} (so 0 \pi^{0}) [410 MeV]
   kHNLDcyMu3Pi_12  = 15,       // N --> \mu^{\mp} + (3\pi)^{\pm}, 1 \pi^{\pm} and 2 \pi^{0} [515 MeV]
   kHNLDcyMu3Pi_30  = 16,       // N --> \mu^{\mp} + (3\pi)^{\pm}, 2 \pi^{\pm} and 1 \pi^{\mp} (so 0 \pi^{0}) [525 MeV]
   kHNLDcyNu4Pi0_04 = 17,       // N --> \nu + (4\pi)^{0}, 4 pi0 [540 MeV]
   kHNLDcyE4Pi_13   = 18,       // N --> e^{\mp} + (4pi)^{\pm}, 3 pi0 [545 MeV]
   kHNLDcyNuEta     = 19,       // N --> \nu + \eta (any kind) [548 MeV]
   kHNLDcyNu4Pi0_22 = 20,       // N --> \nu + (4\pi)^{0}, 2 pi0, 1 each of \pi^{\pm} [550 MeV]
   kHNLDcyE4Pi_31   = 21,       // N --> e^{\mp} + (4\pi)^{\pm}, 1 pi0 [555 MeV]
   kHNLDcyNu4Pi0_40 = 22,       // N --> e^{\mp} + (4\pi)^{0}, 0 pi0 [560 MeV]
   kHNLDcyMu4Pi_13  = 23,       // N --> \mu^{\mp} + (4\pi)^{\pm}, 3 pi0 [649 MeV]
   kHNLDcyMu4Pi_31  = 24,       // N --> \mu^{\mp} + (4\pi)^{\pm}, 1 pi0 [659 MeV]
   kHNLDcyNu5Pi0_05 = 25,       // N --> \nu + (5\pi)^{0}, 5 pi0 [675 MeV]
   kHNLDcyE5Pi_14   = 26,       // N --> e^{\mp} + (5\pi)^{\pm}, 4 pi0 [680 MeV]
   kHNLDcyNu5Pi0_23 = 27,       // N --> \nu + (5\pi)^{0}, 3 pi0 [685 MeV]
   kHNLDcyE5Pi_32   = 28,       // N --> e^{\mp} + (5\pi)^{\pm}, 2 pi0 [690 MeV]
   kHNLDcyNu5Pi0_41 = 29,       // N --> \nu + (5\pi)^{0}, 1 pi0 [695 MeV]
   kHNLDcyE5Pi_50   = 30,       // N --> e^{\mp} + (5\pi)^{\pm}, 0 pi0 [700 MeV]
   kHNLDcyNuomega   = 31,       // N --> \nu (any kind) + \omega [783 MeV]
   kHNLDcyMu5Pi_14  = 32,       // N --> \mu^{\mp} + (5\pi)^{\pm}, 4 pi0 [785 MeV]
   kHNLDcyMu5Pi_32  = 33,       // N --> \mu^{\mp} + (5\pi)^{\pm}, 2 pi0 [795 MeV]
   kHNLDcyMu5Pi_50  = 34,       // N --> \mu^{\mp} + (5\pi)^{\pm}, 0 pi0 [805 MeV]
   kHNLDcyNuEtaP    = 35,       // N --> \nu + \etaPrime [958 MeV]
   kHNLDcyKKNu      = 36,       // N --> K^{\pm} + K^{\mp} + \nu (any kind) [987 MeV]
   kHNLDcyK0K0Nu    = 37,       // N --> K^{0} + K^{0} + \nu (any kind) [994 MeV]
   kHNLDcyNuPhi     = 38,       // N --> \nu + \phi [1019 MeV]
   kHNLDcyNuTauE    = 39,       // N --> \nu_{\tau/e} e^{\mp} \tau^{\pm}. Only W. Summed over nue and nutau [1780 MeV]
   kHNLDcyNuTauMu   = 40,       // N --> \nu_{\tau/\mu} \mu^{\mp} \tau^{\pm}. Only W. Summed over numu and nutau [1880 MeV]
   kHNLDcyEDs       = 41,       // N --> e^{\mp} + D_{s}^{\pm} [1970 MeV]
   kHNLDcyMuDs      = 42,       // N --> \mu^{\mp} + D_{s}^{\pm} [2070 MeV]
   kHNLDcyEDsX      = 43,       // N --> e^{\mp} + D_{s}^{*;\pm} [2110 MeV]
   kHNLDcyMuDsX     = 44,       // N --> \mu^{\mp} + D_{s}^{*;\pm} [2210 MeV]
   kHNLDcyNuEtac    = 45,       // N --> \nu + \eta_{c} [2980 MeV]
   kHNLDcyTEST      = 99, // N --> vv. Test only, not a valid FS.

 } HNLDecayMode_t;

  } // namespace HNL
} // namespace genie
#endif

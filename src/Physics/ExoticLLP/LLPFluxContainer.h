//____________________________________________________________________________
/*!

\class    genie::llp::FluxContainer

\brief    A GENIE flux container specific for LLP.

\author   John Plows

\created  June 24th, 2024

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _LLP_FLUX_CONTAINER_H_
#define _LLP_FLUX_CONTAINER_H_

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <sstream>
#include <cassert>
#include <climits>
#include <algorithm>
#include <iomanip>

#include <TVector3.h>
#include <TLorentzVector.h>

//#include "Framework/EventGen/GFluxI.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/StringUtils.h"

#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TSystem.h>
#include "TRegexp.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/UnitUtils.h"

class TFile;
class TChain;
class TTree;
class TBranch;

using std::string;
using std::ostream;

namespace genie{
  
  namespace llp {
    
    class FluxContainer;
    ostream & operator << (ostream & stream, const FluxContainer & gnmf);
    
    class FluxContainer: public TObject {
      
    public: 
      FluxContainer();
      virtual ~FluxContainer() {};
      
      void ResetCopy() const;
      void Print(const Option_t * /* opt */) const;
      
      friend ostream & operator << (ostream & stream, const FluxContainer & gnmf);

      // members

      mutable int evtno; ///< Index corresponding to flux ntuple entry
 
      mutable TLorentzVector v4;      ///< Production vertex in NEAR frame [m, ns]
      mutable TLorentzVector v4_user; ///< Production vertex in USER frame [m, ns]

      mutable TLorentzVector p4_parent;       ///< Parent momentum in NEAR frame, LLP production [GeV]
      mutable TLorentzVector p4_parent_user ; ///< Parent momentum in USER frame, LLP production [GeV]

      mutable TLorentzVector entry;      ///< Entry point of ray into detector volume, NEAR frame [m, ns]
      mutable TLorentzVector entry_user; ///< Entry point of ray into detector volume, USER frame [m, ns]
      
      mutable TLorentzVector exit;       ///< Exit point of ray from detector volume, NEAR frame [m, ns]
      mutable TLorentzVector exit_user;  ///< Exit point of ray from detector volume, USER frame [m, ns]

      mutable TLorentzVector p4;      ///< LLP momentum in NEAR frame [GeV]
      mutable TLorentzVector p4_user; ///< LLP momentum in USER frame [GeV]
      
      mutable double wgt_xy;
      
    }; // class FluxContainer
    
  } // namespace llp
  
} // namespace genie

#endif // #ifndef _LLP_FLUX_CONTAINER_H_

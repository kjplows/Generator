//----------------------------------------------------------------------------
/*!

  The underlying concrete instance of an exotic LLP.

\namespace  genie::llp

\class      genie::llp::ExoticLLP

\brief      An LLP with definite mass and branching ratios.

\author     John Plows <kjplows \at liverpool.ac.uk>

\created    May 10, 2024

\cpright    Copyright (c) 2003-2024, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _LLP_EXOTIC_LLP_H_
#define _LLP_EXOTIC_LLP_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <TStreamerElement.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

typedef std::vector< std::pair< double, std::vector<int> > > ModeVector;

namespace genie {

  namespace llp {

    struct ModeObject {

      int fIndex;                //! Index of the mode object
      int fMode;                 //! 0 = production ( vis --> LLP ), 1 = decay ( LLP --> vis )
      std::string fName;         //! Unique mode identifier
      double fScore;             //! Conditional probability for this mode to occur
      std::vector<int> fPDGList; //! List of the PDG codes. Parent first, daughters later

      bool IsProductionMode() { return ( fMode == 0 ); }
      bool IsDecayMode() { return ( fMode == 1 ); }

      std::string GetName() { return fName; }
      double GetScore() { return fScore; }
      
      std::vector<int> GetPDGList() { return fPDGList; }

    }; // struct ModeObject

    class ExoticLLP {

    public:

      // C'tors
      ExoticLLP();
      ExoticLLP( double mass, ModeVector  modes );
      ExoticLLP( double mass, ModeVector productionModes, ModeVector decayModes );
      /*virtual*/ ~ExoticLLP();

      // Getters
      inline double GetMass() const { return fMass; }
      inline std::vector< ModeObject > GetProductionModes() const { return fProductionModes; }
      inline std::vector< ModeObject > GetDecayModes() const { return fDecayModes; }

      // Inspect the LLP
      void Print(const Option_t * /* opt */) const;
      friend std::ostream & operator << (std::ostream & stream, const ExoticLLP & LLP);

    private:

      double fMass;                                       //! Mass in MeV
      mutable std::vector< genie::llp::ModeObject > fProductionModes; //! The channels that can produce LLP
      mutable std::vector< genie::llp::ModeObject > fDecayModes;      //! The channels that LLP can produce

      void ConstructModes( ModeVector modes ) const;
      ModeObject MakeModeObject( std::pair< double, std::vector< int > > mode ) const;

      ClassDef(ExoticLLP, 1)
    }; // class ExoticLLP

  } // namespace llp

} // namespace genie

#endif // #ifndef _LLP_EXOTIC_LLP_H_

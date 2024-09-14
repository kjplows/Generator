//____________________________________________________________________________
/*!

\class   genie::llp::ChannelCalculatorI

\brief   Pure abstract base class. Defines the ChannelCalculatorI interface
         to be implemented by LLPConfigurator Algorithm.

\author  John Plows <kjplows \at liverpool.ac.uk>
         University of Liverpool

	 based off AxialFormFactorModelI by
	 Aaron Meyer <asmeyer2012 \at uchicago.edu>

	 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	 University of Liverpool & STFC Rutherford Appleton Laboratory

\created May 10, 2024

\cpright Copyright (c) 2003-2024, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _LLP_CHANNEL_CALCULATOR_I_H_
#define _LLP_CHANNEL_CALCULATOR_I_H_

#include "Framework/Algorithm/Algorithm.h"

namespace genie {

  class Registry;

  namespace llp {

    class ChannelCalculatorI : public Algorithm {

    public:

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      virtual void Configure(const Registry & config) = 0;
      virtual void Configure(string config) = 0;

    protected:

      ChannelCalculatorI();
      ChannelCalculatorI(string name);
      ChannelCalculatorI(string name, string config);

      ClassDef(ChannelCalculatorI, 1)
    }; // class ChannelCalculatorI

  } // namespace llp

} // namespace genie

#endif // #ifndef _LLP_CHANNEL_CALCULATOR_I_H_

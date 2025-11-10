//____________________________________________________________________________
/*!

\class   genie::SplinePostProcessorI

\brief   Interface for the SplinePostProcessor.
         Concrete implementations of this interface use the 'Visitor' design
	 to perform an operation on a Spline.

\author  John Plows <kplows \at liverpool.ac.uk>
         University of Liverpool

\created October 29, 2025

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit 
	 http://genie-mc.github.io/copyright.html
*/
//____________________________________________________________________________

#ifndef _SPLINE_POST_PROCESSOR_I_H_
#define _SPLINE_POST_PROCESSOR_I_H_

#include "Framework/Algorithm/Algorithm.h"

namespace genie {

  class Spline;

  class SplinePostProcessorI : public Algorithm {

  public:

    virtual ~SplinePostProcessorI();

    //-- define the SplinePostProcessorI interface

    virtual std::vector<double> ProcessSpline(const std::vector<double> E,
					      const std::vector<double> xsec) const = 0;

  protected:

    SplinePostProcessorI();
    SplinePostProcessorI(string name);
    SplinePostProcessorI(string name, string config);

  }; // class SplinePostProcessorI

} // namespace genie

#endif // #ifndef _SPLINE_POST_PROCESSOR_I_H_

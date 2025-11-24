//____________________________________________________________________________
/*!

\class   genie::SplinePostProcessor

\brief   Post processor class to smooth out a Spline object.
         This acts on the xsec values (not the abscissae, E) of each Spline
	 by implementing a Gaussian kernel of variable width.

	 The width of this kernel (in GeV) is user-configurable, and is linearly interpolated
	 between stops. 
	 See en.wikipedia.org/wiki/Kernel_smoother#Gaussian_kernel_smoother

\author  John Plows <kplows \at liverpool.ac.uk>
         University of Liverpool

\created November 3, 2025

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit 
	 http://genie-mc.github.io/copyright.html
*/
//____________________________________________________________________________

#ifndef _SPLINE_POST_PROCESSOR_H_
#define _SPLINE_POST_PROCESSOR_H_

#include "Framework/Numerical/SplinePostProcessorI.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgId.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <unordered_set>

namespace genie {

  // Add an exception handler for spline post processing
  // make it part of genie::exceptions
  namespace exceptions {
    class SplineProcessingException {

    public:
      SplineProcessingException();
      SplineProcessingException(const SplineProcessingException & exception);
      SplineProcessingException(std::string reason, bool interrupt = false);
      ~SplineProcessingException();

      void Init(void);
      void Copy(const SplineProcessingException & exception);
      void Print(std::ostream & stream) const;

      void SetReason(std::string reason) { fReason = reason; }
      void SetInterrupt(bool interrupt) { fInterrupt = interrupt; }
      string ShowReason  (void) const { return fReason; }
      bool   DoInterrupt (void) const { return fInterrupt; }

      friend std::ostream & operator << ( std::ostream & stream, 
					  const SplineProcessingException & exception );

    private:

      std::string fReason;
      bool        fInterrupt;

    }; // class genie::exceptions::SplineProcessingException
  } // namespace exceptions

  class Spline;

  class SplinePostProcessor: public SplinePostProcessorI {

  public:
    SplinePostProcessor();
    SplinePostProcessor(std::string name);
    SplinePostProcessor(std::string name, std::string config);
    ~SplinePostProcessor();

    //-- implement the SplinePostProcessorI interface
    std::vector<double> ProcessSpline(const std::vector<double> E,
				      const std::vector<double> xsec) const;

    //-- override the Algorithm::Configure methods to load configuration
    //   data to private data members
    void Configure (const Registry & config);
    void Configure (std::string param_set);

    //-- check if an algorithm is handled
    bool IsHandled (std::string alg_name, std::string alg_config);
    bool IsHandled (const Algorithm * alg);

    //-- Check if we've enabled post-processing
    bool UsePostProcessing() const { return fUsePostProcessing; }

    //-- Workhorse method
    Spline * ProcessSpline( const Spline * spl ) const;

  private:

    void Init       (void);
    void LoadConfig (void);

    void SetThresholds( std::vector<double> & stops, std::vector<double> & widths,
			double xmin, double xmax ) const;

    //-- private data members
    std::unordered_set<std::string> fAlgs; ///< Algorithms handled by the SplinePostProcessor. 

    std::vector<double> fStops;  ///< Abscissae for the interpolation of the Gaussian kernel smoother
    std::vector<double> fWidths; ///< Widths of the Gaussian kernel smoother at the stops

    bool fUsePostProcessing; ///< If true, do post processing. If false, do nothing.

  }; // class SplinePostProcessor

} // namespace genie

#endif // #ifndef _SPLINE_POST_PROCESSOR_H_

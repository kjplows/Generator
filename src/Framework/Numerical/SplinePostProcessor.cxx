//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org 

 John Plows <kplows \at liverpool.ac.uk>
 University of Liverpool 
*/
//____________________________________________________________________________

#include "Framework/Numerical/SplinePostProcessor.h"
#include "Framework/Messenger/Messenger.h"

using std::string;
using std::ostream;
using std::vector;

using namespace genie;

//___________________________________________________________________________
namespace genie::exceptions {
  ostream & operator << (ostream & stream, 
			 const SplineProcessingException & exception ) {
    exception.Print(stream);
    return stream;
  }
//___________________________________________________________________________
  SplineProcessingException::SplineProcessingException() {
    this->Init();
  }
//___________________________________________________________________________
  SplineProcessingException::SplineProcessingException( const SplineProcessingException & exception ) {
    this->Copy(exception);
  }
//___________________________________________________________________________
  SplineProcessingException::SplineProcessingException( string reason, bool interrupt ) {
    this->Init();
    this->SetReason(reason);
    this->SetInterrupt(interrupt);
  }
//___________________________________________________________________________
  SplineProcessingException::~SplineProcessingException() {
  
  }
//___________________________________________________________________________
  void SplineProcessingException::Init(void) {
    fReason = "";
    fInterrupt = false;
  }
//___________________________________________________________________________
  void SplineProcessingException::Copy(const SplineProcessingException & exception) {
    fReason    = exception.fReason;
    fInterrupt = exception.fInterrupt;
  }
//___________________________________________________________________________
  void SplineProcessingException::Print(ostream & stream) const {
    stream << "** EXCEPTION in SplinePostProcessor! Reason : " << this->ShowReason();
    if( this->DoInterrupt() ) { 
      stream << " *~*~*~* INTERRUPTING EXECUTION!! *~*~*~* \n";
    }
  }
//___________________________________________________________________________
} // namespace exceptions

using namespace genie::exceptions;

//___________________________________________________________________________
SplinePostProcessor::SplinePostProcessor() : SplinePostProcessorI("genie::SplinePostProcessor") {
  this->Init();
}
//___________________________________________________________________________
SplinePostProcessor::SplinePostProcessor(string config) : 
  SplinePostProcessorI("genie::SplinePostProcessor", config) {
  this->Init();
}
//___________________________________________________________________________
SplinePostProcessor::SplinePostProcessor(string name, string config) :
  SplinePostProcessorI(name, config) {
  this->Init();
}
//___________________________________________________________________________
SplinePostProcessor::~SplinePostProcessor() {
  fStops.clear();
  fWidths.clear();
  fAlgs.clear();
}
//___________________________________________________________________________
void SplinePostProcessor::Init(void) {
  LOG("SplinePostProcessor", pDEBUG) << "Initialising SplinePostProcessor...";
  this->Configure("Default");
}
//___________________________________________________________________________
void SplinePostProcessor::LoadConfig(void) {
  fStops.clear(); fWidths.clear(); fAlgs.clear();
  LOG("SplinePostProcessor", pDEBUG) << "Loading the SplinePostProcessor configuration...";

  vector<string> alg_vect;
  this->GetParamVect( "HandledAlgorithms", alg_vect );
  // Drop these into the set for quick lookup
  for( string alg : alg_vect ) { 
    fAlgs.insert(alg);
    LOG("SplinePostProcessor", pDEBUG) << "Will handle algorithm named " << alg;
  }

  this->GetParamVect( "Stops", fStops );
  this->GetParamVect( "Widths", fWidths );
  if( fStops.size() != fWidths.size() ) {
    // The Algorithm has been misconfigured. It's early, so complain and ask user to fix.
    LOG("SplinePostProcessor", pFATAL) << "Different number of stops and widths! "
				       << " (" << fStops.size() << " vs " << fWidths.size() << ")";
    throw SplineProcessingException("N stops != N widths", true);
  }
  
}
//___________________________________________________________________________
void SplinePostProcessor::Configure(const Registry & config) {
  Algorithm::Configure(config);
  try {
    this->LoadConfig();
  } catch (SplineProcessingException exception) {
    LOG("SplinePostProcessor", pFATAL) << exception;
    exit(1);
  }
}
//___________________________________________________________________________
void SplinePostProcessor::Configure(std::string param_set) {
  Algorithm::Configure(param_set);
  try {
    this->LoadConfig();
  } catch (SplineProcessingException exception) {
    LOG("SplinePostProcessor", pFATAL) << exception;
    exit(1);
  }
}
//___________________________________________________________________________
bool SplinePostProcessor::IsHandled(std::string alg_name, std::string alg_config) {
  std::string alg_key = alg_name + "/" + alg_config;
  return (fAlgs.count(alg_key) != 0);
}
//___________________________________________________________________________
bool SplinePostProcessor::IsHandled(const Algorithm * alg) {
  std::string alg_key = alg->Id().Key();
  return (fAlgs.count(alg_key) != 0);
}
//___________________________________________________________________________
std::vector<double> SplinePostProcessor::ProcessSpline( const std::vector<double> E,
							const std::vector<double> xsec ) const {

  // Are the thresholds sensible?
  std::vector<double> stops;
  std::vector<double> widths;
  try {
    this->SetThresholds( stops, widths, E.front(), E.back() );
  } catch (SplineProcessingException exception) {
    if(exception.DoInterrupt()) {
      LOG("SplinePostProcessor", pFATAL) << exception;
      exit(1);
    } else {
      LOG("SplinePostProcessor", pERROR) << exception;
    }
  } // Set thresholds

  // Read out the knots from the spline and put them into arrays
  const int NKnots = E.size();
  std::vector<double> new_xsec;

  // Apply Gaussian kernel smoother here, with a linear activation function for widths
  std::vector<double>::iterator it_stop  = stops.begin()+1;
  std::vector<double>::iterator it_width = widths.begin()+1;

  double w1 = widths.front(); double w2 = *it_width;
  for( int i = 0; i < NKnots; i++ ) {
    // Linear activation function between s1 and s2
    auto activation = [&](const double s1, const double s2) {
      if(E[i] <= s1) return 0.0;
      if(E[i] >= s2) return 1.0;
      return (E[i] - s1)/(s2 - s1);
    };

    if(E[i] > *it_stop && E[i] < stops.back()) { // Update iterators
      it_stop++; it_width++;
      w1 = w2; w2 = *(it_width);
    }
    
    double rat = activation(*(it_stop-1), *(it_stop));
    double width = w1 + (w2 - w1) * rat; // between w1 and w2
    double num = 0.0; double den = 0.0;
    for(int j = 0; j < NKnots; j++) {
      double w = std::exp(-0.5 * std::pow((E[i] - E[j]) / width, 2.0)); // Gauss
      num += w * xsec[j]; // weight * value
      den += w; // total weight
    }

    // Knot value out
    new_xsec.emplace_back( num / den );
    
  } // Loop over knots of spline; smooth each one by configurable width

  return new_xsec;
}
//___________________________________________________________________________
void SplinePostProcessor::SetThresholds( std::vector<double> & stops,
					 std::vector<double> & widths,
					 double xmin, double xmax ) const {
  // Copy fStops and fWidths to the temp objects, and then start working on them
  stops = fStops; widths = fWidths;
  if( stops.front() > xmax || stops.back() < xmin ) { // Bad configuration, this will be nonsense.
    throw SplineProcessingException("Thresholds out of bounds", true);
    return;
  }
  
  int n_rejected_front = 0; int n_rejected_back = 0;
  while( stops.front() < xmin ) { // Throw out elements until we're in the range
    stops.erase(stops.begin());
    widths.erase(widths.begin());
    n_rejected_front++;
  }
  while( stops.back() > xmax ) { // Throw out elements until we're in the range
    stops.erase(stops.end()-1);
    widths.erase(widths.end()-1);
    n_rejected_back++;
  }

  if( stops.size() == 1 ) {  // Too few elements. Stop execution
    throw SplineProcessingException("Not enough elements left to interpolate", true);
  }

  // If the xmin and/or xmax are outside first / last stop, add some with copies of width
  if( xmin < stops.front() ) {
    stops.insert(stops.begin(), xmin);
    widths.insert(widths.begin(), widths.front());
  }
  if( xmax > stops.back() ) {
    stops.emplace_back(xmax);
    widths.emplace_back(widths.back());
  }
  
  if(n_rejected_front > 0 || n_rejected_back > 0) { // Complain, but don't stop execution
    LOG("SplinePostProcessor", pERROR) << "Removed " << n_rejected_front << " elements from the front"
				       << " of the stops / widths vectors, and " << n_rejected_back
				       << " from the back. Take note!";
    throw SplineProcessingException("Using different stops to those configured");
  }
  return;
}

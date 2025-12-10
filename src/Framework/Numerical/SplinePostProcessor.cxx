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
SplinePostProcessor * SplinePostProcessor::fInstance = 0;
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
  fStopMap.clear();
  fWidthMap.clear();
  fAlgs.clear();
}
//___________________________________________________________________________
SplinePostProcessor * SplinePostProcessor::Instance()
{
  if(fInstance == 0) fInstance = new SplinePostProcessor;
  return fInstance;
}
//___________________________________________________________________________
void SplinePostProcessor::Init(void) {
  if(fAlgs.size() > 0) return;
  LOG("SplinePostProcessor", pDEBUG) << "Initialising SplinePostProcessor...";
  this->Configure("Default");
}
//___________________________________________________________________________
void SplinePostProcessor::LoadConfig(void) {
  fStops.clear(); fWidths.clear(); fAlgs.clear(); fStopMap.clear(); fWidthMap.clear();
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
  fStopMap[0]  = fStops;
  fWidthMap[0] = fWidths;

  // If there are any more keys, typically delineated with "@Pdg=...", add these to the map
  // Similar architecture to NuclearModelMap's implementation
  RgIMap entries = (this->GetConfig()).GetItemMap();
  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it){
    const std::string  key  = it->first; // Notice this prepends "N" and appends "s"
    const std::string  sens = "@Pdg=";
    if( key.find(sens.c_str()) != std::string::npos
	&& strcmp(key.substr(0, 1).c_str(), "N") == 0 ) {
      // Figure out the size of this thing
      int n;
      LOG("SplinePostProcessor", pDEBUG) << "Getting N stops from here " << key;
      this->GetParam( key, n );
      
      // Decide if this is stops or widths
      // See https://stackoverflow.com/questions/14265581
      std::string comm_name(key);
      comm_name.erase(0, 1);
      comm_name.erase(comm_name.size() - 1);
      std::string dup_name(comm_name); dup_name.append(",");
      LOG("SplinePostProcessor", pDEBUG) << "Key = " << comm_name << " with size "
					 << comm_name.size();
      if( key.find("Stops") != std::string::npos ) {
	std::string subs = dup_name.substr(10); // 5 "Stops" + 5 "@Pdg="
	std::vector<int> pdgs;
	size_t pos = 0;
	while( (pos = subs.find(",", 0)) != std::string::npos ) {
	  int pdg = std::stoi(subs.substr(0, pos));
	  LOG("SplinePostProcessor", pDEBUG) << "Adding key " << pdg << " to stops";
	  subs.erase(0, pos+1); // "," has length 1
	  std::vector<double> stops; stops.resize(n);
	  for( int i = 0; i < n; ++i ) {
	    RgKey tmp_key = Algorithm::BuildParamVectKey(comm_name, i);
	    this->GetParam( tmp_key, stops[i] );
	  }
	  fStopMap[pdg] = stops;
	} // split strings
	  
      } else if ( key.find("Widths") != std::string::npos ) {
	std::string subs = dup_name.substr(11); // 6 "Widths" + 5 "@Pdg="
	std::vector<int> pdgs;
	size_t pos = 0;
	while( (pos = subs.find(",", 0)) != std::string::npos ) {
	  int pdg = std::stoi(subs.substr(0, pos));
	  LOG("SplinePostProcessor", pDEBUG) << "Adding key " << pdg << " to widths";
	  subs.erase(0, pos+1); // "," has length 1
	  std::vector<double> widths; widths.resize(n);
	  for( int i = 0; i < n; ++i ) {
	    RgKey tmp_key = Algorithm::BuildParamVectKey(comm_name, i);
	    this->GetParam( tmp_key, widths[i] );
	  }
	  fWidthMap[pdg] = widths;
	} // split strings
      } // stops or widths
    } // extract PDG entries
  } // loop over registry entries of this Alg

  // Check for misconfigurations
  for ( const auto & mapkey : fStopMap ) {
    if( fStopMap[mapkey.first].size() != fWidthMap[mapkey.first].size() ) {
      // The Algorithm has been misconfigured. It's early, so complain and ask user to fix.
      LOG("SplinePostProcessor", pFATAL) << "Different number of stops and widths! "
					 << " (" << fStopMap[mapkey.first].size()
					 << " vs " << fWidthMap[mapkey.first].size() << ")"
					 << " at PDG key = " << mapkey.first
					 << " (key 0 means default)";
      throw SplineProcessingException("N stops != N widths", true);
    } // misconfiguration check
  } // over each key

  this->GetParam( "UsePostProcessing", fUsePostProcessing );
  
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
							const std::vector<double> xsec,
							const int pdg ) const {

  // If we don't want post processing, return out now.
  if( !fUsePostProcessing ) {
    return xsec;
  } // no op if no post processing

  // Set the stops and widths to be what we asked for
  int elem = (fStopMap.find(pdg) != fStopMap.end()) ? pdg : 0;
  this->SetStops(fStopMap.at(elem)); this->SetWidths(fWidthMap.at(elem));

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
      LOG("SplinePostProcessor", pWARN) << exception;
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
    LOG("SplinePostProcessor", pWARN)  << "Removed " << n_rejected_front << " elements from the front"
				       << " of the stops / widths vectors, and " << n_rejected_back
				       << " from the back. Take note!";
    throw SplineProcessingException("Using different stops to those configured");
  }
  return;
}

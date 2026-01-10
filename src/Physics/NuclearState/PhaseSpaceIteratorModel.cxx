 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 John Plows <kplows \at liverpool.ac.uk>
 University of Liverpool

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/PhaseSpaceIteratorModel.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________
PhaseSpaceIteratorModel::PhaseSpaceIteratorModel() :
  NuclearModelI("genie::PhaseSpaceIteratorModel"), fX(0), fY(0) {}
//____________________________________________________________________________
PhaseSpaceIteratorModel::PhaseSpaceIteratorModel(string config) :
  NuclearModelI("genie::PhaseSpaceIteratorModel", config), fX(0), fY(0) {}
//____________________________________________________________________________
void PhaseSpaceIteratorModel::Configure(const Registry & config) {
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhaseSpaceIteratorModel::Configure(string param_set) {
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhaseSpaceIteratorModel::LoadConfig(void)
{

  NuclearModelI::LoadConfig(); // sets fKFTable inherited from NuclearModelI. See "FermiMomentumTable"
  
  // Set the minimum and maximum (pN, Eb) and number of bins
  std::vector<double> pN_limits; std::vector<double> Eb_limits;
  int pN_bins; int Eb_bins;
  this->GetParamVect( "NucleonMomentumRange", pN_limits );
  this->GetParamVect( "BindingEnergyRange", Eb_limits );
  this->GetParamDef( "NucleonMomentumBins", pN_bins, 100 );
  this->GetParamDef( "BindingEnergyBins", Eb_bins, 100 );

  fPhaseSpace = std::make_unique<TH2D>( "fPhaseSpace", "Phase space: (pN, Eb);pN [GeV/c];Eb [MeV]",
					pN_bins, pN_limits[0], pN_limits[1],
					Eb_bins, Eb_limits[0], Eb_limits[1] );
}
//____________________________________________________________________________
// Next() iterates along X first, then along Y, and stops when both fX and fY hit NBins
// (i.e. the top right corner of phase space)
bool PhaseSpaceIteratorModel::Next() const {
  if( !fPhaseSpace ) return false;

  const int nx = fPhaseSpace->GetNbinsX();
  const int ny = fPhaseSpace->GetNbinsY();

  if( fX == 0 && fY == 0 ) { fX = 1; fY = 1; return true; } // Start
  if( fX < nx ) { fX++; return true; } // Move along row (X)
  if( fY < ny ) { fX = 1; fY++; return true; } // Move along column (Y) and restart X
  return false; // Finished!
}
//____________________________________________________________________________
double PhaseSpaceIteratorModel::Prob(const double p, const double w, const Target & /* tgt */) const {
  // apply a uniform on the phase space
  if( !fPhaseSpace ) return 0.;
  if( fPhaseSpace->FindBin( p, w ) < 0 ) return 0.; // outside bin range
  
  double xmin = fPhaseSpace->GetXaxis()->GetXmin();
  double xmax = fPhaseSpace->GetXaxis()->GetXmax();
  double ymin = fPhaseSpace->GetYaxis()->GetXmin();
  double ymax = fPhaseSpace->GetYaxis()->GetXmax();

  double area = (xmax - xmin) * (ymax - ymin);
  // If the area is zero, complain now, because that suggests the axes are misconfigured
  if( area == 0.0 ) {
    LOG( "PhaseSpaceIterator", pERROR ) << "Returning zero probability due to misconfigured axes:"
					<< " (pNmin, pNmax) x (Ebmin, Ebmax) = ("
					<< xmin << ", " << xmax << ") x ("
					<< ymin << ", " << ymax << "). Please recheck your configuration of"
					<< " the PhaseSpaceIteratorModel!";
    return 0.0;
  }
  // maybe remove the probability of nucleons outside the bin
  if( p < fPhaseSpace->GetXaxis()->GetBinLowEdge(fX) || p > fPhaseSpace->GetXaxis()->GetBinUpEdge(fX) ) return 0.0;
  if( w < fPhaseSpace->GetYaxis()->GetBinLowEdge(fY) || w > fPhaseSpace->GetYaxis()->GetBinUpEdge(fY) ) return 0.0;

  return (area > 0.) ? 1.0 / area : 0.0;
}
//____________________________________________________________________________
bool PhaseSpaceIteratorModel::GenerateNucleon(const double pN, const double Eb,
					      const Target & /*tgt*/) const {
  // We explicitly ask the nucleon to be generated pointing along +z, and then use
  // NuclearModelI::SetMomentum3(const TVector3 &) to modify.
  fCurrRemovalEnergy = Eb;
  fCurrMomentum = TVector3(0.0, 0.0, pN);
  return true;
}
//____________________________________________________________________________

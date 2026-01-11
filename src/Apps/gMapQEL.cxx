//____________________________________________________________________________
/*!

\program gqelmap

\brief   GENIE utility program constructing the cross section of the QEL model
         as a function of struck nucleon momentum and binding energy.
	 Accepts a single neutrino energy.
	 Constructs one TH2D with two axes:
	   x --> struck nucleon p [GeV/c, lab frame]
	   y --> nuclear binding energy [MeV]

         Syntax :

	   gqelmap -p nupdg
	           -t tgtpdg
		   [-o output root file]
		   -E Enu
		   [-a, --angular-approximation-order]
		      Order of the approximation for angular dependence of the xsec

                   // command line args handled by RunOpt:
                   [--tune tune_name]
		   [--message-thresholds your_messenger.xml]
		   [--seed random_number_seed]

		   -h Show syntax and exit

         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.

         Options :
           -p
               A nu PDG code.
	   -t
	       A nuclear target PDG code.
           -o
               Path to the output ROOT file.
	       Default: `qel_map.root`
	   -E  
	       Neutrino energy [GeV]
	   -a, --angular-approximation-order
	       What order of approximation to use (Gauss-Legendre quadrature) for defining
	       the dependence on the polar angle of the struck nucleon.
	       Default: 10

           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
              [default: "G18_02a_00_000"].
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
	   --seed
	      User-specified random number seed.

        ***  See the User Manual for more details and examples. ***

\author  John Plows <kplows \at liverpool.ac.uk>
 University of Liverpool

\created January 3rd, 2026

\cpright Copyright (c) 2003-2026, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
//#include <sstream>
#include <chrono>

#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
#include <fenv.h> // for `feenableexcept`
#endif

#include <TSystem.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include "TMath.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/EventGeneratorList.h"
#include "Framework/EventGen/EventGeneratorListAssembler.h"
#include "Framework/EventGen/InteractionGeneratorMap.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/InteractionSelectorI.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/PhaseSpaceIteratorModel.h"
#include "Physics/QuasiElastic/XSection/NewQELXSec.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

#ifdef __GENIE_HDF5_ENABLED__
// Having this guard early on in the file is inappropriate.
// This needs to be processed after GBuild tells GENIE about flags.
#include "H5Cpp.h"
using namespace H5;
#endif

using std::string;
using std::vector;
using namespace std::chrono ;

using namespace genie;
using namespace genie::utils;
using namespace genie::utils::gsl;
using namespace genie::utils::math;

// Prototypes
void          GetCommandLineArgs (int argc, char ** argv);
void          PrintSyntax        (void);
#ifdef __GENIE_HDF5_ENABLED__
void          WriteToHDF5        (const TH2D & h, double Enu, const std::string fname);
#endif
// Reimplementing ComputeFullQELPXSec for absolute control over *exactly* what's happening.
// Note Eb is *not a reference*. Also, we always assume the nucleon is bound.
// Because genie::utils::BindHitNucleon() isn't called, we just pass the pN and Eb.
double        ComputeGridQELPXSec(Interaction * interaction,
				  const PhaseSpaceIteratorModel * nucl_model,
				  const XSecAlgorithmI * xsec_model,
				  double cth0, double ph0, 
				  const double pN, const double Eb,
				  GaussLegQuad approximant);
// Helpers
TVector3 COMframe2Lab(const genie::InitialState & initialState);

// Global variables
// -- probe PDG
int gOptProbePDG = 0;
// -- target PDG
int gOptTargetPDG = 0;
// -- output file
string kDefOutMapFile = "qel_map.root"; string gOptOutMapFile = kDefOutMapFile;
// -- Neutrino energy
double gOptEnu = -1.;
// -- Order of the Gauss-Legendre angular approximation
int kDefApproxOrder = 10; int gOptApproxOrder = kDefApproxOrder;
// -- random number seed
long int kDefRanSeed = -1; long int gOptRanSeed = kDefRanSeed;

TDatabasePDG* tb = TDatabasePDG::Instance();

// Class definition for our own implementation of FullQELdXSec, I want more control
// Its name will be: GridQELdXSec
//____________________________________________________________________________
//____________________________________________________________________________
class GridQELdXSec : public ROOT::Math::IBaseFunctionMultiDim {
public:
  
  GridQELdXSec(const XSecAlgorithmI* xsec_model, const Interaction* interaction);
  virtual ~GridQELdXSec();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int NDim(void) const {return 2;}
  double DoEval(const double* xin) const;
  ROOT::Math::IBaseFunctionMultiDim* Clone(void) const {
    return new GridQELdXSec(fXSecModel, fInteraction);
  }
  
  Interaction* GetInteractionPtr() {return fInteraction;}
  const Interaction& GetInteraction() const {return *fInteraction;}

  void SetEb(double Eb) {fEb = Eb;}
  void SetpN(double pN) {fpN = pN;}
  void SetApproximant(GaussLegQuad approx) {fApproximant = approx;}

private:
  
  const XSecAlgorithmI* fXSecModel;
  const PhaseSpaceIteratorModel* fNuclModel;
  Interaction* fInteraction;
  
  GaussLegQuad fApproximant;

  double fEb;
  double fpN;
}; // class GridQELXSec
// END of definition

// Implementation of GridQELXSec

//____________________________________________________________________________
//____________________________________________________________________________
GridQELdXSec::GridQELdXSec(const XSecAlgorithmI * xsec_model, const Interaction * interaction) :
  fXSecModel(xsec_model), fInteraction(new Interaction(*interaction)), fEb(0.0), fpN(0.0) {
  // Return the PhaseSpaceIteratorModel
  fNuclModel = dynamic_cast<
    const PhaseSpaceIteratorModel*>( (fXSecModel->SubAlg("IntegralNuclearModel"))
				     ->SubAlg("NuclearModel"));
  assert(fNuclModel);
}
//____________________________________________________________________________
GridQELdXSec::~GridQELdXSec() {
  delete fInteraction;
}
//____________________________________________________________________________
double GridQELdXSec::DoEval(const double * xin) const {
  /* 
   * Let xin be as in FullQELdXSec. xin[0] --> cos_theta0, xin[1] --> phi_theta0
   * where (theta, phi) are the polar, azimuth of the angle between the COM frame
   * 3-momentum of the outgoing lepton, and the COM frame velocity measured in lab frame
   */
  //return fEb;

  double cth0 = xin[0];
  double ph0 = xin[1];

  double ev = ComputeGridQELPXSec(fInteraction, fNuclModel, fXSecModel, cth0, ph0, fpN, fEb,
				  fApproximant);
  return ev;
}

//____________________________________________________________________________
// Helper
TVector3 COMframe2Lab(const genie::InitialState& initialState) {
  TLorentzVector k4 = *(initialState.GetProbeP4( genie::kRfLab ) );
  TLorentzVector p4 = *(initialState.TgtPtr()->HitNucP4Ptr());
  TLorentzVector totMom = k4 + p4;

  TVector3 beta = totMom.BoostVector();

  return beta;
}
//____________________________________________________________________________

// Our own implementation of ComputeFullQELdXSec, including bits from BindHitNucleon()
double ComputeGridQELPXSec(Interaction * interaction,
			   const PhaseSpaceIteratorModel * nucl_model,
			   const XSecAlgorithmI * xsec_model,
			   double cth0, double ph0, 
			   const double pN, const double Eb,
			   GaussLegQuad approximant) {
  
  // First, ensure we've done everything BindHitNucleon() would have done for a bound nucleon
  Target * tgt = interaction->InitState().TgtPtr();
  TLorentzVector * p4Ni = tgt->HitNucP4Ptr();

  // First, set the nucleon to have the appropriate momentum in the nucleus
  // Don't set 4-momentum yet as the nucleon is off shell
  TVector3 p3N_ref( 0.0, 0.0, pN );
  
  // Look up the (on-shell) mass of the initial nucleon
  double mNi = tb->GetParticle( tgt->HitNucPdg() )->Mass(); 

  // Set the (possibly off-shell) initial nucleon energy based on
  // the selected binding energy mode. Always put the initial nucleon
  // on shell if it is not part of a composite nucleus
  double ENi = 0.;

  // For a nuclear target with a bound initial struck nucleon, take binding
  // energy effects and Pauli blocking into account when computing QE
  // differential cross sections
  interaction->ResetBit( kIAssumeFreeNucleon );

  // Initial nucleus mass
  double Mi = tgt->Mass();

  // Final nucleus mass using GENIE convention for the Bodek/Ritchie Fermi gas model
  // RETHERE: Check L312 of QELUtils.cxx for how SpectralFunc does this calc differently!
  double Mf = Mi + Eb - mNi;

  // The (lab-frame) off-shell initial nucleon energy is the difference
  // between the lab frame total energies of the initial and remnant nuclei
  ENi = Mi - std::sqrt( Mf*Mf + pN*pN );
  if( ENi < 0.0 ) { return 0.0; }

  // Update the initial nucleon lab-frame 4-momentum
  p4Ni->SetVect( p3N_ref );
  p4Ni->SetE( ENi );

  // END section from BindHitNucleon()

  // A very high-momentum bound nucleon (which is far off the mass shell)
  // can have a momentum greater than its total energy. This leads to numerical
  // issues (NaNs) since the invariant mass of the nucleon becomes imaginary.
  // In such cases, just return zero to avoid trouble.
  if ( interaction->InitState().Tgt().HitNucP4().M() <= 0. ) return 0.;

  // Mass of the outgoing lepton
  double lepMass = interaction->FSPrimLepton()->Mass();

  // Look up the (on-shell) mass of the final nucleon
  double mNf = tb->GetParticle( interaction->RecoilNucleonPdg() )->Mass();

  // In here, we will have to do one evaluation for each choice of the nucleon cos_theta,
  // i.e. the direction of p4Ni.Vect(). This will influence the kinematics
  // Sample from the approximant
  auto const & weights = approximant.weights;
  auto const & cthetas = approximant.nodes;

  double xsec = 0.;
  for( int i = 0; i < weights.size() ; i++ ) {

    double xsec_cth = 0.;
    // Due to phi invariance about nucleon azimuthal angle we just pick +'ve root for sin theta
    double wgt = weights[i]; double cth = cthetas[i]; double sth = std::sqrt(1.0 - cth*cth);
    TVector3 p3Ni( 0.0, pN * sth, pN * cth );
    p4Ni->SetVect( p3Ni );

    // Mandelstam s for the probe/hit nucleon system
    double s = std::pow( interaction->InitState().CMEnergy(), 2 );

    // Return a differential cross section of zero if we're below threshold (and
    // therefore need to sample a new event)
    if ( std::sqrt(s) < lepMass + mNf ) return 0.;
    
    double outLeptonEnergy = ( s - mNf*mNf + lepMass*lepMass ) / (2 * std::sqrt(s));
    
    if (outLeptonEnergy*outLeptonEnergy - lepMass*lepMass < 0.) return 0.;
    double outMomentum = TMath::Sqrt(outLeptonEnergy*outLeptonEnergy - lepMass*lepMass);

    // Compute the boost vector for moving from the COM frame to the
    // lab frame, i.e., the velocity of the COM frame as measured
    // in the lab frame.
    TVector3 beta = COMframe2Lab( interaction->InitState() );

    // FullDifferentialXSec depends on theta_0 and ph0, the lepton COM
    // frame angles with respect to the direction of the COM frame velocity
    // as measured in the lab frame. To generate the correct dependence
    // here, first set the lepton COM frame angles with respect to +z
    // (via TVector3::SetTheta() and TVector3::SetPhi()).
    TVector3 lepton3Mom(0., 0., outMomentum);
    lepton3Mom.SetTheta( TMath::ACos(cth0) );
    lepton3Mom.SetPhi( ph0 );

    // Then rotate the lepton 3-momentum so that the old +z direction now
    // points along the COM frame velocity (beta)
    TVector3 zvec(0., 0., 1.);
    TVector3 rot = ( zvec.Cross(beta) ).Unit();
    double angle = beta.Angle( zvec );

    // Handle the edge case where beta is along -z, so the
    // cross product above vanishes
    if ( beta.Perp() == 0. && beta.Z() < 0. ) {
      rot = TVector3(0., 1., 0.);
      angle = genie::constants::kPi;
    }

    // Rotate if the rotation vector is not 0
    if ( rot.Mag() >= genie::controls::kASmallNum ) {
      lepton3Mom.Rotate(angle, rot);
    }

    // Construct the lepton 4-momentum in the COM frame
    TLorentzVector lepton(lepton3Mom, outLeptonEnergy);

    // The final state nucleon will have an equal and opposite 3-momentum
    // in the COM frame and will be on the mass shell
    TLorentzVector outNucleon(-1*lepton.Px(),-1*lepton.Py(),-1*lepton.Pz(),
			      TMath::Sqrt(outMomentum*outMomentum + mNf*mNf));
    
    // Boost the 4-momenta for both particles into the lab frame
    lepton.Boost(beta);
    outNucleon.Boost(beta);

    TLorentzVector  nuP4 = *(interaction->InitState().GetProbeP4( genie::kRfLab ));
    TLorentzVector qP4 = nuP4 - lepton;
    double Q2 = -1 * qP4.Mag2();

    interaction->KinePtr()->SetFSLeptonP4( lepton );
    interaction->KinePtr()->SetHadSystP4( outNucleon );
    interaction->KinePtr()->SetQ2( Q2 );

    // Check the Q2 range. If we're outside of it, don't bother
    // with the rest of the calculation.
    Range1D_t Q2lim = interaction->PhaseSpace().Q2Lim();
    if (Q2 < Q2lim.min || Q2 > Q2lim.max) return 0.;

    // Compute the QE cross section for the current kinematics
    xsec_cth = xsec_model->XSec(interaction, genie::kPSQELEvGen);
    // update the xsec
    xsec += wgt * xsec_cth;
  }
 
  return xsec;
}

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gqelmap", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // throw on NaNs and Infs...
#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  // Init
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);

  GaussLegendreQuadrature * GLQuadLib = GaussLegendreQuadrature::Instance();
  GaussLegQuad approximant = GLQuadLib->GetGLQuad(gOptApproxOrder);

  // Load in the QEL model from the tune
  // Because we explicitly iterate over the entire (pN, Eb) space, we don't want a nuclear
  // model to constrain us. So we have to build what GEVGDriver does from scratch...
  LOG("gqelmap", pNOTICE)
       << "Setting event generator list: " << RunOpt::Instance()->EventGeneratorList();
  EventGeneratorListAssembler evglist_assembler(RunOpt::Instance()->EventGeneratorList());
  EventGeneratorList * evglist = evglist_assembler.AssembleGeneratorList();
  InteractionGeneratorMap * intgen_map = new InteractionGeneratorMap;
  intgen_map->UseGeneratorList(evglist);

  // Configure the interaction
  intgen_map->BuildMap(InitialState(gOptTargetPDG, gOptProbePDG)); 

  // Now we can get the AlgorithmFactory involved. The interaction selector points to config!
  AlgFactory * algf = AlgFactory::Instance();
  InteractionSelectorI * int_selector = dynamic_cast<
    InteractionSelectorI *>( algf->AdoptAlgorithm("genie::PhysInteractionSelector", "Default") );

  // We now work similarly to GEVGDriver::CreateSplines().
  EventGeneratorList::const_iterator evgliter = evglist->begin();

  // Loop over all EventGenerator objects used in the current job
  for( ; evgliter != evglist->end() ; ++evgliter ) {
    // current event generator
    const EventGeneratorI * evgen = *evgliter;
    
    /*
    // ask the event generator to produce a list of all interaction it can
    // generate for the input initial state
    const InteractionListGeneratorI * ilstgen = evgen->IntListGenerator();
    InteractionList * ilst = ilstgen->CreateInteractionList(free_init);
    if(!ilst) continue;
    */

    // total cross section algorithm used by the current EventGenerator
    const XSecAlgorithmI * alg = evgen->CrossSectionAlg();
    
    // get the energy range from the EventGenerator validity context
    double Emin = evgen->ValidityContext().Emin() ;
    double Emax = std::min(evgen->ValidityContext().Emax(), gOptEnu) ;

    // We now have an algorithm we would normally ask for its Integral().
    // Check it really is QEL
    std::string alg_key = alg->Id().Key();
    if( alg_key.find("QEL") == std::string::npos ) { continue; }
    if( alg_key.find("CC") == std::string::npos ) {
      // Okay. This is a different QEL cross section, which I will not deal with for now.
      // But let the user know this is being skipped.
      LOG("gqelmap", pNOTICE) << "Skipping map for algorithm with key " << alg_key 
			      << " (not yet implemented)";
    }
    LOG("gqelmap", pNOTICE) << "Performing map for algorithm with key " << alg_key;
    
    /*
    // Load up the integrator for the model
    const XSecIntegratorI * integrator = dynamic_cast<
      const XSecIntegratorI*>( alg->SubAlg("XSec-Integrator") );
    LOG("gqelmap", pNOTICE) << "The xsec integrator is " << integrator->Id().Key();
    */

    // Now, access the nuclear model
    // Note, this is done by going "via" the NuclearModelMap <-- IntegralNuclearModel
    const NuclearModelI * nucl_model = dynamic_cast<
      const NuclearModelI*>( (alg->SubAlg("IntegralNuclearModel"))->SubAlg("NuclearModel") );
    //LOG("gqelmap", pNOTICE) << "The nuclear model is " << nucl_model->Id().Key();
    // If the model is not the PhaseSpaceIterator, complain and exit
    std::string nucl_model_name = nucl_model->Id().Name();
    LOG("gqelmap", pFATAL) << "name is " << nucl_model_name;
    if( strcmp( nucl_model_name.c_str(), "genie::PhaseSpaceIteratorModel" ) != 0 ) {
      LOG("gqelmap", pFATAL) << "You have configured a nuclear model that is not the "
			     << "PhaseSpaceIteratorModel. Please switch to that model to run this"
			     << " utility. Exiting now to save you time.";
      return 2;
    }
    const PhaseSpaceIteratorModel * ps_model = dynamic_cast<
				const PhaseSpaceIteratorModel*>( nucl_model );

    // Integration internals I don't want to define over and over again
    ROOT::Math::IntegrationMultiDim::Type ig_type = 
      utils::gsl::IntegrationNDimTypeFromString( "adaptive" );
    // For the purposes of this exercise, we always use the nuclear model
    std::string bind_mode_str = "UseNuclearModel";
    QELEvGen_BindingMode_t bind_mode = StringToQELBindingMode( bind_mode_str );

    // Integration ranges for the lepton COM frame scattering angles (in the
    // kPSQELEvGen phase space, these are measured with respect to the COM
    // velocity as observed in the lab frame)
    Range1D_t cth0_lim( -1., 1. );
    Range1D_t ph0_lim( 0., 2.*constants::kPi );
    
    double kine_min[2] = { cth0_lim.min, ph0_lim.min };
    double kine_max[2] = { cth0_lim.max, ph0_lim.max };

    double abstol = 1e-16; // We mostly care about relative tolerance. Hardcoded to 1.0e-3, 100 evals
    double reltol = 1e-3;
    int maxeval = 100;
    
    // Get all the interactions that can be generated by this driver (n or p, depending on probe PDG)
    const InteractionList & ilst = intgen_map->GetInteractionList();
    // Loop over all interactions & compute cross sections
    InteractionList::const_iterator intliter = ilst.begin();
    for( ; intliter != ilst.end(); ++intliter ) {
      // get current interaction. Raw pointer, ugh.
      Interaction * interaction = new Interaction(**intliter);

      double Ethr = interaction->PhaseSpace().Threshold();
      LOG("gqelmap", pNOTICE) << "Threshold = " << Ethr << " GeV";

      if (Ethr>Emax) {
	SLOG("gqelmap", pFATAL) << "Energy threshold higher than requested energy.";
	SLOG("gqelmap", pFATAL) << "Energy threshold = " << Ethr << " GeV";
	SLOG("gqelmap", pFATAL) << "Energy requested = " << Emax << " GeV";
	return 4;
      }

      // Get a reference to the histogram of the model. 
      // We may not modify the reference in this App, the model takes care of that.
      const TH2D & phase_space = ps_model->PhaseSpace();
      LOG("gqelmap", pFATAL) << "phase_space has " << phase_space.GetNbinsX() << " x "
			     << phase_space.GetNbinsY() << " bins";
      
      TLorentzVector nup4(0., 0., gOptEnu, gOptEnu);
      interaction->InitStatePtr()->SetProbeP4(nup4);
      
      Target* tgt = interaction->InitState().TgtPtr();      
      // Construct the integrating function.
      
      // Using raw pointer here, that's what ROOT wants
      // last 0. is min angle for EM scattering, normally configurable but I'll hard-code it
      // also, I'll hard-code the integration type to be adaptive
      GridQELdXSec * func = 
	new GridQELdXSec( alg, interaction );
      func->SetApproximant( approximant );
      
      // Switch to using the copy of the interaction in the integrator rather than
      // the copy that we made in this function
      delete interaction;
      interaction = func->GetInteractionPtr();
      
      // Also update the pointer to the Target
      tgt = interaction->InitState().TgtPtr();
      
      // The actual integrating function.
      ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, reltol, maxeval);
      
      // Now setup the custom loop.
      /*
       * While the nuclear model still has (pN, Eb) to look over,
       * generate a nucleon at that (pN, Eb) bin (centre of bin). (Note Eb needs to be in GeV)
       * Then make as many evaluations as the order of your angular approximant,
       * given that you need to evaluate at costh = <that abscissa>.
       * 
       * Once done, sum up the xsec from each iteration of the nuclear target and give the answer
       */
      
      steady_clock::time_point start = steady_clock::now();
	
      int iy_prev = -1;
      while( ps_model->Next() ) {
	
	// sample from the centre of the bin
	int ix = ps_model->X(); int iy = ps_model->Y();
	double pN = phase_space.GetXaxis()->GetBinCenter(ix);
	double Eb = phase_space.GetYaxis()->GetBinCenter(iy);
	Eb *= units::MeV / units::GeV; // to GeV
	func->SetpN(pN); func->SetEb(Eb);
	if( iy != iy_prev ) {
	  LOG("gqelmap", pDEBUG) << "iy --> " << iy << " , Eb = " << 1000. * Eb << " MeV" ;
	  iy_prev = iy;
	}
	

	//LOG("gqelmap", pDEBUG) << "trying pN = " << pN << " GeV/c, Eb = "
	//		       << Eb << " GeV";
	// all the machinery is dealt with internally by GridQELXSec
	double xsec = ig.Integral(kine_min, kine_max);
	
	//LOG("gqelmap", pINFO) << "total xsec = " << xsec;
	// error?
	ps_model->SetXSec(xsec, 0.);
      } // loop over (pN, Eb) bins
      
      delete func;
      //ps_model->Reset();
      
      steady_clock::time_point end = steady_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(end - start);
      LOG("gqelmap", pNOTICE) << "Evaluated map at Enu = " << gOptEnu << " GeV in " 
			      << time_span.count() << " s";
      
      // Write this out as a ROOT file
      TFile fout(gOptOutMapFile.c_str(), "RECREATE");
      phase_space.Write();
      fout.Close();

      // and if we have HDF5 write an HDF5 too
#ifdef __GENIE_HDF5_ENABLED__
      std::string hdf5_output = gOptOutMapFile.substr(0, gOptOutMapFile.find_last_of("."));
      hdf5_output.append(".h5"); // replace ".root" with ".h5"
      WriteToHDF5(phase_space, gOptEnu, hdf5_output);
#endif
      
    } // loop over all interactions
  } // loop over all EventGenerator objects
  
  // Manually take care of cleaning up all these pointers
  if( evglist ){ delete evglist; evglist = 0; }
  if( intgen_map ){ delete intgen_map; intgen_map = 0; }
  if( int_selector ){ delete int_selector; int_selector = 0; }

  return 0;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gqelmap", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gqelmap -p nupdg"
    << "\n    -t tgtpdg"
    << "\n    [-o output_file.root]"
    << "\n    -E Enu"
    << "\n    [-a, --angular-approximation-order n]"
    << "\n    [--tune tune_name]"
    << "\n    [--message-thresholds your-messenger.xml]"
    << "\n    [--seed random_number_seed]"
    << RunOpt::RunOptSyntaxString(false)
    << "\n";

}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gqelmap", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->EnableBareXSecPreCalc(true);
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  if( parser.OptionExists('h') ) {
    PrintSyntax();
    exit(0);
  }

  // output XML file name
  if( parser.OptionExists('o') ) {
    LOG("gqelmap", pINFO) << "Reading output filename";
    gOptOutMapFile = parser.ArgAsString('o');
  } else { LOG("gqelmap", pINFO) << "Unspecified filename - Using default"; }

  // Neutrino PDG
  if( parser.OptionExists('p') ) {
    LOG("gqelmap", pINFO) << "Reading probe PDG";
    gOptProbePDG = parser.ArgAsInt('p');
  } else {
    LOG("gqelmap", pFATAL) << "Please specify the probe PDG code. Exiting now.";
    PrintSyntax();
    exit(2);
  }

  // Target PDG
  if( parser.OptionExists('t') ) {
    LOG("gqelmap", pINFO) << "Reading target PDG";
    gOptTargetPDG = parser.ArgAsInt('t');
  } else {
    LOG("gqelmap", pFATAL) << "Please specify the target PDG code. Exiting now.";
    PrintSyntax();
    exit(2);
  }

  // Max Enu
  if( parser.OptionExists('E') ) {
    LOG("gqelmap", pINFO) << "Reading neutrino energy";
    gOptEnu = parser.ArgAsDouble('E');
  } else {
    LOG("gqelmap", pFATAL) << "Please specify the neutrino energy. Exiting now.";
    PrintSyntax();
    exit(2);
  }

  // Gauss-Legendre approximant order
  if( parser.OptionExists('a') || parser.OptionExists("angular-approximation-order") ) {
    LOG("qelmap", pINFO) << "Reading angular approximation order";
    gOptApproxOrder = (parser.OptionExists('a')) ? 
      parser.ArgAsInt('a') : parser.ArgAsInt("angular-approximation-order");
  } else { LOG("qelmap", pINFO) << "Unspecified angular approximation order - Using default"; }

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gqelmap", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else { LOG("gqelmap", pINFO) << "Unspecified random number seed - Using default"; }

  //
  // print the command-line options
  //
  LOG("gqelmap", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gqelmap job configuration")
     << "\n Neutrino PDG code : " << gOptProbePDG
     << "\n Target PDG code : " << gOptTargetPDG
     << "\n Output map file : " << gOptOutMapFile
     << "\n Neutrino energy : " << gOptEnu << " GeV"
     << "\n Order of approximation for nucleon angular direction: " << gOptApproxOrder
     << "\n Random number seed : " << gOptRanSeed
     << "\n";

  LOG("gqelmap", pNOTICE) << *RunOpt::Instance();
}
//____________________________________________________________________________
// Writes a TH2D of (nx, ny) bins into a dataset, "xsec", of shape (1, nx, ny)
// and an "enu" dataset of shape (1,)
// TODO: add metadata showing the QEL model, nucleus, target...
#ifdef __GENIE_HDF5_ENABLED__
void WriteToHDF5(const TH2D & h, double Enu, const std::string fname) {
  const int nx = h.GetNbinsX();
  const int ny = h.GetNbinsY();

  // flatten into row-major contiguous buffer
  std::vector<double> buffer(nx*ny, 0.0);
  for( int ix = 1; ix <= nx; ix++ ) {
    for( int iy = 1; iy <= ny; iy++ ) {
      buffer[(ix-1)*ny + (iy-1)] = h.GetBinContent(ix, iy);
    }
  }

  H5File fhout(fname, H5F_ACC_TRUNC);

  hsize_t dims[3] = {1, static_cast<hsize_t>(nx), static_cast<hsize_t>(ny)};
  DataSpace space(3, dims);
  
  DataSet dset = fhout.createDataSet("xsec", PredType::NATIVE_DOUBLE, space);
  dset.write(buffer.data(), PredType::NATIVE_DOUBLE);

  // and an enu dataset
  hsize_t edims[1] = {1};
  DataSpace espace(1, edims);

  DataSet eset = fhout.createDataSet("Enu", PredType::NATIVE_DOUBLE, espace);
  eset.write(&Enu, PredType::NATIVE_DOUBLE);
}
#endif

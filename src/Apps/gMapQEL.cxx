//____________________________________________________________________________
/*!

\program gqelmap

\brief   GENIE utility program constructing the cross section of the QEL model
         as a function of struck nucleon momentum and binding energy.
	 Constructs one TH3D with three axes:
	   x --> struck nucleon p [GeV/c, lab frame]
	   y --> nuclear binding energy [MeV]
	   z --> probe energy [GeV]

         Syntax :

	   gqelmap -p nupdg
		   [-o output root file]
		   [-e, --probe-E-max    max Enu ]
		       [--probe-E-bins   N bins  ]
		   [-q, --nucleon-p-max  max pN  ]
		       [--nucleon-p-bins N bins  ]
		   [-b, --binding-E-max  max Eb  ]
		       [--binding-E-bins N bins  ]

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
           -o
               Path to the output ROOT file.
	       Default: `qel_map.root`
	   -e, --probe-E-max
	       Maximum neutrino energy [GeV].
	       Default: 10 GeV
           --probe-E-bins
	       How many bins in neutrino energy to generate.
	       Default: 100 bins
	   -q, --nucleon-p-max
	       Maximum struck nucleon momentum [lab frame, GeV/c].
	       Default: 2 GeV/c
	   --nucleon-p-bins
	       How many bins in struck nucleon momentum to generate.
	       Default: 400 bins
	   -b, --binding-E-max
	       Maximum nuclear binding energy [MeV].
	       Default: 500 MeV
	   --binding-E-bins
	       How many bins in binding energy to generate.
	       Default: 200 bins
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
#include <sstream>

#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
#include <fenv.h> // for `feenableexcept`
#endif

#include <TSystem.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"

using std::string;
using std::vector;

using namespace genie;
using namespace genie::utils::math;

// Prototypes
void          GetCommandLineArgs (int argc, char ** argv);
void          PrintSyntax        (void);

// Global variables
// -- probe PDG
int gOptProbePDG = 0;
// -- output file
string kDefOutMapFile = "qel_map.root"; string gOptOutMapFile = kDefOutMapFile;
// -- probe Emax
double kDefProbeEmax = 10.0; double gOptProbeEmax = kDefProbeEmax;
int kDefProbeEmaxNbins = 100; int gOptProbeEmaxNbins = kDefProbeEmaxNbins;
// -- struck nucleon pmax
double kDefNucleonPmax = 2.0; double gOptNucleonPmax = kDefNucleonPmax;
int kDefNucleonPmaxNbins = 400; int gOptNucleonPmaxNbins = kDefNucleonPmaxNbins;
// -- binding Emax. Note I use GeV here, and will convert to/from MeV on I/O.
double kDefBindingEmax = 0.5; double gOptBindingEmax = kDefBindingEmax;
int kDefBindingEmaxNbins = 200; int gOptBindingEmaxNbins = kDefBindingEmaxNbins;
// -- Order of the Gauss-Legendre angular approximation
int kDefApproxOrder = 10; int gOptApproxOrder = kDefApproxOrder;
// -- random number seed
long int kDefRanSeed = -1; long int gOptRanSeed = kDefRanSeed;

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

  std::ostringstream asts;
  asts << "Here is my approximant: "
       << "\nOrder: " << approximant.n
       << "\nNodes: { ";
  for( const double & node : approximant.nodes ) { asts << node << " "; }
  asts << "}\nWeights: { ";
  for( const double & wgt : approximant.weights ) { asts << wgt << " "; }
  asts << "}\nError coefficient: " << approximant.err;
  LOG("gqelmap", pFATAL) << asts.str();

  return 0;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gqelmap", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gqelmap -p nupdg"
    << "\n    [-o output_file.root]"
    << "\n    [-e, --probe-E-max max_Enu (GeV)]"
    << "\n    [--probe-E-bins N_Enu_bins]"
    << "\n    [-q, --nucleon-p-max max_pN (GeV/c)]"
    << "\n    [--nucleon-p-bins N_pN_bins]"
    << "\n    [-b, --binding-E-max max_Eb (MeV)]"
    << "\n    [--binding-E-bins N_Eb_bins]"
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

  // probe Emax
  if( parser.OptionExists('e') || parser.OptionExists("probe-E-max") ) {
    LOG("qelmap", pINFO) << "Reading maximum neutrino energy";
    gOptProbeEmax = (parser.OptionExists('e')) ? 
      parser.ArgAsDouble('e') : parser.ArgAsDouble("probe-E-max");
  } else { LOG("qelmap", pINFO) << "Unspecified maximum neutrino energy - Using default"; }

  // probe E Nbins
  if( parser.OptionExists("probe-E-bins") ) {
    LOG("qelmap", pINFO) << "Reading neutrino energy N bins";
    gOptProbeEmaxNbins = parser.ArgAsInt("probe-E-bins");
  } else { LOG("qelmap", pINFO) << "Unspecified neutrino energy N bins - Using default"; }

  // struck nucleon Pmax
  if( parser.OptionExists('q') || parser.OptionExists("nucleon-p-max") ) {
    LOG("qelmap", pINFO) << "Reading maximum struck nucleon momentum";
    gOptNucleonPmax = (parser.OptionExists('q')) ? 
      parser.ArgAsDouble('q') : parser.ArgAsDouble("nucleon-p-max");
  } else { LOG("qelmap", pINFO) << "Unspecified maximum struck nucleon momentum - Using default"; }

  // struck nucleon P Nbins
  if( parser.OptionExists("nucleon-p-bins") ) {
    LOG("qelmap", pINFO) << "Reading struck nucleon momentum N bins";
    gOptNucleonPmaxNbins = parser.ArgAsInt("probe-p-bins");
  } else { LOG("qelmap", pINFO) << "Unspecified struck nucleon momentum N bins - Using default"; }

  // binding Emax
  if( parser.OptionExists('b') || parser.OptionExists("binding-E-max") ) {
    LOG("qelmap", pINFO) << "Reading maximum binding energy";
    gOptBindingEmax = (parser.OptionExists('b')) ? 
      parser.ArgAsDouble('b') : parser.ArgAsDouble("binding-E-max");
  } else { LOG("qelmap", pINFO) << "Unspecified maximum binding energy - Using default"; }

  // binding E Nbins
  if( parser.OptionExists("binding-E-bins") ) {
    LOG("qelmap", pINFO) << "Reading binding energy N bins";
    gOptBindingEmaxNbins = parser.ArgAsInt("binding-E-bins");
  } else { LOG("qelmap", pINFO) << "Unspecified binding energy N bins - Using default"; }

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
     << "\n Output map file : " << gOptOutMapFile
     << "\n Maximum neutrino energy [GeV] : " << gOptProbeEmax
     << "\n N bins in neutrino energy : " << gOptProbeEmaxNbins
     << "\n Maximum struck nucleon momentum [GeV/c] : " << gOptNucleonPmax
     << "\n N bins in struck nucleon momentum : " << gOptNucleonPmaxNbins
     << "\n Maximum binding energy [MeV] : " << gOptBindingEmax * units::GeV / units::MeV
     << "\n N bins in binding energy : " << gOptBindingEmaxNbins
     << "\n Order of approximation for nucleon angular direction: " << gOptApproxOrder
     << "\n Random number seed : " << gOptRanSeed
     << "\n";

  LOG("gqelmap", pNOTICE) << *RunOpt::Instance();
}

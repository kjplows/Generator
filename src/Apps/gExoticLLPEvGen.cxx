//________________________________________________________________________________________
/*!

\author  John Plows <kjplows \at liverpool.ac.uk>
         University of Liverpool

\created May 8, 2024

\cpright Copyright (c) 2003-2024, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>

#include <TSystem.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/TuneId.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/CmdLnArgParser.h"

#include "Physics/ExoticLLP/ExoticLLP.h"
#include "Physics/ExoticLLP/LLPConfigurator.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;
using namespace genie::llp;

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX__
//#include "Tools/Flux/GNuMIFlux.h"
#include <TH1.h>
#endif // #ifdef __GENIE_FLUX_DRIVERS_ENABLED__

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#define __CAN_USE_ROOT_GEOM__
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

// function prototypes
void   GetCommandLineArgs (int argc, char ** argv);
void   PrintSyntax        (void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
//void     FillFluxNonsense        (FluxContainer &ggn);
//void     FillFlux                (FluxContainer &ggn, FluxContainer &tgn);
#endif // #ifdef __GENIE_FLUX_DRIVERS_ENABLED__

//TLorentzVector GeneratePosition( GHepRecord * event );
#ifdef __CAN_USE_ROOT_GEOM__
void  InitBoundingBox    (void);
#endif // #ifdef __CAN_USE_ROOT_GEOM__

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format
string          kDefOptEvFilePrefix = "gntp";
string          kDefOptFluxFilePath = "./input-flux.root";

string          gOptEvFilePrefix = kDefOptEvFilePrefix; // event file prefix

string          kDefOptSName   = "genie::EventGenerator";
string          kDefOptSConfig = "ExoticLLP"; // just copy from the HNL config
string          kDefOptSTune   = "GLLP24_00a_00_000";

//
Long_t           gOptRunNu        = 1000;                // run number
int              gOptNev          = 10;                  // number of events to generate

double           gOptEnergyLLP    = -1;                  // LLP energy [ GeV ]
double           gOptMassLLP      = -1;                  // LLP mass [ MeV ]

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
string           gOptFluxFilePath = kDefOptFluxFilePath; // where flux files live
map<string,string> gOptFluxShortNames;
bool             gOptIsUsingTrees = false;               // using flat-tree flux files?
int              gOptFirstEvent   = 0;                  // skip to this entry in flux tree
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

bool             gOptUsingRootGeom = false;              // using root geom or target mix?
string           gOptRootGeom;                           // input ROOT file with realistic detector geometry

#ifdef __CAN_USE_ROOT_GEOM__
TGeoManager *    gOptRootGeoManager = 0;                 // the workhorse geometry manager
TGeoVolume  *    gOptRootGeoVolume  = 0;
#endif // #ifdef __CAN_USE_ROOT_GEOM__

string           gOptRootGeomTopVol = "";                // input geometry top event generation volume
double           gOptGeomLUnits = 0;                     // input geometry length units
long int         gOptRanSeed = -1;                       // random number seed

// Geometry bounding box and origin - read from the input geometry file (if any)
double fdx = 0; // half-length - x
double fdy = 0; // half-length - y
double fdz = 0; // half-length - z
double fox = 0; // origin - x
double foy = 0; // origin - y
double foz = 0; // origin - z

// vector of the allowed decay modes
struct LLPDecayMode {
  int idx; //! index
  std::vector< int > pdg_code_list; //! list of PDG codes produced
  double score; //! conditional probability of this mode
};
std::vector< LLPDecayMode > allowed_decay_modes;

//_________________________________________________________________________________________
int main(int argc, char ** argv)
{
  // suppress ROOT Info messages
  gErrorIgnoreLevel = kWarning;

  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  // Init messenger and random number seed
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);

  __attribute__((unused)) RandomGen * rnd = RandomGen::Instance();

  const Algorithm * algLLPConfigurator = AlgFactory::Instance()->GetAlgorithm("genie::llp::LLPConfigurator", "Default");

  const LLPConfigurator * LLP_configurator = dynamic_cast< const LLPConfigurator * >( algLLPConfigurator );

  // Let's ensure that the LLP Configurator plays ball. Ask it what the LLP mass is
  ExoticLLP llp = LLP_configurator->RetrieveLLP();

  LOG( "gevgen_exotic_llp", pFATAL ) << llp;

  LOG( "gevgen_exotic_llp", pFATAL )
    << "This is a TEST. Goodbye world!";
    
  return 0;
}
//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_exotic_llp", pINFO) << "Parsing command line arguments";

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
    PrintSyntax();
    exit(0);
  }
  
  // force the app to look at LLP tune by default
  // if user passes --tune argument, look at the user input tune instead
  char * expargv[ argc + 2 ];
  bool didFindUserInputTune = false;
  std::string stExtraTuneBit = kDefOptSTune;
  
  if( parser.OptionExists("tune") ){
    didFindUserInputTune = true;
    stExtraTuneBit = parser.ArgAsString("tune");
    LOG( "gevgen_exotic_llp", pWARN )
      << "Using input LLP tune " << parser.ArgAsString("tune");
  } else {
    LOG( "gevgen_exotic_llp", pWARN )
      << "Using default LLP tune " << kDefOptSTune;
  }
  // append this to argv
  for( int iArg = 0; iArg < argc; iArg++ ){
    expargv[iArg] = argv[iArg];
  }
  if( !didFindUserInputTune ){
    char * chBit = const_cast< char * >( stExtraTuneBit.c_str() ); // ugh. Ugly.
    std::string stune("--tune"); char * tBit = const_cast< char * >( stune.c_str() );
    expargv[argc] = tBit;
    expargv[argc+1] = chBit;

    LOG( "gevgen_exotic_llp", pNOTICE ) << stune << " " << stExtraTuneBit;
  }

  // Common run options.
  int expargc = ( didFindUserInputTune ) ? argc : argc+2;
  std::string stnull(""); char * nBit = const_cast< char * >( stnull.c_str() );
  expargv[expargc] = nBit;

  RunOpt::Instance()->ReadFromCommandLine(expargc,expargv);

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevgen_exotic_llp", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_exotic_llp", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_exotic_llp", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_exotic_llp", pFATAL)
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  // get LLP mass 
  const Algorithm * algLLPConf = AlgFactory::Instance()->GetAlgorithm("genie::llp::LLPConfigurator", "Default");
  const LLPConfigurator * llp_conf = dynamic_cast< const LLPConfigurator * >( algLLPConf );
  gOptMassLLP = llp_conf->RetrieveLLP().GetMass();

  /*
  bool isMonoEnergeticFlux = true;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  if( parser.OptionExists('f') ) {
    LOG("gevgen_exotic_llp", pDEBUG)
      << "A flux has been offered. Searching this path: " << parser.ArgAsString('f');
    isMonoEnergeticFlux = false;
    gOptFluxFilePath = parser.ArgAsString('f');
    
    // check if this is valid path (assume these are dk2nu files)
    if( gSystem->OpenDirectory( gOptFluxFilePath.c_str() ) != NULL ){
      gOptIsUsingDk2nu = true;
      LOG("gevgen_exotic_llp", pDEBUG)
	<< "dk2nu flux files detected. Will create flux spectrum dynamically.";
    } else {
      LOG("gevgen_exotic_llp", pFATAL)
	<< "Invalid flux file path " << gOptFluxFilePath;
      exit(1);
    }
  } else {
    // we need the 'E' option! Log it and pass below
    LOG("gevgen_exotic_llp", pINFO)
      << "No flux file offered. Assuming monoenergetic flux.";
  } //-f
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

  // LLP energy (only relevant if we do not have an input flux)
  gOptEnergyLLP = -1;
  if( isMonoEnergeticFlux ){
    if( parser.OptionExists('E') ) {
      LOG("gevgen_exotic_llp", pDEBUG)
        << "Reading LLP energy";
      gOptEnergyLLP = parser.ArgAsDouble('E');
    } else {
      LOG("gevgen_exotic_llp", pFATAL)
        << "You need to specify the LLP energy";
      PrintSyntax();
      exit(0);
    } //-E
    assert(gOptEnergyLLP > gOptMassLLP);
  }

  gOptIsMonoEnFlux = isMonoEnergeticFlux;

  // first flux entry to read
  if( parser.OptionExists("firstEvent") ) {
    gOptFirstEvent = parser.ArgAsInt("firstEvent");
    LOG( "gevgen_exotic_llp", pINFO )
      << "Starting flux readin at first event = " << gOptFirstEvent;
  } // --firstEvent

  // LLP decay mode
  int mode = -1;
  if( parser.OptionExists('m') ) {
    LOG("gevgen_exotic_llp", pDEBUG)
        << "Reading LLP decay mode";
    mode = parser.ArgAsInt('m');
  } else {
    LOG("gevgen_exotic_llp", pINFO)
        << "No decay mode specified - will sample from allowed decay modes";
  } //-m
  gOptDecayMode = (LLPDecayMode_t) mode;

  bool allowed = utils::hnl::IsKinematicallyAllowed(gOptDecayMode, gOptMassLLP);
  if(!allowed) {
    LOG("gevgen_exotic_llp", pFATAL)
      << "Specified decay is not allowed kinematically for the given LLP mass";
    PrintSyntax();
    exit(0);
  }

  //
  // geometry
  //

  string geom = "";
  string lunits;
#ifdef __CAN_USE_ROOT_GEOM__
  // string dunits;
  if( parser.OptionExists('g') ) {
    LOG("gevgen_exotic_llp", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');

    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file =
            ! (gSystem->AccessPathName(geom.c_str()));
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    } else {
      LOG("gevgen_exotic_llp", pFATAL)
	<< "Geometry option is not a ROOT file. Please use ROOT geom.";
      PrintSyntax();
      exit(1);
    }
  } else {
      // LOG("gevgen_exotic_llp", pFATAL)
      //   << "No geometry option specified - Exiting";
      // PrintSyntax();
      // exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // length units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_exotic_llp", pDEBUG)
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_exotic_llp", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // // density units:
     // if( parser.OptionExists('D') ) {
     //    LOG("gevgen_exotic_llp", pDEBUG)
     //       << "Checking for input geometry density units";
     //    dunits = parser.ArgAsString('D');
     // } else {
     //    LOG("gevgen_exotic_llp", pDEBUG) << "Using default geometry density units";
     //    dunits = kDefOptGeomDUnits;
     // } // -D
     gOptGeomLUnits = utils::units::UnitFromString(lunits);
     // gOptGeomDUnits = utils::units::UnitFromString(dunits);

  } // using root geom?
#endif // #ifdef __CAN_USE_ROOT_GEOM__
  */

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_exotic_llp", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_exotic_llp", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_exotic_llp", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_exotic_llp", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  //
  // >>> print the command line options
  //

  ostringstream gminfo;
  if (gOptUsingRootGeom) {
    gminfo << "Using ROOT geometry - file: " << gOptRootGeom
           << ", top volume: "
           << ((gOptRootGeomTopVol.size()==0) ? "<master volume>" : gOptRootGeomTopVol);
    // << ", length  units: " << lunits;
    // << ", density units: " << dunits;
  } else gminfo << "No ROOT geometry loaded";

  LOG("gevgen_exotic_llp", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_exotic_llp job configuration");

  LOG("gevgen_exotic_llp", pNOTICE)
     << "\n @@ Run number    : " << gOptRunNu
     << "\n @@ Random seed   : " << gOptRanSeed
     << "\n @@ LLP mass      : " << gOptMassLLP << " GeV"
     << "\n @@ Geometry      : " << gminfo.str()
     << "\n @@ Statistics    : " << gOptNev << " events";
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_exotic_llp", pFATAL)
   << "\n **Syntax**"
   << "\n gevgen_exotic_llp [-h] "
   << "\n            [-r run#]"
   << "\n             -n n_of_events"
   << "\n             -f path/to/flux/files"
   << "\n            [--firstEvent first_event_for_flux_readin]"  
   << "\n            [-g geometry (ROOT file)]"
   << "\n            [-L length_units_at_geom]"
   << "\n            [-o output_event_file_prefix]"
   << "\n            [--seed random_number_seed]"
   << "\n            [--message-thresholds xml_file]"
   << "\n            [--event-record-print-level level]"
   << "\n            [--mc-job-status-refresh-rate  rate]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << " or look at the source code: $GENIE/src/Apps/gExoticLLPEvGen.cxx"
   << "\n";
}
//_________________________________________________________________________________________

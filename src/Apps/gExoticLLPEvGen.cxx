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
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
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
#include "Physics/ExoticLLP/Decayer.h"
#include "Physics/ExoticLLP/LLPConfigurator.h"
#include "Physics/ExoticLLP/LLPFluxContainer.h"
#include "Physics/ExoticLLP/LLPFluxCreator.h"
#include "Physics/ExoticLLP/LLPVertexGenerator.h"
#include "Physics/ExoticLLP/VolumeSeeker.h"

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
void     SetFluxBranchAddresses  (TTree * tree);
void     CreateFluxBranches      (void);
void     CopyFromPointerFlux     (void);
#endif // #ifdef __GENIE_FLUX_DRIVERS_ENABLED__

//TLorentzVector GeneratePosition( GHepRecord * event );
#ifdef __CAN_USE_ROOT_GEOM__
void  InitBoundingBox    (void);
TLorentzVector * CastX4ToUser( VolumeSeeker * sk, TLorentzVector vc );
TLorentzVector * CastP4ToUser( VolumeSeeker * sk, TLorentzVector vc );
#endif // #ifdef __CAN_USE_ROOT_GEOM__

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format
string          kDefOptEvFilePrefix = "gntp";
string          kDefOptFluxFilePath = "./input-flux.root";
string          kDefOptTopVolName   = "TOP"; // default top volume name 

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
int              gOptFirstEvent   = 0;                   // skip to this entry in flux tree
FluxContainer    gOptFluxInfo;                           // custom LLP flux driver

// a FluxContainer-like struct to hold addresses to pointers
struct FluxPointerContainer {
  int evtno;
  int pdg;
  TLorentzVector * v4;
  TLorentzVector * v4_user;
  TLorentzVector * p4_parent; 
  TLorentzVector * p4_parent_user;
  TLorentzVector * entry;
  TLorentzVector * entry_user;
  TLorentzVector * exit;
  TLorentzVector * exit_user;
  TLorentzVector * p4;
  TLorentzVector * p4_user;
  double wgt_xy;
  double boost_factor;
  double wgt_collimation;

  FluxPointerContainer() {
    evtno = 0;
    pdg = 0;
    v4 = new TLorentzVector(); v4_user = new TLorentzVector();
    p4_parent = new TLorentzVector(); p4_parent_user = new TLorentzVector();
    entry = new TLorentzVector(); entry_user = new TLorentzVector();
    exit = new TLorentzVector(); exit_user = new TLorentzVector();
    p4 = new TLorentzVector(); p4_user = new TLorentzVector();
    wgt_xy = 0.0;
    boost_factor = 0.0;
    wgt_collimation = 0.0;
  }

  ~FluxPointerContainer() {
    evtno = 0;
    pdg = 0;
    delete v4; delete v4_user;
    delete p4_parent; delete p4_parent_user;
    delete entry; delete entry_user;
    delete exit; delete exit_user;
    delete p4; delete p4_user;
    wgt_xy = 0.0;
    boost_factor = 0.0;
    wgt_collimation = 0.0;
  }
};

FluxPointerContainer gOptFluxPtInfo;

bool             gOptWriteAngularAcceptance = false;     // write an angular acceptance branch to tree
TFile *          gOutputFluxFile = 0;                    // output flux file (clone of input)
TTree *          gOutputFluxTree = 0;                    // output flux tree
AngularRegion    gAngularRegion;                         // Event-by-event angular region, for event generation purposes
double           gAngularRegionSize = 0;                 // Size of this region / 4\pi
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

bool             gOptUsingRootGeom = false;              // using root geom or target mix?
string           gOptRootGeom;                           // input ROOT file with realistic detector geometry

#ifdef __CAN_USE_ROOT_GEOM__
TGeoManager *    gOptRootGeoManager = 0;                 // the workhorse geometry manager
TGeoVolume  *    gOptRootGeoVolume  = 0;
#endif // #ifdef __CAN_USE_ROOT_GEOM__

bool             gOptTopVolSelected = false;             // did the user ask for a specific top volume?
string           gOptTopVolName = kDefOptTopVolName;     // input geometry top event generation volume
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

// LLP lifetime in rest frame
double CoMLifetime = -1.0; // GeV^{-1}

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

  [[maybe_unused]] RandomGen * rnd = RandomGen::Instance();
  VolumeSeeker * vsek = VolumeSeeker::Instance();

  // Load the LLP and sub-algorithms
  const Algorithm * algLLPConfigurator = AlgFactory::Instance()->GetAlgorithm("genie::llp::LLPConfigurator", "Default");

  const LLPConfigurator * LLP_configurator = dynamic_cast< const LLPConfigurator * >( algLLPConfigurator );

  ExoticLLP llp = LLP_configurator->RetrieveLLP();
  // inspect the LLP
  LOG( "gevgen_exotic_llp", pFATAL ) << llp;

  const Algorithm * algFluxCreator = AlgFactory::Instance()->GetAlgorithm("genie::llp::FluxCreator", "Default");
  const Algorithm * algVtxGen = AlgFactory::Instance()->GetAlgorithm("genie::llp::VertexGenerator", "Default");

  const FluxCreator * fluxCreator = dynamic_cast< const FluxCreator * >( algFluxCreator );
  const VertexGenerator * vtxGen = dynamic_cast< const VertexGenerator * >( algVtxGen );

  // check the config is correct
  vsek->PrintConfig();

  // pass the geometry information
  //vsek->SetGeomFile( gOptRootGeom, gOptTopVolName );
  // initialise the elements of VolumeSeeker
  vsek->ClearEvent();

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu, gOptRanSeed);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  LOG("gevgen_exotic_llp", pNOTICE)
    << "Initialised Ntuple Writer";

  // add flux info to the tree
  if( gOptWriteAngularAcceptance ) {
    // fill the flux object with nonsense to start with
    FluxContainer * ptGnmf = new FluxContainer();
    gOptFluxInfo = *ptGnmf;
    delete ptGnmf;
  }

  // Create a MC job monitor for a periodically updated status file
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  LOG("gevgen_hnl", pNOTICE)
    << "Initialised MC job monitor";

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

#ifdef __CAN_USE_ROOT_GEOM__
  // Read geometry bounding box - for vertex position generation
  if( gOptUsingRootGeom ){
    InitBoundingBox();
  }
#endif // #ifdef __CAN_USE_ROOT_GEOM__

  // Event loop
  int iflux = (gOptFirstEvent < 0) ? 0 : gOptFirstEvent; int ievent = iflux;
  int maxFluxEntries = -1;
  //vtxGen->SetGeomFile( gOptRootGeom, gOptTopVolName );
  vsek->SetGeomFile( gOptRootGeom, gOptTopVolName );

  // An extremely, extremely weird bug. What.
  std::cout.flush(); // why? ROOT? what do you do?
  std::cerr.flush(); // why? ROOT? what do you do?
  LOG("gevgen_exotic_llp", pDEBUG) << "FLUSHED COUT AND FLUSHING AGAIN.";
  std::cerr << "FLUSHED COUT AND FLUSHING AGAIN." << std::endl;
  std::cout.flush(); // why? ROOT? what do you do?
  std::cerr.flush(); // why? ROOT? what do you do?

  TChain * flux_tree = 0; 
  if( ! gSystem->AccessPathName( gOptFluxFilePath.c_str() ) ) {
    flux_tree = new TChain("flux"); flux_tree->Add( gOptFluxFilePath.c_str() );

    if( gOptWriteAngularAcceptance ){
      gOutputFluxFile = TFile::Open( "./exotic_llp_output_flux.root", "RECREATE" );
      //gOutputFluxTree = dynamic_cast<TTree *>( flux_tree->CloneTree(0) );
      gOutputFluxTree = new TTree( flux_tree->GetName(), flux_tree->GetName() );

      CreateFluxBranches();

      // Also set a couple of addresses...
      flux_tree->SetBranchAddress( "v4", &(gOptFluxPtInfo.v4) );
      flux_tree->SetBranchAddress( "p4_parent", &(gOptFluxPtInfo.p4_parent) );
      flux_tree->SetBranchAddress( "parent_pdg", &gOptFluxPtInfo.pdg );
    } // if write angular acceptance
    
    else {
      SetFluxBranchAddresses( flux_tree );
    } // if use pre-calculated angular acceptance and flux

  } else {
    LOG( "gevgen_exotic_llp", pFATAL ) << "Could not find input flux file. Exiting.";
    exit(1);
  }

  bool tooManyEntries = false;
  if( gOptNev < 0 ) gOptNev = flux_tree->GetEntries();
  while (1) {
    if( tooManyEntries ){
      if( gOptNev >= 10000 ){
	if( (ievent-gOptFirstEvent) % (gOptNev / 1000) == 0 ){
	  int irat = (iflux-gOptFirstEvent) / ( gOptNev / 1000 );
	  std::cerr << 0.1 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev  << " processed ) \r" << std::flush;
	}
      } else if( gOptNev >= 10 ) {
	if( (ievent-gOptFirstEvent) % (gOptNev / 10) == 0 ){
	  int irat = (iflux-gOptFirstEvent) / ( gOptNev / 10 );
	  std::cerr << 10.0 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev  << " processed ) \r" << std::flush;
	}
      }
    } else {
      if( gOptNev >= 10000 ){
	if( (ievent-gOptFirstEvent) % (gOptNev / 1000) == 0 ){
	  int irat = (ievent-gOptFirstEvent) / ( gOptNev / 1000 );
	  std::cerr << 0.1 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev <<  " processed ) \r" << std::flush;
	}
      } else if( gOptNev >= 10 ) {
	if( (ievent-gOptFirstEvent) % (gOptNev / 10) == 0 ){
	  int irat = (ievent-gOptFirstEvent) / ( gOptNev / 10 );
	  std::cerr << 10.0 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev  << " processed ) \r" << std::flush;
	}
      }
    }

    LOG( "gevgen_exotic_llp", pDEBUG ) << "ievent = " << ievent
				       << ", gOptFirstEvent = " << gOptFirstEvent
				       << ", gOptNev = " << gOptNev;
    
    if( tooManyEntries && ((iflux-gOptFirstEvent) == gOptNev) ) break;
    else if( (ievent-gOptFirstEvent) == gOptNev ) break;
    
    if( ievent < gOptFirstEvent ){ ievent++; continue; }
    
    assert( ievent >= gOptFirstEvent && gOptFirstEvent >= 0 && "First event >= 0" );

    //LOG("gevgen_exotic_llp", pNOTICE)
    //  << " *** Getting entry for event............ " << (ievent-gOptFirstEvent);

    flux_tree->GetEntry(ievent);

    //LOG("gevgen_exotic_llp", pNOTICE)
    //  << " *** GOT ENTRY for event............ " << (ievent-gOptFirstEvent);

    // If we want to calculate angular acceptance, don't actually process events but just populate branches.

    if( gOptWriteAngularAcceptance ){
      LOG("gevgen_exotic_llp", pNOTICE)
	<< " *** Filling flux for event............ " << (ievent-gOptFirstEvent);
      
      TLorentzVector * v4 = gOptFluxPtInfo.v4;
      TLorentzVector * p4_parent = gOptFluxPtInfo.p4_parent;

      LOG("gevgen_exotic_llp", pNOTICE)
	<< " *** Origin v4 (to populate event) is " << utils::print::X4AsString( v4 );

      TVector3 origin_point = v4->Vect();
      TVector3 momentum = p4_parent->Vect();
      vsek->PopulateEvent( origin_point, momentum );
      //[[maybe_unused]] bool result = vsek->RaytraceDetector();
      gAngularRegion = vsek->AngularAcceptance();
      gAngularRegionSize = vsek->AngularSize( gAngularRegion );

      LOG("gevgen_exotic_llp", pDEBUG) << "Calculating angular acceptance for event " << ievent
				       << " at vz = " << v4->Z()
				       << " as wgt_xy = " << gAngularRegionSize / (4.0 * constants::kPi); 

      //gOptFluxPtInfo = *(new FluxPointerContainer()); NONONONO

      // Here we set all the branches' values. It is inelegant, but has to be done.

      gOptFluxPtInfo.evtno = ievent;

      //TLorentzVector v4( 0.0, -60.0, vz, 0.0 );
      gOptFluxPtInfo.v4_user = CastX4ToUser( vsek, *(gOptFluxPtInfo.v4) );
      
      gOptFluxPtInfo.p4_parent_user = CastP4ToUser( vsek, *p4_parent );

      gOptFluxPtInfo.wgt_xy = gAngularRegionSize / (4.0 * constants::kPi);

      gOutputFluxTree->Fill();

      vsek->ClearEvent();
    } else {
      LOG("gevgen_exotic_llp", pNOTICE)
	<< " *** Generating event............ " << (ievent-gOptFirstEvent);
      
      EventRecord * event = new EventRecord;
      event->SetWeight(1.0);
      event->SetProbability( 1.0 );

      // Assume that there is a flux and that it contains all the information. 
      // So read it from the input file.

      CopyFromPointerFlux(); // because ROOT wants addresses to pointers, which is annoying...

      LOG( "gevgen_exotic_llp", pDEBUG )
	<< "\nTesting: read v4 " << utils::print::X4AsString( &gOptFluxInfo.v4 )
	<< "\nAnd momentum = " << utils::print::P4AsString( &gOptFluxInfo.p4_parent )
	<< "\nwith wgt_xy = " << gOptFluxInfo.wgt_xy
	<< "\nPointers: read v4 " << utils::print::X4AsString( gOptFluxPtInfo.v4 )
	<< "\nAnd momentum = " << utils::print::P4AsString( gOptFluxPtInfo.p4_parent )
	<< "\nwith wgt_xy = " << gOptFluxPtInfo.wgt_xy;

      // First, we need to construct the input LLP. Based on the parent kinematics
      Decayer * decayer = Decayer::Instance();
      decayer->ClearEvent();
      fluxCreator->UpdateFluxInfo( gOptFluxInfo );
      fluxCreator->ProcessEventRecord(event);
      // check if we should move on -- might be the flux was unable to produce an entry.
      if( event->Probability() < 0 ){

	LOG( "gevgen_exotic_llp", pWARN )
	  << "Unable to produce event " << (ievent-gOptFirstEvent);

      } else { // Could produce event from flux, continue on

	gOptFluxInfo  = fluxCreator->RetrieveFluxInfo();
	gOptEnergyLLP = (gOptFluxInfo.p4).E();

	LOG( "gevgen_exotic_llp", pDEBUG ) << "Sanity: boost factor is " << gOptFluxInfo.boost_factor;

	int decay = 0;
	int typeMod = ( gOptFluxInfo.pdg >= 0 ) ? 1 : -1; 
	Interaction * interaction = Interaction::LLP(typeMod * genie::kPdgLLP, gOptEnergyLLP, decay);
	
	event->AttachSummary(interaction);
	
	LOG("gevgen_exotic_llp", pINFO)
	  << "Generated event: " << *event;
	
	// Add event at the output ntuple, refresh the mc job monitor & clean-up
	ntpw.AddEventRecord(ievent, event);
	mcjmonitor.Update(ievent,event);

      }
      
      delete event;
    } // main event generation
    
    ievent++;
    
  } // end event loop

  if( gOptWriteAngularAcceptance ){ 
    gOutputFluxTree->Write(); 
    gOutputFluxFile->Close(); 
  } else { ntpw.Save(); }

  //delete v4; delete p4_parent;

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

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  // Flag for incoming flux files
  if( parser.OptionExists('f') ) {
    LOG("gevgen_exotic_llp", pDEBUG)
      << "A flux has been offered. Searching this path: " << parser.ArgAsString('f');
    gOptFluxFilePath = parser.ArgAsString('f');
    
    // check if this is valid path
    if( 
       ( gSystem->OpenDirectory( gOptFluxFilePath.c_str() ) != NULL ) ||
       ( ! gSystem->AccessPathName( gOptFluxFilePath.c_str() ) ) 
	){
      LOG("gevgen_exotic_llp", pDEBUG)
	<< "Flux files detected. Good!";
    } else {
      LOG("gevgen_exotic_llp", pFATAL)
	<< "Invalid flux file path " << gOptFluxFilePath;
      exit(1);
    }
  } //-f

#ifdef __CAN_USE_ROOT_GEOM__
  // Flag to perform and save angular acceptance calculations
  if( parser.OptionExists("calculate_angular_acceptance") ){
    LOG("gevgen_exotic_llp", pNOTICE)
      << "\nWill calculate angular acceptance for each event in the flux."
      << "\nWill NOT generate events. Use the output of this stage to generate events later.";

    gOptWriteAngularAcceptance = true;
  } // --calculate_angular_acceptance
#endif // #ifdef __CAN_USE_ROOT_GEOM__
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

  /*
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
  */

  //
  // geometry
  //

  std::string geom = "";
  std::string lunits;
#ifdef __CAN_USE_ROOT_GEOM__
  // std::string dunits;
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

     // check for top volume selection
     if( parser.OptionExists("top_volume") ) {
       gOptTopVolSelected = true;
       gOptTopVolName = parser.ArgAsString("top_volume");
       LOG("gevgen_hnl", pINFO)
	 << "Using the following volume as top: " << gOptTopVolName;
     } else {
       LOG("gevgen_hnl", pINFO)
	 << "Using default top_volume";
     } // --top_volume

  } // using root geom?
#endif // #ifdef __CAN_USE_ROOT_GEOM__

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
           << ((gOptTopVolName.size()==0) ? "<master volume>" : gOptTopVolName);
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
void SetFluxBranchAddresses( TTree * tree )
{
  tree->SetBranchAddress( "evtno",            &(gOptFluxPtInfo.evtno)            );
  tree->SetBranchAddress( "parent_pdg",       &(gOptFluxPtInfo.pdg)              );
  tree->SetBranchAddress( "v4",               &(gOptFluxPtInfo.v4)               );
  tree->SetBranchAddress( "v4_user",          &(gOptFluxPtInfo.v4_user)          );
  tree->SetBranchAddress( "p4_parent",        &(gOptFluxPtInfo.p4_parent)        );
  tree->SetBranchAddress( "p4_parent_user",   &(gOptFluxPtInfo.p4_parent_user)   );
  tree->SetBranchAddress( "entry",            &(gOptFluxPtInfo.entry)            );
  tree->SetBranchAddress( "entry_user",       &(gOptFluxPtInfo.entry_user)       );
  tree->SetBranchAddress( "exit",             &(gOptFluxPtInfo.exit)             );
  tree->SetBranchAddress( "exit_user",        &(gOptFluxPtInfo.exit_user)        );
  tree->SetBranchAddress( "p4",               &(gOptFluxPtInfo.p4)               );
  tree->SetBranchAddress( "p4_user",          &(gOptFluxPtInfo.p4_user)          );
  tree->SetBranchAddress( "wgt_xy",           &(gOptFluxPtInfo.wgt_xy)           );
  tree->SetBranchAddress( "boost_factor",     &(gOptFluxPtInfo.boost_factor)     );
  tree->SetBranchAddress( "wgt_collimaution", &(gOptFluxPtInfo.wgt_collimation)  );
}
//_________________________________________________________________________________________
void CreateFluxBranches(void)
{
  gOutputFluxTree->Branch( "evtno",          &(gOptFluxPtInfo.evtno),         "evtno/I"  );
  gOutputFluxTree->Branch( "parent_pdg",     &(gOptFluxPtInfo.pdg),           "parent_pdg/I"  );
  gOutputFluxTree->Branch( "v4",             &(gOptFluxPtInfo.v4)                        );
  gOutputFluxTree->Branch( "v4_user",        &(gOptFluxPtInfo.v4_user)                   );
  gOutputFluxTree->Branch( "p4_parent",      &(gOptFluxPtInfo.p4_parent)                 );
  gOutputFluxTree->Branch( "p4_parent_user", &(gOptFluxPtInfo.p4_parent_user)            );
  gOutputFluxTree->Branch( "entry",          &(gOptFluxPtInfo.entry)                     );
  gOutputFluxTree->Branch( "entry_user",     &(gOptFluxPtInfo.entry_user)                );
  gOutputFluxTree->Branch( "exit",           &(gOptFluxPtInfo.exit)                      );
  gOutputFluxTree->Branch( "exit_user",      &(gOptFluxPtInfo.exit_user)                 );
  gOutputFluxTree->Branch( "p4",             &(gOptFluxPtInfo.p4)                        );
  gOutputFluxTree->Branch( "p4_user",        &(gOptFluxPtInfo.p4_user)                   );
  gOutputFluxTree->Branch( "wgt_xy",         &(gOptFluxPtInfo.wgt_xy),        "wgt_xy/D" );
  gOutputFluxTree->Branch( "boost_factor",   &(gOptFluxPtInfo.boost_factor),  "boost_factor/D" );
  gOutputFluxTree->Branch( "wgt_collimation", &(gOptFluxPtInfo.wgt_collimation), "wgt_collimation/D" );
}
//_________________________________________________________________________________________
void CopyFromPointerFlux(void)
{
  gOptFluxInfo.evtno            = gOptFluxPtInfo.evtno;
  gOptFluxInfo.pdg              = gOptFluxPtInfo.pdg;
  gOptFluxInfo.v4               = *(gOptFluxPtInfo.v4);
  gOptFluxInfo.v4_user          = *(gOptFluxPtInfo.v4_user);
  gOptFluxInfo.p4_parent        = *(gOptFluxPtInfo.p4_parent);
  gOptFluxInfo.p4_parent_user   = *(gOptFluxPtInfo.p4_parent_user);
  gOptFluxInfo.entry            = *(gOptFluxPtInfo.entry);
  gOptFluxInfo.entry_user       = *(gOptFluxPtInfo.entry_user);
  gOptFluxInfo.exit             = *(gOptFluxPtInfo.exit);
  gOptFluxInfo.exit_user        = *(gOptFluxPtInfo.exit_user);
  gOptFluxInfo.wgt_xy           = gOptFluxPtInfo.wgt_xy;
  gOptFluxInfo.boost_factor     = gOptFluxPtInfo.boost_factor;
  gOptFluxInfo.wgt_collimation  = gOptFluxPtInfo.wgt_collimation;
}
//_________________________________________________________________________________________
TLorentzVector * CastX4ToUser( VolumeSeeker * sk, TLorentzVector vc )
{
  TVector3 vc_space = vc.Vect();
  TVector3 vc_user_space = sk->TranslateToUser( vc_space );
  vc_user_space = sk->RotateToUser( vc_user_space );
  return new TLorentzVector( vc_user_space.X(), vc_user_space.Y(), vc_user_space.Z(), vc.T() );
}
//_________________________________________________________________________________________
TLorentzVector * CastP4ToUser( VolumeSeeker * sk, TLorentzVector vc )
{
  TVector3 vc_space = vc.Vect();
  TVector3 vc_user_space = sk->RotateToUser( vc_space );
  return new TLorentzVector( vc_user_space.Px(), vc_user_space.Py(), vc_user_space.Pz(), vc.E() );
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
   << "\n            [--calculate_angular_acceptance]"  
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
//............................................................................
#ifdef __CAN_USE_ROOT_GEOM__
void InitBoundingBox(void)
{
// Initialise geometry bounding box, used for generating HNL vertex positions

  LOG("gevgen_exotic_llp", pINFO)
    << "Initialising geometry bounding box.";

  fdx = 0; // half-length - x
  fdy = 0; // half-length - y
  fdz = 0; // half-length - z
  fox = 0; // origin - x
  foy = 0; // origin - y
  foz = 0; // origin - z

  if(!gOptUsingRootGeom){ // make a unit-m sided box
    LOG("gevgen_exotic_llp", pINFO)
      << "No geometry file input detected, making a unit-m side box volume.";

    TGeoManager * geom = new TGeoManager( "box1", "A simple box detector" );

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

    //--- make the top object volume
    //const double boxSideX = 2.5, boxSideY = 2.5, boxSideZ = 2.5; // m
    //const double bigBoxSide = 2.0 * std::max( boxSideX, std::max( boxSideY, boxSideZ ) ); // m
    //const double worldLen = 1.01 * bigBoxSide; // m

    TGeoVolume * topvol = geom->MakeBox( "TOP", Vacuum, 101.0, 101.0, 101.0 );
    geom->SetTopVolume( topvol );

    //--- make the detector box container
    TGeoVolume * boxvol = geom->MakeBox( "VOL", Vacuum, 100.5, 100.5, 100.5 );
    boxvol->SetVisibility(kFALSE);

    //--- origin is at centre of the box
    TGeoVolume * box = geom->MakeBox( "BOX", Al, 100.0, 100.0, 100.0 );
    //TGeoTranslation * tr0 = new TGeoTranslation( 0.0, 0.0, 0.0 );
    TGeoRotation * rot0 = new TGeoRotation( "rot0", 90.0, 0.0, 90.0, 90.0, 0.0, 0.0 );

    //--- add directly to top volume
    topvol->AddNode( box, 1, rot0 );
    
    gOptRootGeoManager = geom;

    return;
  } 

  bool geom_is_accessible = ! (gSystem->AccessPathName(gOptRootGeom.c_str()));
  if (!geom_is_accessible) {
    LOG("gevgen_exotic_llp", pFATAL)
      << "The specified ROOT geometry doesn't exist! Initialization failed!";
    exit(1);
  }

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 
  if( !gOptTopVolSelected ){
    TGeoVolume * main_volume = gOptRootGeoManager->GetTopVolume();
    gOptTopVolName = main_volume->GetName();
    LOG("gevgen_exotic_llp", pINFO) << "Using top volume name " << gOptTopVolName;
  }

  TGeoVolume * top_volume = gOptRootGeoManager->GetVolume(gOptTopVolName.c_str());
  assert( top_volume && "Top volume exists" );
  TGeoShape * ts  = top_volume->GetShape();

  TGeoBBox *  box = (TGeoBBox *)ts;

  //get box origin and dimensions (in the same units as the geometry)
  fdx = box->GetDX();
  fdy = box->GetDY();
  fdz = box->GetDZ();
  fox = (box->GetOrigin())[0];
  foy = (box->GetOrigin())[1];
  foz = (box->GetOrigin())[2];

  LOG("gevgen_exotic_llp", pDEBUG)
    << "Before conversion the bounding box has:"
    << "\nOrigin = ( " << fox << " , " << foy << " , " << foz << " )"
    << "\nDimensions = " << fdx << " x " << fdy << " x " << fdz
    << "\n1cm = 1.0 unit";

  // Convert from local to SI units
  fdx *= gOptGeomLUnits;
  fdy *= gOptGeomLUnits;
  fdz *= gOptGeomLUnits;
  fox *= gOptGeomLUnits;
  foy *= gOptGeomLUnits;
  foz *= gOptGeomLUnits;

  LOG("gevgen_exotic_llp", pINFO)
    << "Initialised bounding box successfully.";

}
#endif // #ifdef __CAN_USE_ROOT_GEOM__
//............................................................................

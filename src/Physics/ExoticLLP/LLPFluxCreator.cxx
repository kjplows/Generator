//____________________________________________________________________________
/*
  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  
  Author: John Plows <kplows \at liverpool.ac.uk>
          University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/LLPFluxCreator.h"

using namespace genie;
using namespace genie::llp;

//____________________________________________________________________________
FluxCreator::FluxCreator() : FluxRecordVisitorI("genie::llp::FluxCreator")
{

}
//____________________________________________________________________________
FluxCreator::FluxCreator(string name) : FluxRecordVisitorI(name)
{

}
//____________________________________________________________________________
FluxCreator::FluxCreator(string name, string config) : FluxRecordVisitorI(name, config)
{

}
//____________________________________________________________________________
FluxCreator::~FluxCreator()
{
  
}
//____________________________________________________________________________
void FluxCreator::UpdateFluxInfo( FluxContainer info ) const
{
  fFluxInfo = info;
}
//____________________________________________________________________________
FluxContainer FluxCreator::RetrieveFluxInfo() const
{
  return fFluxInfo;
}
//____________________________________________________________________________
void FluxCreator::ProcessEventRecord(GHepRecord * evrec) const
{
  // Adds in the initial state LLP, and nothing else. 
  LOG( "ExoticLLP", pDEBUG ) << "Processing flux info with parent PDG = " << fFluxInfo.pdg;
  
  // First, run some quick calculations on the flux container.
  //this->AddInfoToFlux();
  
  // The calculation in this module will always be in NEAR coordinates
  TLorentzVector parent_p4 = fFluxInfo.p4_parent; // GeV/c
  TLorentzVector parent_v4 = fFluxInfo.v4; // m

  int parent_pdg = fFluxInfo.pdg;

  // The information that goes into the EventRecord will always be in USER coordinates 

  // Get a random point in the top volume to evaluate the flux at
  VolumeSeeker * vsek = VolumeSeeker::Instance();
  TVector3 rand_point = vsek->GetRandomPointInTopVolNEAR();
  LOG( "ExoticLLP", pWARN ) << "Evaluating flux at NEAR point " 
			    << utils::print::Vec3AsString( &rand_point );

  // First, apply the constraint of kinematics.
  TVector3 separation = rand_point - parent_v4.Vect(); // m
  
  // Calculate the opening angle between this ray and the parent momentum
  double cos_zeta = separation.Dot( rand_point ) / ( separation.Mag() * rand_point.Mag() );
  double zeta = std::acos( cos_zeta ) * 180.0 / constants::kPi;

  // Find a decay channel for the LLP.

  // First, make a decay mode. Read in the scores from each decay mode, and get the winning decay mode

  // Poll only from the object of modes that corresponds to the parent!
  std::unordered_map< int, std::vector<ModeObject> >::iterator active_modes =
    fGroupedModes.find( std::abs( parent_pdg ) );
  if( active_modes == fGroupedModes.end() ){
    LOG( "ExoticLLP", pFATAL )
      << "\nNo production modes found for parent with PDG code " << parent_pdg << " !!!"
      << "\nUpdate the config file in $GENIE/config/GLLP* and retry running...";
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("No production mode associated with parent particle");
    exception.SwitchOnFastForward();
    throw exception;
  }

  std::vector< ModeObject > llp_production_modes = (*active_modes).second;

  LOG( "ExoticLLP", pDEBUG ) << "For the parent with PDG code "
			     << parent_pdg << ", there are " << llp_production_modes.size()
			     << " production modes";

  std::vector< ModeObject >::iterator it_modes = llp_production_modes.begin();

  RandomGen * rnd = RandomGen::Instance();
  double production_score = rnd->RndGen().Rndm();

  double score_seen = (*it_modes).GetScore();
  
  while( it_modes != llp_production_modes.end() && score_seen < production_score ) {
    ++it_modes; score_seen += (*it_modes).GetScore(); 
  }
  if( it_modes == llp_production_modes.end() ) --it_modes;
  ModeObject chosen_production = *(it_modes);
  LOG( "ExoticLLP", pDEBUG ) 
    << "With thrown score " << production_score << " we picked the channel with name " 
    << chosen_production.GetName();

  int llp_pdg = kPdgLLP;
  if( fFluxInfo.pdg < 0 ) llp_pdg = kPdgAntiLLP;
  int type_mod = (llp_pdg > 0) ? 1 : -1;

  // Get the PDG product list from the decay mode
  std::vector<int> chosen_PDGVector = chosen_production.GetPDGList();
  // turn this into a PDGCodeList
  // RETHERE: how do we deal with two LLP?
  bool allow_duplicate = true;
  PDGCodeList chosen_PDGList(allow_duplicate);
  for( std::vector<int>::iterator it_vec = chosen_PDGVector.begin();
       it_vec != chosen_PDGVector.end(); ++it_vec ) {
    // first check that typeMod * pdg code exists. If not, it's something like -1 * pi0
    // (which is its own antiparticle), and we drop the typeMod
    if( chosen_PDGList.ExistsInPDGLibrary( type_mod * (*it_vec) ) )
      chosen_PDGList.push_back( type_mod * (*it_vec) );
    else
      chosen_PDGList.push_back( *it_vec );
  }

  LOG( "ExoticLLP", pDEBUG ) << "\nPDGList is " << chosen_PDGList;
  
  // Tell the Decayer instance about this PDGList
  Decayer * decayer = Decayer::Instance();
  decayer->SetProducts( chosen_PDGList );

  // Perform a phase space decay
  // RETHERE: Add in angular correlations of final state products?
  [[maybe_unused]] bool decay_ok = decayer->UnpolarisedDecay();
  
  // Read in the results, these will be REST FRAME
  std::vector< GHepParticle > decayed_results = decayer->GetRestFrameResults();
  GHepParticle parent_particle = *(decayed_results.begin());
  TLorentzVector parent_p4_rest = *(parent_particle.GetP4());
  LOG( "ExoticLLP", pDEBUG ) << "The 4-momentum of the 0th entry is "
			     << utils::print::P4AsString( &parent_p4_rest );

  TLorentzVector probe_p4(0.0, 0.0, 1.0, 1.0);
  TLorentzVector probe_v4(0.0, 0.0, 0.0, 0.0);

  GHepParticle ptLLP( llp_pdg, kIStInitialState, -1, -1, -1, -1, probe_p4, probe_v4 );
  evrec->AddParticle( ptLLP );
}
//____________________________________________________________________________
void FluxCreator::AddInfoToFlux() const
{
  // At this point the flux container already has some information, and needs more.
  // We will populate with the LLP 4-momentum, the boost correction, and the acceptance correction.
  // Throughout this section, we will work in NEAR coordinates.

  this->SetCurrentEntry( fFluxInfo.evtno );

  // The LLP energy is the same along the ray along the detector. Only the angles matter.
  
  TLorentzVector p4par  = fFluxInfo.p4_parent;
  TLorentzVector delta4 = fFluxInfo.entry - fFluxInfo.v4;

  TVector3 p3par = p4par.Vect();
  TVector3 delta = delta4.Vect();

  // Get the opening angle by projecting out the longitudinal delta component
  TVector3 deltaT = delta - (p3par.Unit()).Dot(delta) * p3par.Unit();
  
  double cos_zeta = (p3par.Unit()).Dot((deltaT.Unit()));
  double zeta_rad = std::acos(cos_zeta);
  double zeta     = zeta_rad * 180.0 / constants::kPi;

  LOG( "ExoticLLP", pDEBUG )
    << "\np3par  = " << utils::print::Vec3AsString(&p3par)
    << "\ndelta  = " << utils::print::Vec3AsString(&delta)
    << "\ndeltaT = " << utils::print::Vec3AsString(&deltaT)
    << "\nzeta   = " << zeta << " [deg]";

  // First, make a decay mode. Read in the scores from each decay mode, and get the PDG product list
  std::vector< ModeObject > llp_production_modes = fExoticLLP.GetProductionModes();
  [[maybe_unused]] RandomGen * rnd = RandomGen::Instance();
  double production_score = rnd->RndGen().Rndm();

  double score_seen = 0.0;
  std::vector< ModeObject >::iterator it_modes = llp_production_modes.begin();
  while( it_modes != llp_production_modes.end() && score_seen < production_score ) {
    score_seen += (*it_modes).GetScore(); ++it_modes;
  }
  if( it_modes == llp_production_modes.end() ) --it_modes;
  ModeObject chosen_production = *(it_modes);
  LOG( "ExoticLLP", pDEBUG ) 
    << "With thrown score " << production_score << " we picked the channel with name " 
    << chosen_production.GetName();

  // RETHERE implement
  return;
}
//____________________________________________________________________________
TLorentzVector FluxCreator::LLPEnergy() const
{
  // RETHERE implement by reading in p4_parent from fFluxInfo
  return TLorentzVector( 0.0, 0.0, 0.0, fMass );
}
//____________________________________________________________________________
std::pair< double, double > 
FluxCreator::CalculateAcceptanceCorrection( TLorentzVector p4_LLP, double zeta ) const
{
  // RETHERE implement. First is Lab frame LLP energy, second is accCorr
  return std::pair<double, double>( 0.0, 0.0 );
}
//____________________________________________________________________________
double FluxCreator::labangle( double * x, double * par )
{
  double xrad = x[0] * TMath::DegToRad();
  double Ehad = par[0], pxhad = par[1], pyhad = par[2], pzhad = par[3];
  double pllp = par[4], Ellp = par[5];

  TLorentzVector p4had( pxhad, pyhad, pzhad, Ehad );
  TVector3 boost_vec = p4had.BoostVector(); // beta of parent in lab frame

  // assume phi invariance so create LLP rest-frame momentum along y'z' plane
  // This assumption is valid even in the full 2D treatment, as phi is uniquely determined
  // by a single quaternion rotation
  TLorentzVector pncm( 0.0, pllp * TMath::Sin( xrad ), pllp * TMath::Cos( xrad ), Ellp );

  // boost into lab frame
  pncm.Boost( boost_vec );
  
  // return lab frame theta wrt parent momentum in deg
  double num = pxhad * pncm.X() + pyhad * pncm.Y() + pzhad * pncm.Z();
  double den = p4had.P() * pncm.P();
  double theta = TMath::ACos( num / den ) * 180.0 / constants::kPi;
  return theta;
}
//____________________________________________________________________________
void FluxCreator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FluxCreator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FluxCreator::LoadConfig(void)
{
  if( fIsConfigLoaded ) return;

  LOG( "ExoticLLP", pDEBUG ) << "Loading flux-creation parameters from config...";
  
  this->GetParam( "Mass", fMass );

  this->GetParam( "FluxTree", ctree_name );
  this->GetParam( "MetaTree", cmeta_name );

  // RETHERE add tree names so that we can read into the AliasedBranch objects
  this->GetParam( "FluxBranch_parent_pdg", m_ct_parent_pdg_alias );
  
  this->InitialiseTree(); this->InitialiseMeta();

  // Get the LLP configurator and from that the ExoticLLP
  const Algorithm * algLLPConfigurator = AlgFactory::Instance()->GetAlgorithm("genie::llp::LLPConfigurator", "Default");
  const LLPConfigurator * LLP_configurator = dynamic_cast< const LLPConfigurator * >( algLLPConfigurator );

  fExoticLLP = LLP_configurator->RetrieveLLP();

  // Now get the production modes from the LLP, and group them by parent
  std::vector< ModeObject > llp_production_modes = fExoticLLP.GetProductionModes();
  std::vector< ModeObject >::iterator it_modes = llp_production_modes.begin();
  while( it_modes != llp_production_modes.end() ) {
    ModeObject mobj = *it_modes; 
    int parent_pdg = std::abs( (mobj.GetPDGList()).at(0) ); // treat pi- as pi+ etc
    std::unordered_map< int, 
			std::vector< ModeObject > >::iterator mmap = fGroupedModes.find( parent_pdg );
    if( mmap == fGroupedModes.end() ) {
      std::vector< ModeObject > tmp_mobj = { mobj };
      std::pair< int, std::vector< ModeObject > > tmp_pair = { parent_pdg, tmp_mobj };
      fGroupedModes.insert( tmp_pair );
      LOG( "ExoticLLP", pDEBUG ) << "Inserted new parent with PDG " << parent_pdg;
    } else {
      std::vector< ModeObject > mobj_vec = (*mmap).second;
      mobj_vec.emplace_back( mobj );
      fGroupedModes[ parent_pdg ] = mobj_vec; // update the vector
      LOG( "ExoticLLP", pDEBUG ) << "Added mode to parent with PDG " << parent_pdg;
    }

    //LOG( "ExoticLLP", pDEBUG ) << "Parent with PDG " << parent_pdg << " now has "
    //			       << (*mmap).second.size() << " production modes";
    
    ++it_modes;
  }

  LOG( "ExoticLLP", pDEBUG ) << "For all parents together, there are " << llp_production_modes.size()
			     << " production modes";

  //delete LLP_configurator; LLP_configurator = 0;
  //delete algLLPConfigurator; algLLPConfigurator = 0;

  fIsConfigLoaded = true;

  LOG( "ExoticLLP", pDEBUG ) << "FluxCreator configured.";
}
//____________________________________________________________________________
void FluxCreator::SetInputFluxPath(std::string finpath) const
{
  LOG( "ExoticLLP", pNOTICE ) << "Reading fluxes from " << finpath;
  fCurrPath = finpath;
}
//____________________________________________________________________________
int FluxCreator::GetNFluxEntries() const
{
  if( fNEntries <= 0 ) this->OpenFluxInput( fCurrPath );
  return fNEntries;
}
//____________________________________________________________________________
void FluxCreator::SetCurrentEntry( int iCurr ) const
{
  iCurrEntry = iCurr;
}
//____________________________________________________________________________
void FluxCreator::SetFirstFluxEntry( int iFirst ) const
{
  fFirstEntry = iFirst;
}
//____________________________________________________________________________
void FluxCreator::OpenFluxInput( std::string finpath ) const
{
  if( fPathLoaded ) return;

  iCurrEntry = fFirstEntry;
  fCurrPath = finpath;

  LOG( "ExoticLLP", pDEBUG ) << "Getting flux input from finpath = " << finpath;
  if( !ctree ){ ctree = new TChain( ctree_name.c_str() ); cmeta = new TChain( cmeta_name.c_str() ); }

  // recurse over files in this directory and add to chain
  if( this->ends_with( finpath, std::string(".root") ) ){ // passed a single file, add to chain
    LOG( "ExoticLLP", pNOTICE ) << "Adding flux file " << finpath;
    ctree->Add( finpath.c_str() ); cmeta->Add( finpath.c_str() );
  } else {
    if( !this->ends_with( finpath, std::string("/") ) ) finpath.append("/");
    std::list<TString> files = this->RecurseOverDir( finpath );
    std::list<TString>::iterator itFiles = files.begin();
    while( itFiles != files.end() ){
      TString fullpath = *itFiles;
      LOG( "ExoticLLP", pNOTICE ) << "Adding flux file " << fullpath;
      ctree->Add( fullpath ); cmeta->Add( fullpath );
      ++itFiles;
    } // add files to chain
  }

  if( !ctree ) LOG( "ExoticLLP", pFATAL ) << "Could not open flux tree with name " << ctree_name << "!";
  if( !cmeta ) LOG( "ExoticLLP", pFATAL ) << "Could not open flux meta with name " << cmeta_name << "!";
  assert( ctree && cmeta && "Could open flux and meta trees" );

  fNEntries = ctree->GetEntries();

  // Now set branch addresses!
  //ctree->SetBranchAddress( m_ct_parent_pdg.Alias, &(m_ct_parent_pdg.Value) );

  fPathLoaded = true;
}
//----------------------------------------------------------------------------
std::list<TString> FluxCreator::RecurseOverDir( std::string finpath ) const
{
  // Cascades down in the directory and finds all files ending in .root

  TSystemDirectory topDir( finpath.c_str(), finpath.c_str() );
  std::list<TString> files; int nFiles = 0;
  std::list<TString> dirNames;
  std::list<TObject *> dirs; // this will take all directories that have not been opened yet.
  dirs.emplace_front( &topDir );
  dirNames.emplace_front( topDir.GetName() );

  LOG( "ExoticLLP", pDEBUG )
    << "Starting to add files to input. Current size is " << dirs.size();

  while( dirs.size() > 0 ){ // there is still stuff we haven't looked at.
    int nNow = dirs.size();
    LOG( "ExoticLLP", pDEBUG ) 
      << "Scanning directory " << dirNames.front() << " with " << nNow << " elements...";

    // go to first object and get the structure next level down
    TSystemDirectory * currDir = dynamic_cast<TSystemDirectory *>( dirs.front() );
    TString dirPath = dirNames.front();

    // first, strip the first two elements . and ..
    TList * rootElements = currDir->GetListOfFiles(); rootElements->Sort();
    rootElements->ls();
    rootElements->Remove( rootElements->First() ); // .
    rootElements->Remove( rootElements->First() ); // ..

    if( rootElements->GetEntries() == 0 ) continue;
    else rootElements->ls();

    //TSystemFile * elem;
    TIter next(rootElements);
    TObject * sFile;
    TIter sNext( rootElements );

    // put all the files in the list, and all the directories in the dirs list.
    // for names, add the full path (== dirPath + "/" + name of next dir)
    // TSystemDirectory inherits from TSystemFile
    while( (sFile = sNext()) != NULL ){ // you read that right, this is one =.
      TString fullPath = dirPath + "/" + sFile->GetName() ;
      if( dynamic_cast< TSystemDirectory * >( sFile ) ) {
	dirs.emplace_back( sFile );
	dirNames.emplace_back( fullPath );
	LOG( "ExoticLLP", pDEBUG ) 
	  << "Adding directory " << fullPath.Data() << " to linked list...";
      } else
	files.emplace_back( fullPath );
    }

    dirs.pop_front();
    dirNames.pop_front();
  } // while dirs.size() > 0

  nFiles = files.size();
  LOG( "ExoticLLP", pDEBUG )
    << "Found " << nFiles << " files in total.";
  return files;
}
//----------------------------------------------------------------------------
void FluxCreator::InitialiseTree() const
{
  m_ct_parent_pdg     = AliasedBranch<int>( 0, m_ct_parent_pdg_alias );
  m_ct_parent_p4      = AliasedBranch<TLorentzVector>( TLorentzVector( 0.0, 0.0, 0.0, 0.0 ) );
  m_ct_production_v4  = AliasedBranch<TLorentzVector>( TLorentzVector( 0.0, 0.0, 0.0, 0.0 ) );
}
//----------------------------------------------------------------------------
void FluxCreator::InitialiseMeta() const
{
  m_cm_run_no = AliasedBranch<int>( 0 );
}
// A whole bunch of auxiliary functions here made to help with the acceptance correction
//----------------------------------------------------------------------------
double FluxCreator::AccCorr_Sqrt( double thetalab, double mass, 
				  double EPar, double MPar, double ENu ) const
{
  double theta = thetalab * TMath::DegToRad();
  double tanTheta = TMath::Tan( theta );
  
  double pPar = std::sqrt( EPar * EPar - MPar * MPar );
  double bPar = pPar / EPar;
  double gPar = EPar / MPar;

  double arg1 = ( ENu * ENu - mass * mass ) * 
    ( 1.0 + tanTheta * tanTheta * gPar * gPar );
  double arg2 = ( tanTheta * bPar * gPar * ENu ) * ( tanTheta * bPar * gPar * ENu );

  return std::sqrt( std::max( 0.0, arg1 - arg2 ) );
}
//----------------------------------------------------------------------------
double FluxCreator::AccCorr_Denom( double thetalab, double mass, 
				   double EPar, double MPar, double ENu ) const
{
  double theta = thetalab * TMath::DegToRad();
  double tanTheta = TMath::Tan( theta );
  
  double pPar = std::sqrt( EPar * EPar - MPar * MPar );
  double bPar = pPar / EPar;
  double gPar = EPar / MPar;

  double qNu = std::sqrt( ENu * ENu - mass * mass );

  double oth = tanTheta * bPar * gPar * ENu;
  
  return ( ( qNu - oth ) * ( qNu + oth ) );
}
//----------------------------------------------------------------------------
double FluxCreator::AccCorr_SolnArgs( double thetalab, double mass, 
				      double EPar, double MPar, double ENu, bool isPos ) const
{
  double theta = thetalab * TMath::DegToRad();
  double tanTheta = TMath::Tan( theta );

  double pPar = std::sqrt( EPar * EPar - MPar * MPar );
  double bPar = pPar / EPar;
  double gPar = EPar / MPar;

  double denom = this->AccCorr_Denom( thetalab, mass, EPar, MPar, ENu );

  double qNu = std::sqrt( ENu * ENu - mass * mass );
  double sqt = this->AccCorr_Sqrt( thetalab, mass, EPar, MPar, ENu );

  int coeff = ( isPos ) ? 1 : -1;
  double numer = qNu * qNu + coeff * bPar * ENu * sqt;

  return tanTheta * gPar * numer / denom;
}
//----------------------------------------------------------------------------
double FluxCreator::AccCorr_Solution( double thetalab, double mass, 
				      double EPar, double MPar, double ENu, bool isPos ) const
{
  double pPar = std::sqrt( EPar * EPar - MPar * MPar );
  double bPar = pPar / EPar;
  double gPar = EPar / MPar;

  double qNu = std::sqrt( ENu * ENu - mass * mass );
  double bNu = qNu / ENu;
  //double gNu = (bNu < 1.0) ? ENu / mass : -1.0;

  double tanMaxTheta = ( bNu >= bPar ) ? 180.0 : 1.0 / ( gPar * std::sqrt( ( bPar / bNu ) * ( bPar / bNu ) - 1.0 ) );
  double maxTheta = TMath::ATan( tanMaxTheta ) * TMath::RadToDeg();
  if( maxTheta < 0.0 ) maxTheta += 180.0;

  double arg = 0.0;
  if( isPos ){ // positive solution
    if( bNu >= bPar ){ // labangle is surjective
      arg = this->AccCorr_SolnArgs( thetalab, mass, EPar, MPar, ENu, isPos );
      double tsol = TMath::ATan( arg ) * TMath::RadToDeg();
      if( tsol < 0.0 ) tsol += 180.0;
      return tsol;
    } else { // labangle is not surjective, check if preimage can be found
      if( thetalab > maxTheta ) return 0.0;

      arg = this->AccCorr_SolnArgs( thetalab, mass, EPar, MPar, ENu, isPos );
      double tsol = TMath::ATan( arg ) * TMath::RadToDeg();
      if( tsol < 0.0 ) tsol += 180.0;
      return tsol;
    }
  } else { // negative solution
    if( bNu >= bPar ){ // LLP is too fast, so labangle is monotonic and no negative solution
      return 0.0;
    } else { // check if preimage can be found
      if( thetalab > maxTheta ) return 0.0;

      arg = this->AccCorr_SolnArgs( thetalab, mass, EPar, MPar, ENu, isPos );
      double tsol = TMath::ATan( arg ) * TMath::RadToDeg();
      if( tsol < 0.0 ) tsol += 180.0;
      return tsol;
    }
  }

  return -1.0; // you should never see this.
}
//----------------------------------------------------------------------------
double FluxCreator::Forwards_Fcn( double Theta, TLorentzVector p4par, TLorentzVector p4LLP ) const
{
  /*
   * This method calculates the lab-frame emission angle from a rest-frame emission angle.
   *
   * In the following y = tan\theta
   */

  double ThetaRad = Theta * TMath::DegToRad();
  double ct = TMath::Cos( ThetaRad ); double st = TMath::Sin( ThetaRad );

  double q = p4LLP.P(); // rest frame momentum
  double E = p4LLP.E(); // rest frame energy
  
  double gamma = p4par.E() / p4par.M();
  double beta  = p4par.P() / p4par.E();

  double thetaRad = TMath::ATan( q * st / ( gamma * ( beta * E + q * ct ) ) );
  double theta = thetaRad * TMath::RadToDeg();
  if( theta < 0.0 ) theta = 180.0 - theta;
  
  return theta;
}
//----------------------------------------------------------------------------
double FluxCreator::Inverted_Fcn( double theta, TLorentzVector p4par, TLorentzVector p4LLP, bool backwards ) const
{
  /*
   * This method calculates the rest-frame emission angle that yields a lab-frame emission angle.
   * If the rest --> lab formula is not monotonic, the restriction on either the [0, max] (increasing)
   * or [max, 180] (decreasing) is understood. By default the increasing is chosen, unless
   * backwards == true ==> decreasing
   *
   * In the following y = tan\theta, x = cos\Theta ( as sin\Theta >= 0 always for \Theta <= 180 deg )
   */

  double inv = 0.0;

  double q = p4LLP.P(); // rest frame momentum
  double E = p4LLP.E(); // rest frame energy
  
  double gamma = p4par.E() / p4par.M();
  double beta  = p4par.P() / p4par.E();

  if( theta == 90.0 ) return TMath::ACos( -beta * E / q ) * TMath::RadToDeg();

  double y = TMath::Tan( theta * TMath::DegToRad() );
  
  // Answer as calculated by Mathematica. This is Cos[Theta]
  int constMod = backwards ? 1 : -1;
  double constTerm = constMod * beta * E * std::pow( gamma * y, 2.0 );
  
  double sqrt_1  = std::pow( q, 2.0 );
  double sqrt_2a = std::pow( y * gamma, 2.0 );
  double sqrt_2b = q * q - std::pow( beta * E, 2.0 );
  double sqrt   = sqrt_1 + sqrt_2a * sqrt_2b;
  if( sqrt < 0.0 ) sqrt = 0.0;
  sqrt = std::sqrt( sqrt );

  double denom = q * (1.0 + y * y * gamma * gamma);

  double c1 = (constTerm + sqrt) / denom;
  double c2 = (constTerm - sqrt) / denom; if( c2 == c1 ) c2 = -1.0;

  inv = ( theta <= 90.0 ) ? TMath::ACos( c1 ) : TMath::ACos( c2 );
  inv = backwards ? 180.0 - inv * TMath::RadToDeg() : inv * TMath::RadToDeg() ;
  
  return inv;
}
//----------------------------------------------------------------------------
double FluxCreator::Derivative( double theta, TLorentzVector p4par, TLorentzVector p4LLP ) const
{
  // Gives the derivative of the collimation function at theta
  double der = 0.0;

  double beta  = p4par.P() / p4par.E();
  double gamma = p4par.E() / p4par.M();
  
  double q = p4LLP.P();
  double E = p4LLP.E();

  double cx = TMath::Cos( TMath::DegToRad() * theta );
  double sx = TMath::Sin( TMath::DegToRad() * theta );

  double term1 = q * cx / ( gamma * ( beta * E + q * cx ) );
  double term2 = q * q * sx * sx / ( gamma * std::pow( beta * E + q * cx , 2.0 ) );
  double numer = term1 + term2;
  double denom = 1.0 + std::pow( q * sx / ( gamma * ( beta * E + q * cx ) ) , 2.0 );

  der = numer / denom;

  return der;
}

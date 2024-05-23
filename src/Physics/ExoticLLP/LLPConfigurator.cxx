//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kjplows \at liverpool.ac.uk>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/LLPConfigurator.h"

using namespace genie;
using namespace genie::llp;

//____________________________________________________________________________
LLPConfigurator::LLPConfigurator() : 
  ChannelCalculatorI("genie::llp::LLPConfigurator")
{
  fIsConfigLoaded = false;
}
//____________________________________________________________________________
LLPConfigurator::LLPConfigurator(string name) : ChannelCalculatorI(name)
{
  fIsConfigLoaded = false;
}
//____________________________________________________________________________
LLPConfigurator::LLPConfigurator(string name, string config) : 
  ChannelCalculatorI(name, config)
{
  fIsConfigLoaded = false;
}
//____________________________________________________________________________
LLPConfigurator::~LLPConfigurator()
{

}
//____________________________________________________________________________
void LLPConfigurator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LLPConfigurator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LLPConfigurator::LoadConfig(void)
{
  if( fIsConfigLoaded ) return;

  this->GetParam( "Mass", fMass );

  // RETHERE Majorana vs Dirac?
  // RETHERE spin?

  // Path to the XML files
  std::string relative_path;
  this->GetParam( "InputChannels", relative_path );
  
  // Get the tune directory from RunOpt
  TuneId * tune = RunOpt::Instance()->Tune();
  // ask it where the TuneDirectory lives
  std::string tune_dir = tune->TuneDirectory();

  fInputPath = tune_dir;
  fInputPath.append("/"); fInputPath.append( relative_path.c_str() );

  LOG( "ExoticLLP", pFATAL )
    << "\nTune directory is " << tune_dir
    << "\nRelative path is " << relative_path
    << "\n==> full path is " << fInputPath;

  // Now that we have the path, read in the XML file
  this->ParseInputFile();
  
  // construc the LLP!
  ModeVector all_modes;
  for( ModeVector::iterator mvit = fProdChannels.begin() ; 
       mvit != fProdChannels.end() ; ++mvit ) all_modes.emplace_back( *mvit );
  for( ModeVector::iterator mvit = fDecayChannels.begin() ; 
       mvit != fDecayChannels.end() ; ++mvit ) all_modes.emplace_back( *mvit );

  ExoticLLP * ptLLP = new ExoticLLP( fMass, all_modes );

  // For now, just make a single production mode K+ --> X + e+ and decay mode X --> e- e+
  /*
  std::vector<int> prodDummy = { kPdgKP, kPdgLLP, kPdgPositron };
  std::vector<int> prod2Dummy = { kPdgKP, kPdgLLP, kPdgAntiMuon };
  std::vector<int> prod3Dummy = { kPdgDP, kPdgLLP, kPdgAntiTau };
  std::vector<int> decDummy  = { kPdgLLP, kPdgElectron, kPdgPositron };
  std::vector<int> dec2Dummy  = { kPdgLLP, kPdgMuon, kPdgAntiMuon };

  std::pair< double, std::vector<int> > dprod_pair = std::make_pair( 0.5, prodDummy );
  std::pair< double, std::vector<int> > dprod2_pair = std::make_pair( 0.5, prod2Dummy );
  std::pair< double, std::vector<int> > dprod3_pair = std::make_pair( 0.5, prod3Dummy );
  std::pair< double, std::vector<int> > ddec_pair  = std::make_pair( 0.5, decDummy  );
  std::pair< double, std::vector<int> > ddec2_pair  = std::make_pair( 0.5, dec2Dummy  );

  ModeVector dmodes = { dprod_pair, dprod2_pair, dprod3_pair, ddec_pair, ddec2_pair };

  ExoticLLP * ptLLP = new ExoticLLP( fMass, dmodes );
  */
  
  fLLP = *ptLLP;
}
//____________________________________________________________________________
void LLPConfigurator::ParseInputFile() const 
{
  // from the input path, construct an xmlDocPtr and have the thing parse files
  xmlDocPtr input_xml = xmlParseFile( fInputPath.c_str() );

  // xmlDocGetRootElement points to <alg_conf>
  xmlNodePtr root_element = xmlDocGetRootElement( input_xml );
  xmlNodePtr loaded_set = root_element->xmlChildrenNode;

  // make sure we grab the param_set
  for( ; loaded_set != NULL ; loaded_set = loaded_set->next ) {
    std::string pname = utils::xml::TrimSpaces( const_cast<xmlChar*>( loaded_set->name ) );
    if( pname == "param_set" ) break;
  }

  LOG( "ExoticLLP", pFATAL ) << utils::xml::TrimSpaces( const_cast<xmlChar*>( loaded_set->name ) );

  //fXmlChannels = this->ParseParamSet( input_xml, loaded_set );
  tie( fProdChannels, fDecayChannels ) = this->ParseParamSet( input_xml, loaded_set );
}
//____________________________________________________________________________
std::tuple< 
  ModeVector, ModeVector 
  > LLPConfigurator::ParseParamSet( xmlDocPtr & doc, xmlNodePtr pset ) const
{

  ModeVector prod_channels;
  ModeVector decay_channels;

  // First, print some details about this

  xmlNodePtr child = pset->xmlChildrenNode; 
  
  for( ; child != NULL ; child = child->next ) {

    std::string pname = utils::xml::TrimSpaces( const_cast<xmlChar*>( child->name ) );

    if( pname == "param" ) {
      std::string ppname = utils::str::TrimSpaces( utils::xml::GetAttribute( child, "name" ) );
      if( ppname == "Details" ) {
	std::string details = std::string( (const char*) (xmlNodeListGetString( doc, child->xmlChildrenNode, 1 )) );
	LOG( "ExoticLLP", pFATAL )
	  << "\nFound a param_set from input XML at path " << fInputPath
	  << "\nHere are some details: " << details;
      } // found details
    } // <param> under <param_set>

    // now look at the <*Modes>
    
    if( pname == "productionModes" )
      prod_channels = this->ParseModes( doc, child );
    else if( pname == "decayModes" )
      decay_channels = this->ParseModes( doc, child );

  } // children of <param_set>

  return std::make_tuple( prod_channels, decay_channels );
}
//____________________________________________________________________________
ModeVector LLPConfigurator::ParseModes( xmlDocPtr & doc, xmlNodePtr node ) const 
{
  ModeVector return_vector;

  xmlNodePtr child = node->xmlChildrenNode;
  std::vector<double> knot_masses;
  ModeKnotVector knot_modes;
  double interp = 0.0;
  
  for( ; child != NULL ; child = child->next ) {
    
    std::string gpname = utils::xml::TrimSpaces( const_cast<xmlChar*>( child->name ) );

    // first get the masses at which this is evaluated
    if( gpname == "param" ) {
	     
      if( utils::str::TrimSpaces( utils::xml::GetAttribute( child, "type" ) ) == "vec-double" &&
	  utils::str::TrimSpaces( utils::xml::GetAttribute( child, "name" ) ) == "Masses" ) {
	std::string val = utils::xml::TrimSpaces( xmlNodeListGetString( doc, child->xmlChildrenNode, 1 ) );
	knot_masses = this->GetDoubleVector(val);
      } // populate knot_masses

    } // <param> under <productionModes>

    // Now get the full objects of modes
    else if( gpname == "modeObject" ) { 

      std::vector< double > knot_scores;
      std::vector< int > mode_pdgList;

      xmlNodePtr gdChild = child->xmlChildrenNode;
	  
      for( ; gdChild != NULL ; gdChild = gdChild->next ) {

	std::string ggname = utils::xml::TrimSpaces( const_cast<xmlChar*>( gdChild->name ) );
	if( ggname == "param" ) {
	    
	  if( utils::str::TrimSpaces( utils::xml::GetAttribute( gdChild, "type" ) ) == "vec-double" &&
	      utils::str::TrimSpaces( utils::xml::GetAttribute( gdChild, "name" ) ) == "Scores" ) {
	    std::string val = utils::xml::TrimSpaces( xmlNodeListGetString( doc, gdChild->xmlChildrenNode, 1 ) );
	    knot_scores = this->GetDoubleVector(val);
	  }
	  else if( utils::str::TrimSpaces( utils::xml::GetAttribute( gdChild, "type" ) ) == "vec-int" &&
		   utils::str::TrimSpaces( utils::xml::GetAttribute( gdChild, "name" ) ) == "PDGList" ) {
	    std::string val = utils::xml::TrimSpaces( xmlNodeListGetString( doc, gdChild->xmlChildrenNode, 1 ) );
	    mode_pdgList = this->GetIntVector(val);
	      
	  } // either scores or pdg list of mode

	} // <param> under <modeObject>

      } // children of <modeObject>
	  
      // Got both the scores and the pdgList from the modeObject
      knot_modes.emplace_back( std::make_pair< 
			       std::vector<double>, std::vector<int> 
			       >( std::move(knot_scores), std::move(mode_pdgList) ) );
	  
    } // <modeObject> under <productionModes>

  } // children of <productionModes>


  // Now that we parsed the masses, let's figure out which bit to use
  // first, assert that the mass is sensible, otherwise ask the user to fix
  double mass_small = *(knot_masses.begin());
  double mass_large = *(knot_masses.end()-1);
  if( fMass < mass_small || fMass > mass_large ){
    LOG( "ExoticLLP", pFATAL )
      << "Fatal error: the LLP mass is incompatible with the given XML inputs!"
      << "\nInput mass: " << fMass << " GeV, mass range [ " << mass_small
      << ", " << mass_large << " ] GeV";
    assert( fMass >= mass_small && fMass <= mass_large && "LLP mass in range specified from XML inputs at production" );
  }
      
  std::vector<double>::iterator kit = knot_masses.begin();
  while( kit != knot_masses.end() && *kit <= fMass )
    ++kit;

      
  double nearest_mass_small = *(knot_masses.begin());
  double nearest_mass_large = *(knot_masses.end()-1);
  if( kit != knot_masses.begin() ) {
    nearest_mass_small = *(kit-1);
    nearest_mass_large = *kit;
	
    interp = ( fMass - nearest_mass_small ) / ( nearest_mass_large - nearest_mass_small );
  }
      
  LOG( "ExoticLLP", pFATAL ) << "Interpolating between [ "
			     << nearest_mass_small << ", " << nearest_mass_large
			     << " ] ==> interp = " << interp;

  // Okay, so now let's construct for each mode the scores
  ModeKnotVector::iterator kmit = knot_modes.begin();
      
  for( ; kmit != knot_modes.end(); ++kmit ) {

    std::vector< int > tmpPDGList = (*kmit).second;
    std::vector< double > tmpScores = (*kmit).first;
	
    int i1 = std::max( static_cast<int>( (kit-1) - knot_masses.begin() ), 0 );
    int i2 = std::min( 1, static_cast<int>( kit - knot_masses.begin() ) );
    double s1 = tmpScores.at(i1), s2 = tmpScores.at(i2);

    double win_score = s1 + ( s2 - s1 ) * interp;

    return_vector.emplace_back( std::make_pair< 
				double, std::vector<int> 
				> ( std::move( win_score ), std::move( tmpPDGList ) ) );
	
  } // loop over all input modes and determine scores

  return return_vector;
}
//____________________________________________________________________________
std::vector<double> LLPConfigurator::GetDoubleVector( std::string str ) const
{
  // copied from GNuMIFlux.cxx
  // turn string into vector<double>
  // be liberal about separators, users might punctuate for clarity
  std::vector<std::string> strtokens = genie::utils::str::Split(str, " ,;:()[]=\t\n");
  std::vector<double> vect;
  size_t ntok = strtokens.size();

  for (size_t i=0; i<ntok; ++i) {
    std::string trimmed = utils::str::TrimSpaces(strtokens[i]);
    if( " " == trimmed || "" == trimmed ) continue; // skip empty strings
    double val = strtod(trimmed.c_str(), (char**)NULL);
    vect.push_back(val);
  }

  return vect;
}
//____________________________________________________________________________
std::vector<int> LLPConfigurator::GetIntVector( std::string str ) const
{
  // copied from GNuMIFlux.cxx
  // turn string into vector<int>
  // be liberal about separators, users might punctuate for clarity
  std::vector<std::string> strtokens = genie::utils::str::Split(str, " ,;:()[]=\t\n");
  std::vector<int> vect;
  size_t ntok = strtokens.size();

  for (size_t i=0; i<ntok; ++i) {
    std::string trimmed = utils::str::TrimSpaces(strtokens[i]);
    if( " " == trimmed || "" == trimmed ) continue; // skip empty strings
    int val = strtod(trimmed.c_str(), (char**)NULL);
    vect.push_back(val);
  }

  return vect;
}
//____________________________________________________________________________

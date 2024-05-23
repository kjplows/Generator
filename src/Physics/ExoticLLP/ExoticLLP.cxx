//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kjplows \at liverpool.ac.uk>
 University of Liverpool

*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Physics/ExoticLLP/ExoticLLP.h"

using namespace genie;
using namespace genie::llp;

//____________________________________________________________________________
ExoticLLP::ExoticLLP() 
{

}
//____________________________________________________________________________
ExoticLLP::ExoticLLP( double mass, 
		      std::vector< std::pair< double, std::vector< int > > > modes ) :
  fMass( mass )
{
  ConstructModes( modes );
}
//____________________________________________________________________________
ExoticLLP::~ExoticLLP()
{

}
//____________________________________________________________________________
void ExoticLLP::ConstructModes( std::vector< std::pair< double, 
				std::vector< int > > > modes ) const 
{
  // The PDG database GENIE uses
  TDatabasePDG * dbase = PDGLibrary::Instance()->DBase();
  
  std::vector< std::pair< double, std::vector< int > > >::iterator mit = modes.begin();

  for( ; mit != modes.end(); ++mit ) {

    std::vector< int > pdgList = (*mit).second;
    std::vector< int > nuPdgList = pdgList; // to reorder products later

    LOG( "ExoticLLP", pFATAL ) << "Testing mode " << (mit - modes.begin());
    std::ostringstream asts; asts << "Mode " << (mit - modes.begin()) << ": ";
    for( std::vector<int>::iterator pdgit = nuPdgList.begin() ;
	 pdgit != nuPdgList.end() ; ++pdgit ) asts << *pdgit << " ";
    LOG( "ExoticLLP", pFATAL ) << asts.str();
    
    // Production or decay?
    bool isProduction = false; bool isDecay = false;
    if( std::abs(pdgList.at(0)) == kPdgLLP ) { 
      isDecay = true;
    } else { // Search for the LLP in the daughter list

      bool llp_found = false;
      std::vector< int >::iterator pit = pdgList.begin(); pit++;
      for( ; pit != pdgList.end(); ++pit ) {
	if( std::abs(*pit) == kPdgLLP ) { llp_found = true; break; }
      }

      if( !llp_found ){ // Huh. We didn't find the LLP. Warn and exit out
	LOG( "ExoticLLP", pFATAL )
	  << "\nThere is a mode in the input XML file that does not contain the LLP PDG code: "
	  << kPdgLLP << " or its CP conjugate."
	  << "\nThis indicates a problem with the input XML file. Exiting now to prevent problems downstream.";
      }
      assert( llp_found && "LLP PDG code in each given mode in the XML" );

      isProduction = true;
      /*
      // Reorder the list so that LLP is the first daughter
      std::vector< int > tmpList; tmpList.emplace_back( pdgList.at(0) );
      tmpList.emplace_back( kPdgLLP );
      for( std::vector<int>::iterator ppit = pdgList.begin(); ppit != pdgList.end(); ++ppit ) {
	LOG( "ExoticLLP", pFATAL ) << "Checking  " << *ppit << " at position " << tmpList.size();
	if( *ppit != pdgList.at(0) && *pit != kPdgLLP ) { 
	  tmpList.emplace_back( *ppit );
	  LOG( "ExoticLLP", pFATAL ) << "Emplacing back " << *ppit << " at position " << tmpList.size();
	}
      }
      nuPdgList = tmpList; // now contains LLP in position 1
      */
    } // production modes

    LOG( "ExoticLLP", pFATAL ) << "Testing mode " << (mit - modes.begin()) << " again!";
    std::ostringstream bsts; bsts << "Mode " << (mit - modes.begin()) << ": ";
    for( std::vector<int>::iterator pdgit = nuPdgList.begin() ;
	 pdgit != nuPdgList.end() ; ++pdgit ) bsts << *pdgit << " ";
    LOG( "ExoticLLP", pFATAL ) << bsts.str();

    assert( ( isProduction || isDecay ) && "Each mode is either an LLP production or LLP decay mode." );

    // Build a ModeObject
    double score = (*mit).first;
    if( isDecay ) {
      int idx = static_cast<int>(fDecayModes.size());
      // build the name
      std::vector<int>::iterator pit = nuPdgList.begin();
      std::string name = "";
      if( *pit < 0 ) name.append("LLPBarTo:");
      else name.append("LLPTo:");
      pit++;

      for( ; pit != nuPdgList.end(); ++pit ) {
	TParticlePDG * tmp_particle = dbase->GetParticle(*pit);
	name.append( tmp_particle->GetName() );
	if( pit < nuPdgList.end() - 1 ) name.append( ":" );
      } // name construction

      ModeObject mobj;
      mobj.fIndex = idx;
      mobj.fMode = 1;
      mobj.fName = name;
      mobj.fScore = score;
      mobj.fPDGList = nuPdgList;

      // and append to the modes
      fDecayModes.emplace_back( mobj );
    } // decay modes
    else {
      int idx = static_cast<int>(fProductionModes.size());
      // build the name
      std::vector<int>::iterator pit = nuPdgList.begin();
      TParticlePDG * parent_particle = dbase->GetParticle(*pit);
      std::string name = "";
      name.append( parent_particle->GetName() );

      pit++;
      if( *pit < 0 ) name.append("To:LLPBar:");
      else name.append("To:LLP:");

      pit++;
      for( ; pit != nuPdgList.end(); ++pit ) {
	TParticlePDG * tmp_particle = dbase->GetParticle(*pit);
	name.append( tmp_particle->GetName() );
	if( pit < nuPdgList.end() - 1 ) name.append( ":" );
      } // name construction

      ModeObject mobj;
      mobj.fIndex = idx;
      mobj.fMode = 0;
      mobj.fName = name;
      mobj.fScore = score;
      mobj.fPDGList = nuPdgList;

      // and append to the modes
      fProductionModes.emplace_back( mobj );
    } // production modes
    
  } // mit = modes.begin() --> modes.end()
}
//____________________________________________________________________________
void ExoticLLP::Print(const Option_t * /* opt */) const
{
  std::cout << *this << std::endl;
}
//____________________________________________________________________________

namespace genie {
  namespace llp {
    std::ostream & operator << ( std::ostream & stream, const ExoticLLP & LLP ) 
    {
      stream << "Test print"
	     << "\nLLP mass: " << LLP.GetMass() << " MeV / c^2";

      std::vector< ModeObject > prodModes = LLP.GetProductionModes();
      std::vector< ModeObject > decayModes = LLP.GetDecayModes();

      stream << "\nProduction modes:";
      for( std::vector<ModeObject>::iterator vit = prodModes.begin(); vit != prodModes.end() ; ++vit )
	stream << "\n" << (*vit).GetName();
      stream << "\nDecay modes:";
      for( std::vector<ModeObject>::iterator vit = decayModes.begin(); vit != decayModes.end() ; ++vit )
	stream << "\n" << (*vit).GetName();
      return stream;
    }
  } 
}
//____________________________________________________________________________

//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kplows \at liverpool.ac.uk>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/LLPFluxContainer.h"

using namespace genie;
using namespace genie::llp;

//___________________________________________________________________________
FluxContainer::FluxContainer() : TObject()
{
  this->ResetCopy();
}
//___________________________________________________________________________
void FluxContainer::ResetCopy() const
{
  evtno = 0;

  v4.SetXYZT( 0.0, 0.0, 0.0, 0.0 );
  v4_user.SetXYZT( 0.0, 0.0, 0.0, 0.0 );

  p4_parent.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );
  p4_parent_user.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );

  entry.SetXYZT( 0.0, 0.0, 0.0, 0.0 );
  entry_user.SetXYZT( 0.0, 0.0, 0.0, 0.0 );

  exit.SetXYZT( 0.0, 0.0, 0.0, 0.0 );
  exit_user.SetXYZT( 0.0, 0.0, 0.0, 0.0 );

  p4.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );
  p4_user.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );
  
  wgt_xy = 0.0;
}
//___________________________________________________________________________
void FluxContainer::Print(const Option_t* /* opt */ ) const
{
  LOG("ExoticLLP", pNOTICE) << *this;
}
//___________________________________________________________________________

namespace genie {
namespace llp  {
  ostream & operator << (
    ostream & stream, const FluxContainer & info)
    {
      stream << "Printing information on this ExoticLLP...";
      stream << "\nevtno = " << info.evtno;
      stream << "\nProduction vertex [NEAR, m, ns] = " << utils::print::X4AsString( &(info.v4) )
	     << "\nProduction vertex [USER, m, ns] = " << utils::print::X4AsString( &(info.v4_user) )
	     << "\nLLP momentum    [NEAR, GeV] = " << utils::print::P4AsString( &(info.p4) )
	     << "\nLLP momentum    [USER, GeV] = " << utils::print::P4AsString( &(info.p4_user) )
	     << "\nParent momentum [NEAR, GeV] = " << utils::print::P4AsString( &(info.p4_parent) )
	     << "\nParent momentum [USER, GeV] = " << utils::print::P4AsString( &(info.p4_parent_user) )
	     << "\nEntry point     [NEAR, m, ns]   = " << utils::print::X4AsString( &(info.entry) )
	     << "\nEntry point     [USER, m, ns]   = " << utils::print::X4AsString( &(info.entry_user) )
	     << "\nExit point      [NEAR, m, ns]   = " << utils::print::X4AsString( &(info.exit) )
	     << "\nExit point      [USER, m, ns]   = " << utils::print::X4AsString( &(info.exit_user) )
	     << "\nwgt_xy = " << info.wgt_xy;
      return stream;
    }
}
}

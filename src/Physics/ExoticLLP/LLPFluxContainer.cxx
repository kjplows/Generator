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

ClassImp( FluxContainer )

//___________________________________________________________________________
FluxContainer::FluxContainer() : TObject()
{
  this->ResetCopy();
}
//___________________________________________________________________________
FluxContainer::FluxContainer(const FluxContainer & flc) : TObject(flc)
{
  mass = flc.mass;
  lifetime = flc.lifetime;

  evtno = flc.evtno;
  pdg   = flc.pdg;
  
  v4 = flc.v4;
  v4_user = flc.v4_user;

  p4_parent = flc.p4_parent;
  p4_parent_user = flc.p4_parent_user;

  entry = flc.entry;
  entry_user = flc.entry_user;
  
  exit = flc.exit;
  exit_user = flc.exit_user;

  p4 = flc.p4;
  p4_user = flc.p4_user;

  decay = flc.decay;
  decay_user = flc.decay_user;

  wgt_xy = flc.wgt_xy;
  boost_factor = flc.boost_factor;
  wgt_collimation = flc.wgt_collimation;

  wgt_survival = flc.wgt_survival;
  wgt_detdecay = flc.wgt_detdecay;

  vtx_rng = flc.vtx_rng;
}
//___________________________________________________________________________
FluxContainer & FluxContainer::operator = (const FluxContainer & flc)
{
  this->mass = flc.mass;
  this->lifetime = flc.lifetime;

  this->evtno = flc.evtno;
  this->pdg   = flc.pdg;

  this->v4 = flc.v4;
  this->v4_user = flc.v4_user;

  this->p4_parent = flc.p4_parent;
  this->p4_parent_user = flc.p4_parent_user;

  this->entry = flc.entry;
  this->entry_user = flc.entry_user;
  
  this->exit = flc.exit;
  this->exit_user = flc.exit_user;

  this->p4 = flc.p4;
  this->p4_user = flc.p4_user;

  this->decay = flc.decay;
  this->decay_user = flc.decay_user;

  this->wgt_xy = flc.wgt_xy;
  this->boost_factor = flc.boost_factor;
  this->wgt_collimation = flc.wgt_collimation;

  this->wgt_survival = flc.wgt_survival;
  this->wgt_detdecay = flc.wgt_detdecay;

  this->vtx_rng = flc.vtx_rng;

  return *this;
}
//___________________________________________________________________________
void FluxContainer::ResetCopy() const
{
  mass = 0.0;
  lifetime = 0.0;

  evtno = 0;
  pdg   = 0;

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

  decay.SetXYZT( 0.0, 0.0, 0.0, 0.0 );
  decay_user.SetXYZT( 0.0, 0.0, 0.0, 0.0 );
  
  wgt_xy = 0.0;
  boost_factor = 0.0;
  wgt_collimation = 0.0;

  wgt_survival = 0.0;
  wgt_detdecay = 0.0;

  vtx_rng = 0.0;
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
      stream << "\nmass = " << info.mass << " [GeV]"
	     << "\nlifetime c*tau = " << info.lifetime << " [m]";
      stream << "\nevtno = " << info.evtno
	     << "\nParent PDG = " << info.pdg;
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
	     << "\nDecay point     [NEAR, m, ns]   = " << utils::print::X4AsString( &(info.decay) )
	     << "\nDecay point     [USER, m, ns]   = " << utils::print::X4AsString( &(info.decay_user) )
	     << "\nwgt_xy          = " << info.wgt_xy
	     << "\nboost_factor    = " << info.boost_factor
	     << "\nwgt_collimation = " << info.wgt_collimation
	     << "\nwgt_survival    = " << info.wgt_survival
	     << "\nwgt_detdecay    = " << info.wgt_detdecay
	     << "\nRandom uniform number used for vertex generation = " << info.vtx_rng;
      
      return stream;
    }
}
}

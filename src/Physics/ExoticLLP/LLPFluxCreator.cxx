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

  // RETHERE: Need to reduce the scores to only those that are kinematically allowed.
  std::vector< ModeObject > llp_vetted_modes;
  std::vector< ModeObject >::iterator it_vetmd = llp_production_modes.begin();
  double total_vetted_score = 0.0;

  while( it_vetmd != llp_production_modes.end() ) {
    // Check if the mode can actually make LLP
    std::vector<int> mode_pdg_list = (*it_vetmd).GetPDGList();
    double mass_sum = 0.0;
    double mass_parent = PDGLibrary::Instance()->Find( mode_pdg_list.at(0) )->Mass();
    for( std::vector<int>::iterator it_pdg = mode_pdg_list.begin()+1; 
	 it_pdg != mode_pdg_list.end(); ++it_pdg ) 
      mass_sum += PDGLibrary::Instance()->Find( *it_pdg )->Mass();

    if( mass_sum <= mass_parent ) {
      llp_vetted_modes.emplace_back( *it_vetmd );
      total_vetted_score += (*it_vetmd).GetScore();
    }
    ++it_vetmd;
  }

  std::vector< ModeObject >::iterator it_modes = llp_vetted_modes.begin();
  // Note that if we have no vetted modes, or if the total score is zero, we must bail.
  if( total_vetted_score == 0.0 || it_modes == llp_vetted_modes.end() ){
    std::vector<int> mode_pdg_list = (*(llp_production_modes).begin()).GetPDGList();
    LOG( "ExoticLLP", pWARN ) << "Cannot make an LLP of mass " 
			      << PDGLibrary::Instance()->Find( kPdgLLP )->Mass() << " GeV from parent "
			      << "with mass = " 
			      << PDGLibrary::Instance()->Find( *(mode_pdg_list.begin()) )->Mass() 
			      << " GeV, so bailing.";
    evrec->SetProbability(0);
    return;
  }

  RandomGen * rnd = RandomGen::Instance();
  double production_score = rnd->RndGen().Rndm();

  double score_seen = (*it_modes).GetScore() / total_vetted_score;
  
  
  while( it_modes != llp_vetted_modes.end() && score_seen < production_score ) {
    ++it_modes; score_seen += (*it_modes).GetScore(); 
  }
  if( it_modes == llp_vetted_modes.end() ) --it_modes;
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
  GHepParticle probe_particle = *(decayed_results.begin());
  TLorentzVector probe_p4_rest = *(probe_particle.GetP4());
  //LOG( "ExoticLLP", pDEBUG ) << "The 4-momentum of the probe in the rest frame is "
  //			     << utils::print::P4AsString( &probe_p4_rest );

  /* 
   * Now we need the acceptance calculations. 
   * First, choose a point inside the top volume
   * Then compare the requested deviation angle with the max deviation
   * If not, reroll the decay and a point and try again
   * If good, get the boost factor and acceptance correction and write it out
   */

  // The information that goes into the EventRecord will always be in USER coordinates 

  // Get a random point in the top volume to evaluate the flux at
  VolumeSeeker * vsek = VolumeSeeker::Instance();
  TVector3 rand_point = vsek->GetRandomPointInTopVolNEAR();
  TVector3 origin     = parent_v4.Vect();
  LOG( "ExoticLLP", pWARN ) << "Evaluating flux at NEAR point " 
			    << utils::print::Vec3AsString( &rand_point );

  // First, apply the constraint of kinematics.
  TVector3 separation = rand_point - origin; // m

  TVector3 parent_p3 = parent_p4.Vect();
  TVector3 parent_p3_unit = parent_p3.Unit();
  
  // Calculate the opening angle between this ray and the parent momentum
  double cos_zeta = separation.Dot( parent_p3_unit ) / ( separation.Mag() * parent_p3_unit.Mag() );
  double zeta = std::acos( cos_zeta ) * 180.0 / constants::kPi; // std::acos has support on [0, \pi]

  /*
  LOG( "ExoticLLP", pDEBUG ) << "\nSeparation between:\n"
			     << utils::print::Vec3AsString( &rand_point )
			     << " and\n"
			     << utils::print::Vec3AsString( &origin )
			     << " is:\n"
			     << utils::print::Vec3AsString( &separation );
  LOG( "ExoticLLP", pDEBUG ) << "\nThe opening angle between that separation and the unit vector,\n"
			     << utils::print::Vec3AsString( &parent_p3_unit )
			     << " is zeta = " << zeta << " [deg]";
  LOG( "ExoticLLP", pDEBUG ) << "\nThe parent 4-momentum is " 
			     << utils::print::P4AsString( &parent_p4 );
  */

  // Just note that if we've started in the decay volume, there is no correction from collimation.
  double boost_factor = 1.0; double acc_corr = 1.0;
  TLorentzVector probe_p4(0.0, 0.0, 1.0, 1.0);
  if( vsek->IsInTop( origin ) ) {

    LOG( "ExoticLLP", pWARN ) << "Decay of parent particle in the top volume!";

    // Just boost probe_p4_rest to the lab frame
    probe_p4 = probe_p4_rest;
    probe_p4.Boost( parent_p4.BoostVector() );
  } else {
    // First, condition for monotonicity loss
    double beta_LLP_rest   = probe_p4_rest.Beta();
    double beta_parent_lab = parent_p4.Beta();
    bool monotonicity_lost = ( beta_LLP_rest < beta_parent_lab );
    double gamma_parent_lab = parent_p4.Gamma();

    [[maybe_unused]] bool is_ok = true;
    if( ! monotonicity_lost ) { // this can always be accepted
      is_ok = true;
    } else { // we must compare to the maximum achievable angle
      
      double velocity_ratio = beta_LLP_rest / beta_parent_lab; // also the cos of the rest-frame emission polar angle at which the maximum deviation from the parent in the lab frame occurs
      double max_tangent = velocity_ratio / ( gamma_parent_lab * std::sqrt( 1.0 - velocity_ratio ) );
      double max_angle = std::atan( max_tangent ) * 180.0 / constants::kPi; // std::atan has support on [-pi/2, pi/2]
      if( max_angle < 0.0 ) max_angle += 180.0;

      is_ok = ( max_angle >= zeta );
      /*
      if( is_ok )
	LOG( "ExoticLLP", pDEBUG ) << "This event, at max angle = " << max_angle << " and zeta = " << zeta << ", will be accepted. Moving on.";
      else
	LOG( "ExoticLLP", pDEBUG ) << "This event, at max angle = " << max_angle << " and zeta = " << zeta << ", will NOT be accepted. Rerolling.";
      */

      // RETHERE reroll a few times until we move on...

      int itry = 0;
      while( ! is_ok && itry < 10 ) {
	itry++;

	rand_point = vsek->GetRandomPointInTopVolNEAR();
	LOG( "ExoticLLP", pWARN ) << "Evaluating flux at NEAR point " 
				  << utils::print::Vec3AsString( &rand_point );
	separation = rand_point - origin; // m
	cos_zeta = separation.Dot( parent_p3_unit ) / ( separation.Mag() * parent_p3_unit.Mag() );
	zeta = std::acos( cos_zeta ) * 180.0 / constants::kPi; // std::acos has support on [0, \pi]

	is_ok = ( max_angle >= zeta );
      }

      if( !is_ok ) { // bail and continue
	evrec->SetProbability(0);
	return;
      }
      
    }
    // RETHERE if ! is_ok bail

    // Since we are accepted, the LLP will always have a region of forwards emission where it will be accepted, and potentially a backwards region too.
    // So calculate these, and if monotonicity was lost, do a throw to decide which region it was

    double pos_soln = this->AccCorr_Solution( zeta, fMass,
					      parent_p4.E(), parent_p4.M(), 
					      probe_p4_rest.E(), true );
    double neg_soln = this->AccCorr_Solution( zeta, fMass,
					      parent_p4.E(), parent_p4.M(), 
					      probe_p4_rest.E(), false );

    double inv_der_pos = std::abs( 1.0 / this->Derivative( pos_soln, parent_p4, probe_p4_rest ) );
    double inv_der_neg = 0.0;
    
    bool forwards = true;

    /*
    std::ostringstream asts;
    asts << "Dumping stats:"
	 << "\nzeta = " << zeta
	 << "\nfMass = " << fMass
	 << "\nparent_p4 = " << utils::print::P4AsString( &parent_p4 )
	 << "\nprobe_p4_rest = " << utils::print::P4AsString( &probe_p4_rest )
	 << "\npos_soln = " << pos_soln
	 << "\nder_pos = " << this->Derivative( pos_soln, parent_p4, probe_p4_rest );
    */

    double image_score = rnd->RndGen().Rndm();    
    if( neg_soln > 0.0 ) {
      inv_der_neg = std::abs( 1.0 / this->Derivative( neg_soln, parent_p4, probe_p4_rest ) );

      if( image_score > inv_der_pos / ( inv_der_pos + inv_der_neg ) ) forwards = false;

      /*
      asts << "\nneg_soln = " << neg_soln
	   << "\nder_neg = " << -1.0 * this->Derivative( neg_soln, parent_p4, probe_p4_rest )
	   << "\nimage_score = " << image_score;
      */
      
    } // if negative solution allowed
    
    /*
     * We need to compare to the acceptance of a massless LLP.
     * If the production mode is 2-body this can be done immediately analytically.
     * If not, then we have to perform a phase-space decay (again RETHERE: including angular correlations?)
     * such that we get a proxy for the rest-frame energy of the massless LLP.
     */
    
    double energy_massless = 0.0;
    double m_parent = parent_p4.M();
    if( decayed_results.size() == 2 )
      energy_massless = probe_p4_rest.E() - fMass * fMass / (2.0 * m_parent);
    else {
      decay_ok = decayer->UnpolarisedDecay(true); // fudge a decay with massless LLP and same other products
      energy_massless = decayer->GetMasslessEnergy();
    }
    //LOG( "ExoticLLP", pDEBUG ) << "A massless LLP would have energy " << energy_massless << " GeV";

    double massless_soln = this->AccCorr_Solution( zeta, 0.0,
						   parent_p4.E(), parent_p4.M(),
						   energy_massless, true ); // always forward emitted
    TLorentzVector massless_p4( energy_massless * probe_p4_rest.Px() / probe_p4_rest.P(),
				energy_massless * probe_p4_rest.Py() / probe_p4_rest.P(),
				energy_massless * probe_p4_rest.Pz() / probe_p4_rest.P(),
				energy_massless );
    double inv_der_massless = std::abs( 1.0 / this->Derivative( massless_soln, parent_p4, massless_p4 ) );
      
    acc_corr = ( forwards ) ? inv_der_pos / inv_der_massless : inv_der_neg / inv_der_massless;
      
    /*
    asts << "\nmassless_soln = " << massless_soln
	 << "\nder_massless = " << this->Derivative( massless_soln, parent_p4, massless_p4 )
	 << "\n==> acc_corr = " << acc_corr;
    */

    //LOG( "ExoticLLP", pDEBUG ) << asts.str();

    double rest_frame_angle = this->Inverted_Fcn( zeta, parent_p4, probe_p4_rest, !forwards );

    // Fantastic. From the rest-frame emission angle, one can calculate the energy.
    double lab_frame_energy = this->CalculateLabFrameEnergy( rest_frame_angle, parent_p4, probe_p4_rest );
    double lab_frame_momentum = std::sqrt( lab_frame_energy * lab_frame_energy - fMass * fMass );

    /*
    LOG( "ExoticLLP", pDEBUG ) << "Dumping stats:"
			       << "\nzeta = " << zeta
			       << "\nfMass = " << fMass
			       << "\nparent_p4 = " << utils::print::P4AsString( &parent_p4 )
			       << "\nprobe_p4_rest = " << utils::print::P4AsString( &probe_p4_rest )
			       << "\npos_soln = " << pos_soln
			       << "\nneg_soln = " << neg_soln
			       << "\nder_pos = " << this->Derivative( pos_soln, parent_p4, probe_p4_rest )
			       << "\nder_neg = " << -1.0 * this->Derivative( neg_soln, parent_p4, probe_p4_rest )
      			       << "\nimage_score = " << image_score
			       << "\nrest_frame_angle = " << rest_frame_angle
			       << "\nlab_frame_energy = " << lab_frame_energy
			       << "\nlab_frame_momentum = " << lab_frame_momentum;
    */

    // and force the momentum to point along the separation unit vector
    TVector3 sep_unit = separation.Unit();
    probe_p4.SetPxPyPzE( lab_frame_momentum * sep_unit.X(),
			 lab_frame_momentum * sep_unit.Y(),
			 lab_frame_momentum * sep_unit.Z(),
			 lab_frame_energy );
    
  } // LLP not produced in top volume

  // Update the boost factor and collimation weights in the FluxContainer
  boost_factor = probe_p4.E() / probe_p4_rest.E();
  fFluxInfo.boost_factor = boost_factor;
  fFluxInfo.wgt_collimation = acc_corr;

  /*
  LOG( "ExoticLLP", pDEBUG ) << "\nBoost factor = " << fFluxInfo.boost_factor
			     << "\nWgt_collimation = " << fFluxInfo.wgt_collimation;
  */

  // Finally, keep track of all the rest-frame kinematics consistent with the flux we just calculated

  TLorentzVector llp_rest_frame_p4_final = probe_p4;
  llp_rest_frame_p4_final.Boost( -parent_p4.BoostVector() ); // Back into parent rest frame
  TVector3 desired_direction = llp_rest_frame_p4_final.Vect().Unit();

  // Feed the randomly generated particles into the system and spit them out
  std::vector< GHepParticle > constrained_rest_stack = 
    this->FindStackKinematics( decayed_results, desired_direction );

  // Sanity check: the zeroth element of this stack should boost back to the p4 of the LLP!
  //TLorentzVector * constrained_llp_p4 = (*constrained_rest_stack.begin()).P4();
  //TLorentzVector * lab_llp_p4 = constrained_llp_p4; lab_llp_p4->Boost( parent_p4.BoostVector() );

  // update the coproduced stuff
  fFluxInfo.cop_pdgs.emplace_back( llp_pdg );
  fFluxInfo.cop_p4xs.emplace_back( llp_rest_frame_p4_final.Px() );
  fFluxInfo.cop_p4ys.emplace_back( llp_rest_frame_p4_final.Py() );
  fFluxInfo.cop_p4zs.emplace_back( llp_rest_frame_p4_final.Pz() );
  fFluxInfo.cop_p4Es.emplace_back( llp_rest_frame_p4_final.E()  );
  for( std::vector<GHepParticle>::iterator it_fss = constrained_rest_stack.begin();
       it_fss != constrained_rest_stack.end(); ++it_fss ) {
    GHepParticle part = *it_fss;
    fFluxInfo.cop_pdgs.emplace_back( part.Pdg() );
    fFluxInfo.cop_p4xs.emplace_back( part.P4()->Px() );
    fFluxInfo.cop_p4ys.emplace_back( part.P4()->Py() );
    fFluxInfo.cop_p4zs.emplace_back( part.P4()->Pz() );
    fFluxInfo.cop_p4Es.emplace_back( part.P4()->E()  );
  }
  
  TLorentzVector probe_v4(0.0, 0.0, 0.0, 0.0); // probe v4 will be updated to event vertex

  GHepParticle ptLLP( llp_pdg, kIStInitialState, -1, -1, -1, -1, probe_p4, probe_v4 );
  evrec->AddParticle( ptLLP );
}
//____________________________________________________________________________
std::vector<GHepParticle> FluxCreator::FindStackKinematics( std::vector<GHepParticle> random_system,
							    TVector3 desired_direction ) const
{
  // The 0th element of the vector is the LLP
  GHepParticle llp_particle = *random_system.begin();
  TVector3 llp_start_momentum = llp_particle.P4()->Vect();
  TVector3 uvec = llp_start_momentum.Unit();
  
  // Assume desired_direction is in the rest frame; the boost should be handled before this method
  TVector3 vvec = desired_direction.Unit();

  // And construct the "half-vector" to do a quaternion rotation;
  // see https://stackoverflow.com/questions/1171849/
  /*
   * In the axis-angle representation, the quaternion is (C, XS, YS, ZS) where C, S = cos(sin)(th/2);
   * and (X, Y, Z) define the axis of the rotation by which one rotates.
   */
  TVector3 wvec = (uvec + vvec).Unit();

  double   dot   = uvec.Dot(wvec);
  TVector3 cross = uvec.Cross(wvec);
  double   crx   = cross.X();
  double   cry   = cross.Y();
  double   crz   = cross.Z();

  // and from that we can construct the rotation matrix!
  double mat_xx = dot*dot + crx*crx - cry*cry - crz*crz;
  double mat_xy = 2.0 * (crx*cry - dot*crz);
  double mat_xz = 2.0 * (crx*crz + dot*cry);
  double mat_yx = 2.0 * (crx*cry + dot*crz);
  double mat_yy = dot*dot - crx*crx + cry*cry - crz*crz;
  double mat_yz = 2.0 * (cry*crz - dot*crx);
  double mat_zx = 2.0 * (crx*crz - dot*cry);
  double mat_zy = 2.0 * (cry*crz + dot*crx);
  double mat_zz = dot*dot - crx*crx - cry*cry + crz*crz;
  
  // now, for every particle in the input stack, apply the rotation and get the output.
  std::vector<GHepParticle> final_system = random_system;
  for( std::vector<GHepParticle>::iterator it_fss = final_system.begin();
       it_fss != final_system.end(); ++it_fss ) {
    TLorentzVector * this_p4 = (*it_fss).P4();
    TVector3 this_direction = (this_p4->Vect()).Unit(); double mom = this_p4->P();
    
    double ox = this_direction.X(), oy = this_direction.Y(), oz = this_direction.Z();
    double tx = mat_xx * ox + mat_xy * oy + mat_xz * oz;
    double ty = mat_yx * ox + mat_yy * oy + mat_yz * oz;
    double tz = mat_zx * ox + mat_zy * oy + mat_zz * oz;

    tx *= mom; ty *= mom; tz *= mom;
    
    (*it_fss).SetMomentum( tx, ty, tz, this_p4->E() ); // Note this is a REST FRAME momentum.
  }

  return final_system;
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

  //LOG( "ExoticLLP", pDEBUG ) << "About to retrieve LLP from LLPConfigurator";
  fExoticLLP = LLP_configurator->RetrieveLLP();
  //LOG( "ExoticLLP", pDEBUG ) << "Successfully retrieved LLP from LLPConfigurator";
  fMass = fExoticLLP.GetMass();

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
    } else {
      std::vector< ModeObject > mobj_vec = (*mmap).second;
      mobj_vec.emplace_back( mobj );
      fGroupedModes[ parent_pdg ] = mobj_vec; // update the vector
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

  //LOG( "ExoticLLP", pDEBUG ) << "FluxCreator configured.";
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

  /*
  LOG( "ExoticLLP", pDEBUG )
    << "Starting to add files to input. Current size is " << dirs.size();
  */

  while( dirs.size() > 0 ){ // there is still stuff we haven't looked at.
    int nNow = dirs.size();
    /*
    LOG( "ExoticLLP", pDEBUG ) 
      << "Scanning directory " << dirNames.front() << " with " << nNow << " elements...";
    */

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

  double velocity_ratio = bNu / bPar;
  //double tanMaxTheta = ( bNu >= bPar ) ? 180.0 : 1.0 / ( gPar * std::sqrt( ( bPar / bNu ) * ( bPar / bNu ) - 1.0 ) );
  double tanMaxTheta = ( bNu >= bPar ) ? 10000.0 : velocity_ratio / ( gPar * std::sqrt( 1.0 - velocity_ratio ) );
  double maxTheta = ( bNu >= bPar ) ? 180.0 : TMath::ATan( tanMaxTheta ) * TMath::RadToDeg();
  if( maxTheta < 0.0 ) maxTheta += 180.0;

  /*
  LOG( "ExoticLLP", pDEBUG ) 
    << "\nthetalab = " << thetalab << " [deg]"
    << "\nLLP mass = " << mass << " [GeV]"
    << "\nEPar, MPar = " << EPar << ", " << MPar << " [GeV]"
    << "\nENu = " << ENu
    << "\nisPos = " << isPos
    << "\n~*~*~*~*~*~*~"
    << "\npPar, bPar, gPar = " << pPar << " [GeV], " << bPar << ", " << gPar
    << "\nqNu, bNu = " << qNu << " [GeV], " << bNu
    << "\n~*~*~*~*~*~*~"
    << "\ntanMaxTheta = " << tanMaxTheta
    << "\nmaxTheta = " << maxTheta;
  */

  double arg = 0.0;
  if( isPos ){ // positive solution
    if( bNu >= bPar ){ // labangle is surjective
      arg = this->AccCorr_SolnArgs( thetalab, mass, EPar, MPar, ENu, isPos );
      double tsol = TMath::ATan( arg ) * TMath::RadToDeg();
      if( tsol < 0.0 ) tsol += 180.0;
      return tsol;
    } else { // labangle is not surjective, check if preimage can be found
      if( thetalab > maxTheta ) {
	LOG( "ExoticLLP", pDEBUG ) << "Too much deviation! Can't be accepted";
	return -1.0;
      }

      arg = this->AccCorr_SolnArgs( thetalab, mass, EPar, MPar, ENu, isPos );
      double tsol = TMath::ATan( arg ) * TMath::RadToDeg();
      if( tsol < 0.0 ) tsol += 180.0;
      return tsol;
    }
  } else { // negative solution
    if( bNu >= bPar ){ // LLP is too fast, so labangle is monotonic and no negative solution
      return 0.0;
    } else { // check if preimage can be found
      if( thetalab > maxTheta ) {
	LOG( "ExoticLLP", pDEBUG ) << "Too much deviation! Can't be accepted";
	return -1.0;
      }

      arg = this->AccCorr_SolnArgs( thetalab, mass, EPar, MPar, ENu, isPos );
      double tsol = TMath::ATan( arg ) * TMath::RadToDeg();
      if( tsol < 0.0 ) tsol += 180.0;
      return tsol;
    }
  }

  return -1.0; // you should never see this.
}
//----------------------------------------------------------------------------
double FluxCreator::CalculateLabFrameEnergy( double Theta, 
					     TLorentzVector p4par, TLorentzVector p4LLP ) const
{
  /*
   * Given a rest-frame emission angle, the rest-frame 4-momentum of the LLP, and the lab-frame 
   * 4-momentum of the parent, one can calculate the lab-frame energy of the LLP.
   */

  double beta  = p4par.Beta();
  double gamma = p4par.Gamma();
  
  double E = p4LLP.E();
  double p = p4LLP.P();

  double ctheta = std::cos( Theta * TMath::DegToRad() );
  
  return gamma * ( E + beta * p * ctheta );
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

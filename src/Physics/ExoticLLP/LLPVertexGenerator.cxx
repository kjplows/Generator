//____________________________________________________________________________
/*
  Copyright (c) 2003-2023, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  
  Author: John Plows <kplows \at liverpool.ac.uk>
          University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/LLPVertexGenerator.h"

using namespace genie;
using namespace genie::llp;
using namespace genie::units;

//____________________________________________________________________________
VertexGenerator::VertexGenerator() :
  Algorithm("genie::llp::VertexGenerator")
{

}
//____________________________________________________________________________
VertexGenerator::VertexGenerator(string name) :
  Algorithm(name)
{

}
//____________________________________________________________________________
VertexGenerator::VertexGenerator(string name, string config) :
  Algorithm(name, config)
{

}
//____________________________________________________________________________
VertexGenerator::~VertexGenerator()
{

}
//____________________________________________________________________________
void VertexGenerator::ReadFluxContainer( const FluxContainer flc ) const
{
  fFluxContainer = flc;
}
//____________________________________________________________________________
void VertexGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
  LOG( "ExoticLLP", pDEBUG ) << "Processing vertex positioning...";
  std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
  fDecayPoint = TLorentzVector(0.0, 0.0, 0.0, 0.0);

  // get the entry and exit points from the volume, for a ray
  // first, make sure the VolumeSeeker is appropriately populated
  VolumeSeeker * vsek = VolumeSeeker::Instance();
  vsek->ClearEvent();
  vsek->PopulateEvent( fFluxContainer.v4.Vect(), fFluxContainer.p4.Vect() ); // only in NEAR 
  
  // Raytrace the detector and get out the entry and exit points
  double reached_raytrace = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count();
  //vsek->RaytraceDetector(); // this is slow
  vsek->IntersectDetector(); // this is less slow and as correct
  double finished_raytrace = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count() - reached_raytrace;
  TVector3 entry_point = vsek->GetEntryPoint(); // USER, m
  TVector3 exit_point  = vsek->GetExitPoint(); // USER, m

  /*
  LOG( "ExoticLLP", pDEBUG ) << "Here are the entry and exit points..."
			     << "\nUser origin: " << utils::print::X4AsString( &fFluxContainer.v4_user )
			     << " [m]"
			     << "\nUser mom:    " << utils::print::P4AsString( &fFluxContainer.p4_user )
			     << " [GeV]"
			     << "\nEntry point: " << utils::print::Vec3AsString( &entry_point )
			     << " [m]"
			     << "\nExit  point: " << utils::print::Vec3AsString( &exit_point )
			     << " [m]";
  */

  // Yay! Now we can use a lifetime to get a vertex
  // Practically, it is a Uniform number.
  // RETHERE: Perhaps we should save the random number too, to cast from tau1 |--> tau2...

  TVector3 ray_in_detector = exit_point - entry_point; // m
  double max_length = ray_in_detector.Mag() * units::m / lunits; // lunits

  // Get the LLP configurator and from that the ExoticLLP
  const Algorithm * algLLPConfigurator = AlgFactory::Instance()->GetAlgorithm("genie::llp::LLPConfigurator", "Default");
  const LLPConfigurator * LLP_configurator = dynamic_cast< const LLPConfigurator * >( algLLPConfigurator );

  fLifetime = LLP_configurator->RetrieveLLP().GetLifetime() * units::m / lunits; // lunits
  
  // Map the uniform number to an exponential
  
  double beta  = fFluxContainer.p4_user.Beta();
  double gamma = fFluxContainer.p4_user.Gamma();
  double kappa = 1.0 / ( beta * gamma * fLifetime ); // (lunits)^{-1}

  double uniform = RandomGen::Instance()->RndGen().Rndm();

  //double elapsed_length = uniform * max_length;

  double elapsed_length = -std::log( 1.0 - uniform * (1.0 - std::exp( -kappa * max_length )) ) / kappa;
  elapsed_length *= lunits / units::m; // m

  // now place the decay vertex
  TVector3 vtx = entry_point + elapsed_length * ray_in_detector.Unit(); // m
  
  /*
  LOG( "ExoticLLP", pDEBUG ) << "\nWith max_length = " << max_length 
			     << " [ " << lunitString << " ]"
			     << "\nand uniform number = " << uniform
			     << "\nand beta, gamma, kappa = " << beta << ", " << gamma 
			     << ", " << kappa 
			     << " [ 1 / " << lunitString << " ]"
			     << "\nwe got an elapsed length = " << elapsed_length 
			     << " [ " << lunitString << " ]"
			     << "\nand therefore put the vertex at "
			     << utils::print::Vec3AsString( &vtx )
			     << " [ m ]";
  */

  double reached_tof = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count() - finished_raytrace;

  // We can also calculate the time of arrival of this LLP
  // We'll use the parent decay vertex time as an input, and add the time-of-flight on top

  TLorentzVector v4_origin = fFluxContainer.v4_user; // m, ns
  TVector3 ray_to_detector = entry_point - v4_origin.Vect(); // m
  double t0 = v4_origin.T() * units::ns / tunits; // tunits

  double velocity = beta * kNewSpeedOfLight; // lunits / tunits
  double full_distance = ray_to_detector.Mag() * units::m / lunits + elapsed_length; // get to detector, and travel in it a bit. lunits
  double tof = full_distance / velocity; // tunits

  double full_time = (t0 + tof) * tunits / units::ns; // ns
  // also update the times for the entry and exit points
  double entry_time = ( t0 + ray_to_detector.Mag() * (units::m / lunits) / 
			velocity ) * tunits / units::ns;
  double exit_time  = ( t0 + ( ray_to_detector.Mag() + max_length ) * (units::m / lunits) /
			velocity ) * tunits / units::ns;
  
  TLorentzVector entry_v4_user( entry_point.X(), entry_point.Y(), entry_point.Z(), entry_time );
  TVector3 entry_point_near = vsek->RotateToNear( entry_point );
  entry_point_near = vsek->TranslateToNear( entry_point_near );
  TLorentzVector entry_v4( entry_point_near.X(), entry_point_near.Y(), entry_point_near.Z(), entry_time );

  TLorentzVector exit_v4_user( exit_point.X(), exit_point.Y(), exit_point.Z(), exit_time );
  TVector3 exit_point_near = vsek->RotateToNear( exit_point );
  exit_point_near = vsek->TranslateToNear( exit_point_near );
  TLorentzVector exit_v4( exit_point_near.X(), exit_point_near.Y(), exit_point_near.Z(), exit_time );

  fFluxContainer.entry = entry_v4; fFluxContainer.entry_user = entry_v4_user;
  fFluxContainer.exit  = exit_v4;  fFluxContainer.exit_user  = exit_v4_user;

  // now set the final decay point
  fDecayPoint.SetXYZT( vtx.X(), vtx.Y(), vtx.Z(), full_time ); // m, ns

  double finished_decayPoint = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count() - reached_tof;

  /*
  LOG( "ExoticLLP", pDEBUG ) << "\nThe velocity v = beta * c = " << velocity 
			     << " [ " << lunitString << " / " << tunitString << " ]";
  LOG( "ExoticLLP", pDEBUG ) << "\nFrom a full distance of "
			     << "\nto detector: " << ray_to_detector.Mag()
			     << " [ " << lunitString << " ]"
			     << "\nin detector: " << elapsed_length
			     << " [ " << lunitString << " ]"
			     << "\nwe get a tof = " << tof
			     << " [ " << tunitString << " ]"
			     << "\nand from an origin time " << t0
			     << " [ " << tunitString << " ]"
			     << "\nso the full time is " << full_time
			     << " [ ns ]";
  LOG( "ExoticLLP", pDEBUG ) << "\nThe final decay point is " << utils::print::X4AsString( &fDecayPoint );
  */

  // also get the decay point in NEAR
  TVector3 decay_point_near = vsek->RotateToNear( vtx );
  decay_point_near = vsek->TranslateToNear( decay_point_near );
  TLorentzVector decay_v4_near( decay_point_near.X(), decay_point_near.Y(), decay_point_near.Z(), full_time );
  fFluxContainer.decay = decay_v4_near; fFluxContainer.decay_user = fDecayPoint;
  fFluxContainer.vtx_rng = uniform;

  event_rec->SetVertex( fDecayPoint );
  fFluxContainer.wgt_survival = std::exp( -ray_to_detector.Mag() / ( beta * gamma * fLifetime ) );
  fFluxContainer.wgt_detdecay = 1.0 - std::exp( -max_length / ( beta * gamma * fLifetime ) );

  double end_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count() - finished_decayPoint;

  //LOG( "ExoticLLP", pDEBUG ) << fFluxContainer;

  LOG( "ExoticLLP", pDEBUG ) 
    << "\nTook " << reached_raytrace << " ms to reach raytrace"
    << "\nTook " << finished_raytrace << " ms to finish raytrace"
    << "\nTook " << reached_tof << " ms to reach ToF"
    << "\nTook " << finished_decayPoint << " ms to finish decayPoint"
    << "\nTook " << end_time << " ms to wrap up";
}
//____________________________________________________________________________
double VertexGenerator::CalcTravelLength( double betaMag, double CoMLifetime, double maxLength ) const
{
  // decay probability P0(t) = 1 - exp( -t/tau ) where:
  // t   = time-of-flight (in rest frame)
  // tau = CoMLifetime

  assert( betaMag > 0.0 && betaMag < 1.0 && "LLP is massive and moving" ); // massive moving particle
  double maxLabTime = maxLength / ( betaMag * kNewSpeedOfLight );
  double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  double maxRestTime = maxLabTime / gamma ; // this is how "wide" the detector looks

  // if P(DL=0) = 1, P(DL = LMax) = exp( - LMax / c * 1/( beta * gamma ) * 1 / CoMLifetime )
  double PExit = std::exp( - maxRestTime / CoMLifetime );

  // from [0,1] we'd reroll anything in [0, PExit] and keep (PExit, 1]. That's expensive.
  // Instead, let 1 ==> maxRestTime, 0 ==> 0, exponential decay
  
  RandomGen * rnd = RandomGen::Instance();
  double ranthrow = rnd->RndGen().Uniform();

  double S0 = (1.0 - PExit) * ranthrow + PExit; 
  double rest_time = CoMLifetime * std::log( 1.0 / S0 );
  double elapsed_time = rest_time * gamma;
  double elapsed_length = elapsed_time * betaMag * kNewSpeedOfLight;

  return elapsed_length;
}
//____________________________________________________________________________
double VertexGenerator::GetMaxLength() const
{
  TLorentzVector entryPoint = fFluxContainer.entry_user;
  TLorentzVector exitPoint  = fFluxContainer.exit_user;

  double ex = entryPoint.X(); double ey = entryPoint.Y(); double ez = entryPoint.Z();
  double xx = exitPoint.X(); double xy = exitPoint.Y(); double xz = exitPoint.Z();

  return std::sqrt( (ex-xx)*(ex-xx) + (ey-xy)*(ey-xy) + (ez-xz)*(ez-xz) );
}
//____________________________________________________________________________
TLorentzVector VertexGenerator::GetDecayPoint( double travelLength ) const
{

  TLorentzVector entryPoint = fFluxContainer.entry_user;
  TLorentzVector momentum = fFluxContainer.p4;

  double ex = entryPoint.X(); double ey = entryPoint.Y(); double ez = entryPoint.Z();
  double et = entryPoint.T();
  double px = momentum.X(); double py = momentum.Y(); double pz = momentum.Z();
  double p2 = px*px + py*py + pz*pz; double p = std::sqrt(p2);
  px *= 1./p; py *= 1./p; pz *= 1./p;

  double dx = ex + travelLength * px;
  double dy = ey + travelLength * py;
  double dz = ez + travelLength * pz;

  double dt = et + travelLength / ( kNewSpeedOfLight * momentum.P() / momentum.E() );

  TLorentzVector decayPoint( dx, dy, dz, dt );
  return decayPoint;
}
//____________________________________________________________________________
void VertexGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void VertexGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void VertexGenerator::LoadConfig()
{
  if( fIsConfigLoaded ) return;

  LOG( "ExoticLLP", pDEBUG )
    << "Loading geometry parameters from file. . .";

  fIsConfigLoaded = true;
}
//____________________________________________________________________________
void VertexGenerator::GetInterestingPoints( TLorentzVector & entryPoint, TLorentzVector & exitPoint, 
					    TLorentzVector & decayPoint ) const
{
  entryPoint = fFluxContainer.entry_user;
  exitPoint  = fFluxContainer.exit_user;
  decayPoint = fDecayPoint;
}
/*
//____________________________________________________________________________
void VertexGenerator::SetGeomFile( string geomfile, string topVolume ) const
{
  fGeomFile = geomfile;
  fTopVolume = topVolume;
}
*/

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
void VertexGenerator::ReadFluxContainer( FluxContainer flc ) const
{
  FluxContainer * pt_flc = new FluxContainer( flc );
  fFluxContainer = *pt_flc;
  delete pt_flc;
}
//____________________________________________________________________________
void VertexGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
  /*
  // before anything else: find the geometry!
  if( !fGeoManager ){
    LOG( "ExoticLLP", pINFO )
      << "Getting geometry information from " << fGeomFile;

    fGeoManager = TGeoManager::Import(fGeomFile.c_str());
    
    TGeoVolume * main_volume = fGeoManager->GetTopVolume();
    TGeoVolume * top_volume = fGeoManager->GetVolume( fTopVolume.c_str() );
    assert( top_volume && "Top volume exists" );
    // now get the translation of the top volume
    if( main_volume != top_volume ) {
      //main_volume->FindMatrixOfDaughterVolume(top_volume);
      //TGeoHMatrix * hmat = fGeoManager->GetHMatrix();
      TGeoMatrix * hmat = this->FindFullTransformation( main_volume, top_volume );
      const Double_t * tran = hmat->GetTranslation();
      fTx = tran[0] * units::cm / units::m;
      fTy = tran[1] * units::cm / units::m;
      fTz = tran[2] * units::cm / units::m;
      LOG( "ExoticLLP", pDEBUG )
	<< "Got translation of volume with name " << top_volume->GetName() << " which is ( " 
	<< fTx << ", " << fTy << ", " << fTz << " ) [m]";
    }
    fGeoManager->SetTopVolume(top_volume);
    TGeoShape * ts = top_volume->GetShape();
    TGeoBBox * box = (TGeoBBox *) ts;
    
    this->ImportBoundingBox( box );
  }

  this->SetStartingParameters( event_rec );
  LOG( "ExoticLLP", pDEBUG ) << "Starting parameters SET.";

  double weight = 1.0; // pure geom weight

  TVector3 startPoint, momentum, entryPoint, exitPoint; // USER mm, GeV/GeV
  startPoint.SetXYZ( fSx, fSy, fSz );
  momentum.SetXYZ( fPx, fPy, fPz );
  
  bool didIntersectDet = this->VolumeEntryAndExitPoints( startPoint, momentum, entryPoint, exitPoint, fGeoManager, fGeoVolume );

  if( !isParticleGun && isUsingDk2nu ) assert( didIntersectDet && "Forced to hit detector somewhere" ); // forced to hit detector somewhere!
  else {
    std::vector< double > * newProdVtx = new std::vector< double >();
    newProdVtx->emplace_back( startPoint.X() );
    newProdVtx->emplace_back( startPoint.Y() );
    newProdVtx->emplace_back( startPoint.Z() );

  }
  if( !didIntersectDet ){ // bail
    LOG( "ExoticLLP", pERROR )
      << "Bailing...";
    TLorentzVector v4dummy( -999.9, -999.9, -999.9, -999.9 );
    event_rec->SetVertex( v4dummy );
    return;
  }

  this->EnforceUnits( "mm", "rad", "ns" );

  // move fCoMLifetime to ns from GeV^{-1}
  fCoMLifetime *= 1.0 / ( units::ns * units::GeV );

  double maxDx = exitPoint.X() - entryPoint.X();
  double maxDy = exitPoint.Y() - entryPoint.Y();
  double maxDz = exitPoint.Z() - entryPoint.Z(); // mm

  double maxLength = std::sqrt( std::pow( maxDx , 2.0 ) +
				std::pow( maxDy , 2.0 ) +
				std::pow( maxDz , 2.0 ) ); // mm
  
  TLorentzVector * p4HNL = event_rec->Particle(0)->GetP4();
  double betaMag = p4HNL->P() / p4HNL->E();
  double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  
  double elapsed_length = this->CalcTravelLength( betaMag, fCoMLifetime, maxLength ); //mm
  __attribute__((unused)) double ratio_length = elapsed_length / maxLength;

  TVector3 HNLOriginPoint( fNx * units::m / units::mm, fNy * units::m / units::mm, fNz * units::m / units::mm ); // USER mm
  
  // from these we can also make the weight. It's P( survival ) * P( decay in detector | survival )
  double distanceBeforeDet = std::sqrt( std::pow( (HNLOriginPoint.X() - entryPoint.X()), 2.0 ) + 
					std::pow( (HNLOriginPoint.Y() - entryPoint.Y()), 2.0 ) + 
					std::pow( (HNLOriginPoint.Z() - entryPoint.Z()), 2.0 ) ); // mm
  
  double timeBeforeDet = distanceBeforeDet / ( betaMag * kNewSpeedOfLight ); // ns lab
  double timeInsideDet = maxLength / ( betaMag * kNewSpeedOfLight ); // ns lab
  
  double LabToRestTime = 1.0 / ( gamma );
  timeBeforeDet *= LabToRestTime; // ns rest
  timeInsideDet *= LabToRestTime; // ns rest
  
  double survProb = std::exp( - timeBeforeDet / fCoMLifetime );
  weight *= 1.0 / survProb;
  double decayProb = 1.0 - std::exp( - timeInsideDet / fCoMLifetime );
  weight *= 1.0 / decayProb;

  LOG( "ExoticLLP", pDEBUG ) 
    << "\nOutput of lifetime calc:"
    << "\nLifetime = " << fCoMLifetime << " [CoM ns]"
    << "\nDistance to detector = " << distanceBeforeDet * units::mm / units::m
    << " : inside detector = " << maxLength * units::mm / units::m << " [m]"
    << "\nTime before detector = " << timeBeforeDet / LabToRestTime
    << " : inside detector = " << timeInsideDet / LabToRestTime << " [LAB ns]"
    << "\nTime before detector = " << timeBeforeDet
    << " : inside detector = " << timeInsideDet << " [CoM ns]"
    << "\n"
    << "PSurv = " << survProb
    << " : PDec = " << decayProb;

  // save the survival and decay probabilities
  if( event_rec->Particle(1) && event_rec->Particle(2) ){
    event_rec->Particle(1)->SetPosition( 0.0, 0.0, 0.0, survProb );
    event_rec->Particle(2)->SetPosition( 0.0, 0.0, 0.0, decayProb );
  }

  // update the weight
  event_rec->SetWeight( event_rec->Weight() * weight );

  TVector3 decayPoint = this->GetDecayPoint( elapsed_length, entryPoint, momentum ); // USER, mm

  // write out vtx in [m, ns]
  TLorentzVector x4( decayPoint.X() * units::mm / units::m,
		     decayPoint.Y() * units::mm / units::m,
		     decayPoint.Z() * units::mm / units::m,
		     event_rec->Vertex()->T() );

  event_rec->SetVertex(x4);

  // the validation app doesn't run the Decayer. So we will insert two neutrinos (not a valid
  // decay mode), to store entry and exit point
  if( !isUsingDk2nu ){
    LOG( "ExoticLLP", pDEBUG ) << "About to insert two neutrinos into a NULL event record";
    assert( !event_rec->Particle(1) && "Event record only has HNL if gevald_hnl -M 3" );
    
    TLorentzVector tmpp4( 0.0, 0.0, 0.0, 0.5 );
    TLorentzVector ex4( 0.0, 0.0, 0.0, 0.0 );
    ex4.SetXYZT( entryPoint.X(), entryPoint.Y(), entryPoint.Z(), 0.0 );
    TLorentzVector xx4( 0.0, 0.0, 0.0, 0.0 );
    xx4.SetXYZT( exitPoint.X(), exitPoint.Y(), exitPoint.Z(), 0.0 );

    GHepParticle nu1( genie::kPdgNuMu, kIStStableFinalState, -1, -1, -1, -1, tmpp4, ex4 );
    GHepParticle nu2( genie::kPdgAntiNuMu, kIStStableFinalState, -1, -1, -1, -1, tmpp4, xx4 );

    event_rec->AddParticle( nu1 ); event_rec->AddParticle( nu2 );

    // save the survival and decay probabilities
    // event_rec->Particle(1)->SetPolarization( survProb, decayProb );
    event_rec->Particle(1)->SetPosition( 0.0, 0.0, 0.0, survProb );
    event_rec->Particle(2)->SetPosition( 0.0, 0.0, 0.0, decayProb );
    event_rec->SetWeight(weight);
  }

  // also set entry and exit points. Do this in x4 of Particles(1,2)
  if( event_rec->Particle(1) && event_rec->Particle(2) ){
    //(event_rec->Particle(1))->SetPosition( entryPoint.X(), entryPoint.Y(), entryPoint.Z(), event_rec->Particle(1)->Vt() );
    //(event_rec->Particle(2))->SetPosition( exitPoint.X(), exitPoint.Y(), exitPoint.Z(), event_rec->Particle(2)->Vt() );
    (event_rec->Particle(1))->SetPosition( entryPoint.X(), entryPoint.Y(), entryPoint.Z(), survProb );
    (event_rec->Particle(2))->SetPosition( exitPoint.X(), exitPoint.Y(), exitPoint.Z(), decayProb );
  }

  delete p4HNL;
  */
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

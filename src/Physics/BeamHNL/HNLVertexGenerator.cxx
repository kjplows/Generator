//____________________________________________________________________________
/*
  Copyright (c) 2003-2023, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  
  Author: John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford
*/
//____________________________________________________________________________

#include "Physics/BeamHNL/HNLVertexGenerator.h"

using namespace genie;
using namespace genie::hnl;
using namespace genie::units;

//____________________________________________________________________________
VertexGenerator::VertexGenerator() :
  GeomRecordVisitorI("genie::hnl::VertexGenerator")
{

}
//____________________________________________________________________________
VertexGenerator::VertexGenerator(string name) :
  GeomRecordVisitorI(name)
{

}
//____________________________________________________________________________
VertexGenerator::VertexGenerator(string name, string config) :
  GeomRecordVisitorI(name, config)
{

}
//____________________________________________________________________________
VertexGenerator::~VertexGenerator()
{

}
//____________________________________________________________________________
void VertexGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
  /*!
   *  Uses ROOT's TGeoManager to find out where the intersections with the detector volume live
   *  Label them as entry and exit point. Then use them to determine:
   *  1) A decay vertex within the detector
   *  2) A time-of-decay (== delay of HNL to reach the decay vertex wrt a massless SM v)
   *  3) Geom weight: Survival to detector * decay within detector.
   */

  // before anything else: find the geometry!
  if( !fGeoManager ){
    LOG( "HNL", pINFO )
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
      LOG( "HNL", pDEBUG )
	<< "Got translation of volume with name " << top_volume->GetName() << " which is ( " 
	<< fTx << ", " << fTy << ", " << fTz << " ) [m]";
    }
    fGeoManager->SetTopVolume(top_volume);
    TGeoShape * ts = top_volume->GetShape();
    TGeoBBox * box = (TGeoBBox *) ts;
    
    this->ImportBoundingBox( box );
    
    //this->ImportBoundingBox(box);
  }

  this->SetStartingParameters( event_rec );

  double weight = 1.0; // pure geom weight

  TVector3 startPoint, momentum, entryPoint, exitPoint;
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
    LOG( "HNL", pERROR )
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
				std::pow( maxDz , 2.0 ) );
  
  TLorentzVector * p4HNL = event_rec->Particle(0)->GetP4();
  double betaMag = p4HNL->P() / p4HNL->E();
  double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  
  double elapsed_length = this->CalcTravelLength( betaMag, fCoMLifetime, maxLength ); //mm
  __attribute__((unused)) double ratio_length = elapsed_length / maxLength;
  
  // from these we can also make the weight. It's P( survival ) * P( decay in detector | survival )
  double distanceBeforeDet = std::sqrt( std::pow( (entryPoint.X() - startPoint.X()), 2.0 ) + 
					std::pow( (entryPoint.Y() - startPoint.Y()), 2.0 ) + 
					std::pow( (entryPoint.Y() - startPoint.Z()), 2.0 ) ); // mm
  
  double timeBeforeDet = distanceBeforeDet / ( betaMag * kNewSpeedOfLight ); // ns lab
  double timeInsideDet = maxLength / ( betaMag * kNewSpeedOfLight ); // ns lab
  
  double LabToRestTime = 1.0 / ( gamma );
  timeBeforeDet *= LabToRestTime; // ns rest
  timeInsideDet *= LabToRestTime; // ns rest
  
  double survProb = std::exp( - timeBeforeDet / fCoMLifetime );
  weight *= 1.0 / survProb;
  double decayProb = 1.0 - std::exp( - timeInsideDet / fCoMLifetime );
  weight *= 1.0 / decayProb;

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
    (event_rec->Particle(1))->SetPosition( entryPoint.X(), entryPoint.Y(), entryPoint.Z(), event_rec->Particle(1)->Vt() );
    (event_rec->Particle(2))->SetPosition( exitPoint.X(), exitPoint.Y(), exitPoint.Z(), event_rec->Particle(2)->Vt() );
  }

  delete p4HNL;
}
//____________________________________________________________________________
void VertexGenerator::EnforceUnits( std::string length_units, std::string angle_units, std::string time_units ) const{

  double old_lunits = lunits;
  __attribute__((unused)) double old_aunits = aunits;
  double old_tunits = tunits;

  lunits = utils::units::UnitFromString( length_units ); lunitString = length_units;
  aunits = utils::units::UnitFromString( angle_units );
  tunits = utils::units::UnitFromString( time_units ); tunitString = time_units;

  LOG( "HNL", pWARN )
    << "Switching units:"
    << "\nTo:   " << length_units.c_str() << " , " << angle_units.c_str() << " , " << time_units.c_str()
    << "\nConversion factors: " << lunits/old_lunits << ", " << aunits / old_aunits << ", " << tunits / old_tunits;

  // convert to new units
  fSx /= lunits/old_lunits; fSy /= lunits/old_lunits; fSz /= lunits/old_lunits;
  fPx /= lunits/old_lunits; fPy /= lunits/old_lunits; fPz /= lunits/old_lunits;
  fEx /= lunits/old_lunits; fEy /= lunits/old_lunits; fEz /= lunits/old_lunits;
  fXx /= lunits/old_lunits; fXy /= lunits/old_lunits; fXz /= lunits/old_lunits;
  fLx /= lunits/old_lunits; fLy /= lunits/old_lunits; fLz /= lunits/old_lunits;

  fDx /= lunits/old_lunits; fDy /= lunits/old_lunits; fDz /= lunits/old_lunits;
  fOx /= lunits/old_lunits; fOy /= lunits/old_lunits; fOz /= lunits/old_lunits;

  kNewSpeedOfLight /= (lunits / old_lunits) / (tunits / old_tunits);

  LOG( "HNL", pDEBUG )
    << "kNewSpeedOfLight = " << kNewSpeedOfLight << " [" << lunitString.c_str() << "/"
    << tunitString.c_str() << "]";
}
//____________________________________________________________________________
double VertexGenerator::CalcTravelLength( double betaMag, double CoMLifetime, double maxLength ) const
{
  // decay probability P0(t) = 1 - exp( -t/tau ) where:
  // t   = time-of-flight (in rest frame)
  // tau = CoMLifetime

  assert( betaMag > 0.0 && betaMag < 1.0 && "HNL is massive and moving" ); // massive moving particle
  double maxLabTime = maxLength / ( betaMag * kNewSpeedOfLight );
  double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
  double maxRestTime = maxLabTime / gamma ; // this is how "wide" the detector looks like

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

  /*
  LOG( "HNL", pDEBUG )
    << "\nbetaMag, maxLength, CoMLifetime = " << betaMag << ", " << maxLength << ", " << CoMLifetime
    << "\nbetaMag = " << betaMag << " ==> gamma = " << gamma
    << "\n==> maxLength [" << tunitString.c_str()
    << "] = " << maxRestTime << " (rest frame) = " << maxLabTime << " (lab frame)"
    << "\nranthrow = " << ranthrow << ", PExit = " << PExit
    << "\n==> S0 = " << S0 << " ==> rest_time [" << lunitString.c_str() << "] = " << rest_time
    << " ==> elapsed_time [" << tunitString.c_str()
    << "] = " << elapsed_time << " ==> elapsed_length [" << lunitString.c_str()
    << "] = " << elapsed_length;
  */

  return elapsed_length;
}
//____________________________________________________________________________
TVector3 VertexGenerator::GetDecayPoint( double travelLength, TVector3 & entryPoint, TVector3 & momentum ) const
{
  double ex = entryPoint.X(); double ey = entryPoint.Y(); double ez = entryPoint.Z();
  double px = momentum.X(); double py = momentum.Y(); double pz = momentum.Z();
  double p2 = px*px + py*py + pz*pz; double p = std::sqrt(p2);
  px *= 1./p; py *= 1./p; pz *= 1./p;

  double dx = ex + travelLength * px; fDx = dx;
  double dy = ey + travelLength * py; fDy = dy;
  double dz = ez + travelLength * pz; fDz = dz;

  fDxROOT = fDx * lunits / units::cm;
  fDyROOT = fDy * lunits / units::cm;
  fDzROOT = fDz * lunits / units::cm;

  TVector3 decayPoint( dx, dy, dz );
  return decayPoint;
}
//____________________________________________________________________________
double VertexGenerator::GetMaxLength( TVector3 & entryPoint, TVector3 & exitPoint ) const
{
  double ex = entryPoint.X(); double ey = entryPoint.Y(); double ez = entryPoint.Z();
  double xx = exitPoint.X(); double xy = exitPoint.Y(); double xz = exitPoint.Z();

  return std::sqrt( (ex-xx)*(ex-xx) + (ey-xy)*(ey-xy) + (ez-xz)*(ez-xz) );
}
//____________________________________________________________________________
void VertexGenerator::MakeSDV() const
{
  fOx = 0.0; fOy = 0.0; fOz = 0.0;
  fLx = 1.0; fLy = 1.0; fLz = 1.0; // m

  lunits = utils::units::UnitFromString( "m" );
  aunits = utils::units::UnitFromString( "rad" );
  tunits = utils::units::UnitFromString( "ns" );
  
  kNewSpeedOfLight = genie::units::kSpeedOfLight 
  * (genie::units::m / lunits)
  / (genie::units::s / tunits);

  LOG("HNL", pDEBUG)
    << "Setting simple decay volume with unit-m side."
    << "\nSetting units to \"mm\", \"rad\", \"ns\"";

  EnforceUnits("mm","rad","ns");
}
//____________________________________________________________________________
// if entry and exit points, populate TVector3's with their coords. If not, return false
bool VertexGenerator::SDVEntryAndExitPoints( TVector3 & startPoint, TVector3 momentum,
					    TVector3 & entryPoint, TVector3 & exitPoint ) const
{
  assert( fOx == 0.0 && fOy == 0.0 && fOz == 0.0 && 
	  fLx == 1000.0 && fLy == 1000.0 && fLz == 1000.0 &&
	  "Decay volume is unit-m side, centred at origin"); // SDV, mm
  fSx = startPoint.X(); fSy = startPoint.Y(); fSz = startPoint.Z(); // mm
  fPx = momentum.X(); fPy = momentum.Y(); fPz = momentum.Z(); // GeV
  double fP2 = fPx*fPx + fPy*fPy + fPz*fPz; double fP = std::sqrt(fP2); // GeV
  fPx *= 1.0/fP; fPy *= 1.0/fP; fPz *= 1.0/fP; // GeV / GeV

  // calc parameter for line at each face [mm]
  double txP = (  fLx - fSx ) / fPx;
  double txM = ( -fLx - fSx ) / fPx;
  double tyP = (  fLy - fSy ) / fPy;
  double tyM = ( -fLy - fSy ) / fPy;
  double tzP = (  fLz - fSz ) / fPz;
  double tzM = ( -fLz - fSz ) / fPz;

  // do we have an entry or exit anywhere?
  // entry from face Q = const <==> {  pr(momentum, Q) points to origin && within bounding square }
  // exit  from face Q = const <==> { -pr(momentum, Q) points to origin && within bounding square }
  double q1t = 0.0, q2t = 0.0;
  bool pointsOnX = false, pointsOnY = false, pointsOnZ = false;

  // case x = +fLx
  q1t = fSy + txP * fPy; q2t = fSz + txP * fPz;
  if( std::abs( q1t ) <= fLy && std::abs( q2t ) <= fLz ){ // within bounding square
    pointsOnX = true;
    if( fSx * fPx < 0 ){ // pointing towards origin
      fEx = fLx; fEy = q1t; fEz = q2t;
    } else if( fSx * fPx > 0 ){ // pointing away from origin
      fXx = fLx; fXy = q1t; fXz = q2t;
    } else return false; // treat tangent as no entry
  }
  // case x = -fLx
  q1t = fSy + txM * fPy; q2t = fSz + txM * fPz;
  if( std::abs( q1t ) <= fLy && std::abs( q2t ) <= fLz ){ // within bounding square
    pointsOnX = true;
    if( fSx * fPx < 0 ){ // pointing towards origin
      fEx = -fLx; fEy = q1t; fEz = q2t;
    } else if( fSx * fPx > 0 ){ // pointing away from origin
      fXx = -fLx; fXy = q1t; fXz = q2t;
    } else return false; // treat tangent as no entry
  }

  // case y = +fLy
  q1t = fSz + tyP * fPz; q2t = fSx + tyP * fPx;
  if( std::abs( q1t ) <= fLz && std::abs( q2t ) <= fLx ){ // within bounding square
    pointsOnY = true;
    if( fSy * fPy < 0 ){ // pointing towards origin
      fEx = q2t; fEy = fLy; fEz = q1t;
    } else if( fSy * fPy > 0 ){ // pointing away from origin
      fXx = q2t; fXy = fLy; fXz = q1t;
    } else return false; // treat tangent as no entry
  }
  // case y = -fLy
  q1t = fSz + tyM * fPz; q2t = fSx + tyM * fPx;
  if( std::abs( q1t ) <= fLz && std::abs( q2t ) <= fLx ){ // within bounding square
    pointsOnY = true;
    if( fSy * fPy < 0 ){ // pointing towards origin
      fEx = q2t; fEy = -fLy; fEz = q1t;
    } else if( fSy * fPy > 0 ){ // pointing away from origin
      fXx = q2t; fXy = -fLy; fXz = q1t;
    } else return false; // treat tangent as no entry
  }

  // case z = +fLz
  q1t = fSx + tzP * fPx; q2t = fSy + tzP * fPy;
  if( std::abs( q1t ) <= fLx && std::abs( q2t ) <= fLy ){ // within bounding square
    pointsOnZ = true;
    if( fSz * fPz < 0 ){ // pointing towards origin
      fEx = q1t; fEy = q2t; fEz = fLz;
    } else if( fSz * fPz > 0 ){ // pointing away from origin
      fXx = q1t; fXy = q2t; fXz = fLz;
    } else return false; // treat tangent as no entry
  }
  // case z = -fLz
  q1t = fSx + tzM * fPx; q2t = fSy + tzM * fPy;
  if( std::abs( q1t ) <= fLx && std::abs( q2t ) <= fLy ){ // within bounding square
    pointsOnZ = true;
    if( fSz * fPz < 0 ){ // pointing towards origin
      fEx = q1t; fEy = q2t; fEz = -fLz;
    } else if( fSz * fPz > 0 ){ // pointing away from origin
      fXx = q1t; fXy = q2t; fXz = -fLz;
    } else return false; // treat tangent as no entry
  }

  bool finalPoints = ( pointsOnX || pointsOnY || pointsOnZ );
  if( finalPoints ){
    entryPoint.SetXYZ( fEx, fEy, fEz );
    exitPoint.SetXYZ( fXx, fXy, fXz );
    return true;
  }
  
  // missed detector
  return false;
}
//____________________________________________________________________________
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
void VertexGenerator::ImportBoundingBox( TGeoBBox * box ) const
{
  fLx = box->GetDX() * units::cm / lunits;
  fLy = box->GetDY() * units::cm / lunits;
  fLz = box->GetDZ() * units::cm / lunits;
  fOx = (box->GetOrigin())[0] * units::cm / lunits;
  fOy = (box->GetOrigin())[1] * units::cm / lunits;
  fOz = (box->GetOrigin())[2] * units::cm / lunits;

  fLxROOT = box->GetDX();
  fLyROOT = box->GetDY();
  fLzROOT = box->GetDZ();

  fOxROOT = (box->GetOrigin())[0];
  fOyROOT = (box->GetOrigin())[1];
  fOzROOT = (box->GetOrigin())[2];
}
//____________________________________________________________________________
void VertexGenerator::SetStartingParameters( GHepRecord * event_rec ) const
{
  isUsingDk2nu = (event_rec->Particle(1) != NULL); // validation App doesn't run Decayer
  isParticleGun = (event_rec->Particle(0)->FirstMother() < -1); // hack
  isUsingRootGeom = true;

  /*
  LOG( "HNL", pDEBUG ) << "isParticleGun = " << (int) isParticleGun
		       << " with first mother " << event_rec->Particle(0)->FirstMother();
  */

  //uMult = ( isUsingDk2nu ) ? units::m / units::mm : units::cm / units::mm;
  //uMult = units::m / units::mm;
  //xMult = ( isUsingDk2nu ) ? units::cm / units::mm : 1.0;
  xMult = ( !isParticleGun ) ? units::cm / units::mm : 1.0;

  fCoMLifetime = event_rec->Probability();

  assert( event_rec->Particle(0) && "Event record has HNL" );

  TVector3 dumori(0.0, 0.0, 0.0); // tgt-hall frame origin is 0
  /*
  TVector3 detori( (fCx + fDetTranslation.at(0)) * units::m / units::cm,
		   (fCy + fDetTranslation.at(1)) * units::m / units::cm,
		   (fCz + fDetTranslation.at(2)) * units::m / units::cm ); // for rotations of the detector
  if( isParticleGun ){
    detori.SetXYZ( fDetTranslation.at(0) * units::m / units::mm,
		   fDetTranslation.at(1) * units::m / units::mm,
		   fDetTranslation.at(2) * units::m / units::mm );
  }
  */
  TVector3 detori( (fCx) * units::m / units::cm,
		   (fCy) * units::m / units::cm,
		   (fCz) * units::m / units::cm ); // for rotations of the detector
  if( isParticleGun ){
    detori.SetXYZ( 0.0, 0.0, 0.0 );
  }

  TLorentzVector * x4HNL = event_rec->Particle(0)->GetX4(); // NEAR, cm ns
  TVector3 xHNL_near = x4HNL->Vect();
  TVector3 xHNL_user = ( isParticleGun ) ? this->ApplyUserRotation( xHNL_near, detori, fDetRotation, false ) : xHNL_near; // tgt-hall --> user
  TLorentzVector * x4HNL_user = new TLorentzVector();
  if( !isParticleGun ){
    x4HNL_user->SetXYZT( xHNL_user.X() - (fCx) * units::m / units::cm, 
			 xHNL_user.Y() - (fCy) * units::m / units::cm,
			 xHNL_user.Z() - (fCz) * units::m / units::cm,
			 x4HNL->T() ); // USER, cm ns
  } else {
    /*
    x4HNL_user->SetXYZT( xHNL_user.X() - fDetTranslation.at(0) * units::m / units::mm,
			 xHNL_user.Y() - fDetTranslation.at(1) * units::m / units::mm,
			 xHNL_user.Z() - fDetTranslation.at(2) * units::m / units::mm,
			 x4HNL->T() ); // USER, mm ns
    */
    x4HNL_user->SetXYZT( xHNL_user.X() * units::m / units::mm,
			 xHNL_user.Y() * units::m / units::mm,
			 xHNL_user.Z() * units::m / units::mm,
			 x4HNL->T() ); // USER, mm ns
  }

  // set starting point for calculations. Where would the flux go?
  TLorentzVector * x4Flux = 0;
  if( !isParticleGun ){
    x4Flux = event_rec->Particle(1)->GetX4(); // USER, m
    x4Flux->SetXYZT( x4Flux->X() * units::m / units::cm,
		     x4Flux->Y() * units::m / units::cm,
		     x4Flux->Z() * units::m / units::cm, 0.0 ); // USER, cm
  } else {
    x4Flux = x4HNL_user;
  }
  TVector3 startPoint( xMult * x4Flux->X(), xMult * x4Flux->Y(), xMult * x4Flux->Z() ); // USER mm

  LOG( "HNL", pDEBUG )
    << "\nx4HNL_user = " << utils::print::X4AsString( x4HNL_user ) << " [mm]"
    << "\nstartPoint = " << utils::print::Vec3AsString( &startPoint ) << " [mm]";

  //double mtomm = units::m / units::mm;
  
  TLorentzVector * p4HNL = event_rec->Particle(0)->GetP4();
  TVector3 momentum( p4HNL->Px(), p4HNL->Py(), p4HNL->Pz() );


  fSx = startPoint.X(); fSy = startPoint.Y(); fSz = startPoint.Z();
  fSxROOT = fSx * units::mm / units::cm;
  fSyROOT = fSy * units::mm / units::cm;
  fSzROOT = fSz * units::mm / units::cm;
  fPx = momentum.X(); fPy = momentum.Y(); fPz = momentum.Z();

  delete p4HNL;
  delete x4HNL;
  delete x4HNL_user;
  if( x4Flux ) delete x4Flux;
}
//____________________________________________________________________________
bool VertexGenerator::VolumeEntryAndExitPoints( TVector3 & startPoint, TVector3 & momentum,
						TVector3 & entryPoint, TVector3 & exitPoint,
						TGeoManager * /* gm */, TGeoVolume * /* vol */ ) const
{
  const double mmtolunits = units::mm / lunits;

  // setup initial quantities
  double sx = startPoint.X(); double sy = startPoint.Y(); double sz = startPoint.Z();
  sx *= mmtolunits; sy *= mmtolunits; sz *= mmtolunits;
  double px = momentum.X(); double py = momentum.Y(); double pz = momentum.Z();
  double p2 = px*px + py*py + pz*pz; double p = std::sqrt(p2);
  px *= 1./p; py *= 1./p; pz *= 1./p;

  fSx = sx; fSy = sy; fSz = sz;
  fSxROOT = fSx * lunits / units::cm; fSyROOT = fSy * lunits / units::cm; fSzROOT = fSz * lunits / units::cm;
  fPx = px; fPy = py; fPz = pz;

  // offset the volume translation. Turn this back once we're done.
  double firstXROOT = fSxROOT - fTx * units::m / units::cm, 
    firstYROOT = fSyROOT - fTy * units::m / units::cm, 
    firstZROOT = fSzROOT - fTz * units::m / units::cm;
  // also offset the detector offset and reapply once we're done
  firstXROOT -= fDetTranslation.at(0) * units::m / units::cm;
  firstYROOT -= fDetTranslation.at(1) * units::m / units::cm;
  firstZROOT -= fDetTranslation.at(2) * units::m / units::cm;

  TVector3 dumori(0.0, 0.0, 0.0);

  // we will demand that this is inside the required geometry.
  gGeoManager->SetCurrentPoint( firstXROOT, firstYROOT, firstZROOT );
  gGeoManager->SetCurrentDirection( px, py, pz );
  
  LOG( "HNL", pINFO )
    << "Starting to figure out entrances:"
    << "\nStarting point is ( " << fSx << ", " << fSy << ", " << fSz << " ) [ " << lunitString.c_str() << " ]"
    << "\nIn our volume this becomes ( " 
    << fSx - (fTx + fDetTranslation.at(0)) * units::m / lunits << ", " 
    << fSy - (fTy + fDetTranslation.at(1)) * units::m / lunits << ", " 
    << fSz - (fTz + fDetTranslation.at(2)) * units::m / lunits << " ) [ " 
    << lunitString.c_str() << " ]"
    << "\nStarting dirn  is ( " << fPx << ", " << fPy << ", " << fPz << " ) ";

  std::string pathString = this->CheckGeomPoint( firstXROOT, firstYROOT, firstZROOT );
  LOG( "HNL", pDEBUG ) << "Here is the pathString: " << pathString;

  LOG( "HNL", pDEBUG ) << "Starting to search for intersections...";

  if( isParticleGun ) {
    // We start outside the detector. First, write out the starting point in top_volume coordinates
    //RandomGen * rnd = RandomGen::Instance();
    LOG( "HNL" , pDEBUG )
      << "isParticleGun!!"
      << "\nstartPoint = ( " << (gGeoManager->GetCurrentPoint())[0]
      << ", " << (gGeoManager->GetCurrentPoint())[1] << ", " << (gGeoManager->GetCurrentPoint())[2]
      << " ) [top_volume, cm]"
      << "\nDirection = ( " << (gGeoManager->GetCurrentDirection())[0]
      << ", " << (gGeoManager->GetCurrentDirection())[1] 
      << ", " << (gGeoManager->GetCurrentDirection())[2] << " ) [GeV/GeV]";

    assert( (gGeoManager->GetCurrentDirection())[2] != 0.0 && "HNL propagates along USER z" );

    // Check along the detector for entering the detector geometry

    //double xBack = -fLxROOT/2.0; double xFront = fLxROOT/2.0;
    //double yBack = -fLyROOT/2.0; double yFront = fLyROOT/2.0;
    double zBack = -fLzROOT/2.0; double zFront = fLzROOT/2.0; // USER cm

    //double zTest = (rnd->RndGen()).Uniform( zBack, zFront ); // USER cm
    double zTest = zBack + 0.01 * (zFront - zBack);

    double dz = zTest - (gGeoManager->GetCurrentPoint())[2]; // cm
    double dxdz = (gGeoManager->GetCurrentDirection())[0]/(gGeoManager->GetCurrentDirection())[2];
    double dydz = (gGeoManager->GetCurrentDirection())[1]/(gGeoManager->GetCurrentDirection())[2];
    double dx = dxdz * dz; // cm
    double dy = dydz * dz; // cm
    double xTest = (gGeoManager->GetCurrentPoint())[0] + dx;
    double yTest = (gGeoManager->GetCurrentPoint())[1] + dy;

    // now check to see if this is within the geometry
    pathString = this->CheckGeomPoint( xTest, yTest, zTest );
    while( pathString.find( fTopVolume.c_str() ) == string::npos &&
	   zTest < zFront ){
      dz = zTest - (gGeoManager->GetCurrentPoint())[2]; // cm
      dxdz = (gGeoManager->GetCurrentDirection())[0]/(gGeoManager->GetCurrentDirection())[2];
      dydz = (gGeoManager->GetCurrentDirection())[1]/(gGeoManager->GetCurrentDirection())[2];
      dx = dxdz * dz; // cm
      dy = dydz * dz; // cm
      xTest = (gGeoManager->GetCurrentPoint())[0] + dx;
      yTest = (gGeoManager->GetCurrentPoint())[1] + dy;

      pathString = this->CheckGeomPoint( xTest, yTest, zTest );

      zTest += 0.01 * (zFront - zBack);
    }
    if( pathString.find( fTopVolume.c_str() ) == string::npos ){
      LOG( "HNL", pDEBUG )
	<< "This trajectory does NOT intersect the detector. Bailing...";
      return false;
    } else { // we are ok! Set the new starting point and continue to find entry and exit points.
      LOG( "HNL", pDEBUG )
	<< "This trajectory DOES intersect the detector. Good!";
      gGeoManager->SetCurrentPoint( xTest, yTest, zTest );
      firstXROOT = xTest; firstYROOT = yTest; firstZROOT = zTest;
      LOG( "HNL", pDEBUG )
	<< "New starting point is at ( " << xTest << ", " << yTest << ", " << zTest << " ) [top_volume, cm]";
    }
  } // if( isParticleGun )

  // we are inside the top volume. We want to exit it twice, once going backwards and once forwards.
  // The track enters in the former point and exits in the latter.

  // -- Entry point
  gGeoManager->SetCurrentPoint( firstXROOT, firstYROOT, firstZROOT );
  gGeoManager->SetCurrentDirection( -px, -py, -pz );

  pathString = this->CheckGeomPoint( (gGeoManager->GetCurrentPoint())[0],
				     (gGeoManager->GetCurrentPoint())[1],
				     (gGeoManager->GetCurrentPoint())[2] );

  entryPoint.SetXYZ( firstXROOT, firstYROOT, firstZROOT );
  const double smallStep = 0.1; // cm --> have mm precision
  //std::max( 0.1, 1.0e-3 * std::min( fLxROOT, std::min( fLyROOT, fLzROOT ) ) ); // cm --> have mm precision
  while( pathString.find( fTopVolume.c_str() ) != string::npos ){
    double newX = entryPoint.X() + (gGeoManager->GetCurrentDirection())[0] * smallStep;
    double newY = entryPoint.Y() + (gGeoManager->GetCurrentDirection())[1] * smallStep;
    double newZ = entryPoint.Z() + (gGeoManager->GetCurrentDirection())[2] * smallStep;

    pathString = this->CheckGeomPoint( newX, newY, newZ );
    if( pathString.find( fTopVolume.c_str() ) != string::npos )
      entryPoint.SetXYZ( newX, newY, newZ );
  } // exit out
  
  // Let's save this point
  fEx = ( gGeoManager->GetCurrentPoint() )[0] * genie::units::cm / lunits;
  fEy = ( gGeoManager->GetCurrentPoint() )[1] * genie::units::cm / lunits;
  fEz = ( gGeoManager->GetCurrentPoint() )[2] * genie::units::cm / lunits;

  fExROOT = ( gGeoManager->GetCurrentPoint() )[0];
  fEyROOT = ( gGeoManager->GetCurrentPoint() )[1];
  fEzROOT = ( gGeoManager->GetCurrentPoint() )[2];

  fEx += (fTx + fDetTranslation.at(0)) * units::m / lunits; 
  fEy += (fTy + fDetTranslation.at(1)) * units::m / lunits; 
  fEz += (fTz + fDetTranslation.at(2)) * units::m / lunits;
  fExROOT += (fTx + fDetTranslation.at(0)) * units::m / units::cm; 
  fEyROOT += (fTy + fDetTranslation.at(1)) * units::m / units::cm; 
  fEzROOT += (fTz + fDetTranslation.at(2)) * units::m / units::cm;

  entryPoint.SetXYZ( fEx, fEy, fEz ); // ensure correct units

  TVector3 entryPoint_user( fExROOT * units::cm / units::m,
			    fEyROOT * units::cm / units::m,
			    fEzROOT * units::cm / units::m ); // USER, m

  TVector3 entryPoint_near = this->ApplyUserRotation( entryPoint_user, dumori, fDetRotation, true );
  entryPoint_near.SetXYZ( entryPoint_near.X() + (fCx + fDetTranslation.at(0)),
			  entryPoint_near.Y() + (fCy + fDetTranslation.at(1)),
			  entryPoint_near.Z() + (fCz + fDetTranslation.at(2)) );

  LOG( "HNL", pDEBUG )
    << "\nEntry point found at ( " << fEx << ", " << fEy << ", " << fEz << " ) [" << lunitString.c_str() << "]"
    << "\nIn ROOT, entry at    ( " << fExROOT << ", " << fEyROOT << ", " << fEzROOT << " ) [cm]"; 

  // now go back to the beginning, reset the direction to go forwards, and try again.
  gGeoManager->SetCurrentPoint( firstXROOT, firstYROOT, firstZROOT );
  gGeoManager->SetCurrentDirection( px, py, pz );

  pathString = this->CheckGeomPoint( (gGeoManager->GetCurrentPoint())[0],
				     (gGeoManager->GetCurrentPoint())[1],
				     (gGeoManager->GetCurrentPoint())[2] );

  // -- Exit point

  exitPoint.SetXYZ( firstXROOT, firstYROOT, firstZROOT );
  pathString = this->CheckGeomPoint( (gGeoManager->GetCurrentPoint())[0],
				     (gGeoManager->GetCurrentPoint())[1],
				     (gGeoManager->GetCurrentPoint())[2] );
  while( pathString.find( fTopVolume.c_str() ) != string::npos ){
    double newX = exitPoint.X() + (gGeoManager->GetCurrentDirection())[0] * smallStep;
    double newY = exitPoint.Y() + (gGeoManager->GetCurrentDirection())[1] * smallStep;
    double newZ = exitPoint.Z() + (gGeoManager->GetCurrentDirection())[2] * smallStep;

    pathString = this->CheckGeomPoint( newX, newY, newZ );
    if( pathString.find( fTopVolume.c_str() ) != string::npos )
      exitPoint.SetXYZ( newX, newY, newZ );
  } // exit out

    // Let's save this point
  fXx = ( gGeoManager->GetCurrentPoint() )[0] * genie::units::cm / lunits;
  fXy = ( gGeoManager->GetCurrentPoint() )[1] * genie::units::cm / lunits;
  fXz = ( gGeoManager->GetCurrentPoint() )[2] * genie::units::cm / lunits;

  fXxROOT = ( gGeoManager->GetCurrentPoint() )[0];
  fXyROOT = ( gGeoManager->GetCurrentPoint() )[1];
  fXzROOT = ( gGeoManager->GetCurrentPoint() )[2];

  fXx += (fTx + fDetTranslation.at(0)) * units::m / lunits; 
  fXy += (fTy + fDetTranslation.at(1)) * units::m / lunits; 
  fXz += (fTz + fDetTranslation.at(2)) * units::m / lunits;
  fXxROOT += (fTx + fDetTranslation.at(0)) * units::m / units::cm; 
  fXyROOT += (fTy + fDetTranslation.at(1)) * units::m / units::cm; 
  fXzROOT += (fTz + fDetTranslation.at(2)) * units::m / units::cm;

  exitPoint.SetXYZ( fXx, fXy, fXz ); // ensure correct units

  TVector3 exitPoint_user( fXxROOT * units::cm / units::m,
			   fXyROOT * units::cm / units::m,
			   fXzROOT * units::cm / units::m ); // USER, m

  TVector3 exitPoint_near = this->ApplyUserRotation( exitPoint_user, dumori, fDetRotation, true );
  exitPoint_near.SetXYZ( exitPoint_near.X() + (fCx + fDetTranslation.at(0)),
			 exitPoint_near.Y() + (fCy + fDetTranslation.at(1)),
			 exitPoint_near.Z() + (fCz + fDetTranslation.at(2)) );

  LOG( "HNL", pDEBUG )
    << "\nExit point found at ( " << fXx << ", " << fXy << ", " << fXz << " ) [" << lunitString.c_str() << "]"
    << "\nIn ROOT, exit at    ( " << fXxROOT << ", " << fXyROOT << ", " << fXzROOT << " ) [cm]"; 

  return true;
  
}
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__
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

  LOG( "HNL", pDEBUG )
    << "Loading geometry parameters from file. . .";

  this->GetParamVect( "Near2User_T", fB2UTranslation );
  this->GetParamVect( "Near2User_R", fDetRotation );
  this->GetParamVect( "Near2Beam_R", fB2URotation );
  this->GetParamVect( "DetCentre_User", fDetTranslation );
  fCx = fB2UTranslation.at(0); fCy = fB2UTranslation.at(1); fCz = fB2UTranslation.at(2);
  fUx = fDetTranslation.at(0); fUy = fDetTranslation.at(1); fUz = fDetTranslation.at(2);
  fAx1 = fB2URotation.at(0); fAz = fB2URotation.at(1); fAx2 = fB2URotation.at(2);
  fBx1 = fDetRotation.at(0); fBz = fDetRotation.at(1); fBx2 = fDetRotation.at(2);

  fTx = 0.0; fTy = 0.0; fTz = 0.0;

  fIsConfigLoaded = true;
}
//____________________________________________________________________________
void VertexGenerator::GetInterestingPoints( TVector3 & entryPoint, TVector3 & exitPoint, TVector3 & decayPoint ) const
{
  entryPoint.SetXYZ( fEx, fEy, fEz );
  exitPoint.SetXYZ( fXx, fXy, fXz );
  decayPoint.SetXYZ( fDx, fDy, fDz );
}
//____________________________________________________________________________
TVector3 VertexGenerator::ApplyUserRotation( TVector3 vec, bool doBackwards ) const
{
  double vx = vec.X(), vy = vec.Y(), vz = vec.Z();

  double Ax2 = ( doBackwards ) ? -fAx2 : fAx2;
  double Az  = ( doBackwards ) ? -fAz  : fAz;
  double Ax1 = ( doBackwards ) ? -fAx1 : fAx1;

  // Ax2 first
  double x = vx, y = vy, z = vz;
  vy = y * std::cos( Ax2 ) - z * std::sin( Ax2 );
  vz = y * std::sin( Ax2 ) + z * std::cos( Ax2 );
  y = vy; z = vz;
  // then Az
  vx = x * std::cos( Az )  - y * std::sin( Az );
  vy = x * std::sin( Az )  + y * std::cos( Az );
  x = vx; y = vy;
  // Ax1 last
  vy = y * std::cos( Ax1 ) - z * std::sin( Ax1 );
  vz = y * std::sin( Ax1 ) + z * std::cos( Ax1 );

  TVector3 nvec( vx, vy, vz );
  return nvec;
}
//____________________________________________________________________________
TVector3 VertexGenerator::ApplyUserRotation( TVector3 vec, TVector3 oriVec, std::vector<double> rotVec, bool doBackwards ) const
{
  double vx = vec.X(), vy = vec.Y(), vz = vec.Z();
  double ox = oriVec.X(), oy = oriVec.Y(), oz = oriVec.Z();
  
  vx -= ox; vy -= oy; vz -= oz; // make this rotation about detector origin

  assert( rotVec.size() == 3 && "3 Euler angles" ); // want 3 Euler angles, otherwise this is unphysical.
  double Ax2 = ( doBackwards ) ? -rotVec.at(2) : rotVec.at(2);
  double Az  = ( doBackwards ) ? -rotVec.at(1) : rotVec.at(1);
  double Ax1 = ( doBackwards ) ? -rotVec.at(0) : rotVec.at(0);

  // Ax2 first
  double x = vx, y = vy, z = vz;
  vy = y * std::cos( Ax2 ) - z * std::sin( Ax2 );
  vz = y * std::sin( Ax2 ) + z * std::cos( Ax2 );
  y = vy; z = vz;
  // then Az
  vx = x * std::cos( Az )  - y * std::sin( Az );
  vy = x * std::sin( Az )  + y * std::cos( Az );
  x = vx; y = vy;
  // Ax1 last
  vy = y * std::cos( Ax1 ) - z * std::sin( Ax1 );
  vz = y * std::sin( Ax1 ) + z * std::cos( Ax1 );

  // back to beam frame
  vx += ox; vy += oy; vz += oz;
  TVector3 nvec( vx, vy, vz );
  return nvec;
}
//____________________________________________________________________________
void VertexGenerator::SetGeomFile( string geomfile, string topVolume ) const
{
  fGeomFile = geomfile;
  fTopVolume = topVolume;
}
//____________________________________________________________________________
std::string VertexGenerator::CheckGeomPoint( Double_t x, Double_t y, Double_t z ) const
{
  Double_t point[3];
  Double_t local[3];
  point[0] = x;
  point[1] = y;
  point[2] = z;
  TGeoVolume *vol = gGeoManager->GetVolume(fTopVolume.c_str()); //gGeoManager->GetTopVolume();
  TGeoNode *node = gGeoManager->FindNode(point[0], point[1], point[2]);
  gGeoManager->MasterToLocal(point, local);
  return gGeoManager->GetPath();
}
//____________________________________________________________________________
TGeoMatrix * VertexGenerator::FindFullTransformation( TGeoVolume * top_vol, TGeoVolume * tar_vol ) const
{
  // Recurses over the ROOT file geometry structure to find the target volume tar_vol.
  // Returns the full transformation matrix of the daughter volume as a composition of matrices

  std::list<TGeoNode *> nodes; // to parse hierarchy (i.e. daughters)
  std::list<std::string> paths; // store path of nodes checked
  std::list<TGeoMatrix *> mats; // compositions of matrices here!

  assert( top_vol && tar_vol && "Top and target volumes both accessible" );

  std::string targetPath( tar_vol->GetName() );

  // Start by grabbing the daughter structure of the top volume and parse until found
  TGeoNode * top_node = top_vol->GetNode(0); nodes.emplace_back( top_node );
  std::string top_path( top_node->GetName() ); paths.emplace_back( top_path );
  TGeoMatrix * top_mat = top_node->GetMatrix(); mats.emplace_back( top_mat );

  std::string test = paths.front();
  // strip all slashes from test
  while( test.find("/") != string::npos ){
    int idx = test.find("/");
    test = test.substr(idx+1);
  }
  // and strip tailing underscore
  int ididx = test.rfind("_");
  test = test.substr( 0, ididx );

  LOG( "HNL", pNOTICE )
    << "Looking for this targetPath: " << targetPath;

  // could be we hit the top volume, in which case we skip the loop
  bool foundPath = (strcmp( test.c_str(), targetPath.c_str() ) == 0);

  //while( test.find( targetPath.c_str() ) == string::npos ){ // still looking for the path.
  while( strcmp( test.c_str(), targetPath.c_str() ) != 0 && !foundPath ){ // still looking
    TGeoNode * node = nodes.front();
    std::string path = paths.front();
    TGeoMatrix * mat = mats.front();

    assert( node  && "Node is not null" );
    assert( mat && "Matrix is not null" );

    int nDaughters = node->GetNdaughters();
    LOG( "HNL", pDEBUG ) << "Node with name " << path << " has " << nDaughters << " daughters...";
    for( int iDaughter = 0; iDaughter < nDaughters; iDaughter++ ){
      TGeoNode * dNode = node->GetDaughter(iDaughter);
      assert( dNode && "Daughter node not null" );
      std::string dPath( path );
      dPath.append( "/" ); dPath.append( dNode->GetName() );
      TGeoMatrix * nodeMat = dNode->GetMatrix();

      LOG( "HNL", pDEBUG ) << "Got node, path, and matrix for daughter node "
			   << iDaughter << " / " << nDaughters-1 << "...";

      // construct the full updated matrix from multiplying dMat on the left of mat
      const Double_t * nodeRot = nodeMat->GetRotationMatrix();
      const Double_t * nodeTra = nodeMat->GetTranslation();
      
      const Double_t * baseRot = mat->GetRotationMatrix();
      const Double_t * baseTra = mat->GetTranslation();

      const Double_t compTra[3] = { baseTra[0] + nodeTra[0], 
				    baseTra[1] + nodeTra[1],
				    baseTra[2] + nodeTra[2] }; // this was easy.
      const Double_t compRot[9] = { nodeRot[0] * baseRot[0] + nodeRot[1] * baseRot[3] + nodeRot[2] * baseRot[6],
				    nodeRot[0] * baseRot[1] + nodeRot[1] * baseRot[4] + nodeRot[2] * baseRot[7],
				    nodeRot[0] * baseRot[2] + nodeRot[1] * baseRot[5] + nodeRot[2] * baseRot[8],
				    nodeRot[3] * baseRot[0] + nodeRot[4] * baseRot[3] + nodeRot[5] * baseRot[6],
				    nodeRot[3] * baseRot[1] + nodeRot[4] * baseRot[4] + nodeRot[5] * baseRot[7],
				    nodeRot[3] * baseRot[2] + nodeRot[4] * baseRot[5] + nodeRot[5] * baseRot[8],
				    nodeRot[6] * baseRot[0] + nodeRot[7] * baseRot[3] + nodeRot[8] * baseRot[6],
				    nodeRot[6] * baseRot[1] + nodeRot[7] * baseRot[4] + nodeRot[8] * baseRot[7],
				    nodeRot[6] * baseRot[2] + nodeRot[7] * baseRot[5] + nodeRot[8] * baseRot[8] }; // less easy but ok.

      // construct a TGeoMatrix * from these quantities...
      TGeoHMatrix * hmat = new TGeoHMatrix( dPath.c_str() );
      hmat->SetTranslation( compTra );
      hmat->SetRotation( compRot );
      TGeoMatrix * dMat = dynamic_cast< TGeoMatrix * >( hmat );

      /*
      LOG( "HNL", pDEBUG )
	<< "\nNode with name " << targetPath << " not yet found."
	<< "\nParsing node with name " << dPath << "..."
	<< "\n\nThis node had the following translations: ( " 
	<< nodeTra[0] << ", " << nodeTra[1] << ", " << nodeTra[2] << " )"
	<< " composed onto ( " << baseTra[0] << ", " << baseTra[1] << ", " << baseTra[2] << " ),"
	<< "\ngiving a final translation ( " << compTra[0] << ", " << compTra[1] << ", " 
	<< compTra[2] << " )."
	<< "\n\nThis node had the following rotation matrix: ( ( " 
	<< nodeRot[0] << ", " << nodeRot[1] << ", " << nodeRot[2] << " ), ( "
	<< nodeRot[3] << ", " << nodeRot[4] << ", " << nodeRot[5] << " ), ( "
	<< nodeRot[6] << ", " << nodeRot[7] << ", " << nodeRot[8] << " ) ),"
	<< "\ncomposed onto ( ( "
	<< baseRot[0] << ", " << baseRot[1] << ", " << baseRot[2] << " ), ( "
	<< baseRot[3] << ", " << baseRot[4] << ", " << baseRot[5] << " ), ( "
	<< baseRot[6] << ", " << baseRot[7] << ", " << baseRot[8] << " ) ),"
	<< "\ngiving a final rotation ( ( "
	<< compRot[0] << ", " << compRot[1] << ", " << compRot[2] << " ), ( "
	<< compRot[3] << ", " << compRot[4] << ", " << compRot[5] << " ), ( "
	<< compRot[6] << ", " << compRot[7] << ", " << compRot[8] << " ) ).";
      */

      // add to list TAIL and strike away list HEAD
      nodes.emplace_back( dNode );
      paths.emplace_back( dPath );
      mats.emplace_back( dMat );

      // break if we found the target path to ensure the TAIL always points to desired node
      //if( dPath.find( targetPath.c_str() ) ) break;
      while( dPath.find("/") != string::npos ){
	int idx = dPath.find("/");
	dPath = dPath.substr(idx+1);
      }
      ididx = dPath.rfind("_");
      dPath = dPath.substr( 0, ididx );
      if( strcmp( dPath.c_str(), targetPath.c_str() ) == 0 ){ foundPath = true; break; }
    } // loop over daughters

    if( !foundPath ){ // prevent popping out of last element!
      nodes.pop_front();
      paths.pop_front();
      mats.pop_front();
      
      test = paths.front();
      while( test.find("/") != string::npos ){
	int idx = test.find("/");
	test = test.substr(idx+1);
      }
      ididx = test.rfind("_");
      test = test.substr( 0, ididx );
    }
  } // while path not found

  std::string final_path = paths.back();
  while( final_path.find("/") != string::npos ){
    int idx = final_path.find("/");
    final_path = final_path.substr(idx+1);
  }
  ididx = final_path.rfind("_");
  final_path = final_path.substr( 0, ididx );
  assert( strcmp( final_path.c_str(), targetPath.c_str() ) == 0 && foundPath &&
	  "Found the target volume's path in the ROOT geometry hierarchy" );
  // found the path! The matrix is at the end.
  TGeoMatrix * final_mat = mats.back();

  const Double_t * final_tra = final_mat->GetTranslation();
  const Double_t * final_rot = final_mat->GetRotationMatrix();

  LOG( "HNL", pINFO )
    << "Found the target volume! Here is its path and full matrix:"
    << "\nPath: " << paths.back()
    << "\nTranslations: ( " << final_tra[0] << ", " << final_tra[1] << ", " << final_tra[2]
    << " ) [cm]"
    << "\nRotation matrix: ( ( " 
    << final_rot[0] << ", " << final_rot[1] << ", " << final_rot[2] << " ), ( "
    << final_rot[3] << ", " << final_rot[4] << ", " << final_rot[5] << " ), ( "
    << final_rot[6] << ", " << final_rot[7] << ", " << final_rot[8] << " ) )";

  return final_mat;
}

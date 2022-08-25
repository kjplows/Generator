//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include "Physics/NeutralHeavyLepton/NHLDecayVolume.h"

using namespace genie;
using namespace genie::NHL;
using namespace genie::units;

//____________________________________________________________________________
NHLDecayVolume::NHLDecayVolume() :
  EventRecordVisitorI("genie::NHL::NHLDecayVolume")
{

}
//____________________________________________________________________________
NHLDecayVolume::NHLDecayVolume(string config) :
  EventRecordVisitorI("genie::NHL::NHLDecayVolume", config)
{

}
//____________________________________________________________________________
NHLDecayVolume::~NHLDecayVolume()
{

}
//____________________________________________________________________________
void NHLDecayVolume::ProcessEventRecord(GHepRecord * event_rec) const
{
  /*!
   *  Uses ROOT's TGeoManager to find out where the intersections with the detector volume live
   *  Label them as entry and exit point. Then use them to determine:
   *  1) A decay vertex within the detector
   *  2) A time-of-decay (== delay of NHL to reach the decay vertex wrt a massless SM v)
   *  3) Geom weight: Survival to detector * decay within detector.
   */

  LOG( "NHL", pDEBUG )
    << "Entering ProcessEventRecord...";

  int trajIdx = 0, trajMax = 20;
  double weight = 1.0; // pure geom weight

  TVector3 startPoint, momentum, entryPoint, exitPoint;
  startPoint.SetXYZ( fSx, fSy, fSz );
  momentum.SetXYZ( fPx, fPy, fPz );
  
  if( !fGeoManager )
    fGeoManager = TGeoManager::Import(fGeomFile.c_str());
  bool didIntersectDet = this->VolumeEntryAndExitPoints( startPoint, momentum, entryPoint, exitPoint, fGeoManager, fGeoVolume );

  LOG( "NHL", pDEBUG )
    << "\n startPoint = " << utils::print::Vec3AsString( &startPoint )
    << "\n momentum = " << utils::print::Vec3AsString( &momentum );

  if( isUsingDk2nu ) assert( didIntersectDet ); // forced to hit detector somewhere!
  else {

    const Algorithm * algNHLGen = AlgFactory::Instance()->GetAlgorithm("genie::NHLPrimaryVtxGenerator", "Default");
    
    const NHLPrimaryVtxGenerator * nhlgen = dynamic_cast< const NHLPrimaryVtxGenerator * >( algNHLGen );
    
    std::vector< double > * newProdVtx = new std::vector< double >();
    newProdVtx->emplace_back( startPoint.X() );
    newProdVtx->emplace_back( startPoint.Y() );
    newProdVtx->emplace_back( startPoint.Z() );

    while( !didIntersectDet && trajIdx < trajMax ){
      // sample prod vtx and momentum... again
      LOG( "NHL", pDEBUG )
	<< "Sampling another trajectory (index = " << trajIdx << ")";
      newProdVtx  = nhlgen->GenerateDecayPosition( event_rec );
      
      startPoint.SetXYZ( newProdVtx->at(0), newProdVtx->at(1), newProdVtx->at(2) );
      LOG( "NHL", pDEBUG )
	<< "Set start point for this trajectory = ( " << startPoint.X() << ", " << startPoint.Y() << ", " << startPoint.Z() << " ) [cm]";
      
      trajIdx++;
      didIntersectDet = this->VolumeEntryAndExitPoints( startPoint, momentum, entryPoint, exitPoint, fGeoManager, fGeoVolume );

      newProdVtx->clear();
    }
    LOG("NHL", pNOTICE) << "Called NHLDecayVolume::VolumeEntryAndExitPoints " << trajIdx + 1 << " times";

  }
  if( trajIdx == trajMax && !didIntersectDet ){ // bail
    LOG( "NHL", pERROR )
      << "Unable to make a single good trajectory that intersects the detector after " << trajIdx << " tries! Bailing...";
    TLorentzVector v4dummy( -999.9, -999.9, -999.9, -999.9 );
    event_rec->SetVertex( v4dummy );
    return;
  }

  LOG( "NHL", pDEBUG )
    << "Intersected detector";
  this->EnforceUnits( "mm", "rad", "ns" );

  // move fCoMLifetime to ns from GeV^{-1}
  fCoMLifetime *= 1.0 / ( units::ns * units::GeV );

  double maxDx = exitPoint.X() - entryPoint.X();
  double maxDy = exitPoint.Y() - entryPoint.Y();
  double maxDz = exitPoint.Z() - entryPoint.Z();

  double maxLength = std::sqrt( std::pow( maxDx , 2.0 ) +
				std::pow( maxDy , 2.0 ) +
				std::pow( maxDz , 2.0 ) );
  
  TLorentzVector * p4NHL = event_rec->Particle(0)->GetP4();
  double betaMag = p4NHL->P() / p4NHL->E();
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

  // update the weight
  event_rec->SetWeight( event_rec->Weight() * weight );

  TVector3 decayPoint = this->GetDecayPoint( elapsed_length, entryPoint, momentum );

  // write out vtx in [m, ns]
  TLorentzVector x4( decayPoint.X() * units::mm / units::m,
		     decayPoint.Y() * units::mm / units::m,
		     decayPoint.Z() * units::mm / units::m,
		     event_rec->Vertex()->T() );

  event_rec->SetVertex(x4);
  
}
//____________________________________________________________________________
void NHLDecayVolume::EnforceUnits( std::string length_units, std::string angle_units, std::string time_units ) const{
  
  LOG( "NHL", pDEBUG )
    << "Switching units to " << length_units.c_str() << " , " << angle_units.c_str() << " , " << time_units.c_str();

  double old_lunits = lunits;
  double old_aunits = aunits;
  double old_tunits = tunits;

  lunits = utils::units::UnitFromString( length_units ); lunitString = length_units;
  aunits = utils::units::UnitFromString( angle_units );
  tunits = utils::units::UnitFromString( time_units ); tunitString = time_units;

  // convert to new units
  fSx /= lunits/old_lunits; fSy /= lunits/old_lunits; fSz /= lunits/old_lunits;
  fPx /= lunits/old_lunits; fPy /= lunits/old_lunits; fPz /= lunits/old_lunits;
  fEx /= lunits/old_lunits; fEy /= lunits/old_lunits; fEz /= lunits/old_lunits;
  fXx /= lunits/old_lunits; fXy /= lunits/old_lunits; fXz /= lunits/old_lunits;
  fLx /= lunits/old_lunits; fLy /= lunits/old_lunits; fLz /= lunits/old_lunits;

  fDx /= lunits/old_lunits; fDy /= lunits/old_lunits; fDz /= lunits/old_lunits;
  fOx /= lunits/old_lunits; fOy /= lunits/old_lunits; fOz /= lunits/old_lunits;

  kNewSpeedOfLight /= (lunits / old_lunits) / (tunits / old_tunits);

  LOG( "NHL", pDEBUG )
    << "kNewSpeedOfLight = " << kNewSpeedOfLight << " [" << lunitString.c_str() << "/"
    << tunitString.c_str() << "]";
}
//____________________________________________________________________________
double NHLDecayVolume::CalcTravelLength( double betaMag, double CoMLifetime, double maxLength ) const
{
  // decay probability P0(t) = 1 - exp( -t/tau ) where:
  // t   = time-of-flight (in rest frame)
  // tau = CoMLifetime

  assert( betaMag > 0.0 && betaMag < 1.0 ); // massive moving particle
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

  LOG( "NHL", pDEBUG )
    << "betaMag, maxLength, CoMLifetime = " << betaMag << ", " << maxLength << ", " << CoMLifetime
    << "\nbetaMag = " << betaMag << " ==> gamma = " << gamma
    << "\n==> maxLength [" << tunitString.c_str()
    << "] = " << maxRestTime << " (rest frame) = " << maxLabTime << " (lab frame)"
    << "\nranthrow = " << ranthrow << ", PExit = " << PExit
    << "\n==> S0 = " << S0 << " ==> rest_time [" << lunitString.c_str() << "] = " << rest_time
    << " ==> elapsed_time [" << tunitString.c_str()
    << "] = " << elapsed_time << " ==> elapsed_length [" << lunitString.c_str()
    << "] = " << elapsed_length;

  return elapsed_length;
}
//____________________________________________________________________________
TVector3 NHLDecayVolume::GetDecayPoint( double travelLength, TVector3 & entryPoint, TVector3 & momentum ) const
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

  LOG( "NHL", pDEBUG )
    << "\ndecayPoint = ( " << dx << ", " << dy << ", " << dz << " ) ["
    << lunitString.c_str() << "]"
    << "\ndecayPoint(ROOT) = ( " << fDxROOT << ", " << fDyROOT << ", " << fDzROOT << " ) [cm]";

  TVector3 decayPoint( dx, dy, dz );
  return decayPoint;
}
//____________________________________________________________________________
double NHLDecayVolume::GetMaxLength( TVector3 & entryPoint, TVector3 & exitPoint ) const
{
  double ex = entryPoint.X(); double ey = entryPoint.Y(); double ez = entryPoint.Z();
  double xx = exitPoint.X(); double xy = exitPoint.Y(); double xz = exitPoint.Z();

  return std::sqrt( (ex-xx)*(ex-xx) + (ey-xy)*(ey-xy) + (ez-xz)*(ez-xz) );
}
//____________________________________________________________________________
void NHLDecayVolume::MakeSDV() const
{
  fOx = 0.0; fOy = 0.0; fOz = 0.0;
  fLx = 1.0; fLy = 1.0; fLz = 1.0; // m

  lunits = utils::units::UnitFromString( "m" );
  aunits = utils::units::UnitFromString( "rad" );
  tunits = utils::units::UnitFromString( "ns" );
  
  kNewSpeedOfLight = genie::units::kSpeedOfLight 
  * (genie::units::m / lunits)
  / (genie::units::s / tunits);

  LOG("NHL", pDEBUG)
    << "Setting simple decay volume with unit-m side."
    << "\nSetting units to \"mm\", \"rad\", \"ns\"";

  EnforceUnits("mm","rad","ns");
}
//____________________________________________________________________________
// if entry and exit points, populate TVector3's with their coords. If not, return false
bool NHLDecayVolume::SDVEntryAndExitPoints( TVector3 & startPoint, TVector3 momentum,
					    TVector3 & entryPoint, TVector3 & exitPoint ) const
{
  assert( fOx == 0.0 && fOy == 0.0 && fOz == 0.0 && 
	  fLx == 1000.0 && fLy == 1000.0 && fLz == 1000.0 ); // SDV, mm
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
void NHLDecayVolume::ImportBoundingBox( TGeoBBox * box ) const
{
  fLx = 2.0 * box->GetDX() * units::cm / lunits;
  fLy = 2.0 * box->GetDY() * units::cm / lunits;
  fLz = 2.0 * box->GetDZ() * units::cm / lunits;
  fOx = (box->GetOrigin())[0] * units::cm / lunits;
  fOy = (box->GetOrigin())[1] * units::cm / lunits;
  fOz = (box->GetOrigin())[2] * units::cm / lunits;

  fLxROOT = 2.0 * box->GetDX();
  fLyROOT = 2.0 * box->GetDY();
  fLzROOT = 2.0 * box->GetDZ();

  fOxROOT = (box->GetOrigin())[0];
  fOyROOT = (box->GetOrigin())[1];
  fOzROOT = (box->GetOrigin())[2];

  LOG( "NHL", pDEBUG )
    << "\nImported bounding box with origin at ( " << fOx << ", " << fOy << ", " << fOz << " ) and sides " << fLx << " x " << fLy << " x " << fLz << " [units: " << lunitString.c_str() << "]"
    << "\nIn ROOT units this is origin at ( " << fOxROOT << ", " << fOyROOT << ", " << fOzROOT << " ) and sides " << fLxROOT << " x " << fLyROOT << " x " << fLzROOT << " [cm]";
}
//____________________________________________________________________________
void NHLDecayVolume::SetStartingParameters( GHepRecord * event_rec, double NHLCoMTau, bool usingDk2nu, bool usingRootGeom, string geomfile ) const
{
  isUsingDk2nu = usingDk2nu;
  uMult = ( isUsingDk2nu ) ? units::m / units::mm : units::cm / units::mm;
  xMult = ( isUsingDk2nu ) ? units::cm / units::mm : 1.0;

  isUsingRootGeom = usingRootGeom;

  fCoMLifetime = NHLCoMTau;

  assert( event_rec->Particle(0) );

  TLorentzVector * x4NHL = event_rec->Particle(0)->GetX4();
  TVector3 startPoint( xMult * x4NHL->X(), xMult * x4NHL->Y(), xMult * x4NHL->Z() ); // mm
  if( fUseBeamMomentum ){
    double mtomm = units::m / units::mm;
    // passive transformation. First return to tgt-hall frame, then to detector frame
    TVector3 beamOrigin( 0.0, 0.0, 0.0 ), detOrigin( fUx * mtomm, fUy * mtomm, fUz * mtomm );
    TVector3 startUnrotated = startPoint;
    startPoint = this->ApplyUserRotation( startPoint, beamOrigin, fB2URotation, true ); // beam --> tgt-hall
    TVector3 startTgt = startPoint;
    startPoint = this->ApplyUserRotation( startPoint, detOrigin, fDetRotation, true ); // tgt-hall --> det
    
    LOG( "NHL", pDEBUG )
      << "\n\n Unrotated startPoint: " << utils::print::Vec3AsString( &startUnrotated )
      << "\n Tgt-hall startPoint: " << utils::print::Vec3AsString( &startTgt )
      << "\n Final startPoint: " << utils::print::Vec3AsString( &startPoint );
  }
  TLorentzVector * p4NHL = event_rec->Particle(0)->GetP4();
  TVector3 momentum( p4NHL->Px(), p4NHL->Py(), p4NHL->Pz() );

  fSx = startPoint.X(); fSy = startPoint.Y(); fSz = startPoint.Z();
  fSxROOT = fSx * units::mm / units::cm;
  fSyROOT = fSy * units::mm / units::cm;
  fSzROOT = fSz * units::mm / units::cm;
  fPx = momentum.X(); fPy = momentum.Y(); fPz = momentum.Z();

  fGeomFile = geomfile;
  if( !fGeoManager )
    fGeoManager = TGeoManager::Import(geomfile.c_str());
}
//____________________________________________________________________________
bool NHLDecayVolume::VolumeEntryAndExitPoints( TVector3 & startPoint, TVector3 & momentum,
					       TVector3 & entryPoint, TVector3 & exitPoint,
					       TGeoManager * gm, TGeoVolume * /* vol */ ) const
{
  const double mmtolunits = units::mm / lunits;

  double sx = startPoint.X(); double sy = startPoint.Y(); double sz = startPoint.Z();
  sx *= mmtolunits; sy *= mmtolunits; sz *= mmtolunits;
  double px = momentum.X(); double py = momentum.Y(); double pz = momentum.Z();
  double p2 = px*px + py*py + pz*pz; double p = std::sqrt(p2);
  px *= 1./p; py *= 1./p; pz *= 1./p;

  fSx = sx; fSy = sy; fSz = sz;
  fSxROOT = fSx * lunits / units::cm; fSyROOT = fSy * lunits / units::cm; fSzROOT = fSz * lunits / units::cm;
  fPx = px; fPy = py; fPz = pz;

  // put first point slightly inside the bounding box
  double firstZOffset = -0.1; // m
  firstZOffset *= units::m / lunits;

  double firstZ = fOz - fLz/2.0 - firstZOffset;

  LOG( "NHL", pDEBUG )
    << "\nfirstZ = " << firstZ << " [" << lunitString.c_str() << "]";

  // now find which point the line would hit this z at
  double dz = firstZ - sz;
  double tz = dz / pz;
  double dx = tz * px;
  double dy = tz * py;
  double firstX = sx + dx;
  double firstY = sy + dy;

  // now we gotta return everything to cm for ROOT to work its magic.
  double firstXROOT = firstX * lunits / units::cm,
    firstYROOT = firstY * lunits / units::cm, firstZROOT = firstZ * lunits / units::cm;
  
  //assert( gm );
  gm = TGeoManager::Import(fGeomFile.c_str());
  gm->SetCurrentPoint( firstXROOT, firstYROOT, firstZROOT );
  gm->SetCurrentDirection( px, py, pz );

  LOG( "NHL", pINFO )
    << "\nCurrent point     is: ( " << firstX << ", " << firstY << ", " << firstZ << " ) [" << lunitString.c_str() << "]"
    << "\nFrom start point    : ( " << sx << ", " << sy << ", " << sz << " ) [" << lunitString.c_str() << "]"
    << "\nIn ROOT, current is : ( " << firstXROOT << ", " << firstYROOT << ", " << firstZROOT << " ) [cm]"
    << "\nIn ROOT, start is   : ( " << fSxROOT << ", " << fSxROOT << ", " << fSzROOT << " ) [cm]"
    << "\nCurrent direction is: ( " << px << ", " << py << ", " << pz << " ) [GeV/GeV]";

  assert( gm->FindNode() == NULL || gm->FindNode() == gm->GetTopNode() ); // need to be outside volume!

  double stepmax = 1.0e+6; // cm 
  stepmax *= genie::units::cm / lunits;

  // int ibound = 0;
  // const int imax = 10;

  LOG( "NHL", pDEBUG )
    << "Starting to search for intersections...";
  
  // enter the volume.
  TGeoNode * nextNode = gm->FindNextBoundaryAndStep( stepmax );
  // sometimes the TGeoManager likes to hit the BBox and call this an entry point. Step forward again.
  const double * tmpPoint = gm->GetCurrentPoint();
  if( std::abs(tmpPoint[0]) == fLx/2.0 * lunits / units::cm ||
      std::abs(tmpPoint[1]) == fLy/2.0 * lunits / units::cm ||
      std::abs(tmpPoint[2]) == fLz/2.0 * lunits / units::cm )
    nextNode = gm->FindNextBoundaryAndStep();

  if( nextNode == NULL ) return false;

  // entered the detector, let's save this point
  fEx = ( gm->GetCurrentPoint() )[0] * genie::units::cm / lunits;
  fEy = ( gm->GetCurrentPoint() )[1] * genie::units::cm / lunits;
  fEz = ( gm->GetCurrentPoint() )[2] * genie::units::cm / lunits;
  entryPoint.SetXYZ( fEx, fEy, fEz );

  fExROOT = ( gm->GetCurrentPoint() )[0];
  fEyROOT = ( gm->GetCurrentPoint() )[1];
  fEzROOT = ( gm->GetCurrentPoint() )[2];

  LOG( "NHL", pDEBUG )
    << "\nEntry point found at ( " << fEx << ", " << fEy << ", " << fEz << " ) [" << lunitString.c_str() << "]"
    << "\nIn ROOT, entry at    ( " << fExROOT << ", " << fEyROOT << ", " << fEzROOT << " ) [cm]"; 

  // now propagate until we exit again
  
  int bdIdx = 0;
  const int bdIdxMax = 1e+4;

  double sfx = 0.0, sfy = 0.0, sfz = 0.0; // coords of the "safe" points in user units
  double sfxROOT = 0.0, sfyROOT = 0.0, sfzROOT = 0.0; // same, in cm

  // do one big step first
  // then if not outside yet, step by ever smaller steps until some threshold
  //Double_t sNext = std::max( fLx, std::max( fLy, fLz ) ) / 2.0;
  Double_t sNext = std::min( std::max( fLx, std::max( fLy, fLz ) ), 100.0 * lunits / units::cm ) / 2.0;
  Double_t sNextROOT = sNext * lunits / units::cm;
  gm->SetStep( sNextROOT );
  LOG( "NHL", pINFO )
    << "fLx, fLy, fLz = " << fLx << ", " << fLy << ", " << fLz << " ==> sNextROOT = " << sNextROOT;
  gm->Step();
  
  // FindNextBoundaryAndStep() sets step size to distance to next boundary and executes that step
  // so one "step" here is actually one big step + one small step
  while( gm->FindNextBoundaryAndStep() && bdIdx < bdIdxMax ){
    const Double_t * currPoint = gm->GetCurrentPoint();
    if( bdIdx % 100 == 0 ){
      LOG( "NHL", pDEBUG )
	<< "Step " << bdIdx << " : ( " << currPoint[0] << ", " << currPoint[1] << ", " << currPoint[2] << " ) [cm]";
    }
    sfxROOT = currPoint[0]; sfyROOT = currPoint[1]; sfzROOT = currPoint[2];
    if( sNextROOT >= 2.0 * lunits / units::cm ) sNextROOT *= 0.5;
    gm->SetStep( sNextROOT );
    gm->Step();
    bdIdx++;
  }
  if( bdIdx == bdIdxMax ){
    LOG( "NHL", pWARN )
      << "Failed to exit this volume. Dropping this trajectory.";
    return false;
  }

  // Always step back one step
  /*
  const Double_t * ffPoint = gm->GetCurrentPoint();
  if( std::abs(ffPoint[0] - fOxROOT) > fLxROOT/2.0 || 
      std::abs(ffPoint[1] - fOyROOT) > fLyROOT/2.0 || 
      std::abs(ffPoint[2] - fOzROOT) > fLzROOT/2.0 ){
    LOG( "NHL", pDEBUG )
      << "Overstepped bounding box: we're at ( " << ffPoint[0] << ", " << ffPoint[1] << ", " << ffPoint[2] << " ) [cm]";
    const Double_t * sfDir = gm->GetCurrentDirection();
    gm->SetCurrentDirection( -sfDir[0], -sfDir[1], -sfDir[2] );
    __attribute__((unused)) TGeoNode * tmpNode = gm->FindNextBoundaryAndStep();
    LOG( "NHL", pDEBUG )
      << "We turned back with new step = " << gm->GetStep();
    // and set direction back to normal
    gm->SetCurrentDirection( sfDir[0], sfDir[1], sfDir[2] );
  }
  const Double_t * sfPoint = gm->GetCurrentPoint();
  sfxROOT = sfPoint[0]; sfyROOT = sfPoint[1]; sfzROOT = sfPoint[2];
  */
  sfx = sfxROOT * units::cm / lunits; sfy = sfyROOT * units::cm / lunits; sfz = sfzROOT * units::cm / lunits;

  // exited the detector, let's save this point
  fXx = sfx; fXxROOT = sfxROOT;
  fXy = sfy; fXyROOT = sfyROOT;
  fXz = sfz; fXzROOT = sfzROOT;
  exitPoint.SetXYZ( fXx, fXy, fXz );

  LOG( "NHL", pINFO )
    << "\nExit point found at ( " << fXx << ", " << fXy << ", " << fXz << " ) ["
    << lunitString.c_str() << "]"
    << "\nIn ROOT, exit at    ( " << fXxROOT << ", " << fXyROOT << ", " << fXzROOT << " ) [cm]"; 

  return true;
  
}
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__
//____________________________________________________________________________
void NHLDecayVolume::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NHLDecayVolume::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NHLDecayVolume::LoadConfig()
{
  if( fIsConfigLoaded ) return;

  LOG( "NHL", pDEBUG )
    << "Loading geometry parameters from file. . .";

  this->GetParam( "UseBeamMomentum", fUseBeamMomentum );
  this->GetParamVect( "Beam2User_T", fB2UTranslation );
  this->GetParamVect( "Beam2User_R", fB2URotation );
  this->GetParamVect( "Beam2Det_R", fDetRotation );
  this->GetParamVect( "DetCentre_User", fDetTranslation );
  fCx = fB2UTranslation.at(0); fCy = fB2UTranslation.at(1); fCz = fB2UTranslation.at(2);
  fUx = fDetTranslation.at(0); fUy = fDetTranslation.at(1); fUz = fDetTranslation.at(2);
  fAx1 = fB2URotation.at(0); fAz = fB2URotation.at(1); fAx2 = fB2URotation.at(2);
  fBx1 = fDetRotation.at(0); fBz = fDetRotation.at(1); fBx2 = fDetRotation.at(2);

  fIsConfigLoaded = true;
}
//____________________________________________________________________________
void NHLDecayVolume::GetInterestingPoints( TVector3 & entryPoint, TVector3 & exitPoint, TVector3 & decayPoint ) const
{
  LOG( "NHL", pDEBUG ) << "Getting interesting points...";
  entryPoint.SetXYZ( fEx, fEy, fEz );
  exitPoint.SetXYZ( fXx, fXy, fXz );
  decayPoint.SetXYZ( fDx, fDy, fDz );
}
//____________________________________________________________________________
TVector3 NHLDecayVolume::ApplyUserRotation( TVector3 vec, bool doBackwards ) const
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
TVector3 NHLDecayVolume::ApplyUserRotation( TVector3 vec, TVector3 oriVec, std::vector<double> rotVec, bool doBackwards ) const
{
  double vx = vec.X(), vy = vec.Y(), vz = vec.Z();
  double ox = oriVec.X(), oy = oriVec.Y(), oz = oriVec.Z();

  LOG( "NHL", pDEBUG )
    << "\t original vec [mm] : " << utils::print::Vec3AsString( &vec )
    << "\t origin vec [mm] : " << utils::print::Vec3AsString( &oriVec );
  
  vx -= ox; vy -= oy; vz -= oz; // make this rotation about detector origin

  assert( rotVec.size() == 3 ); // want 3 Euler angles, otherwise this is unphysical.
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
  LOG( "NHL", pDEBUG )
    << "\nVector is now " << utils::print::Vec3AsString( &nvec );
  return nvec;
}
//____________________________________________________________________________
//____________________________________________________________________________
/*!

\class     genie::llp::VertexGenerator

\brief     LLP vertex generator
           Takes a point (x0, y0, z0) and a momentum (px, py, pz) as inputs
	   Returns two points (x1, y1, z1) and (x2, y2, z2), 
	   and an event vertex (X, Y, Z, T).

	   (x1/2, y1/2, z1/2) are entry/exit points of the ray 
	   passing by (x0, y0, z0) parallel to (px, py, pz)

	   T is the quantity ( 1 / (beta * c) - 1 ) * ( || ( X - x0, Y - y0, Z - z0 ) || ),
	   i.e. the delay of a massive particle of velocity beta wrt a speed-of-light particle

\author    John Plows <kplows \at liverpool.ac.uk>
           University of Liverpool
          
\created   May 23rd, 2024

\cpright   Copyright (c) 2003-2023, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _LLP_VERTEX_GENERATOR_H_
#define _LLP_VERTEX_GENERATOR_H_

#include <cmath>
#include <cassert>
#include <list>

#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>

//#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/PrintUtils.h"

#include "Physics/ExoticLLP/LLPGeomRecordVisitorI.h"

namespace genie {
namespace llp {

  class ExoticLLP;

  class VertexGenerator : public GeomRecordVisitorI {

  public:

    VertexGenerator();
    VertexGenerator(string name);
    VertexGenerator(string name, string config);
    ~VertexGenerator();

    //-- implement the EventRecordVisitorI interface
    void ProcessEventRecord(GHepRecord * event_rec) const;

    // overload the Algorithm::Configure() methods to load private data
    // members from configuration options
    void Configure(const Registry & config);
    void Configure(string config);

    void SetGeomFile( std::string geomfile, std::string topVolume ) const;

  private:

    void LoadConfig();

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
    // use bounding box origin & sides
    void ImportBoundingBox( TGeoBBox * box ) const;
    TGeoMatrix * FindFullTransformation( TGeoVolume * top_vol, TGeoVolume * tar_vol ) const;
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

    void SetStartingParameters( GHepRecord * event_rec ) const;
    
    // --------------------------------------------------
    // Utilities
    // --------------------------------------------------

    // simple getter
    void GetInterestingPoints( TVector3 & entryPoint, TVector3 & exitPoint, TVector3 & decayPoint ) const;

    // Build a simple 1x1x1 m3 box around (0,0,0).
    void MakeSDV() const;

    // enforce chosen units
    void EnforceUnits( std::string length_units, std::string angle_units, std::string time_units ) const;

    // calculate travel length in detector
    double CalcTravelLength( double betaMag, double CoMLifetime, double maxLength ) const;

    // assign decay point given length
    TVector3 GetDecayPoint( double travelLength, TVector3 & entryPoint, TVector3 & momentum ) const;

    // get max length in detector
    double GetMaxLength( TVector3 & entryPoint, TVector3 & exitPoint ) const;

    // --------------------------------------------------
    // Simple Decay Volume - fallback if no geometry
    // --------------------------------------------------

    // in case of SDV, calculate entry and exit points (if any) of some trajectory
    // if no entry/exit points return false
    bool SDVEntryAndExitPoints( TVector3 & startPoint, TVector3 momentum,
				TVector3 & entryPoint,
				TVector3 & exitPoint ) const;

    // --------------------------------------------------
    // ROOT geometry stuff
    // --------------------------------------------------

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
    
    // get entry & exit points directly from volume
    bool VolumeEntryAndExitPoints( TVector3 & startPoint, TVector3 & momentum,
				   TVector3 & entryPoint, TVector3 & exitPoint, 
				   TGeoManager * gm, TGeoVolume * vol ) const;
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

    TVector3 ApplyUserRotation( TVector3 vec, bool doBackwards ) const;
    TVector3 ApplyUserRotation( TVector3 vec, TVector3 oriVec, std::vector<double> rotVec, bool doBackwards ) const;

    std::string CheckGeomPoint( Double_t x, Double_t y, Double_t z ) const;

    // --------------------------------------------------

    mutable bool fIsConfigLoaded = false;

    mutable double lunits = genie::units::mm; mutable std::string lunitString = "mm";
    mutable double aunits = genie::units::rad;
    mutable double tunits = genie::units::ns; mutable std::string tunitString = "ns";

    mutable double fSx = 0.0, fSy = 0.0, fSz = 0.0; //start point
    mutable double fNx = 0.0, fNy = 0.0, fNz = 0.0; //origin point of HNL
    mutable double fPx = 0.0, fPy = 0.0, fPz = 0.0; //momentum
    mutable double fEx = 0.0, fEy = 0.0, fEz = 0.0; //entry point
    mutable double fXx = 0.0, fXy = 0.0, fXz = 0.0; //exit  point

    mutable double fSxROOT = 0.0, fSyROOT = 0.0, fSzROOT = 0.0; // start point in cm
    mutable double fNxROOT = 0.0, fNyROOT = 0.0, fNzROOT = 0.0; // origin point of HNL in cm
    mutable double fExROOT = 0.0, fEyROOT = 0.0, fEzROOT = 0.0; // entry point in cm
    mutable double fXxROOT = 0.0, fXyROOT = 0.0, fXzROOT = 0.0; // exit  point in cm

    mutable double fDx = 0.0, fDy = 0.0, fDz = 0.0; //decay point
    mutable double fOx = 0.0, fOy = 0.0, fOz = 0.0; //origin
    mutable double fLx = 0.0, fLy = 0.0, fLz = 0.0; //dimensions

    mutable double fDxROOT = 0.0, fDyROOT = 0.0, fDzROOT = 0.0; // decay point in cm
    mutable double fOxROOT = 0.0, fOyROOT = 0.0, fOzROOT = 0.0; // origin in cm
    mutable double fLxROOT = 0.0, fLyROOT = 0.0, fLzROOT = 0.0; // dimensions in cm

    mutable double fCoMLifetime = 0.0; // LLP lifetime in ns
    
    mutable double kNewSpeedOfLight = genie::units::kSpeedOfLight * ( genie::units::m / lunits ) / ( genie::units::s / tunits );

    mutable string fGeomFile = "";
    mutable string fTopVolume = "";
    mutable TGeoManager * fGeoManager = 0;
    mutable TGeoVolume * fGeoVolume = 0;
    
    mutable bool isUsingDk2nu = false;
    mutable bool isParticleGun = false;
    mutable bool isUsingRootGeom = false;
    mutable double uMult = 1.0, xMult = 1.0; // these need to be different.

    mutable double fCx, fCy, fCz;   // translation: from NEAR origin to USER origin [m]
    mutable double fTx, fTy, fTz;   // If user --top_volume passed, where's its origin in USER frame? [m]
    mutable double fUx, fUy, fUz;   // translation: from USER origin to detector centre [m]
    mutable double fAx1, fAz, fAx2; // rotation: from NEAR frame to USER frame
    mutable std::vector< double > fUserOrigin, fDetTranslation, fDetRotation;

    ClassDef(VertexGenerator, 1)
  }; // class VertexGenerator

} // namespace llp
} // namespace genie

#endif // #ifndef _LLP_VERTEX_GENERATOR_H_

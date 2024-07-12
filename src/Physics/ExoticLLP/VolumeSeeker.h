//----------------------------------------------------------------------------
/*!

\class    genie::llp::VolumeSeeker

\brief    Responsible for performing calculations in 3D space and transporting
          an LLP to a detector, calculating lab-frame solid angles,
          getting entry and exit points of ray in 3D space.

	  Due to its nature it will be a singleton.

	  The NEAR (read, global) coordinate system will be indexed by X, Y, Z
	  The USER (read, detector) coordinate system will be indexed by x, y, z

\author   John Plows <kjplows \at liverpool.ac.uk>
          University of Liverpool

\created  May 27th, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _LLP_VOLUMESEEKER_H_
#define _LLP_VOLUMESEEKER_H_

#include <cmath>
#include <cassert>
#include <list>
#include <tuple>

#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>

#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/PrintUtils.h"

typedef std::pair< double, double > Point; //! To store (theta, phi) values
typedef std::pair< Point, Point > PointRaster; //! To associate (thetamin, thetamax) for each phi
typedef std::vector< PointRaster > AngularRegion; //! The bounding shape, in (theta, phi) space

namespace genie {

  namespace llp {

    class VolumeSeeker {

    public:

      //! Access instance
      static VolumeSeeker * Instance();

      //! Let the VolumeSeeker know about which geometry file you're using
      void SetGeomFile( std::string geomfile, std::string topVolume ) const;

      //! Have some module tell you the configuration and DO NOT allow more changes
      void SetConfig( TVector3 user_origin, TVector3 user_rotation );
      //! Output the current config
      void PrintConfig();

      //! Populate all the current members. Accepts input in NEAR coordinates only
      void PopulateEvent( TVector3 origin_point, TVector3 momentum ) const;
      //! Clear all the current members 
      void ClearEvent() const;

      //! Workhorse methods for NEAR <--> USER transformations
      TVector3 Translate( TVector3 input, bool direction ) const;
      TVector3 Rotate( TVector3 input, bool direction ) const;
      TVector3 TranslateToUser( TVector3 input ) const { return VolumeSeeker::Translate( input, true ); }
      TVector3 TranslateToNear( TVector3 input ) const { return VolumeSeeker::Translate( input, false ); }
      TVector3 RotateToUser( TVector3 input ) const { return VolumeSeeker::Rotate( input, true ); }
      TVector3 RotateToNear( TVector3 input ) const { return VolumeSeeker::Rotate( input, false ); }

      //! Given an origin point and a momentum, find the entry and exit points to the detector
      // RETHERE make private
      bool RaytraceDetector( bool grace = false ) const;
      
      //! Define a region in (theta, phi) space that an HNL can be accepted in...
      //! Note that (theta, phi) are angles with respect to the parent momentum
      AngularRegion AngularAcceptance() const;
      //! And calculate its size
      double AngularSize( AngularRegion alpha ) const;
      double Trapezoid( std::vector<Point> up_vec, std::vector<Point> dn_vec ) const;
      double Simpson( std::vector<Point> pt_vec ) const;

      //! A method that allows VolumeSeeker to read in configurable controls without needing to
      //  inherit from Algorithm. This means it can be an instanced class (which is good)
      void AdoptControls( double ct, double cp, double ft, double fp, double gr ) const;

      //! Getter methods to interface with FluxContainer
      TVector3 GetEntryPoint( bool near = false ) const 
      { return ( near ) ? fEntryPointNEAR : fEntryPoint; }
      TVector3 GetExitPoint( bool near = false ) const 
      { return ( near ) ? fExitPointNEAR : fExitPoint; }

    private:

      VolumeSeeker();
      VolumeSeeker( const VolumeSeeker & vsek );
      virtual ~VolumeSeeker();

      //! Build the bounding box and find the top volume
      void ImportBoundingBox( TGeoBBox * box ) const;
      TGeoMatrix * FindFullTransformation( TGeoVolume * top_vol, TGeoVolume * tar_vol ) const;

      //! Obtain the node (in the ROOT sense) where a point exists
      std::string CheckGeomPoint( TVector3 chkpoint ) const;

      //! Given an origin point and a momentum, find the entry and exit points to the detector
      //bool RaytraceDetector() const;

      //! Some controls
      mutable double m_coarse_theta_deflection = 2.0; // modifier
      mutable double m_fine_theta_deflection = 2.0; // modifier
      mutable double m_coarse_phi_deflection = 5.0; // modifier
      mutable double m_fine_phi_deflection = 5.0; // modifier
      mutable double m_grace_decrement = 0.5e-2; 

      //! And utility functions for calling Raytrace() a lot of times
      void Deflect( double & deflection, bool goUp ) const; // calls Raytrace() with set theta, phi
      void Rasterise( AngularRegion & alpha, bool goRight ) const; // calls Deflect() with set phi
      // Note: alpha --> measures the deflections from parent momentum, 
      //        beta --> measures the directional cosines of the geometry profile in USER
      void ConvertToUserAngles( const TVector3 booked_momentum, AngularRegion deflection_region,
				AngularRegion & cosines_region ) const;

      static VolumeSeeker * fInstance;

      bool fInitialized; //! Done initialising singleton?

      mutable double fLunits = genie::units::m;  mutable std::string fLunitString = "m";
      mutable double fAunits = genie::units::rad; mutable std::string fAunitString = "rad";
      mutable double fTunits = genie::units::ns;  mutable std::string fTunitString = "ns";

      mutable double fToROOTUnits = fLunits / genie::units::cm;
      mutable double fToLUnits = 1.0 / fToROOTUnits;
      
      //! Config options for NEAR <--> USER transformations
      mutable TVector3 fUserOrigin;    //! USER (0, 0, 0) in NEAR coords [m]
      mutable TVector3 fUserRotation;  //! USER Euler angles (a, b, c) in NEAR coords. Extrinsic X-Z-X from NEAR --> USER. I.e. (x, y, z) = R_X(c) * R_Z(b) * R_X(a) * (X, Y, Z)
      mutable bool fIsConfigLoaded; 
      mutable bool fIsGeomFileSet;
      
      //! All vectors with ROOT in the name use the ROOT default units [cm] and refer to the local coordinate system of the top_volume node
      mutable TVector3 fZeroPoint, fZeroPointROOT;     //! Ray intercept at z = 0 plane. [m]
      mutable TVector3 fOriginPoint, fOriginPointROOT; //! (x, y, z) of LLP production vertex [m]
      mutable TVector3 fEntryPoint, fEntryPointROOT;   //! (x, y, z) of ray entry into detector [m]
      mutable TVector3 fExitPoint, fExitPointROOT;     //! (x, y, z) of ray exit from detector [m]

      //! Copies of the above members, but in the global frame
      mutable TVector3 fZeroPointNEAR;                 //! Ray intercept at z = 0 plane. [m]
      mutable TVector3 fOriginPointNEAR;               //! (X, Y, Z) of LLP production vertex [m]
      mutable TVector3 fEntryPointNEAR;                //! (X, Y, Z) of ray entry into detector [m]
      mutable TVector3 fExitPointNEAR;                 //! (X, Y, Z) of ray exit from detector [m]

      mutable TVector3 fMomentum, fMomentumNEAR;       //! (px, py, pz) of particle. [GeV/c]

      /*
	The following 3 vectors describe a RH coordinate system in USER space.
	Axis is the momentum axis (unit vector)
	ThetaAxis is the vector that is perpendicular to Axis, and along with Axis defines a plane.
	This plane contains both the OriginPoint as well as the TopVolumeCentre.
	PhiAxis completes the coordinate system (== Axis.Cross(ThetaAxis))
       */
      mutable TVector3 fAxis, fThetaAxis, fPhiAxis;    //! A right handed coordinate system in USER space.
      mutable double fThetaSeed, fPhiSeed;      //! The (theta, phi) deflection needed to get from fMomentum to fTopVolumeOrigin from fOriginPoint

      mutable std::string fGeomFile;            //! Path to the geometry file
      mutable std::string fTopVolume;           //! Name of the top volume to be used
      mutable TGeoManager * fGeoManager;        //! ROOT TGeoManager
      mutable TGeoVolume * fGeoVolume;          //! ROOT TGeoVolume

      mutable double fLx, fLy, fLz;             //! Bounding box dimensions [m]
      mutable double fLxROOT, fLyROOT, fLzROOT; //! Bounding box dimensions [cm]
      mutable double fOx, fOy, fOz;             //! Bounding box origin [m]
      mutable double fOxROOT, fOyROOT, fOzROOT; //! Bounding box origin [cm]

      mutable TVector3 fTopVolumeOrigin;        //! Origin of top_volume in USER coords [m]
      mutable TVector3 fTopVolumeOriginNEAR;    //! Origin of top_volume in NEAR coords [m]

      struct Cleaner {
	void DummyMethodAndSilentCompiler() {} // see RandomGen.h
	~Cleaner() {
	  if( VolumeSeeker::fInstance != 0 ) {
	    delete VolumeSeeker::fInstance;
	    VolumeSeeker::fInstance = 0;
	  }
	}
      }; // struct Cleaner

      friend struct Cleaner;

      ClassDef(VolumeSeeker, 1)
    }; // class VolumeSeeker

  } // namespace llp

} // namespace genie

#endif // #ifndef _LLP_VOLUMESEEKER_H_

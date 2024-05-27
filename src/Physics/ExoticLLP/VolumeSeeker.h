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
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/PrintUtils.h"

namespace genie {

  namespace llp {

    class VolumeSeeker {

    public:

      //! Access instance
      static VolumeSeeker * Instance();

      //! Have some module tell you the configuration and DO NOT allow more changes
      void SetConfig( TVector3 user_origin, TVector3 user_rotation );
      //! Output the current config
      void PrintConfig();

      //! Clear all the current members 
      

    private:

      VolumeSeeker();
      VolumeSeeker( const VolumeSeeker & vsek );
      virtual ~VolumeSeeker();

      static VolumeSeeker * fInstance;

      bool fInitialized; //! Done initialising singleton?

      mutable double fLunits = genie::units::mm;  mutable std::string fLunitString = "mm";
      mutable double fAunits = genie::units::rad; mutable std::string fAunitString = "rad";
      mutable double fTunits = genie::units::ns;  mutable std::string fTunitString = "ns";
      
      //! Config options for NEAR <--> USER transformations
      mutable TVector3 fUserOrigin;    //! USER (0, 0, 0) in NEAR coords [m]
      mutable TVector3 fUserRotation;  //! USER Euler angles (a, b, c) in NEAR coords. Extrinsic X-Z-X from NEAR --> USER. I.e. (x, y, z) = R_X(c) * R_Z(b) * R_X(a) * (X, Y, Z)
      bool fIsConfigLoaded; 
      
      //! All vectors with ROOT in the name use the ROOT default units [cm].
      mutable TVector3 fZeroPoint, fZeroPointROOT;     //! Ray intercept at z = 0 plane. [m]
      mutable TVector3 fOriginPoint, fOriginPointROOT; //! (x, y, z) of LLP production vertex [m]
      mutable TVector3 fEntryPoint, fEntryPointROOT;   //! (x, y, z) of ray entry into detector [m]
      mutable TVector3 fExitPoint, fExitPointROOT;     //! (x, y, z) of ray exit from detector [m]

      //! Copies of the above members, but in the global frame
      mutable TVector3 fZeroPointNEAR;                 //! Ray intercept at z = 0 plane. [m]
      mutable TVector3 fOriginPointNEAR;               //! (X, Y, Z) of LLP production vertex [m]
      mutable TVector3 fEntryPointNEAR;                //! (X, Y, Z) of ray entry into detector [m]
      mutable TVector3 fExitPointNEAR;                 //! (X, Y, Z) of ray exit from detector [m]

      string fGeomFile;           //! Path to the geometry file
      string fTopVolume;          //! Name of the top volume to be used
      TGeoManager * fGeoManager;   //! ROOT TGeoManager
      TGeoVolume * fGeoVolume;     //! ROOT TGeoVolume

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

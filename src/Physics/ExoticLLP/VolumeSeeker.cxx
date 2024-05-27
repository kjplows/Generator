//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <kjplows \at liverpool.ac.uk>
 University of Liverpool

*/
//____________________________________________________________________________

#include "Physics/ExoticLLP/VolumeSeeker.h"

using namespace genie;
using namespace genie::llp;

//____________________________________________________________________________
VolumeSeeker * VolumeSeeker::fInstance = 0;
bool fIsConfigLoaded = false;
string fGeomFile = "";
string fTopVolume = "";
TGeoManager * fGeoManager = 0;
TGeoVolume * fGeoVolume = 0;
//____________________________________________________________________________
VolumeSeeker::VolumeSeeker()
{
  LOG("ExoticLLP", pINFO) << "VolumeSeeker late initialization";

  fInitialized = false;
  fInstance = 0;
  fIsConfigLoaded = false;
  fGeomFile = "";
  fTopVolume = "";
  fGeoManager = 0;
  fGeoVolume = 0;
}
//____________________________________________________________________________
VolumeSeeker::~VolumeSeeker()
{
  fInstance = 0;
  fIsConfigLoaded = false;

  if( fGeoManager ) delete fGeoManager;
  if( fGeoVolume ) delete fGeoVolume;
}
//____________________________________________________________________________
VolumeSeeker * VolumeSeeker::Instance()
{
  if( fInstance == 0 ) {
    static VolumeSeeker::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new VolumeSeeker;
  }
  return fInstance;
}
//____________________________________________________________________________
void VolumeSeeker::SetConfig( TVector3 user_origin, TVector3 user_rotation )
{
  if( fIsConfigLoaded ){
    LOG( "ExoticLLP", pWARN ) << "NEAR --> USER transformation data already loaded. Doing nothing."
			      << " Something in the code is attempting to re-initialise VolumeSeeker!";
    return;
  }

  fUserOrigin = user_origin;
  fUserRotation = user_rotation;
  fIsConfigLoaded = true;
}
//____________________________________________________________________________
void VolumeSeeker::PrintConfig()
{
  LOG( "ExoticLLP", pDEBUG ) << "\nNEAR --> USER translation: " << utils::print::Vec3AsString( &fUserOrigin ) << " [m]"
			     << "\nNEAR --> USER rotation: " << utils::print::Vec3AsString( &fUserRotation ) << " [rad]";
}

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
  fInitialized = false;
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
void VolumeSeeker::SetGeomFile( std::string geomfile, std::string topVolume ) const
{
  fGeomFile = geomfile;
  fTopVolume = topVolume;
  
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
    TGeoMatrix * hmat = VolumeSeeker::FindFullTransformation( main_volume, top_volume );
    LOG( "ExoticLLP", pDEBUG )
      << "Got translation of volume with name " << top_volume->GetName() << " which is "
      << utils::print::Vec3AsString( &fTopVolumeOrigin ) << " [m]";
  }
  fGeoManager->SetTopVolume(top_volume);
  TGeoShape * ts = top_volume->GetShape();
  TGeoBBox * box = (TGeoBBox *) ts;
    
  VolumeSeeker::ImportBoundingBox( box );
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
//____________________________________________________________________________
TVector3 VolumeSeeker::Translate( TVector3 input, bool direction ) const
{
  // true: NEAR --> USER | false : USER --> NEAR
  TVector3 tr_vec = ( direction == true ) ? fUserOrigin : -fUserOrigin;
  return input - tr_vec;
}
//____________________________________________________________________________
TVector3 VolumeSeeker::Rotate( TVector3 input, bool direction ) const
{
  // true: NEAR --> USER | false : USER --> NEAR
  // inverse of matrix multiplication is (CBA)^-1 = A^-1 B^-1 C^-1
  TVector3 rt_vec = ( direction == true ) ? fUserRotation : 
    TVector3( -fUserRotation.Z(), -fUserRotation.Y(), -fUserRotation.X() );

  // Apply Ax1 then Az then Ax2, i.e. R_X(Ax2) R_Z(Az) R_X(Ax1)
  double Ax1 = rt_vec.X(); double Az = rt_vec.Y(); double Ax2 = rt_vec.Z();

  // Initialise
  double x = input.X(); double y = input.Y(); double z = input.Z();
  double vx = x; double vy = y; double vz = z;
  
  // Rotate: R_X (Ax1)
  vy = y * std::cos( Ax1 ) - z * std::sin( Ax1 );
  vz = y * std::sin( Ax1 ) + z * std::cos( Ax1 );
  y = vy; z = vz; // update
  
  // Rotate: R_Z (Az)
  vx = x * std::cos( Az ) - y * std::sin( Az );
  vy = x * std::sin( Az ) + y * std::cos( Az );
  x = vx; y = vy; // update

  // Rotate: R_X (Ax2)
  vy = y * std::cos( Ax2 ) - z * std::sin( Ax2 );
  vz = y * std::sin( Ax2 ) + z * std::cos( Ax2 );
  y = vy; z = vz; // update

  return TVector3( x, y, z );
}
//____________________________________________________________________________
void VolumeSeeker::PopulateEvent( TVector3 origin_point, TVector3 momentum ) const
{
  fOriginPointNEAR = origin_point;
  fMomentumNEAR = momentum.Unit();

  // The momentum just needs to be rotated to USER
  fMomentum = VolumeSeeker::RotateToUser( fMomentumNEAR );

  // The origin needs to be translated to USER origin first and then rotated
  TVector3 translated_origin = VolumeSeeker::TranslateToUser( fOriginPointNEAR );
  fOriginPoint = VolumeSeeker::RotateToUser( translated_origin );

  fOriginPointROOT = fToROOTUnits * fOriginPoint;
}
//____________________________________________________________________________
void VolumeSeeker::ClearEvent() const
{
  // clear the current members... Not the detector as that's run-invariant
  fZeroPoint.SetXYZ( 0.0, 0.0, 0.0 ); fZeroPointROOT.SetXYZ( 0.0, 0.0, 0.0 );
  fOriginPoint.SetXYZ( 0.0, 0.0, 0.0 ); fOriginPointROOT.SetXYZ( 0.0, 0.0, 0.0 );
  fEntryPoint.SetXYZ( 0.0, 0.0, 0.0 ); fEntryPointROOT.SetXYZ( 0.0, 0.0, 0.0 );
  fExitPoint.SetXYZ( 0.0, 0.0, 0.0 ); fExitPointROOT.SetXYZ( 0.0, 0.0, 0.0 );

  fZeroPointNEAR.SetXYZ( 0.0, 0.0, 0.0 );
  fOriginPointNEAR.SetXYZ( 0.0, 0.0, 0.0 );
  fEntryPointNEAR.SetXYZ( 0.0, 0.0, 0.0 );
  fExitPointNEAR.SetXYZ( 0.0, 0.0, 0.0 );

  fMomentum.SetXYZ( 0.0, 0.0, 0.0 ); fMomentumNEAR.SetXYZ( 0.0, 0.0, 0.0 );
}
//____________________________________________________________________________
void VolumeSeeker::ImportBoundingBox( TGeoBBox * box ) const 
{
  fLx = box->GetDX() * fToLUnits;
  fLy = box->GetDY() * fToLUnits;
  fLz = box->GetDZ() * fToLUnits;
  fOx = (box->GetOrigin())[0] * fToLUnits;
  fOy = (box->GetOrigin())[1] * fToLUnits;
  fOz = (box->GetOrigin())[2] * fToLUnits;

  fLxROOT = box->GetDX();
  fLyROOT = box->GetDY();
  fLzROOT = box->GetDZ();
  fOxROOT = (box->GetOrigin())[0];
  fOyROOT = (box->GetOrigin())[1];
  fOzROOT = (box->GetOrigin())[2];
}
//____________________________________________________________________________
TGeoMatrix * VolumeSeeker::FindFullTransformation( TGeoVolume * top_vol, TGeoVolume * tar_vol ) const
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

  LOG( "ExoticLLP", pNOTICE )
    << "Looking for this targetPath: " << targetPath;

  // could be we hit the top volume, in which case we skip the loop
  bool foundPath = (strcmp( test.c_str(), targetPath.c_str() ) == 0);

  while( strcmp( test.c_str(), targetPath.c_str() ) != 0 && !foundPath ){ // still looking
    TGeoNode * node = nodes.front();
    std::string path = paths.front();
    TGeoMatrix * mat = mats.front();

    assert( node  && "Node is not null" );
    assert( mat && "Matrix is not null" );

    int nDaughters = node->GetNdaughters();
    LOG( "ExoticLLP", pDEBUG ) << "Node with name " << path << " has " << nDaughters << " daughters...";
    for( int iDaughter = 0; iDaughter < nDaughters; iDaughter++ ){
      TGeoNode * dNode = node->GetDaughter(iDaughter);
      assert( dNode && "Daughter node not null" );
      std::string dPath( path );
      dPath.append( "/" ); dPath.append( dNode->GetName() );
      TGeoMatrix * nodeMat = dNode->GetMatrix();

      LOG( "ExoticLLP", pDEBUG ) << "Got node, path, and matrix for daughter node "
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
      LOG( "ExoticLLP", pDEBUG )
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

  LOG( "ExoticLLP", pINFO )
    << "Found the target volume! Here is its path and full matrix:"
    << "\nPath: " << paths.back()
    << "\nTranslations: ( " << final_tra[0] << ", " << final_tra[1] << ", " << final_tra[2]
    << " ) [cm]"
    << "\nRotation matrix: ( ( " 
    << final_rot[0] << ", " << final_rot[1] << ", " << final_rot[2] << " ), ( "
    << final_rot[3] << ", " << final_rot[4] << ", " << final_rot[5] << " ), ( "
    << final_rot[6] << ", " << final_rot[7] << ", " << final_rot[8] << " ) )";

  // Also set the member variables at this stage
  //fTopVolumeOriginROOT.SetXYZ( final_tra[0], final_tra[1], final_tra[2] );
  //fTopVolumeOrigin = fTopVolumeOriginROOT * fToLUnits;
  fTopVolumeOrigin.SetXYZ( final_tra[0] * fToLUnits, 
			   final_tra[1] * fToLUnits, 
			   final_tra[2] * fToLUnits );
  fTopVolumeOriginNEAR = VolumeSeeker::RotateToNear( fTopVolumeOrigin );
  fTopVolumeOriginNEAR = VolumeSeeker::TranslateToNear( fTopVolumeOriginNEAR );

  return final_mat;
}
//____________________________________________________________________________
std::string VolumeSeeker::CheckGeomPoint( TVector3 chkpoint ) const
{
  Double_t point[3] = { chkpoint.X(), chkpoint.Y(), chkpoint.Z() };
  Double_t local[3];
  TGeoVolume * vol = fGeoManager->GetVolume( fTopVolume.c_str() );
  TGeoNode * node = fGeoManager->FindNode( point[0], point[1], point[2] );
  fGeoManager->MasterToLocal( point, local ); // don't know why but just do it
  return fGeoManager->GetPath();
}
//____________________________________________________________________________
bool VolumeSeeker::RaytraceDetector( bool grace ) const
{
  // Our point starts out at fOriginPoint and has directional cosines fMomentum
  
  // Important subtlety! The ROOT coordinate system for a top_volume does not know about the transformation matrix
  // Practically, this means that you must subtract fTopVolumeOriginROOT from all your ROOT calcs!!!!
  double t_param = 0.0; TVector3 dev_vec = fTopVolumeOrigin - fOriginPoint;
  if( dev_vec.Mag() > 0.0 ) {
    // First, calculate the starting point: intercept of ray with the z = 0 plane
    // Important caveat: z = 0 is z = 0 of *detector*. This is not the same as z = 0 USER
    double delta = dev_vec.Z(); // dz
    double mom_cos = fMomentum.Z(); // pz,hat
    if( delta == 0.0 ){ delta = dev_vec.Y(); mom_cos = fMomentum.Y(); } // make it dy instead
    if( delta == 0.0 ){ delta = dev_vec.X(); mom_cos = fMomentum.X(); } // OK, it *has* to be dx

    t_param = delta / mom_cos;
  } // get to z = 0 (or y = 0, or x = 0) plane

  // transport to starting point (hopefully close enough to the top volume to be meaningful)
  fZeroPoint.SetXYZ( fOriginPoint.X() + t_param * fMomentum.X(),
		     fOriginPoint.Y() + t_param * fMomentum.Y(),
		     fOriginPoint.Z() + t_param * fMomentum.Z() );
  fZeroPointNEAR = VolumeSeeker::RotateToNear( fZeroPoint );
  fZeroPointNEAR = VolumeSeeker::TranslateToNear( fZeroPointNEAR );
  fZeroPointROOT = (fZeroPoint - fTopVolumeOrigin) * fToROOTUnits; // subtract translation subtlety

  // check that this point lies in the geometry.
  std::string pathString = VolumeSeeker::CheckGeomPoint( fZeroPointROOT );

  //LOG( "ExoticLLP", pDEBUG ) << "Checking point " << utils::print::Vec3AsString( &fZeroPointROOT ) << " [ROOT]";

  // if allowing for grace, check a little bit further in. 
  // distance, in 1% increments
  if( grace ) {
    double grace_modifier = 1.0;
    double t_original = t_param;
    bool inside_bbox = true;
    while( pathString.find( fTopVolume.c_str() ) == string::npos &&
	   grace_modifier > 0.0 && inside_bbox ) {
      grace_modifier -= 0.01;
      t_param = t_original * grace_modifier;
      fZeroPoint.SetXYZ( fOriginPoint.X() + t_param * fMomentum.X(),
			 fOriginPoint.Y() + t_param * fMomentum.Y(),
			 fOriginPoint.Z() + t_param * fMomentum.Z() );
      fZeroPointNEAR = VolumeSeeker::RotateToNear( fZeroPoint );
      fZeroPointNEAR = VolumeSeeker::TranslateToNear( fZeroPointNEAR );
      fZeroPointROOT = (fZeroPoint - fTopVolumeOrigin) * fToROOTUnits; // subtract translation subtlety

      inside_bbox = ( ( ( std::abs( fZeroPointROOT.Z() - fOzROOT ) < fLzROOT ) &&
			dev_vec.Z() != 0.0 ) ||
		      ( ( std::abs( fZeroPointROOT.Y() - fOyROOT ) < fLyROOT ) &&
			dev_vec.Y() != 0.0 ) ||
		      ( ( std::abs( fZeroPointROOT.X() - fOxROOT ) < fLxROOT ) &&
			dev_vec.X() != 0.0 ) );

      //LOG( "ExoticLLP", pDEBUG ) << "Checking point " << utils::print::Vec3AsString( &fZeroPointROOT ) << " [ROOT] with grace modifier = " << grace_modifier << ", inside_bbox = " << (int) inside_bbox;

      pathString = VolumeSeeker::CheckGeomPoint( fZeroPointROOT );
    }
  }

  //LOG( "ExoticLLP", pDEBUG ) << "Here is the pathString: " << pathString;
  //LOG( "ExoticLLP", pDEBUG ) << "Starting to search for intersections...";

  if( pathString.find( fTopVolume.c_str() ) == string::npos ) return false; // No luck.

  // we are inside the top volume. We want to exit it twice, once going backwards and once forwards.
  // The track enters in the former point and exits in the latter.

  // -- Entry point
  fGeoManager->SetCurrentPoint( fZeroPointROOT.X(), fZeroPointROOT.Y(), fZeroPointROOT.Z() );
  fGeoManager->SetCurrentDirection( -fMomentum.X(), -fMomentum.Y(), -fMomentum.Z() );

  TVector3 curr_point( (fGeoManager->GetCurrentPoint())[0],
		       (fGeoManager->GetCurrentPoint())[1], 
		       (fGeoManager->GetCurrentPoint())[2] );
  pathString = VolumeSeeker::CheckGeomPoint( curr_point );

  fEntryPoint = fZeroPoint;
  fEntryPointROOT = fZeroPointROOT;
  const double smallStep = 0.1; // cm --> have mm precision

  while( pathString.find( fTopVolume.c_str() ) != string::npos ){ // still in the top_volume
    double newX = fEntryPointROOT.X() + (fGeoManager->GetCurrentDirection())[0] * smallStep;
    double newY = fEntryPointROOT.Y() + (fGeoManager->GetCurrentDirection())[1] * smallStep;
    double newZ = fEntryPointROOT.Z() + (fGeoManager->GetCurrentDirection())[2] * smallStep;

    curr_point.SetXYZ( newX, newY, newZ );
    pathString = VolumeSeeker::CheckGeomPoint( curr_point );
    if( pathString.find( fTopVolume.c_str() ) != string::npos ) {
      fEntryPointROOT.SetXYZ( newX, newY, newZ );
      fEntryPoint = fEntryPointROOT * fToLUnits;
    }
  } // exit out

  // Now we have to return this to the USER system...
  fEntryPointROOT.SetXYZ( (fGeoManager->GetCurrentPoint())[0],
			  (fGeoManager->GetCurrentPoint())[1],
			  (fGeoManager->GetCurrentPoint())[2] ); 
  //fEntryPointROOT += fTopVolumeOriginROOT;
  // And save the member variables
  fEntryPoint = fEntryPointROOT * fToLUnits;
  fEntryPoint += fTopVolumeOrigin;
  fEntryPointNEAR = VolumeSeeker::RotateToNear( fEntryPoint );
  fEntryPointNEAR = VolumeSeeker::TranslateToNear( fEntryPointNEAR );
    
  // -- Exit point
  fGeoManager->SetCurrentPoint( fZeroPointROOT.X(), fZeroPointROOT.Y(), fZeroPointROOT.Z() );
  fGeoManager->SetCurrentDirection( fMomentum.X(), fMomentum.Y(), fMomentum.Z() );

  curr_point.SetXYZ( (fGeoManager->GetCurrentPoint())[0],
		     (fGeoManager->GetCurrentPoint())[1], 
		     (fGeoManager->GetCurrentPoint())[2] );
  pathString = VolumeSeeker::CheckGeomPoint( curr_point );

  fExitPoint = fZeroPoint;
  fExitPointROOT = fZeroPointROOT;

  while( pathString.find( fTopVolume.c_str() ) != string::npos ){ // still in the top_volume
    double newX = fExitPointROOT.X() + (fGeoManager->GetCurrentDirection())[0] * smallStep;
    double newY = fExitPointROOT.Y() + (fGeoManager->GetCurrentDirection())[1] * smallStep;
    double newZ = fExitPointROOT.Z() + (fGeoManager->GetCurrentDirection())[2] * smallStep;

    curr_point.SetXYZ( newX, newY, newZ );
    pathString = VolumeSeeker::CheckGeomPoint( curr_point );
    if( pathString.find( fTopVolume.c_str() ) != string::npos ) {
      fExitPointROOT.SetXYZ( newX, newY, newZ );
      fExitPoint = fExitPointROOT * fToLUnits;
    }
  } // exit out

  // Now we have to return this to the USER system...
  fExitPointROOT.SetXYZ( (fGeoManager->GetCurrentPoint())[0],
			 (fGeoManager->GetCurrentPoint())[1],
			 (fGeoManager->GetCurrentPoint())[2] ); 
  //fExitPointROOT += fTopVolumeOriginROOT;
  // And save the member variables
  fExitPoint = fExitPointROOT * fToLUnits;
  fExitPoint += fTopVolumeOrigin;
  fExitPointNEAR = VolumeSeeker::RotateToNear( fExitPoint );
  fExitPointNEAR = VolumeSeeker::TranslateToNear( fExitPointNEAR );

  LOG( "ExoticLLP", pDEBUG ) << "Entry: " << utils::print::Vec3AsString( &fEntryPoint )
			     << " - Exit: " << utils::print::Vec3AsString( &fExitPoint );
    
  return true;
}
//____________________________________________________________________________
AngularRegion VolumeSeeker::AngularAcceptance() const
{
  AngularRegion alpha;

  // Bookkeep the origin point and momentum, just in case
  const TVector3 booked_origin_point = fOriginPoint;
  const TVector3 booked_origin_point_ROOT = fOriginPointROOT;
  const TVector3 booked_momentum = fMomentum;

  // First, define the looking direction
  const TVector3 axis = fMomentum.Unit();

  // We make a potentially strong assumption here, that the top_volume is simply connected.
  // The calculation will be garbage if not, and we'll crash out rather than give garbage.
  assert( VolumeSeeker::RaytraceDetector() && "The origin point of the top_volume lies inside the top volume" );

  // Get the separation between top volume origin and start point
  const TVector3 seed_vector = fTopVolumeOrigin - fOriginPoint;

  // Calculate the (theta, phi) of that point
  // ROOT defines theta in [0, pi] and phi in [-pi, pi]
  double seed_theta = seed_vector.Theta();
  double seed_phi   = seed_vector.Phi();

  // Rasterise on phi : calculate the range of theta that gives raytraces
  double thetaMax = 0.0; double thetaMin = 0.0;
  // Starting coarsely, 
  double deflection = 0.0; // radians
  VolumeSeeker::Deflect( seed_theta, seed_phi, deflection, true );
  thetaMax = seed_theta + deflection;
  deflection = 0.0;
  VolumeSeeker::Deflect( seed_theta, seed_phi, deflection, false );
  thetaMin = seed_theta + deflection;

  // Add this to the angular region
  Point min_point = std::pair< double, double >( thetaMin, seed_phi );
  Point max_point = std::pair< double, double >( thetaMax, seed_phi );
  PointRaster seed_raster = std::pair< Point, Point >( min_point, max_point );
  alpha.emplace_back( seed_raster );

  LOG( "ExoticLLP", pDEBUG )
    << "Deflections: (th0, ph0) = ( " << seed_theta * 180.0 / constants::kPi 
    << ", " << seed_phi * 180.0 / constants::kPi
    << " ) -- theta_min, max = " << thetaMin * 180.0 / constants::kPi
    << ", " << thetaMax * 180.0 / constants::kPi << " [deg]";

  // just in case, restore the original member variables
  fOriginPoint = booked_origin_point;
  fOriginPointROOT = booked_origin_point_ROOT;
  fMomentum = booked_momentum;
  
  return alpha;
}
//____________________________________________________________________________
void VolumeSeeker::Deflect( double th0, double ph0, 
			    double & deflection, bool goUp ) const
{
  double delta = (goUp) ? m_coarse_theta_deflection * constants::kPi / 180.0 : 
    -1 * m_coarse_theta_deflection * constants::kPi / 180.0;
  double new_theta = th0;
  while( VolumeSeeker::RaytraceDetector( true ) && new_theta >= 0.0 && new_theta <= constants::kPi ){
    deflection += delta;
    LOG( "ExoticLLP", pDEBUG ) << "Trying deflection " 
			       << deflection * 180.0 / constants::kPi << " deg...";
    new_theta = th0 + deflection;
    TVector3 new_momentum( std::sin( new_theta ) * std::cos( ph0 ),
			   std::sin( new_theta ) * std::sin( ph0 ), std::cos( new_theta ) );
    fMomentum = new_momentum;
  } // coarse loop
  deflection -= delta;
  new_theta = th0 + deflection;
  fMomentum.SetXYZ( std::sin( new_theta ) * std::cos( ph0 ),
		    std::sin( new_theta ) * std::sin( ph0 ), std::cos( new_theta ) );
  double epsilon = (goUp) ? m_fine_theta_deflection * constants::kPi / 180.0 : 
    -1 * m_fine_theta_deflection * constants::kPi / 180.0;
  while( VolumeSeeker::RaytraceDetector( true ) && new_theta >= 0.0 && new_theta <= constants::kPi ){
    deflection += epsilon;
    LOG( "ExoticLLP", pDEBUG ) << "Trying deflection " << deflection * 180.0 / constants::kPi << " deg...";
    new_theta = th0 + deflection;
    TVector3 new_momentum( std::sin( new_theta ) * std::cos( ph0 ),
			   std::sin( new_theta ) * std::sin( ph0 ), std::cos( new_theta ) );
    fMomentum = new_momentum;
  } // fine loop
  deflection -= epsilon;
  new_theta = th0 + deflection;
  fMomentum.SetXYZ( std::sin( new_theta ) * std::cos( ph0 ),
		    std::sin( new_theta ) * std::sin( ph0 ), std::cos( new_theta ) );

  //LOG( "ExoticLLP", pDEBUG ) << "DONE with final deflection " << deflection * 180.0 / constants::kPi
  //			     << " deg.";
}
//____________________________________________________________________________

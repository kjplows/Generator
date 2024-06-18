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
  LOG( "ExoticLLP", pINFO )
    << "Successfully imported geometry";

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
  //fMomentumNEAR = momentum.Unit();
  fMomentumNEAR = momentum;

  // The momentum just needs to be rotated to USER
  fMomentum = VolumeSeeker::RotateToUser( fMomentumNEAR );
  fAxis = fMomentum.Unit();

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

  LOG( "ExoticLLP", pDEBUG ) << "Successfully imported bounding box with dimensions "
			     << fLx << " x " << fLy << " x " << fLz << " [ " << fLunitString << " ]";
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

  LOG( "ExoticLLP", pDEBUG )
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
  __attribute__((unused)) TGeoVolume * vol = fGeoManager->GetVolume( fTopVolume.c_str() );
  __attribute__((unused)) TGeoNode * node = fGeoManager->FindNode( point[0], point[1], point[2] );
  fGeoManager->MasterToLocal( point, local ); // don't know why but just do it
  return fGeoManager->GetPath();
}
//____________________________________________________________________________
bool VolumeSeeker::RaytraceDetector( bool grace ) const
{
  // Our point starts out at fOriginPoint and has directional cosines fMomentum

  // Have to find the appropriate deviation vector that will get you onto the correct plane
  // that is orthogonal to fAxis. 
  // Important subtlety! The ROOT coordinate system for a top_volume does not know about the transformation matrix -- take care of this at the beginning
  double t_param = 0.0; TVector3 dev_vec = fTopVolumeOrigin - fOriginPoint;
  
  if( dev_vec.Mag() > 0.0 ) {
    /*
      Moving along fMomentum, find the plane which is orthogonal to fAxis and intersects fTopVolumeOrigin

      The line is parametrised as (x, y, z)(t) = (x0, y0, z0) + t * (px, py, pz)
      with (x0, y0, z0) == fOriginPoint and (px, py, pz) == fMomentum

      The plane is parametrised as nx * (x-x1) + ny * (y-y1) + nz * (z-z1) = 0
      with (x1, y1, z1) == fTopVolumeOrigin and (nx, ny, nz) == fAxis

      Substitute in the parametrisations to get
      
      t = ( fAxis Dot dev_vec ) / ( fAxis Dot fMomentum )
      The denominator is just fMomentum.Mag() by fAxis's construction
     */

    t_param = fAxis.Dot( dev_vec ) / fAxis.Dot( fMomentum ); // guaranteed nonzero
    /*
    double delta = dev_vec.Z(); // dz
    double mom_cos = fMomentum.Z(); // pz,hat
    if( delta == 0.0 ){ delta = dev_vec.Y(); mom_cos = fMomentum.Y(); } // make it dy instead
    if( delta == 0.0 ){ delta = dev_vec.X(); mom_cos = fMomentum.X(); } // OK, it *has* to be dx

    t_param = delta / mom_cos;
    */
  } // get to z = 0 (or y = 0, or x = 0) plane

  // transport to starting point (hopefully close enough to the top volume to be meaningful)
  fZeroPoint.SetXYZ( fOriginPoint.X() + t_param * fMomentum.X(),
		     fOriginPoint.Y() + t_param * fMomentum.Y(),
		     fOriginPoint.Z() + t_param * fMomentum.Z() );
  fZeroPointNEAR = VolumeSeeker::RotateToNear( fZeroPoint );
  fZeroPointNEAR = VolumeSeeker::TranslateToNear( fZeroPointNEAR );
  fZeroPointROOT = (fZeroPoint - fTopVolumeOrigin) * fToROOTUnits; // subtract translation subtlety

  //LOG( "ExoticLLP", pDEBUG ) << "Checking point " << utils::print::Vec3AsString( &fZeroPointROOT ) << " [ROOT]";

  // check that this point lies in the geometry.
  std::string pathString = VolumeSeeker::CheckGeomPoint( fZeroPointROOT );

  // if allowing for grace, check a little bit further in along the fMomentum direction.
  // in increments of the bounding box diagonal
  if( grace ) {
    const double maximum_grace = 1.0;
    double grace_modifier = maximum_grace;

    //double upper_bound = std::sqrt( fLx*fLx + fLy*fLy + fLz*fLz );
    double upper_bound = std::max( fLx, std::max( fLy, fLz ) );
    bool inside_bbox = true;

    while( pathString.find( fTopVolume.c_str() ) == string::npos &&
	   grace_modifier > 0.0 && inside_bbox ) {

      grace_modifier -= m_grace_decrement;
      double axis_mod = (1.0 - grace_modifier) * upper_bound;

      fZeroPoint -= axis_mod * fMomentum.Unit();
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

  //LOG( "ExoticLLP", pDEBUG ) << "Entry: " << utils::print::Vec3AsString( &fEntryPoint )
  //		     << " - Exit: " << utils::print::Vec3AsString( &fExitPoint );
    
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

  // Get the separation between top volume origin and start point
  const TVector3 seed_vector = fTopVolumeOrigin - fOriginPoint;

  // We make a potentially strong assumption here, that the top_volume is simply connected.
  // The calculation will be garbage if not, and we'll crash out rather than give garbage.
  fMomentum = seed_vector.Unit();
  assert( VolumeSeeker::RaytraceDetector() && "The origin point of the top_volume lies inside the top volume" );
  fMomentum = booked_momentum;

  // Of course, the angles are defined with respect to the momentum axis...
  // We need to first project seed_vector onto the theta = 0 plane
  // where the theta = 0 plane is defined by fMomentum and the tranverse component of seed_vector
  TVector3 transverse_seed = seed_vector - fAxis.Dot( seed_vector ) * fAxis;
  if( transverse_seed.Mag() > 1.0e-3 * seed_vector.Mag() ) fThetaAxis = transverse_seed.Unit();
  else{ // just pick a convenient direction
    /*
    if( fAxis.Y() == 0.0 ) fThetaAxis = TVector3( 0.0, 1.0, 0.0 );
    else if( fAxis.X() == 0.0 ) fThetaAxis = TVector3( 1.0, 0.0, 0.0 );
    else if( fAxis.Z() == 0.0 ) fThetaAxis = TVector3( 0.0, 0.0, 1.0 );
    else fThetaAxis = TVector3( fAxis.Y(), -fAxis.X(), 0.0 );
    */
    if( fAxis.X() == 0.0 ) fThetaAxis = TVector3( 0.0, fAxis.Z(), -fAxis.Y() );
    else if( fAxis.Y() == 0.0 ) fThetaAxis = TVector3( fAxis.Z(), 0.0, fAxis.X() );
    else if( fAxis.Z() == 0.0 ) fThetaAxis = TVector3( fAxis.Y(), -fAxis.X(), 0.0 );
    else fThetaAxis = TVector3( fAxis.Y() * fAxis.Z(),
				-0.5 * fAxis.X() * fAxis.Z(),
				-0.5 * fAxis.X() * fAxis.Y() );
  }

  // Now that lets us complete the "momentum" system of coordinates
  fPhiAxis = fAxis.Cross( fThetaAxis );

  LOG( "ExoticLLP", pDEBUG )
    << "\nfAxis      = " << utils::print::Vec3AsString(&fAxis)
    << "\nfThetaAxis = " << utils::print::Vec3AsString(&fThetaAxis)
    << "\nfPhiAxis   = " << utils::print::Vec3AsString(&fPhiAxis);

  // Get the projection of that vector onto the appropriate coordinate system
  const double sep_from_theta = seed_vector.Dot( fPhiAxis );
  const double sep_from_phi   = seed_vector.Dot( fThetaAxis );

  TVector3 projection_theta = seed_vector - sep_from_theta * fPhiAxis;
  TVector3 projection_phi = seed_vector - sep_from_phi * fThetaAxis;
  
  // acos runs from [0, pi]
  /*
  double seed_theta = std::acos( fAxis.Dot( projection_theta ) /
				 std::max( projection_theta.Mag(), 1.0e-10 ) );
  double seed_phi = std::acos( fAxis.Dot( projection_phi ) / 
			       std::max( projection_phi.Mag(), 1.0e-10 ) );
  if( fAxis.Dot( projection_phi ) < 0.0 ) seed_phi *= -1.0;

  if( transverse_seed.Mag() == 0.0 ) { seed_theta = 0.0; seed_phi = 0.0; } // no separation!
  */
  
  // check that Rasterise() does what you want it to
  // RETHERE first make the seed vector and define fThetaSeed, fPhiSeed
  VolumeSeeker::Rasterise( alpha, true );
  VolumeSeeker::Rasterise( alpha, false );

  // sort the angular region in increasing phi
  std::sort( alpha.begin(), alpha.end(), 
	     []( const PointRaster & a, const PointRaster &b ){
	       // first and second points of a PointRaster have same phi
	       Point pta = a.first, ptb = b.first;
	       return pta.second < ptb.second;
	     } );

  LOG( "ExoticLLP", pDEBUG ) << "Angular region has " << alpha.size() << " rasters.";

  // Now we'd like to construct the actual angles on the unit sphere these deflections
  // correspond to. Let's do it!
  //VolumeSeeker::ConvertToUserAngles( booked_momentum, alpha, beta );

  // just in case, restore the original member variables
  fOriginPoint = booked_origin_point;
  fOriginPointROOT = booked_origin_point_ROOT;
  fMomentum = booked_momentum;
  
  //return std::make_tuple( alpha, beta );
  return alpha;
}
//____________________________________________________________________________
double VolumeSeeker::AngularSize( AngularRegion alpha ) const
{
  // The AngularRegion object is a collection of PointRasters, each of which has two points
  // Fundamentally, we split the region into an upper and a lower bound.
  // Get the size of each and return their difference

  std::vector< Point > upper_points, lower_points;
  for( AngularRegion::iterator ait = alpha.begin() ; ait != alpha.end() ; ++ait ) {
    upper_points.emplace_back( (*ait).second );
    lower_points.emplace_back( (*ait).first );
  }

  double size = VolumeSeeker::Trapezoid( upper_points, lower_points );

  return size;
}
//____________________________________________________________________________
double VolumeSeeker::Trapezoid( std::vector<Point> up_vec, std::vector<Point> dn_vec ) const
{
  double total_up = 0.0, total_dn = 0.0;
  double small_height = 0.0, large_height = 0.0;

  Point previous_point_up = *(up_vec.begin());
  Point previous_point_dn = *(dn_vec.begin());
  for( std::vector<Point>::iterator ptit = up_vec.begin() ; ptit != up_vec.end() ; ++ptit ) {

    if( *ptit == previous_point_up ) continue; // skip the first point

    Point current_point_up = *ptit;
    int ttt = ptit - up_vec.begin();
    Point current_point_dn = *( dn_vec.begin() + ttt );

    // check if we have the same phi.. if we do, then bail and move on
    // But always add the last element
    if( current_point_up.second == previous_point_up.second &&
	ptit != up_vec.end() - 1) continue;
    
    double width = std::abs(current_point_up.second - previous_point_up.second);

    /*
     * For each of the two phis, find the upper and the lower theta.
     * Either upper * lower < 0, in which case a deflection of zero would be accepted
     * Or upper * lower > 0, in which case it would not.
     * In the first case, height = hUpper + hLower, h = 1 - cos
     * In the second case, height = hLong - hShort, long,short = max,min(abs(theta))
     
     * Actually there are three cases! Let U1, D1, U2, D2 be the upper, lower points of prev, current.
     * U1D1 < 0 && U2D2 < 0: Zero-deflection OK throughout the raster
     * U1D1 > 0 && U2D2 > 0: Zero-deflection not OK throughout the raster
     * ( U1D1 < 0 && U2D2 > 0 ) || ( U1D1 > 0 && U2D2 < 0 ): Zero-deflection changes through the raster.
     * In that case, you have to split the width....

     * To be even worse, you could have a case which looks like case 2, but D1<U1<0 and 0<D2<U2.
     * I can't think of anyone who'd build a detector that looks like that.
     */

    double theta_up_prev = previous_point_up.first;
    double theta_dn_prev = previous_point_dn.first;
    double theta_up_curr = current_point_up.first;
    double theta_dn_curr = current_point_dn.first;

    double cand_U1 = 1.0 - std::cos( theta_up_prev );
    double cand_D1 = 1.0 - std::cos( theta_dn_prev );
    double cand_U2 = 1.0 - std::cos( theta_up_curr );
    double cand_D2 = 1.0 - std::cos( theta_dn_curr );

    /*
     * CASE 1: Zero-deflection OK for both rasters.
     */

    if( theta_up_prev * theta_dn_prev <= 0 && theta_up_curr * theta_dn_curr <= 0 ){
      double h1 = cand_U1 + cand_D1;
      double h2 = cand_U2 + cand_D2;

      total_up += width * (h1+h2)/2.0;
      total_dn += 0.0;
    }

    /*
     * CASE 2: Zero-deflection not OK for either raster.
     */

    if( theta_up_prev * theta_dn_prev > 0 && theta_up_curr * theta_dn_curr > 0 ){
      if( theta_up_prev < 0.0 ) { 
	// D1 < U1 < 0 --> 1-cos(D1) > 1-cos(U1) > 0
	// Also D2 < U2 < 0 as there's no crossing of zero

	total_up += width * (cand_D1 + cand_D2)/2.0;
	total_dn += width * (cand_U1 + cand_U2)/2.0;
      } else if( theta_dn_prev >= 0.0 ) {
	// 0 < D1 < U1 --> 1-cos(U1) > 1-cos(D1) > 0
	// Also 0 < D2 < U2 as there's not crossing of zero

	total_up += width * (cand_U1 + cand_U2)/2.0;
	total_dn += width * (cand_D1 + cand_D2)/2.0;
      }
    }

    /*
     * CASE 3: Zero-deflection OK in one raster, but not in the other. 
     * There is a zero-crossing which we must estimate by interpolation
     * and reduce to a simultaneous evaluation of Case 1 and Case 2.

     * RETHERE need to code this in when / if needed
     */

    /*
    // Get the size of the upper region
    double cand_one = 1.0 - std::cos( previous_point_up.first );
    double cand_two = 1.0 - std::cos( current_point_up.first );
    small_height = std::min( cand_one, cand_two );
    large_height = std::max( cand_one, cand_two );
    // trapezoid area is w * (h + H)/2
    total_up += width * ( small_height + large_height ) / 2.0;

    // Get the size of the lower region
    cand_one = 1.0 - std::cos( previous_point_dn.first );
    cand_two = 1.0 - std::cos( current_point_dn.first );
    small_height, large_height = std::min( cand_one, cand_two ), std::max( cand_one, cand_two );
    // trapezoid area is w * (h + H)/2
    total_dn += width * ( small_height + large_height ) / 2.0;
    */

    previous_point_up = current_point_up; // update
    previous_point_dn = current_point_dn; // update
  }

  LOG( "ExoticLLP", pDEBUG ) << "total_up, total_dn = " << total_up << ", " << total_dn;

  return total_up + total_dn;
}
//____________________________________________________________________________
double VolumeSeeker::Simpson( std::vector<Point> pt_vec ) const
{
  // we have a collection of points == pairs of (theta, phi).
  // Need to calculate an integral that goes like (sin\theta d\theta) d\phi
  // So apply Simpson's rule on the theta part (thanks Fubini) and obtain the answer

  double total = 0.0;
  
  Point previous_point = *(pt_vec.begin());
  for( std::vector<Point>::iterator ptit = pt_vec.begin() ; ptit != pt_vec.end() ; ++ptit ) {

    if( *ptit == previous_point ) continue; // skip the first point
    
    Point current_point = *ptit;
    // Check if there are duplicates in theta, which means there was no meaningful change. Same term
    if( current_point.first == previous_point.first ) continue;

    // Simple width on phi
    double phi_bit = current_point.second - previous_point.second;
    // Simpson's 1/3 rule on theta
    double theta_bit = 1.0 / 6.0;
    
    double prev_sin = 1.0 - std::cos( previous_point.first );
    double curr_sin = 1.0 - std::cos( current_point.first );
    double imed_sin = 1.0 - std::cos( ( previous_point.first + current_point.first ) / 2.0 );

    theta_bit *= ( prev_sin + 4.0 * imed_sin + curr_sin );

    total += phi_bit * theta_bit;

    previous_point = current_point; // update
  } // add terms to the integral

  return total;
}
//____________________________________________________________________________
void VolumeSeeker::Rasterise( AngularRegion & alpha, bool goRight ) const
{
  // obligatory momentum protection
  const TVector3 booked_momentum = fMomentum;

  double sweep = 0.0;
  double thetaMax = 0.0, thetaMin = 0.0;
  double deflection_up = 0.0, deflection_down = 0.0;

  // first, check that the middle point isn't in. If it isn't, insert it.
  if( alpha.size() == 0.0 ){
    // One loop to find the maximum possible angles
    VolumeSeeker::Deflect( deflection_up, true );
    VolumeSeeker::Deflect( deflection_down, false );
    
    //bool zero_okay = ( deflection_up * deflection_down <= 0.0 );

    thetaMax = deflection_up; thetaMin = deflection_down;
    
    //thetaMax = std::max( std::abs(deflection_up), std::abs(deflection_down) );
    //thetaMin = zero_okay ? 0.0 : std::min( std::abs(deflection_up), std::abs(deflection_down) );
    
    LOG( "ExoticLLP", pDEBUG )
      << "Deflections: phi = " << 0.0
      << ", theta_min, max = " << thetaMin * 180.0 / constants::kPi
      << ", " << thetaMax * 180.0 / constants::kPi << " [deg]";
    
    // Add this to the angular region
    Point min_point = std::pair< double, double >( thetaMin, sweep );
    Point max_point = std::pair< double, double >( thetaMax, sweep );
    PointRaster raster = std::pair< Point, Point >( min_point, max_point );
    alpha.emplace_back( raster );
    
    fMomentum = booked_momentum;
  } // first element in alpha

  // This is a deflection on phi. So it follows fPhiAxis
  // Each step is a test on phi, and if there is a raytrace anywhere along that phi, for any theta,
  // then start at that raytace and call Deflect.

  // Actually first get the diagonal of the bbox and then build the coarse deflection
  // basically, want to poll N points in each direction

  const double full_diagonal = std::sqrt( fLx*fLx + fLy*fLy + fLz*fLz );
  const TVector3 seed_vector = fTopVolumeOrigin - fOriginPoint;
  const double baseline = seed_vector.Mag();
  const double base_angle = std::atan( full_diagonal / baseline );

  double delta = (goRight) ? base_angle / m_coarse_phi_deflection :
    -1 * base_angle / m_coarse_phi_deflection;
  while( VolumeSeeker::RaytraceDetector( true ) && std::abs(sweep) <= constants::kPi &&
	 ( deflection_up != deflection_down || sweep == 0.0 ) ) {
    deflection_up = 0.0; deflection_down = 0.0; thetaMax = 0.0; thetaMin = 0.0;
    sweep += delta;

    // Calculate the fPhiAxis component
    double scale = booked_momentum.Mag() * std::tan( sweep );
    fMomentum = booked_momentum + scale * fPhiAxis;

    VolumeSeeker::Deflect( deflection_up, true );
    VolumeSeeker::Deflect( deflection_down, false );
    
    bool zero_okay = ( deflection_up * deflection_down <= 0.0 );
    
    thetaMax = std::max( std::abs(deflection_up), std::abs(deflection_down) );
    thetaMin = zero_okay ? 0.0 : std::min( std::abs(deflection_up), std::abs(deflection_down) );

    /*
    LOG( "ExoticLLP", pDEBUG )
      << "Deflections: phi = " << sweep * 180.0 / constants::kPi
      << ", theta_min, max = " << thetaMin * 180.0 / constants::kPi
      << ", " << thetaMax * 180.0 / constants::kPi << " [deg]";
    */

    // and add to the raster
    if( thetaMin != thetaMax ) {
      Point min_point = std::pair< double, double >( thetaMin, sweep );
      Point max_point = std::pair< double, double >( thetaMax, sweep );
      PointRaster raster = std::pair< Point, Point >( min_point, max_point );
      alpha.emplace_back( raster );
    }
  } // coarse loop

  sweep -= delta;
  double coarse_scale = booked_momentum.Mag() * std::tan( sweep );
  fMomentum = booked_momentum + coarse_scale * fPhiAxis;

  // Also remove the last raster point, gets erroneously added
  alpha.pop_back();

  double epsilon = (goRight) ? base_angle / m_fine_phi_deflection :
    -1 * base_angle / m_fine_phi_deflection;
  while( VolumeSeeker::RaytraceDetector( true ) && std::abs(sweep) <= constants::kPi &&
	 ( deflection_up != deflection_down || sweep == 0.0 ) ) {
    deflection_up = 0.0; deflection_down = 0.0; thetaMax = 0.0; thetaMin = 0.0;
    sweep += epsilon;

    // Calculate the fPhiAxis component
    double scale = booked_momentum.Mag() * std::tan( sweep );
    fMomentum = booked_momentum + scale * fPhiAxis;


    VolumeSeeker::Deflect( deflection_up, true );
    VolumeSeeker::Deflect( deflection_down, false );
    
    bool zero_okay = ( deflection_up * deflection_down <= 0.0 );
    
    thetaMax = std::max( std::abs(deflection_up), std::abs(deflection_down) );
    thetaMin = zero_okay ? 0.0 : std::min( std::abs(deflection_up), std::abs(deflection_down) );

    /*
    LOG( "ExoticLLP", pDEBUG )
      << "Deflections: phi = " << sweep * 180.0 / constants::kPi
      << ", theta_min, max = " << thetaMin * 180.0 / constants::kPi
      << ", " << thetaMax * 180.0 / constants::kPi << " [deg]";
    */

    // and add to the raster
    if( thetaMin != thetaMax ) {
      Point min_point = std::pair< double, double >( thetaMin, sweep );
      Point max_point = std::pair< double, double >( thetaMax, sweep );
      PointRaster raster = std::pair< Point, Point >( min_point, max_point );
      alpha.emplace_back( raster );
    }
  } // coarse loop
  sweep -= epsilon;
  
  // all done, reset
  fMomentum = booked_momentum;
}
//____________________________________________________________________________
void VolumeSeeker::Deflect( double & deflection, bool goUp ) const
{
  // bookkeep the momentum that we had
  const TVector3 booked_momentum = fMomentum;

  // This is a deflection on theta. So it follows fThetaAxis
  // Each step sets a deflection angle and adds the appropriately scaled 
  // vector on fThetaAxis to see if there is a raytrace there

  const double full_diagonal = std::sqrt( fLx*fLx + fLy*fLy + fLz*fLz );
  const TVector3 seed_vector = fTopVolumeOrigin - fOriginPoint;
  const double baseline = seed_vector.Mag();
  const double base_angle = std::atan( full_diagonal / baseline );
  
  double delta = (goUp) ? base_angle / m_coarse_theta_deflection : 
    -1 * base_angle / m_coarse_theta_deflection;
  while( VolumeSeeker::RaytraceDetector( true ) && std::abs(deflection) <= constants::kPi ){
    deflection += delta;
    // Now calculate the fThetaAxis component
    double scale = booked_momentum.Mag() * std::tan( deflection );
    fMomentum = booked_momentum + scale * fThetaAxis;
  } // coarse loop
  deflection -= delta;
  double coarse_scale = booked_momentum.Mag() * std::tan( deflection );
  fMomentum = booked_momentum + coarse_scale * fThetaAxis;

  double epsilon = (goUp) ? base_angle / m_fine_theta_deflection : 
    -1 * base_angle / m_fine_theta_deflection;
  while( VolumeSeeker::RaytraceDetector( true ) && std::abs(deflection) <= constants::kPi ){
    deflection += epsilon;
    double scale = booked_momentum.Mag() * std::tan( deflection );
    fMomentum = booked_momentum + scale * fThetaAxis;
  } // fine loop
  deflection -= epsilon;

  fMomentum = booked_momentum;
}
//____________________________________________________________________________
void VolumeSeeker::ConvertToUserAngles( const TVector3 booked_momentum, AngularRegion deflection_region,
					AngularRegion & cosines_region ) const
{
  // Need to calculate the actual angles...
  // Need to decide also which deflections encompass the shape
  
  /* We will check four points, in this order:
   * 1) thetaMax, along + fThetaAxis
   * 2) thetaMax, along - fThetaAxis
   * If both points in detector, then that's the shape. Else:
   * 3) thetaMin, along + fThetaAxis
   * 4) thetaMin, along - fThetaAxis
   * If checked thetaMin, either the ++ or the -- combination will be the correct one.
   */

  AngularRegion::iterator ait = deflection_region.begin();
  
  for( ; ait != deflection_region.end(); ++ait ){

    // obtain the phi and thetaMin, thetaMax from the deflection regions
    PointRaster rst = *ait;
    double phi = (rst.first).second; // a raster has constant phi
    double thetaMin = (rst.first).first;
    double thetaMax = (rst.second).first;

    // we need to calculate the new PointRaster to use.

    //double extracted_phi = mod_p.Phi(); Nope. Wrong formulation.
    /*
      The phi formulation (or "sweep") refers to the deflection along the plane
      spanned by (fAxis, fPhiAxis). It is the angle between those two projections.
      So, calculate the angles directly.

      Using (x, y, z) = r * (sin\theta cos\phi, sin\theta sin\phi, cos\theta)
     */

    double bpx = booked_momentum.X(), bpy = booked_momentum.Y(), bpz = booked_momentum.Z();
    double bp3 = booked_momentum.Mag();

    double base_phi = std::acos( bpx / bp3 ); 
    double base_theta = std::acos( bpz / bp3 );
    if( bpy < 0.0 ) base_phi = 2.0 * constants::kPi - base_phi;

    double extracted_phi = base_phi + phi;
    double extracted_thetaMin = base_theta + thetaMin;
    double extracted_thetaMax = base_theta + thetaMax;

    /*
    LOG( "ExoticLLP", pDEBUG )
      << "\nbooked_momentum  = " << utils::print::Vec3AsString( &booked_momentum )
      << "\nbase_phi         = " << base_phi * 180.0 / constants::kPi
      << "\n==> got phi      = " << extracted_phi * 180.0 / constants::kPi
      << "\nbase_theta       = " << base_theta * 180.0 / constants::kPi
      << "\n==> got thetaMin = " << extracted_thetaMin * 180.0 / constants::kPi
      << "\n==> got thetaMax = " << extracted_thetaMax * 180.0 / constants::kPi;
    */

    // Make the new points and add to the new raster
    Point min_point = std::pair< double, double >( extracted_thetaMin, extracted_phi );
    Point max_point = std::pair< double, double >( extracted_thetaMax, extracted_phi );
    PointRaster raster = std::pair< Point, Point >( min_point, max_point );
    cosines_region.emplace_back( raster );

    fMomentum = booked_momentum;
  } // loop over all the elements in the deflection region
}

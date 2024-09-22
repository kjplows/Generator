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
Vertex::Vertex(): fIndex(0), fX(0.0), fY(0.0), fZ(0.0)
{

}
//____________________________________________________________________________
Vertex::Vertex( int i, double x, double y, double z ) : fIndex(i), fX(x), fY(y), fZ(z)
{

}
//____________________________________________________________________________
Vertex::Vertex( int i, const TVector3 & v ) : fIndex(i), fX(v.X()), fY(v.Y()), fZ(v.Z())
{

}
//____________________________________________________________________________
Vertex::Vertex( int i, const Vertex & v ) : fIndex(i), fX(v.fX), fY(v.fY), fZ(v.fZ)
{

}
//____________________________________________________________________________
Vertex::~Vertex()
{
  fIndex = 0; fX = 0.0; fY = 0.0; fZ = 0.0;
}
//____________________________________________________________________________
double Vertex::Dist( const Vertex v ) const
{
  double dx = fX - v.fX;
  double dy = fY - v.fY;
  double dz = fZ - v.fZ;
  return std::sqrt( dx*dx + dy*dy + dz*dz );
}
//____________________________________________________________________________
TVector3 Vertex::Displacement( const Vertex v ) const
{
  double dx = v.fX - fX;
  double dy = v.fY - fY;
  double dz = v.fZ - fZ;
  return TVector3( dx, dy, dz );
}
//____________________________________________________________________________
double Vertex::Intersection( const Vertex v, const TVector3 unit, const double d ) const
{
  double lambda = -9.9;
  
  TVector3 path = this->Displacement( v );
  double denom = unit.Dot( path );
  if( denom == 0.0 ) return lambda; // no intersection as edge is coplanar with viewing plane

  double numer = d + unit.Dot( this->Displacement( Vertex() ) ); 
  // same as d - unit.Dot( (0,0,0) - *this )
  
  return numer / denom;
}
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
  fIsGeomFileSet = false;

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
void VolumeSeeker::SetOffset( double x, double y, double z ) const
{
  fTopVolumeOffset.SetXYZ( x, y, z );
}
//____________________________________________________________________________
void VolumeSeeker::AdoptControls( bool use_saa, bool use_cmv,
				  double ct, double cp, double ft, double fp, double gr ) const
{
  m_use_saa = use_saa;
  m_use_cmv = use_cmv;
  m_coarse_theta_deflection = ct;
  m_coarse_phi_deflection = cp;
  m_fine_theta_deflection = ft;
  m_fine_phi_deflection = fp;
  m_grace_decrement = gr;
}
//____________________________________________________________________________
void VolumeSeeker::SetGeomFile( std::string geomfile, std::string topVolume ) const
{
  if( fIsGeomFileSet ){ LOG("ExoticLLP", pNOTICE) << "Geom file already set, doing nothing..."; return; }

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

  fIsGeomFileSet = true;
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
  /*
  LOG( "ExoticLLP", pDEBUG ) 
    << "\nInput vec = " << utils::print::Vec3AsString( &input )
    << "\nTranslating by " << utils::print::Vec3AsString( &tr_vec ) << " backwards...";
  */
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

  /*
  LOG( "ExoticLLP", pDEBUG ) << "\norigin_point is " << utils::print::Vec3AsString( &origin_point )
			     << "\nin USER it is   " << utils::print::Vec3AsString( &fOriginPoint )
			     << "\nin ROOT it is   " << utils::print::Vec3AsString( &fOriginPointROOT );
  */

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
    //LOG( "ExoticLLP", pDEBUG ) << "Node with name " << path << " has " << nDaughters << " daughters...";
    for( int iDaughter = 0; iDaughter < nDaughters; iDaughter++ ){
      TGeoNode * dNode = node->GetDaughter(iDaughter);
      assert( dNode && "Daughter node not null" );
      std::string dPath( path );
      dPath.append( "/" ); dPath.append( dNode->GetName() );
      TGeoMatrix * nodeMat = dNode->GetMatrix();

      /*
      LOG( "ExoticLLP", pDEBUG ) << "Got node, path, and matrix for daughter node "
				 << iDaughter << " / " << nDaughters-1 << "...";
      */

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

  /*
  LOG( "ExoticLLP", pINFO )
    << "Found the target volume! Here is its path and full matrix:"
    << "\nPath: " << paths.back()
    << "\nTranslations: ( " << final_tra[0] << ", " << final_tra[1] << ", " << final_tra[2]
    << " ) [cm]"
    << "\nRotation matrix: ( ( " 
    << final_rot[0] << ", " << final_rot[1] << ", " << final_rot[2] << " ), ( "
    << final_rot[3] << ", " << final_rot[4] << ", " << final_rot[5] << " ), ( "
    << final_rot[6] << ", " << final_rot[7] << ", " << final_rot[8] << " ) )";
  */

  // Also set the member variables at this stage
  //fTopVolumeOriginROOT.SetXYZ( final_tra[0], final_tra[1], final_tra[2] );
  //fTopVolumeOrigin = fTopVolumeOriginROOT * fToLUnits;
  fTopVolumeOffset.SetXYZ( final_tra[0] * fToLUnits, 
			   final_tra[1] * fToLUnits, 
			   final_tra[2] * fToLUnits );
  fTopVolumeOrigin = fTopVolumeOffset;
  fTopVolumeOriginNEAR = VolumeSeeker::RotateToNear( fTopVolumeOrigin );
  fTopVolumeOriginNEAR = VolumeSeeker::TranslateToNear( fTopVolumeOriginNEAR );

  return final_mat;
}
//____________________________________________________________________________
TVector3 VolumeSeeker::GetRandomPointInTopVol() const
{
  // This returns in USER coordinates.
  TVector3 outVec(0.0, 0.0, 0.0);
  std::string pathString = "";
  
  // Just in case... sometimes ROOT's IO is funky with this
  std::cout.flush(); std::cerr.flush();
  
  // sample from the bounding box directly
  int nTries = 0; int nMaxTries = 100;
  while( true ) {
    
    RandomGen * rnd = RandomGen::Instance();
    double vxROOT = rnd->RndGen().Uniform( -fLxROOT, fLxROOT );
    double vyROOT = rnd->RndGen().Uniform( -fLyROOT, fLyROOT );
    double vzROOT = rnd->RndGen().Uniform( -fLzROOT, fLzROOT );
    
    outVec.SetXYZ( vxROOT, vyROOT, vzROOT );
    pathString = VolumeSeeker::CheckGeomPoint( outVec );

    if( pathString.find( fTopVolume.c_str() ) != string::npos || nTries > nMaxTries )
      break;
  }
  if( nTries > nMaxTries )
    LOG( "ExoticLLP", pWARN ) << "Could not get a point inside the top volume after "
			      << nMaxTries << " tries. Continuing with last point.";

  outVec.SetXYZ( outVec.X() * fToLUnits + fTopVolumeOrigin.X(), 
		 outVec.Y() * fToLUnits + fTopVolumeOrigin.Y(), 
		 outVec.Z() * fToLUnits + fTopVolumeOrigin.Z() );

  return outVec;
}
//____________________________________________________________________________
TVector3 VolumeSeeker::GetRandomPointInTopVolNEAR() const
{
  // This returns in NEAR coordinates
  TVector3 outVec = VolumeSeeker::GetRandomPointInTopVol();
  //outVec += fTopVolumeOffset; // adds the position of USER origin in NEAR frame
  outVec = VolumeSeeker::RotateToNear( outVec );
  outVec = VolumeSeeker::TranslateToNear( outVec );
  return outVec;
}
//____________________________________________________________________________
bool VolumeSeeker::IsInTop( TVector3 point, bool near ) const
{
  TVector3 chkpoint = point;
  if( near ) { // translate it to USER first
    chkpoint = VolumeSeeker::TranslateToUser( chkpoint );
    chkpoint = VolumeSeeker::RotateToUser( chkpoint );
  }

  // and make into ROOT units...
  chkpoint.SetXYZ( chkpoint.X() * fToROOTUnits, 
		   chkpoint.Y() * fToROOTUnits,
		   chkpoint.Z() * fToROOTUnits );

  return ( VolumeSeeker::CheckGeomPoint(chkpoint) ).find( fTopVolume.c_str() ) != string::npos;
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
void VolumeSeeker::IntersectDetector() const
{
  // Our point starts out at fOriginPoint and has directional cosines fMomentum

  // Have to find the appropriate deviation vector that will get you onto the correct plane
  // that is orthogonal to fAxis. 
  
  //const TVector3 seed_vector = fOriginPoint - fTopVolumeOrigin;
  TVector3 baseline_seed = fMomentum.Unit();

  std::array< TVector3, 3 > aNormalFaces = { TVector3( 1, 0, 0 ),   // +- fLx
					     TVector3( 0, 1, 0 ),   // +- fLy
					     TVector3( 0, 0, 1 ) }; // +- fLz
  std::array< double, 3 > aProductFaces = { baseline_seed.Dot( aNormalFaces[0] ),
					    baseline_seed.Dot( aNormalFaces[1] ),
					    baseline_seed.Dot( aNormalFaces[2] ) };
  std::array< double, 3 > aDimensions = { fLx, fLy, fLz }; // syntactic sugar

  std::vector<TVector3> vIntersections;

  // Check each pair of faces in turn
  for( int iface = 0; iface < 3; iface++ ) {
    if( aProductFaces[iface] == 0.0 ) continue; // no intersection for sure

    TVector3 pos_vec =  aDimensions[iface] * aNormalFaces[iface] + fTopVolumeOffset;
    TVector3 neg_vec = -aDimensions[iface] * aNormalFaces[iface] + fTopVolumeOffset;

    // The line does not go to the origin. It goes to a target location.
    //TVector3 pos_Dvec = fTopVolumeOrigin - pos_vec;
    //TVector3 neg_Dvec = fTopVolumeOrigin - neg_vec;
    TVector3 pos_Dvec = pos_vec - fTarget;
    TVector3 neg_Dvec = neg_vec - fTarget;

    double t_pos = pos_Dvec.Dot( aNormalFaces[iface] ) / aProductFaces[iface];
    double t_neg = neg_Dvec.Dot( aNormalFaces[iface] ) / aProductFaces[iface];

    // evaluate the intersection point at these values
    TVector3 intersection_pos = fTarget + baseline_seed * t_pos;
    TVector3 intersection_neg = fTarget + baseline_seed * t_neg;

    TVector3 dev_pos = intersection_pos - fTopVolumeOrigin;
    TVector3 dev_neg = intersection_neg - fTopVolumeOrigin;

    // if the intersection is not in the face, then scale the vector to be larger than baseline as it is not a valid intersection
    switch( iface ) {
    case 0: // evaluating X, check Y and Z
      if( std::abs( dev_pos.Y() ) <= fLy && std::abs( dev_pos.Z() ) <= fLz ) 
	vIntersections.emplace_back( intersection_pos );
      if( std::abs( dev_neg.Y() ) <= fLy && std::abs( dev_neg.Z() ) <= fLz )
	vIntersections.emplace_back( intersection_neg );
      break;
    case 1: // evaluating Y, check X and Z
      if( std::abs( dev_pos.X() ) <= fLx && std::abs( dev_pos.Z() ) <= fLz ) 
	vIntersections.emplace_back( intersection_pos );
      if( std::abs( dev_neg.X() ) <= fLx && std::abs( dev_neg.Z() ) <= fLz )
	vIntersections.emplace_back( intersection_neg );
      break;
    case 2: // evaluating Z, check X and Y
      if( std::abs( dev_pos.Y() ) <= fLy && std::abs( dev_pos.X() ) <= fLx ) 
	vIntersections.emplace_back( intersection_pos );
      if( std::abs( dev_neg.Y() ) <= fLy && std::abs( dev_neg.X() ) <= fLx )
	vIntersections.emplace_back( intersection_neg );
      break;
    } // argh I could not avoid the switch...
  } // check each pair of faces in turn

  // awesome, now sort the intersections by increasing distance to fOriginPoint
  std::vector< TVector3 > vSortedIntersections;
  while( vIntersections.size() > 0 ) {
    std::vector< TVector3 >::iterator it_int = 
      std::min_element( vIntersections.begin(), vIntersections.end(),
			[this](TVector3 a, TVector3 b){ // capture this to use fOriginPoint
			  return (a-fOriginPoint).Mag() < (b-fOriginPoint).Mag(); } );
    vSortedIntersections.emplace_back( *it_int );
    vIntersections.erase( it_int );
  }

  // and do some accounting for internal variables
  double t_param = 0.0;
  /*
  if( seed_vector.Mag() > 0.0 ) 
    t_param = fAxis.Dot( seed_vector ) / fAxis.Dot( fMomentum ); // guaranteed nonzero
  */

  fZeroPoint.SetXYZ( fOriginPoint.X() + t_param * fMomentum.X(),
		     fOriginPoint.Y() + t_param * fMomentum.Y(),
		     fOriginPoint.Z() + t_param * fMomentum.Z() );
  fZeroPointNEAR = VolumeSeeker::RotateToNear( fZeroPoint );
  fZeroPointNEAR = VolumeSeeker::TranslateToNear( fZeroPointNEAR );
  fZeroPointROOT = (fZeroPoint - fTopVolumeOrigin) * fToROOTUnits; // subtract translation subtlety

  fEntryPoint = vSortedIntersections.at(0);
  fEntryPointROOT = (fEntryPoint - fTopVolumeOrigin) * fToROOTUnits;
  fEntryPointNEAR = VolumeSeeker::RotateToNear( fEntryPoint );
  fEntryPointNEAR = VolumeSeeker::TranslateToNear( fEntryPointNEAR );

  fExitPoint = vSortedIntersections.at(vSortedIntersections.size()-1);
  fExitPointROOT = (fExitPoint - fTopVolumeOrigin) * fToROOTUnits;
  fExitPointNEAR = VolumeSeeker::RotateToNear( fExitPoint );
  fExitPointNEAR = VolumeSeeker::TranslateToNear( fExitPointNEAR );
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

  //LOG( "ExoticLLP", pDEBUG ) << "fMomentum = " << utils::print::Vec3AsString( &fMomentum );
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

  /*
  LOG( "ExoticLLP", pDEBUG ) << "\nbooked_origin_point      = " << utils::print::Vec3AsString( &booked_origin_point )
			     << "\nbooked_origin_point_ROOT = " << utils::print::Vec3AsString( &booked_origin_point_ROOT );
  */

  // First, check if the point is inside the volume. If yes, every emission angle is good!
  std::string original_pathString = this->CheckGeomPoint( booked_origin_point_ROOT );
  //LOG( "ExoticLLP", pDEBUG ) << "original_pathString = " << original_pathString;

  /*
  if( original_pathString.find( fTopVolume.c_str() ) == string::npos )
    LOG("ExoticLLP", pWARN) << "ARGH. NOT IN TOP VOLUME.";
  */

  if( original_pathString.find( fTopVolume.c_str() ) != string::npos ) {
    const double halfpi = constants::kPi / 2.0;
    Point dl_point = std::pair< double, double >( -halfpi, 0.0 );
    Point ul_point = std::pair< double, double >( halfpi, 0.0 );
    Point dr_point = std::pair< double, double >( -halfpi, 4.0 * halfpi );
    Point ur_point = std::pair< double, double >( halfpi, 4.0 * halfpi );
    PointRaster l_raster = std::pair< Point, Point >( dl_point, ul_point );
    PointRaster r_raster = std::pair< Point, Point >( dr_point, ur_point );
    alpha.emplace_back( l_raster );
    alpha.emplace_back( r_raster );
    return alpha;
  }

  // Get the separation between top volume origin and start point
  const TVector3 seed_vector = fTopVolumeOrigin - fOriginPoint;

  //LOG( "ExoticLLP", pDEBUG ) << "seed_vector = " << utils::print::Vec3AsString(&seed_vector);

  // We make a potentially strong assumption here, that the top_volume is simply connected.
  fMomentum = seed_vector.Unit();
  //fMomentum = booked_momentum;

  // Of course, the angles are defined with respect to the momentum axis...
  // We need to first project seed_vector onto the theta = 0 plane
  // where the theta = 0 plane is defined by fMomentum and the tranverse component of seed_vector
  TVector3 transverse_seed = seed_vector - fAxis.Dot( seed_vector ) * fAxis;
  if( transverse_seed.Mag() > 1.0e-3 * seed_vector.Mag() ) fThetaAxis = transverse_seed.Unit();
  else{ // just pick a convenient direction
    if( fAxis.X() == 0.0 ) fThetaAxis = TVector3( 0.0, fAxis.Z(), -fAxis.Y() );
    else if( fAxis.Y() == 0.0 ) fThetaAxis = TVector3( fAxis.Z(), 0.0, fAxis.X() );
    else if( fAxis.Z() == 0.0 ) fThetaAxis = TVector3( fAxis.Y(), -fAxis.X(), 0.0 );
    else fThetaAxis = TVector3( fAxis.Y() * fAxis.Z(),
				-0.5 * fAxis.X() * fAxis.Z(),
				-0.5 * fAxis.X() * fAxis.Y() );
  }

  // Now that lets us complete the "momentum" system of coordinates
  fPhiAxis = fAxis.Cross( fThetaAxis );

  /*
  LOG( "ExoticLLP", pDEBUG )
    << "\nfAxis      = " << utils::print::Vec3AsString(&fAxis)
    << "\nfThetaAxis = " << utils::print::Vec3AsString(&fThetaAxis)
    << "\nfPhiAxis   = " << utils::print::Vec3AsString(&fPhiAxis);
  */

  // Get the projection of that vector onto the appropriate coordinate system
  const double sep_from_theta = seed_vector.Dot( fPhiAxis );
  const double sep_from_phi   = seed_vector.Dot( fThetaAxis );

  TVector3 projection_theta = seed_vector - sep_from_theta * fPhiAxis;
  TVector3 projection_phi = seed_vector - sep_from_phi * fThetaAxis;

  // If using the small angle approximation, go do a FAST calculation no raytracing.
  // Assumes the detector is a square face of size ((the bbox transverse size)) at the baseline.
  // This is an overestimation of the angular acceptance, and it's a configurable option.

  if( m_use_saa ) alpha = VolumeSeeker::SmallAngleRegion();
  else if( m_use_cmv ) alpha = VolumeSeeker::ComputerVision();
  else {
    VolumeSeeker::Rasterise( alpha, true ); // once moving from centre to the right
    VolumeSeeker::Rasterise( alpha, false ); // and once from centre to the left
  }

  // sort the angular region in increasing phi
  std::sort( alpha.begin(), alpha.end(), 
	     []( const PointRaster & a, const PointRaster &b ){
	       // first and second points of a PointRaster have same phi
	       Point pta = a.first, ptb = b.first;
	       return pta.second < ptb.second;
	     } );

  //LOG( "ExoticLLP", pDEBUG ) << "Angular region has " << alpha.size() << " rasters.";

  // just in case, restore the original member variables
  fOriginPoint = booked_origin_point;
  fOriginPointROOT = booked_origin_point_ROOT;
  fMomentum = booked_momentum;
  
  //return std::make_tuple( alpha, beta );
  return alpha;
}
//____________________________________________________________________________
AngularRegion VolumeSeeker::SmallAngleRegion() const
{
  // This AngularRegion will be a square of (-zeta, -zeta) --> (zeta, zeta)
  // Zeta is the diagonal of the BBox over the baseline (or accurately half that: from centre to diag)

  AngularRegion alpha;

  //const double transverse_size = std::sqrt( fLx*fLx + fLy*fLy + fLz*fLz );
  // This syntax is so stupid! I need to write a library that extends __comp to multiple args...
  const double transverse_size = 
    std::max( fLx * std::sqrt( 1.0 - fAxis.X() * fAxis.X() ), 
	      std::max ( fLy * std::sqrt( 1.0 - fAxis.Y() * fAxis.Y() ), 
			 fLz * std::sqrt( 1.0 - fAxis.Z() * fAxis.Z() ) ) );

  /* To get the baseline, calculate the possible intersections of the line from fOriginPoint to fTopVolumeOrigin
   * with each face. Then get the minimum of each.
   * For accounting, we will evaluate the faces at (x = +- fLx, y = +- fLy, z = +- fLz) in that order.
   */
  
  const TVector3 seed_vector = fTopVolumeOrigin - fOriginPoint;
  TVector3 baseline_seed = seed_vector.Unit();
  std::array< TVector3, 3 > aNormalFaces = { TVector3( 1, 0, 0 ),   // +- fLx
					     TVector3( 0, 1, 0 ),   // +- fLy
					     TVector3( 0, 0, 1 ) }; // +- fLz
  std::array< double, 3 > aProductFaces = { baseline_seed.Dot( aNormalFaces[0] ),
					    baseline_seed.Dot( aNormalFaces[1] ),
					    baseline_seed.Dot( aNormalFaces[2] ) };
  std::array< double, 3 > aDimensions = { fLx, fLy, fLz }; // syntactic sugar

  // Check each pair of faces in turn
  double baseline = seed_vector.Mag(); // the distance between origin and topVolumeCentre
  for( int iface = 0; iface < 3; iface++ ) {
    if( aProductFaces[iface] == 0.0 ) continue; // no intersection for sure

    TVector3 pos_vec =  aDimensions[iface] * aNormalFaces[iface] + fTopVolumeOffset;
    TVector3 neg_vec = -aDimensions[iface] * aNormalFaces[iface] + fTopVolumeOffset;

    TVector3 pos_Dvec = pos_vec - fTopVolumeOrigin;
    TVector3 neg_Dvec = neg_vec - fTopVolumeOrigin;

    double t_pos = pos_Dvec.Dot( aNormalFaces[iface] ) / aProductFaces[iface];
    double t_neg = neg_Dvec.Dot( aNormalFaces[iface] ) / aProductFaces[iface];

    // evaluate the intersection point at these values
    TVector3 intersection_pos = fTopVolumeOrigin + baseline_seed * t_pos;
    TVector3 intersection_neg = fTopVolumeOrigin + baseline_seed * t_neg;
    
    TVector3 dev_pos = intersection_pos - fTopVolumeOrigin;
    TVector3 dev_neg = intersection_neg - fTopVolumeOrigin;

    // which is a "seed vector" of...
    TVector3 pos_seed = intersection_pos - fOriginPoint;
    TVector3 neg_seed = intersection_neg - fOriginPoint;

    // if the intersection is not in the face, then scale the vector to be larger than baseline as it is not a valid intersection
    switch( iface ) {
    case 0: // evaluating X, check Y and Z
      if( std::abs( dev_pos.Y() ) > fLy || std::abs( dev_pos.Z() ) > fLz ) pos_seed.SetMag( 2.0 * baseline );
      if( std::abs( dev_neg.Y() ) > fLy || std::abs( dev_neg.Z() ) > fLz ) neg_seed.SetMag( 2.0 * baseline );
      break;
    case 1: // evaluating Y, check X and Z
      if( std::abs( dev_pos.X() ) > fLx || std::abs( dev_pos.Z() ) > fLz ) pos_seed.SetMag( 2.0 * baseline );
      if( std::abs( dev_neg.X() ) > fLx || std::abs( dev_neg.Z() ) > fLz ) neg_seed.SetMag( 2.0 * baseline );
      break;
    case 2: // evaluating Z, check X and Y
      if( std::abs( dev_pos.X() ) > fLx || std::abs( dev_pos.Y() ) > fLy ) pos_seed.SetMag( 2.0 * baseline );
      if( std::abs( dev_neg.X() ) > fLx || std::abs( dev_neg.Y() ) > fLy ) neg_seed.SetMag( 2.0 * baseline );
      break;
    } // argh I could not avoid the switch...

    // only update baseline if the magnitude of either is smaller
    baseline = std::min( baseline, std::min( pos_seed.Mag(), neg_seed.Mag() ) );
    LOG( "ExoticLLP", pDEBUG ) << "Found intersection at \npos = " 
			       << utils::print::Vec3AsString( &intersection_pos )
			       << ", \nneg = " << utils::print::Vec3AsString( &intersection_neg )
			       << "\n==> baseline = " << baseline;

  } // check each pair of faces in turn
  
  const double zeta = std::atan( transverse_size / baseline );

  Point pt_dl = std::pair<double, double>( -zeta, -constants::kPi );
  Point pt_ul = std::pair<double, double>(  zeta, -constants::kPi );
  Point pt_dr = std::pair<double, double>( -zeta,  constants::kPi );
  Point pt_ur = std::pair<double, double>(  zeta,  constants::kPi );

  PointRaster ras_left  = std::pair<Point, Point>( pt_dl, pt_ul );
  PointRaster ras_right = std::pair<Point, Point>( pt_dr, pt_ur );
  
  alpha.emplace_back( ras_left ); alpha.emplace_back( ras_right );

  return alpha;
}
//____________________________________________________________________________
Vertex VolumeSeeker::FindIntersection( std::array< Vertex, 4 > path, 
				       const TVector3 unit, const double d ) const
{
  double lambda = -9.0;
  int active_vertex = -1;
  Vertex v0, v1;

  while( (lambda < 0.0 || lambda > 1.0) && active_vertex < 3 ) {
    active_vertex++;
    v0 = path[active_vertex];
    v1 = path[active_vertex+1];
    lambda = v0.Intersection(v1, unit, d);
  }

  if( lambda < 0.0 || lambda > 1.0 ) // no intersection along this path
    return Vertex( -1, 0.0, 0.0, 0.0 );

  double vx = v0.fX + lambda * ( v0.Displacement(v1) ).X();
  double vy = v0.fY + lambda * ( v0.Displacement(v1) ).Y();
  double vz = v0.fZ + lambda * ( v0.Displacement(v1) ).Z();

  vx += fTopVolumeOffset.X();
  vy += fTopVolumeOffset.Y();
  vz += fTopVolumeOffset.Z();

  return Vertex( 0, vx, vy, vz );
}
//____________________________________________________________________________
AngularRegion VolumeSeeker::ComputerVision() const
{
  // This is a more sophisticated algorithm to get an intersection of the bounding box with 
  // the plane normal to fAxis. 
  // This implements the vertex algorithm from Rezk Salama and Kolb, Proc Vision Modeling and Visualization 2005 115-122.

  AngularRegion alpha;
  
  // Get the unit normal vector of the viewing plane of the box...
  TVector3 norm_vec = -(fOriginPoint.Unit());

  // First, one declares the eight vertices of the bounding box and the starting vertex. 
  Vertex vOrigin( 0, fOriginPoint ); // USER m
  //double origin_dist = vOrigin.Dist( Vertex() ); // Vertex() is just (0,0,0)

  std::vector< Vertex > vVertices;
  vVertices.emplace_back( 1, -fLx, -fLy, -fLz );
  vVertices.emplace_back( 2,  fLx, -fLy, -fLz );
  vVertices.emplace_back( 3,  fLx,  fLy, -fLz );
  vVertices.emplace_back( 4, -fLx,  fLy, -fLz );
  vVertices.emplace_back( 5,  fLx,  fLy,  fLz );
  vVertices.emplace_back( 6, -fLx,  fLy,  fLz );
  vVertices.emplace_back( 7, -fLx, -fLy,  fLz );
  vVertices.emplace_back( 8,  fLx, -fLy,  fLz );

  // Now take the distances of all these vertices from origin
  std::vector< double > distances; // m
  for( std::vector< Vertex >::iterator it_vtx = vVertices.begin();
       it_vtx != vVertices.end(); ++it_vtx ) distances.emplace_back( vOrigin.Dist( *it_vtx ) );
  std::vector< double >::iterator min_distance = std::min_element( distances.begin(), distances.end() );

  int iMinDist = min_distance - distances.begin();
  // Now we construct the paths based on which vertex we picked.
  int first_idx = (vVertices.at(iMinDist)).fIndex;

  std::array< Vertex, 4 > path_1, path_2, path_3;
  std::pair< Vertex, Vertex > phantom_1, phantom_2, phantom_3;
  switch( first_idx ) {
  case 1: // 1 --> 5
    path_1 = { vVertices.at(0), vVertices.at(1), vVertices.at(7), vVertices.at(4) }; // 1, 2, 8, 5
    path_2 = { vVertices.at(0), vVertices.at(3), vVertices.at(2), vVertices.at(4) }; // 1, 4, 3, 5
    path_3 = { vVertices.at(0), vVertices.at(6), vVertices.at(5), vVertices.at(4) }; // 1, 7, 6, 5
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(1), vVertices.at(2) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(3), vVertices.at(5) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(6), vVertices.at(7) );
    break;
  case 5: // 5 --> 1, backwards of 1 --> 5. Note invert paths 2 and 3 to keep the system RH
    path_1 = { vVertices.at(4), vVertices.at(7), vVertices.at(1), vVertices.at(0) }; 
    path_3 = { vVertices.at(4), vVertices.at(2), vVertices.at(3), vVertices.at(0) };
    path_2 = { vVertices.at(4), vVertices.at(5), vVertices.at(6), vVertices.at(0) }; 
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(2), vVertices.at(1) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(5), vVertices.at(3) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(7), vVertices.at(6) );
    break;
  case 2: // 2 --> 6
    path_1 = { vVertices.at(1), vVertices.at(7), vVertices.at(4), vVertices.at(5) }; // 2, 8, 5, 6
    path_2 = { vVertices.at(1), vVertices.at(2), vVertices.at(3), vVertices.at(5) }; // 2, 3, 4, 6
    path_3 = { vVertices.at(1), vVertices.at(0), vVertices.at(6), vVertices.at(5) }; // 2, 1, 7, 6
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(7), vVertices.at(6) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(2), vVertices.at(4) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(0), vVertices.at(3) );
    break;
  case 6: // 6 --> 2
    path_1 = { vVertices.at(5), vVertices.at(4), vVertices.at(7), vVertices.at(1) };
    path_3 = { vVertices.at(5), vVertices.at(3), vVertices.at(2), vVertices.at(1) };
    path_2 = { vVertices.at(5), vVertices.at(6), vVertices.at(0), vVertices.at(1) };
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(6), vVertices.at(7) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(4), vVertices.at(2) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(3), vVertices.at(0) );
    break;
  case 3: // 3 --> 7
    path_1 = { vVertices.at(2), vVertices.at(4), vVertices.at(5), vVertices.at(6) }; // 3, 5, 6, 7
    path_2 = { vVertices.at(2), vVertices.at(3), vVertices.at(0), vVertices.at(6) }; // 3, 4, 1, 7
    path_3 = { vVertices.at(2), vVertices.at(1), vVertices.at(7), vVertices.at(6) }; // 3, 2, 8, 7
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(4), vVertices.at(7) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(3), vVertices.at(5) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(1), vVertices.at(0) );
    break;
  case 7: // 7 --> 3
    path_1 = { vVertices.at(6), vVertices.at(5), vVertices.at(4), vVertices.at(2) };
    path_3 = { vVertices.at(6), vVertices.at(0), vVertices.at(3), vVertices.at(2) };
    path_2 = { vVertices.at(6), vVertices.at(7), vVertices.at(1), vVertices.at(2) };
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(7), vVertices.at(4) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(5), vVertices.at(3) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(0), vVertices.at(1) );
    break;
  case 4: // 4 --> 8
    path_1 = { vVertices.at(3), vVertices.at(2), vVertices.at(4), vVertices.at(7) }; // 4, 3, 5, 8
    path_2 = { vVertices.at(3), vVertices.at(5), vVertices.at(6), vVertices.at(7) }; // 4, 6, 7, 8
    path_3 = { vVertices.at(3), vVertices.at(0), vVertices.at(1), vVertices.at(7) }; // 4, 1, 2, 8
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(2), vVertices.at(1) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(4), vVertices.at(5) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(0), vVertices.at(6) );
    break;
  case 8: // 8 --> 4
    path_1 = { vVertices.at(7), vVertices.at(4), vVertices.at(2), vVertices.at(3) };
    path_3 = { vVertices.at(7), vVertices.at(6), vVertices.at(5), vVertices.at(3) };
    path_2 = { vVertices.at(7), vVertices.at(1), vVertices.at(0), vVertices.at(3) };
    phantom_1 = std::pair< Vertex, Vertex >( vVertices.at(1), vVertices.at(2) );
    phantom_3 = std::pair< Vertex, Vertex >( vVertices.at(5), vVertices.at(4) );
    phantom_2 = std::pair< Vertex, Vertex >( vVertices.at(6), vVertices.at(0) );
    break;
  default: break;
  }

  // Now we have the three independent paths set up with each consisting of four Vertices.
  // Find for each path the intersection
  
  std::vector< Vertex > vIntersections;
  std::array< Vertex, 4 > active_path; 

  // check path 1
  active_path = path_1;
  vIntersections.emplace_back( this->FindIntersection( active_path, norm_vec, 0.0 ) );
  
  // check path 2
  active_path = path_2;
  vIntersections.emplace_back( this->FindIntersection( active_path, norm_vec, 0.0 ) );

  // check path 3
  active_path = path_3;
  vIntersections.emplace_back( this->FindIntersection( active_path, norm_vec, 0.0 ) );

  // Now there are potentially other vertices along the "phantom" edges we didn't pick up. 
  // We need to check them explicitly.
  
  double lambda;
  Vertex v0, v1;
  // check phantom 1
  lambda = -9.0; v0 = phantom_1.first; v1 = phantom_1.second;
  lambda = v0.Intersection( v1, norm_vec, 0.0 );
  if( lambda >= 0.0 && lambda <= 1.0 )
    vIntersections.emplace_back( Vertex( 0, 
					 v0.fX + lambda * ( v0.Displacement(v1) ).X() + fTopVolumeOffset.X(),
					 v0.fY + lambda * ( v0.Displacement(v1) ).Y() + fTopVolumeOffset.Y(),
					 v0.fZ + lambda * ( v0.Displacement(v1) ).Z() + fTopVolumeOffset.Z() ) );
  // check phantom 2
  lambda = -9.0; v0 = phantom_2.first; v1 = phantom_2.second;
  lambda = v0.Intersection( v1, norm_vec, 0.0 );
  if( lambda >= 0.0 && lambda <= 1.0 )
    vIntersections.emplace_back( Vertex( 0, 
					 v0.fX + lambda * ( v0.Displacement(v1) ).X() + fTopVolumeOffset.X(),
					 v0.fY + lambda * ( v0.Displacement(v1) ).Y() + fTopVolumeOffset.Y(),
					 v0.fZ + lambda * ( v0.Displacement(v1) ).Z() + fTopVolumeOffset.Z() ) );
  // check phantom 3
  lambda = -9.0; v0 = phantom_3.first; v1 = phantom_3.second;
  lambda = v0.Intersection( v1, norm_vec, 0.0 );
  if( lambda >= 0.0 && lambda <= 1.0 )
    vIntersections.emplace_back( Vertex( 0, 
					 v0.fX + lambda * ( v0.Displacement(v1) ).X() + fTopVolumeOffset.X(),
					 v0.fY + lambda * ( v0.Displacement(v1) ).Y() + fTopVolumeOffset.Y(),
					 v0.fZ + lambda * ( v0.Displacement(v1) ).Z() + fTopVolumeOffset.Z() ) );

  LOG( "ExoticLLP", pDEBUG ) << "Found " << vIntersections.size() << " intersection vertices";

  /*
   * By now we have from 3 to 6 intersection vertices. Need to sort them, so that they
   * form a *convex* shape! Then we will calculate the area of this volume
   */

  // Need the int to keep track of which vector element we're using. I'll be popping elements to sort
  std::vector< std::pair< int, double > > vAngles; // angle is in (\phi, cos\theta) space wrt centre

  // Define a local coordinate system on the plane
  TVector3 loc_X_axis( norm_vec.Y() * norm_vec.Z(), 
		       norm_vec.X() * norm_vec.Y(),
		       -2.0 * norm_vec.X() * norm_vec.Y() );
  loc_X_axis = loc_X_axis.Unit();
  TVector3 loc_Y_axis = norm_vec.Cross( loc_X_axis );

  LOG( "ExoticLLP", pDEBUG ) 
    << "\nnorm_vec   = " << utils::print::Vec3AsString( &norm_vec )
    << "\nloc_X_axis = " << utils::print::Vec3AsString( &loc_X_axis )
    << "\nloc_Y_axis = " << utils::print::Vec3AsString( &loc_Y_axis );

  // project each vertex onto the localised axes
  std::vector<Point> vPoints;
  std::ostringstream psts;
  for( std::vector< Vertex >::iterator it_vtx = vIntersections.begin();
       it_vtx != vIntersections.end(); ++it_vtx ) {
    double vtx_X = (*it_vtx).fX * loc_X_axis.X() + 
      (*it_vtx).fY * loc_X_axis.Y() + (*it_vtx).fZ * loc_X_axis.Z();
    double vtx_Y = (*it_vtx).fX * loc_Y_axis.X() + 
      (*it_vtx).fY * loc_Y_axis.Y() + (*it_vtx).fZ * loc_Y_axis.Z();
    vPoints.emplace_back( std::pair< double, double >( vtx_X, vtx_Y ) );
    psts << "\nPoint with 3D coordinates ( " << (*it_vtx).fX
	 << ", " << (*it_vtx).fY << ", " << (*it_vtx).fZ << " ) maps to "
	 << "( " << vtx_X << ", " << vtx_Y << " )";
  }
  LOG( "ExoticLLP", pDEBUG ) << psts.str();

  // and from these localised coordinates calculate the angle...
  for( std::vector< Point >::iterator it_cor = vPoints.begin();
       it_cor != vPoints.end(); ++it_cor ) {
    int elem_idx = it_cor - vPoints.begin();
    
    double elem_ang = std::atan( (*it_cor).second / (*it_cor).first );  
    // support on [-pi/2, pi/2] i.e 1st and 4th quadrants
    
    // get the quadrants right
    if( (*it_cor).first >= 0.0 
	&& (*it_cor).second >= 0.0 ) elem_ang *= 1.0; // first 
    else if( (*it_cor).first < 0.0 
	     && (*it_cor).second >= 0.0 ) elem_ang = constants::kPi + elem_ang; // second
    else if( (*it_cor).first < 0.0 
	     && (*it_cor).second < 0.0 ) elem_ang += constants::kPi; // third
    else if( (*it_cor).first >= 0.0 
	     && (*it_cor).second < 0.0 ) elem_ang = 2.0 * constants::kPi + elem_ang; // fourth
    
    vAngles.emplace_back( std::pair< int, double >( elem_idx, elem_ang ) );
  }
  
  // now we have these angles, we can sort the coordinates into a vector such that they're clockwise
  std::vector< int > sorted_intersection_idcs;
  while( vAngles.size() > 0 ) {
    std::vector< std::pair< int, double > >::iterator it_ang = 
      std::min_element( vAngles.begin(), vAngles.end(), 
			[](std::pair<int, double> a, std::pair<int, double> b) 
			{ return a.second < b.second; }); // inject lambda to compare second members

    sorted_intersection_idcs.emplace_back( (*it_ang).first ); // get the index of the smallest angle
    vAngles.erase( it_ang ); // dealt with this
  } // sort the intersection vertices by angle such that they're clockwise

  // and put the full vertices into a sorted vector
  std::vector< Vertex > vSortedVertices;
  std::vector< Point > vSortedPoints;
  for( std::vector<int>::iterator it_sii = sorted_intersection_idcs.begin();
       it_sii != sorted_intersection_idcs.end(); ++it_sii ) {
    vSortedVertices.emplace_back( vIntersections.at( *it_sii ) );
    vSortedPoints.emplace_back( vPoints.at( *it_sii ) );
  }

  std::ostringstream csts;
  csts << "\nHere are the sorted coordinates of each vertex:";
  for( std::vector< Vertex >::iterator it_cor = vSortedVertices.begin();
       it_cor != vSortedVertices.end(); ++it_cor ) {
    int idx = it_cor - vSortedVertices.begin();
    csts << "\nVertex " << idx << ": "
	 << " ( " << (*it_cor).fX << ", " << (*it_cor).fY << ", " << (*it_cor).fZ << " )"
	 << "\n\t(local point = ( " << vSortedPoints.at( idx ).first
	 << ", " << vSortedPoints.at( idx ).second << " ) )";
  }
  LOG( "ExoticLLP", pDEBUG ) << csts.str();

  // Get the area of this polygon
  double area = 0.0;
  std::ostringstream asts;

  // use the half-cross product formula instead
  double ox = 0.0, oy = 0.0, oz = 0.0;
  for( std::vector<Vertex>::iterator it_spt = vSortedVertices.begin();
       it_spt != vSortedVertices.end(); ++it_spt ) {
    ox += (*it_spt).fX / static_cast<double>( vSortedVertices.size() );
    oy += (*it_spt).fY / static_cast<double>( vSortedVertices.size() );
    oz += (*it_spt).fZ / static_cast<double>( vSortedVertices.size() );
  }
  for( int ipt = 0; ipt < vSortedVertices.size(); ipt++ ) {
    int idx = ipt+1; if(idx == vSortedVertices.size()) idx = 0;
    Vertex vA = vSortedVertices.at(ipt);
    Vertex vB = vSortedVertices.at(idx);

    double Ax = vA.fX, Ay = vA.fY, Az = vA.fZ;
    double Bx = vB.fX, By = vB.fY, Bz = vB.fZ;

    TVector3 dOA( Ax - ox, Ay - oy, Az - oz );
    TVector3 dOB( Bx - ox, By - oy, Bz - oz );

    double sq_yz = std::pow( dOA.Y() * dOB.Z() - dOA.Z() * dOB.Y(), 2.0 );
    double sq_zx = std::pow( dOA.Z() * dOB.X() - dOA.X() * dOB.Z(), 2.0 );
    double sq_xy = std::pow( dOA.X() * dOB.Y() - dOA.Y() * dOB.X(), 2.0 );

    area += 0.5 * std::sqrt( sq_yz + sq_zx + sq_xy );

    asts << "\nFrom points O = ( " << ox << ", " << oy << ", " << oz << " ),"
	 << " A = ( " << Ax << ", " << Ay << ", " << Az << " ), and "
	 << "B = ( " << Bx << ", " << By << ", " << Bz << " ) we get an area = "
	 << 0.5 * std::sqrt( sq_yz + sq_zx + sq_xy ); 
  }
  LOG( "ExoticLLP", pDEBUG ) << asts.str();

  // and a circle of this area has radius...
  double rad = std::sqrt( area / constants::kPi );

  // Now also get the baseline!

  /* To get the baseline, calculate the possible intersections of the line from fOriginPoint to fTopVolumeOrigin
   * with each face. Then get the minimum of each.
   * For accounting, we will evaluate the faces at (x = +- fLx, y = +- fLy, z = +- fLz) in that order.
   */
  
  const TVector3 seed_vector = fTopVolumeOrigin - fOriginPoint;
  TVector3 baseline_seed = seed_vector.Unit();
  std::array< TVector3, 3 > aNormalFaces = { TVector3( 1, 0, 0 ),   // +- fLx
					     TVector3( 0, 1, 0 ),   // +- fLy
					     TVector3( 0, 0, 1 ) }; // +- fLz
  std::array< double, 3 > aProductFaces = { baseline_seed.Dot( aNormalFaces[0] ),
					    baseline_seed.Dot( aNormalFaces[1] ),
					    baseline_seed.Dot( aNormalFaces[2] ) };
  std::array< double, 3 > aDimensions = { fLx, fLy, fLz }; // syntactic sugar

  // Check each pair of faces in turn
  double baseline = seed_vector.Mag(); // the distance between origin and topVolumeCentre
  for( int iface = 0; iface < 3; iface++ ) {
    if( aProductFaces[iface] == 0.0 ) continue; // no intersection for sure

    TVector3 pos_vec =  aDimensions[iface] * aNormalFaces[iface] + fTopVolumeOffset;
    TVector3 neg_vec = -aDimensions[iface] * aNormalFaces[iface] + fTopVolumeOffset;

    TVector3 pos_Dvec = pos_vec - fTopVolumeOrigin;
    TVector3 neg_Dvec = neg_vec - fTopVolumeOrigin;

    double t_pos = pos_Dvec.Dot( aNormalFaces[iface] ) / aProductFaces[iface];
    double t_neg = neg_Dvec.Dot( aNormalFaces[iface] ) / aProductFaces[iface];

    // evaluate the intersection point at these values
    TVector3 intersection_pos = fTopVolumeOrigin + baseline_seed * t_pos;
    TVector3 intersection_neg = fTopVolumeOrigin + baseline_seed * t_neg;
    
    TVector3 dev_pos = intersection_pos - fTopVolumeOrigin;
    TVector3 dev_neg = intersection_neg - fTopVolumeOrigin;

    // which is a "seed vector" of...
    TVector3 pos_seed = intersection_pos - fOriginPoint;
    TVector3 neg_seed = intersection_neg - fOriginPoint;

    // if the intersection is not in the face, then scale the vector to be larger than baseline as it is not a valid intersection
    switch( iface ) {
    case 0: // evaluating X, check Y and Z
      if( std::abs( dev_pos.Y() ) > fLy || std::abs( dev_pos.Z() ) > fLz ) pos_seed.SetMag( 2.0 * baseline );
      if( std::abs( dev_neg.Y() ) > fLy || std::abs( dev_neg.Z() ) > fLz ) neg_seed.SetMag( 2.0 * baseline );
      break;
    case 1: // evaluating Y, check X and Z
      if( std::abs( dev_pos.X() ) > fLx || std::abs( dev_pos.Z() ) > fLz ) pos_seed.SetMag( 2.0 * baseline );
      if( std::abs( dev_neg.X() ) > fLx || std::abs( dev_neg.Z() ) > fLz ) neg_seed.SetMag( 2.0 * baseline );
      break;
    case 2: // evaluating Z, check X and Y
      if( std::abs( dev_pos.X() ) > fLx || std::abs( dev_pos.Y() ) > fLy ) pos_seed.SetMag( 2.0 * baseline );
      if( std::abs( dev_neg.X() ) > fLx || std::abs( dev_neg.Y() ) > fLy ) neg_seed.SetMag( 2.0 * baseline );
      break;
    } // argh I could not avoid the switch...

    // only update baseline if the magnitude of either is smaller
    baseline = std::min( baseline, std::min( pos_seed.Mag(), neg_seed.Mag() ) );
    LOG( "ExoticLLP", pDEBUG ) << "Found intersection at \npos = " 
			       << utils::print::Vec3AsString( &intersection_pos )
			       << ", \nneg = " << utils::print::Vec3AsString( &intersection_neg )
			       << "\n==> baseline = " << baseline;

  } // check each pair of faces in turn

  // return an AngularRegion that holds (radius, baseline) in the first point
  Point pt = std::pair< double, double >( rad, baseline );
  PointRaster rs = std::pair<Point, Point>(pt, pt);

  LOG( "ExoticLLP", pDEBUG ) << "rad = " << rad << ", baseline = " << baseline;

  alpha.emplace_back(rs);
    
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

  //double size = VolumeSeeker::Trapezoid( upper_points, lower_points );
  double size = 0.0;

  if( m_use_saa ) {
    double zeta = upper_points.at(0).first;
    size = 2.0 * constants::kPi * ( 1.0 - std::cos(zeta) );
  }

  if( m_use_cmv ) {
    double rad = upper_points.at(0).first;
    double baseline = upper_points.at(0).second;

    double zeta = std::atan( rad / baseline );
    size = 2.0 * constants::kPi * ( 1.0 - std::cos(zeta) ); // This is still valid...
  }

  LOG( "ExoticLLP", pDEBUG ) << "size/4pi = " << size / (4.0 * constants::kPi);
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
     * I can't think of anyone who'd build a detector that looks like that, it's not simply connected.
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
	total_dn -= width * (cand_U1 + cand_U2)/2.0;
      } else if( theta_dn_prev >= 0.0 ) {
	// 0 < D1 < U1 --> 1-cos(U1) > 1-cos(D1) > 0
	// Also 0 < D2 < U2 as there's not crossing of zero

	total_up += width * (cand_U1 + cand_U2)/2.0;
	total_dn -= width * (cand_D1 + cand_D2)/2.0;
      }
    }

    /*
     * CASE 3: Zero-deflection OK in one raster, but not in the other. 
     * There is a zero-crossing which we must estimate by interpolation
     * and reduce to a simultaneous evaluation of Case 1 and Case 2.
     */

    if( ( theta_up_prev * theta_dn_prev <= 0.0 && theta_up_curr * theta_dn_curr > 0.0 ) ||
	( theta_up_prev * theta_dn_prev > 0.0 && theta_up_curr * theta_dn_curr <= 0.0 ) ) {

      // First, find the zero-crossing
      double interp = 0.0;
      double dtheta_dn = theta_dn_curr - theta_dn_prev;
      double dtheta_up = theta_up_curr - theta_up_prev;

      double dphi = current_point_up.second - previous_point_up.second;

      // depending on if crossing is up or dn, interp accordingly
      double dth; double t0 = 0.0;
      if( theta_up_prev * theta_up_curr <= 0.0 ) { // crosses at up
	dth = dtheta_up; t0 = theta_up_curr;
      } else { // crosses at dn
	dth = dtheta_dn; t0 = theta_dn_curr;
      }
      interp = std::abs(t0) / std::abs(dth);

      double theta_zero = 0.0;
      double phi_zero = current_point_up.second + interp * dphi;

      double theta_up_zero, theta_dn_zero;

      // make a new pair of Points
      Point pt_uz, pt_dz; double cand_UZ, cand_DZ;
      if( theta_up_prev * theta_up_curr <= 0.0 ) { // crosses at up
	theta_zero = current_point_dn.first + interp * dth;
	pt_uz = std::pair< double, double >( 0.0, phi_zero ); cand_UZ = 0.0;
	pt_dz = std::pair< double, double >( theta_zero, phi_zero ); cand_DZ = 1.0 - std::cos( theta_zero );
	theta_up_zero = 0.0; theta_dn_zero = theta_zero;
      } else { // crosses at dn
	pt_uz = std::pair< double, double >( theta_zero, phi_zero ); cand_UZ = 1.0 - std::cos( theta_zero );
	pt_dz = std::pair< double, double >( 0.0, phi_zero ); cand_DZ = 0.0;
	theta_up_zero = theta_zero; theta_dn_zero = 0.0;
      }

      // there are two ranges: [curr --> zero] and [zero --> prev]
      // and now we evaluate which is case 1 and which is case 2.

      if( theta_up_curr * theta_dn_curr <= 0.0 ) { // [curr --> zero] 1, [zero --> prev] 2
	// case 1 bit
	double h1_cz = cand_U1 + cand_D1;
	double h2_cz = cand_UZ + cand_DZ;

	total_up += interp * width * (h1_cz+h2_cz)/2.0;
	total_dn += 0.0;

	// case 2 bit
	if( theta_up_prev < 0.0 ) { 
	  // DZ < UZ < 0 --> 1-cos(DZ) > 1-cos(UZ) > 0
	  // Also D2 < U2 < 0 as there's no crossing of zero

	  total_up += (1.0 - interp) * width * (cand_DZ + cand_D2)/2.0;
	  total_dn -= (1.0 - interp) * width * (cand_UZ + cand_U2)/2.0;
	} else if( theta_dn_prev >= 0.0 ) {
	  // 0 < DZ < UZ --> 1-cos(U1) > 1-cos(D1) > 0
	  // Also 0 < D2 < U2 as there's not crossing of zero

	  total_up += (1.0 - interp) * width * (cand_UZ + cand_U2)/2.0;
	  total_dn -= (1.0 - interp) * width * (cand_DZ + cand_D2)/2.0;
	}
      } else { // [curr --> zero] 2, [zero --> prev] 1
	// case 1 bit
	double h1_cz = cand_UZ + cand_DZ;
	double h2_cz = cand_U2 + cand_D2;

	total_up += (1.0 - interp) * width * (h1_cz+h2_cz)/2.0;
	total_dn += 0.0;

	// case 2 bit
	if( theta_up_zero < 0.0 ) { 
	  // DZ < UZ < 0 --> 1-cos(DZ) > 1-cos(UZ) > 0
	  // Also D1 < U1 < 0 as there's no crossing of zero

	  total_up += interp * width * (cand_DZ + cand_D1)/2.0;
	  total_dn -= interp * width * (cand_UZ + cand_U1)/2.0;
	} else if( theta_dn_zero >= 0.0 ) {
	  // 0 < DZ < UZ --> 1-cos(UZ) > 1-cos(DZ) > 0
	  // Also 0 < D1 < U1 as there's not crossing of zero

	  total_up += interp * width * (cand_UZ + cand_U1)/2.0;
	  total_dn -= interp * width * (cand_DZ + cand_D1)/2.0;
	}
      }

    } // case 3


    previous_point_up = current_point_up; // update
    previous_point_dn = current_point_dn; // update
  }

  //LOG( "ExoticLLP", pDEBUG ) << "total_up, total_dn = " << total_up << ", " << total_dn;

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
    
    /*
    LOG( "ExoticLLP", pDEBUG )
      << "Deflections: phi = " << 0.0
      << ", theta_min, max = " << thetaMin * 180.0 / constants::kPi
      << ", " << thetaMax * 180.0 / constants::kPi << " [deg]";
    */
    
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

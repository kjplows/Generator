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

\cpright   Copyright (c) 2003-2024, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _LLP_VERTEX_GENERATOR_H_
#define _LLP_VERTEX_GENERATOR_H_

#include <cmath>
#include <cassert>
#include <list>
#include <tuple>

//#include <TGeoMatrix.h>

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

#include "Physics/ExoticLLP/LLPConfigurator.h"
#include "Physics/ExoticLLP/LLPGeomRecordVisitorI.h"
#include "Physics/ExoticLLP/LLPFluxContainer.h"
#include "Physics/ExoticLLP/VolumeSeeker.h"

namespace genie {
namespace llp {

  class ExoticLLP;
  class VolumeSeeker;
  class LLPConfigurator;

  class VertexGenerator : public Algorithm {

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

    //void SetGeomFile( std::string geomfile, std::string topVolume ) const;

    void ReadFluxContainer( const genie::llp::FluxContainer flc ) const;
    genie::llp::FluxContainer RetrieveFluxContainer() const { return fFluxContainer; }

  private:

    void LoadConfig();
    
    // --------------------------------------------------
    // Utilities
    // --------------------------------------------------

    // simple getter
    void GetInterestingPoints( TLorentzVector & entryPoint, TLorentzVector & exitPoint, TLorentzVector & decayPoint ) const;

    // calculate travel length in detector
    double CalcTravelLength( double betaMag, double CoMLifetime, double maxLength ) const;

    // get max length in detector
    double GetMaxLength() const;

    // assign decay point given length
    TLorentzVector GetDecayPoint( double travelLength ) const;

    // --------------------------------------------------

    //mutable std::string fGeomFile = "";
    //mutable std::string fTopVolume = "";

    mutable bool fIsConfigLoaded = false;

    mutable double lunits = genie::units::m; mutable std::string lunitString = "m";
    mutable double aunits = genie::units::rad;
    mutable double tunits = genie::units::ns; mutable std::string tunitString = "ns";

    mutable genie::llp::FluxContainer fFluxContainer;
    mutable TLorentzVector fDecayPoint; // m, ns
    
    mutable double fLifetime = 0.0; // LLP lifetime c*tau, in m
    
    mutable double kNewSpeedOfLight = genie::units::kSpeedOfLight * ( genie::units::m / lunits ) / ( genie::units::s / tunits );

    ClassDef(VertexGenerator, 2)
  }; // class VertexGenerator

} // namespace llp
} // namespace genie

#endif // #ifndef _LLP_VERTEX_GENERATOR_H_

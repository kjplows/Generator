//____________________________________________________________________________
/*!

\class   genie::llp::GeomRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interfaces
         for the LLP VertexGenerator module.
	 Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Liverpool

	 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	 University of Liverpool & STFC Rutherford Appleton Laboratory

\created May 23rd, 2024

\cpright Copyright (c) 2003-2023, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _LLP_GEOM_RECORD_VISITOR_I_H_
#define _LLP_GEOM_RECORD_VISITOR_I_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"

namespace genie {

  class GHepRecord;

  namespace llp {

    class GeomRecordVisitorI: public EventRecordVisitorI {

    public:
      
      virtual ~GeomRecordVisitorI();

      //-- define the GeomRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      virtual void SetGeomFile( std::string geomfile, std::string topVolume ) const = 0;

    protected:

      GeomRecordVisitorI();
      GeomRecordVisitorI(string name);
      GeomRecordVisitorI(string name, string config);

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
      virtual TGeoMatrix * FindFullTransformation( TGeoVolume * top_vol, TGeoVolume * tar_vol ) const = 0;
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

      ClassDef(GeomRecordVisitorI, 1)
    };

  } // namespace llp

} // namespace genie

#endif // #ifndef _LLP_GEOM_RECORD_VISTOR_I_H_

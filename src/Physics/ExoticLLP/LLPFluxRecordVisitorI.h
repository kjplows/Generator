//____________________________________________________________________________
/*!

\class   genie::llp::FluxRecordVisitorI

\brief   Expands the EventRecordVisitorI interface to include public interfaces
         for the ExoticLLP FluxCreator module.
	 Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Liverpool

	 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
	 University of Liverpool & STFC Rutherford Appleton Laboratory

\created May 23rd, 2024

\cpright Copyright (c) 2003-2024, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _LLP_FLUX_RECORD_VISITOR_I_H_
#define _LLP_FLUX_RECORD_VISITOR_I_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
//#include "Physics/ExoticLLP/LLPFluxContainer.h"

#include <TStreamerElement.h>

namespace genie {

  class GHepRecord;

  namespace llp {

    //class FluxContainer;

    class FluxRecordVisitorI: public EventRecordVisitorI {

    public:

      virtual ~FluxRecordVisitorI();

      //-- define the FluxRecordVisitorI interface

      virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

      //virtual FluxContainer RetrieveFluxInfo() const = 0;
      
      virtual void SetInputFluxPath( std::string finpath ) const = 0;
      virtual int GetNFluxEntries() const = 0;
      virtual void SetFirstFluxEntry( int iFirst ) const = 0;

    protected:

      FluxRecordVisitorI();
      FluxRecordVisitorI(string name);
      FluxRecordVisitorI(string name, string config);

      ClassDef(FluxRecordVisitorI, 2)
    };

  } // namespace llp

} // namespace genie

#endif // #ifndef _LLP_FLUX_RECORD_VISITOR_I_H_

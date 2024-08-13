//----------------------------------------------------------------------------
/*
  Helper class to initalise branches with a user-defined alias

\class   genie::llp::AliasedBranch

\brief   C-struct style: pair of <typename T, std::string>.
         A TTree will pass the T member by reference and the string is used to set address

\author  John Plows <kplows \at liverpool.ac.uk>
         University of Liverpool

\created June 27th, 2024

\cpright    Copyright (c) 2003-2024, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
 */
//----------------------------------------------------------------------------

#ifndef _LLP_ALIASED_BRANCH_H_
#define _LLP_ALIASED_BRANCH_H_

#include <array>
#include <string>

namespace genie {
  namespace llp {
    template<typename T>
      class AliasedBranch {

    public:
      AliasedBranch( T vl = T(), std::string al = "" );

      T Value; // Probably ok to expose this...
      std::string Alias;

    private:
      
    }; // class AliasedBranch
  } // namespace llp
} // namespace genie

// Template definition needs to be included in the header file
template<typename T>
genie::llp::AliasedBranch<T>::AliasedBranch( T vl, std::string al ) :
Value(vl), Alias(al) {}

#endif // #ifndef _LLP_ALIASED_BRANCH_H_

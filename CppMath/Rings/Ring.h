//------------------------------------------------------------------------------
/// \file Ring.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Ring, and ring axioms..
/// \ref Ch. 21 Class Hierarchies, 21.2.Design of Class Hierarchies
///   The C++ Programming Language, 4th Ed., Stroustrup;
/// \details Ring (2 ring operations, ring axioms), following interface
/// implementation.
/// \copyright If you find this code useful, feel free to donate directly
/// (username ernestyalumni or email address above), going directly to:
///
/// paypal.me/ernestyalumni
///
/// which won't go through a 3rd. party like indiegogo, kickstarter, patreon.
/// Otherwise, I receive emails and messages on how all my (free) material on
/// physics, math, and engineering have helped students with their studies, and
/// I know what it's like to not have money as a student, but love physics (or
/// math, sciences, etc.), so I am committed to keeping all my material
/// open-source and free, whether or not sufficiently crowdfunded, under the
/// open-source MIT license: feel free to copy, edit, paste, make your own
/// versions, share, use as you wish.
/// Peace out, never give up! -EY
//------------------------------------------------------------------------------
/// COMPILATION TIPS:
///  g++ -std=c++17 Tuple2_main.cpp -o Tuple2_main
//------------------------------------------------------------------------------
#ifndef _RINGS_RING_H_
#define _RINGS_RING_H_

#include "Groups/AbelianGroup.h"

namespace Rings
{

//------------------------------------------------------------------------------
/// \class Element
/// \details A pure abstract base class for a ring Element. Also, use CRTP
/// pattern, for the return type.
/// \ref https://en.wikipedia.org/wiki/Ring_(mathematics)
/// https://stackoverflow.com/questions/27180342/ \
/// pure-virtual-function-in-abstract-class-with-return-type-of-base-derived-type
/// \tparam R stands for ring R, the ring R that Element belongs to.
//------------------------------------------------------------------------------
template <typename R>
class Element : public Groups::AbelianGroups::Element<R>
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    //--------------------------------------------------------------------------
    /// \details A ring is a monoid under multiplication.
    /// \url https://en.wikipedia.org/wiki/Ring_(mathematics)
    //--------------------------------------------------------------------------

    virtual R operator*(const R& a) const = 0;

    virtual R multiplicative_identity() const = 0; // pure virtual function
};

} // namespace Rings

#endif // _RINGS_RING_H_
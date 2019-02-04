//------------------------------------------------------------------------------
/// \file Field.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Commutative Ring.
/// \ref Ch. 21 Class Hierarchies, 21.2.Design of Class Hierarchies
///   The C++ Programming Language, 4th Ed., Stroustrup;
/// \details Commutative Ring following interface
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
///  g++ --std=c++17 -I ../ RR_main.cpp -o RR_main
//------------------------------------------------------------------------------
#ifndef _RINGS_COMMUTATIVE_RINGS_FIELDS_FIELD_H_
#define _RINGS_COMMUTATIVE_RINGS_FIELDS_FIELD_H_

#include "CommutativeRing.h"

namespace Rings
{

namespace CommutativeRings
{

namespace Fields
{

//------------------------------------------------------------------------------
/// \class Element
/// \details A pure abstract base class for a commutative ring Element. Also,
/// use CRTP pattern, for the return type.
/// Replicated base. If ambiguous, make inheritance public virtual.
/// \ref https://en.wikipedia.org/wiki/Ring_(mathematics)
/// https://stackoverflow.com/questions/27180342/ \
/// pure-virtual-function-in-abstract-class-with-return-type-of-base-derived-type
/// Sec. 21.3.6 Replicated vs. Virtual Bases, Stroustrup.
/// \tparam F stands for field F, the field F that Element belongs to.
//------------------------------------------------------------------------------
template <typename F>
class Element : public Rings::CommutativeRings::Element<F>
{
  public:

    virtual F multiplicative_inverse() const = 0; // pure virtual function
};

} // namespace Fields

} // namespace CommutativeRings

} // namespace Rings

#endif // _RINGS_COMMUTATIVE_RINGS_FIELDS_FIELD_H_
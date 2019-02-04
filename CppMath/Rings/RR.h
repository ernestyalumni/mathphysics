//------------------------------------------------------------------------------
/// \file RR.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Real numbers under addition and multiplication as a field.
/// \ref Ch. 21 Class Hierarchies, 21.2.Design of Class Hierarchies
///   The C++ Programming Language, 4th Ed., Stroustrup;
/// \details Denote real numbers as RR, field of real numbers under addition and
/// multiplication (RR, +, *) denoted throughout.
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
#ifndef _RINGS_COMMUTATIVE_RINGS_FIELDS_RR_H_
#define _RINGS_COMMUTATIVE_RINGS_FIELDS_RR_H_

#include "Field.h"

#include <ostream>
#include <stdexcept> // std::invalid_argument
#include <type_traits>

namespace Rings
{

namespace CommutativeRings
{

namespace Fields
{

template <
  class T, typename = std::enable_if_t<std::is_floating_point<T>::value>>
class RR : public Fields::Element<RR<T>>
{
  public:

    RR():
      a_{0.0}
    {}

    explicit RR(const T a):
      a_{a}
    {}

    RR operator+(const RR& b) const
    {
      return RR{a_ + b.a_};
    }

    RR operator*(const RR& b) const
    {
      return RR{a_ * b.a_};
    }

    //--------------------------------------------------------------------------
    /// \fn identity
    /// \brief Additive identity.
    //--------------------------------------------------------------------------
    RR identity() const
    {
      return RR{0.0};
    }

    //--------------------------------------------------------------------------
    /// \fn inverse
    /// \brief Additive inverse.
    //--------------------------------------------------------------------------
    RR inverse() const
    {
      return RR{-a_};
    }

    RR multiplicative_identity() const
    {
      return RR{1.0};
    }
 
    RR multiplicative_inverse() const
    {
      if (a_ == 0.0)
      {
        throw std::invalid_argument(
          "Only non-zero RRs have a multiplicative inverse");
      }

      return RR{1.0 / a_};
    }
 
    // Accessor
    T data() const
    {
      return a_;
    }

   friend std::ostream& operator<<(std::ostream& os, const RR& a)
    {
      os << a.a_ << ' ';
      return os;      
    }

  private:

    T a_;
}; // class RRs

} // namespace Fields

} // namespace CommutativeRings

} // namespace Rings

#endif // _RINGS_COMMUTATIVE_RINGS_FIELDS_RRS_H_
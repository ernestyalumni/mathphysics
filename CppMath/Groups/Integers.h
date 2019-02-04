//------------------------------------------------------------------------------
/// \file Integers.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Integers under addition as a group.
/// \ref Ch. 21 Class Hierarchies, 21.2.Design of Class Hierarchies
///   The C++ Programming Language, 4th Ed., Stroustrup;
/// \details Denote integers as ZZ, group of integers under addition (ZZ, +)
/// denoted throughout.
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
#ifndef _GROUPS_ABELIAN_GROUPS_INTEGERS_H_
#define _GROUPS_ABELIAN_GROUPS_INTEGERS_H_

#include "AbelianGroup.h"

#include <ostream>
#include <type_traits>

namespace Groups
{

namespace AbelianGroups
{


//------------------------------------------------------------------------------
/// \class Integer
/// \details Use CTRP pattern for the return type.
/// Replicated base. If ambiguous, make inheritance public virtual.
/// https://stackoverflow.com/questions/27180342/ \
/// pure-virtual-function-in-abstract-class-with-return-type-of-base-derived-type
/// Sec. 21.3.6 Replicated vs. Virtual Bases, Stroustrup.
//------------------------------------------------------------------------------
template <class T, typename = std::enable_if_t<std::is_integral<T>::value>>
class Integer : public AbelianGroups::Element<Integer<T>>
{
  public:

    Integer():
      a_{0}
    {}

    explicit Integer(const T a):
      a_{a}
    {}

    Integer operator+(const Integer& b) const
    {
      return Integer{a_ + b.a_};
    }

    Integer identity() const
    {
      return Integer{0};
    }

    Integer inverse() const
    {
      return Integer{-a_};
    }

    T data() const
    {
      return a_;
    }

    friend std::ostream& operator<<(std::ostream& os, const Integer& a)
    {
      os << a.a_ << ' ';
      return os;      
    }

  private:

    T a_;
};

} // namespace AbelianGroups

} // namespace Groups

#endif // _GROUPS_ABELIAN_GROUPS_INTEGERS_H_
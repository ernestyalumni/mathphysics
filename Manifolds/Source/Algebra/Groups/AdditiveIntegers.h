//------------------------------------------------------------------------------
/// \file AdditiveIntegers.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Integers under addition as a group.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details Denote integers as ZZ, group of integers under addition (ZZ, +)
/// denoted throughout.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_GROUPS_ADDITIVE_INTEGERS_H
#define ALGEBRA_GROUPS_ADDITIVE_INTEGERS_H

#include "AbelianGroup.h"

#include <ostream>
#include <type_traits>

namespace Algebra
{
namespace Groups
{
namespace AbelianGroups
{

//------------------------------------------------------------------------------
/// \class Integer
/// \details Use CTRP pattern for the return type.
//------------------------------------------------------------------------------
template <class T, typename = std::enable_if_t<std::is_integral<T>::value>>
class AdditiveInteger : public AbelianGroupElement<AdditiveInteger<T>>
{
  public:

    AdditiveInteger():
      a_{0}
    {}

    explicit AdditiveInteger(const T a):
      a_{a}
    {}

    AdditiveInteger group_law(const AdditiveInteger& b) const
    {
      return operator+(b);
    }

    AdditiveInteger operator+(const AdditiveInteger& b) const
    {
      return AdditiveInteger{a_ + b.a_};
    }

    AdditiveInteger identity() const
    {
      return AdditiveInteger{0};
    }

    AdditiveInteger inverse() const
    {
      return AdditiveInteger{-a_};
    }

    T value() const
    {
      return a_;
    }

    friend std::ostream& operator<<(std::ostream& os, const AdditiveInteger& a)
    {
      os << a.a_ << ' ';
      return os;
    }

  private:

    T a_;
};

} // namespace AbelianGroups
} // namespace Groups
} // namespace Algebra

#endif // ALGEBRA_GROUPS_ADDITIVE_INTEGERS_H
//------------------------------------------------------------------------------
/// \file Matrices2x2.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  2x2 matrices as a ring.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details
//------------------------------------------------------------------------------
#ifndef ALGEBRA_RINGS_MATRICES_2X2_H
#define ALGEBRA_RINGS_MATRICES_2X2_H

#include "Ring.h"

#include <array>
#include <ostream>
#include <type_traits>

namespace Algebra
{
namespace Rings
{

//------------------------------------------------------------------------------
/// \class Matrix2x2
/// \details Use CTRP pattern for the return type.
//------------------------------------------------------------------------------
template <class Field>
class Matrix2x2 : public RingElement<Matrix2x2<Field>>
{
  public:

    Matrix2x2(const Field a, const Field b, const Field c, const Field d):
      elements_{a, b, c, d}
    {}

    Matrix2x2 addition(const Matrix2x2& b) const
    {
      return operator+(b);
    }

    Matrix2x2 operator+(const Matrix2x2& b) const
    {
      return Matrix2x2 {
        elements_.at(0) + b.elements_.at(0),
        elements_.at(1) + b.elements_.at(1),
        elements_.at(2) + b.elements_.at(2),
        elements_.at(3) + b.elements_.at(3)};
    }

    Matrix2x2 additive_identity() const
    {
      return Matrix2x2{0, 0, 0, 0};
    }

    Matrix2x2 additive_inverse() const
    {
      return Matrix2x2{
        -elements_.at(0), -elements_.at(1), -elements_.at(2), -elements_.at(3)};
    }

    std::array<Field, 4> data() const
    {
      return elements_;
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix2x2& a)
    {
      os << a.elements_.at(0) << ' ' << a.elements_.at(1) << "\n";
      os << a.elements_.at(2) << ' ' << a.elements_.at(3);
      return os;
    }

    Matrix2x2 multiplication(const Matrix2x2& b) const
    {
      return operator*(b);
    }

    Matrix2x2 operator*(const Matrix2x2& b) const
    {
      return Matrix2x2 {
        elements_.at(0) * b.elements_.at(0) + 
          elements_.at(1) * b.elements_.at(2),
        elements_.at(0) * b.elements_.at(1) + 
          elements_.at(1) * b.elements_.at(3),
        elements_.at(2) * b.elements_.at(0) + 
          elements_.at(3) * b.elements_.at(2),
        elements_.at(2) * b.elements_.at(1) + 
          elements_.at(3) * b.elements_.at(3)};
    }

  private:

    std::array<Field, 4> elements_;
};

} // namespace Rings
} // namespace Algebra

#endif // ALGEBRA_RINGS_MATRICES_2X2_H
//------------------------------------------------------------------------------
/// \file Field.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Field, as a special case of a commutative ring.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
//------------------------------------------------------------------------------
#ifndef ALGEBRA_RINGS_FIELD_H
#define ALGEBRA_RINGS_FIELD_H

#include "Ring.h"

namespace Algebra
{
namespace Rings
{
namespace Fields
{

//------------------------------------------------------------------------------
/// \class Element
/// \details A base class for a group Element. Also, use CRTP pattern, for the 
/// return type.
/// \ref https://en.wikipedia.org/wiki/Ring_(mathematics)
/// \param F stands for field F, the field F that FieldElement belongs to.
//------------------------------------------------------------------------------
template <typename F>
class FieldElement : public RingElement<FieldElement<F>>
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    // F is an abelian group under addition.

    F addition(const F& b) const
    {
      return object()->addition(b);
    }

    // Always denoted mathematically as 0.
    F additive_identity() const
    {
      return object()->additive_identity();
    }

    F additive_inverse() const
    {
      return object()->additive_inverse();
    }

    /*
    R multiplication(const R& b) const
    {
      return object()->multiplication(b);
    }
    */

  private:

    F& object()
    {
      return static_cast<F&>(*this);
    }
};

} // namespace Fields
} // namespace Rings
} // namespace Algebra

#endif // ALGEBRA_RINGS_FIELD_H



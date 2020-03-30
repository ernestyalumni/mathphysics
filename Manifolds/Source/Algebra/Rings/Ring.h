//------------------------------------------------------------------------------
/// \file Rings.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Ring, and ring axioms.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details Ring (2 ring operations, ring axioms), following CTRP for static
/// polymorphism.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_RINGS_RING_H
#define ALGEBRA_RINGS_RING_H

#include "Algebra/Groups/AbelianGroup.h"
//#include "SemiGroups/"

namespace Algebra
{
namespace Rings
{

//------------------------------------------------------------------------------
/// \class Element
/// \details A base class for a group Element. Also, use CRTP pattern, for the 
/// return type.
/// \ref https://en.wikipedia.org/wiki/Ring_(mathematics)
/// \param R stands for ring R, the ring R that RingElement belongs to.
//------------------------------------------------------------------------------
template <typename R>
class RingElement :
  protected Algebra::Groups::AbelianGroups::AbelianGroupElement<RingElement<R>>
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    // R is an abelian group under addition.

    R addition(const R& b) const
    {
      return object()->addition(b);
    }

    // Always denoted mathematically as 0.
    R additive_identity() const
    {
      return object()->additive_identity();
    }

    R additive_inverse() const
    {
      return object()->additive_inverse();
    }

    /*
    R multiplication(const R& b) const
    {
      return object()->multiplication(b);
    }
    */

  protected:

    R group_law(const R& b) const
    {
      return addition(b);
    }

    R inverse() const
    {
      return additive_inverse();
    }

  private:

    R& object()
    {
      return static_cast<R&>(*this);
    }
};

} // namespace Rings
} // namespace Algebra

#endif // ALGEBRA_RINGS_RING_H



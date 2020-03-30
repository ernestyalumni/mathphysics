//------------------------------------------------------------------------------
/// \file Monoid.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Monoid has single associative binary operation and identity.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details Monoid, following CTRP static polymorphism.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_SEMIGROUPS_MONOID_H
#define ALGEBRA_SEMIGROUPS_MONOID_H

#include "SemiGroup.h"

namespace Algebra
{
namespace SemiGroups
{
namespace Monoids
{

//------------------------------------------------------------------------------
/// \class MonoidElement
/// \details A base class for a monoid Element. Also, use CRTP pattern, for the 
/// return type.
/// https://en.wikipedia.org/wiki/Monoid
/// \param AG stands for abelian group AG, the group AG that Element belongs to.
//------------------------------------------------------------------------------
template <typename M>
class MonoidElement : public SemiGroupElement<MonoidElement<M>>
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    M dot(const M& b) const
    {
      return object()->dot(b);
    }

    M inverse() const
    {
      return object()->inverse();
    }

  private:

    M& object()
    {
      return static_cast<M&>(*this);
    }
};

} // namespace Monoids
} // namespace SemiGroups
} // namespace Algebra

#endif // ALGEBRA_SEMIGROUPS_MONOID_H



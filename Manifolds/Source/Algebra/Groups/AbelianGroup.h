//------------------------------------------------------------------------------
/// \file AbelianGroup.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Abelian group axioms closure, identity, and inverse elements.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details Abelian group (group operations, group axioms), following CTRP
/// static polymorphism. Does operator overload of +.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_GROUPS_ABELIAN_GROUP_H
#define ALGEBRA_GROUPS_ABELIAN_GROUP_H

#include "Group.h"

#include <cassert>

namespace Algebra
{
namespace Groups
{
namespace AbelianGroups
{


//------------------------------------------------------------------------------
/// \class Element
/// \details A base class for a group Element. Also, use CRTP pattern, for the 
/// return type.
/// \ref https://en.wikipedia.org/wiki/Group_(mathematics)
/// \param AG stands for abelian group AG, the group AG that Element belongs to.
//------------------------------------------------------------------------------
template <typename AG>
class AbelianGroupElement : public GroupElement<AbelianGroupElement<AG>>
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    AG group_law(const AG& b) const
    {
      // TODO: Decide if this should be an assert or throw.
      assert(b.group_law(object()) == object()->group_law(b));

      return object()->group_law(b);
    }

    AG identity() const
    {
      return object()->identity();
    }

    AG inverse() const
    {
      return object()->inverse();
    }

  private:

    AG& object()
    {
      return static_cast<AG&>(*this);
    }
};

} // namespace AbelianGroups
} // namespace Groups
} // namespace Algebra

#endif // ALGEBRA_GROUPS_ABELIAN_GROUP_H



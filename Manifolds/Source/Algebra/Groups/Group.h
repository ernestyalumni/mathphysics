//------------------------------------------------------------------------------
/// \file Group.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Group axioms closure, identity, and inverse elements.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details group (group operations, group axioms), following CTRP for static
/// polymorphism.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_GROUPS_GROUP_H
#define ALGEBRA_GROUPS_GROUP_H

namespace Algebra
{
namespace Groups
{

//------------------------------------------------------------------------------
/// \class GroupElement
/// \details A base class for a group Element. Also, use CRTP pattern, for the 
/// return type.
/// \ref https://en.wikipedia.org/wiki/Group_(mathematics)
/// \param G stands for group G, the group G that GroupElement belongs to.
//------------------------------------------------------------------------------
template <typename G>
class GroupElement
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    G group_law(const G& b) const
    {
      return object()->group_law(b);
    }

    G identity() const
    {
      return object()->identity();
    }

    G inverse() const
    {
      return object()->inverse();
    }

  private:

    G& object()
    {
      return static_cast<G&>(*this);
    }
};

} // namespace Groups
} // namespace Algebra

#endif // ALGEBRA_GROUPS_GROUP_H



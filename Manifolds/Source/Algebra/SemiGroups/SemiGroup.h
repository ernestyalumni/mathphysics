//------------------------------------------------------------------------------
/// \file SemiGroup.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Semigroup, a set with associative binary operation.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details group (group operations, group axioms), following CTRP for static
/// polymorphism.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_SEMIGROUPS_SEMIGROUP_H
#define ALGEBRA_SEMIGROUPS_SEMIGROUP_H

namespace Algebra
{
namespace SemiGroups
{

//------------------------------------------------------------------------------
/// \class SemiGroupElement
/// \details A base class for a semigroup Element. Also, use CRTP pattern, for 
/// the return type.
/// \ref https://en.wikipedia.org/wiki/Semigroup
/// \param SG stands for semigroup SG, the semigroup SG that SemiGroupElement
/// belongs to.
//------------------------------------------------------------------------------
template <typename SG>
class SemiGroupElement
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    SG dot(const SG& b) const
    {
      return object()->dot(b);
    }

  private:

    SG& object()
    {
      return static_cast<SG&>(*this);
    }
};

} // namespace SemiGroups
} // namespace Algebra

#endif // ALGEBRA_SEMIGROUPS_SEMIGROUP_H

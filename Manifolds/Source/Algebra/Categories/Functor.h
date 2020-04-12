//------------------------------------------------------------------------------
/// \file Functor.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Functors.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details group (group operations, group axioms), following CTRP for static
/// polymorphism.
//------------------------------------------------------------------------------
#ifndef CATEGORIES_FUNCTORS_FUNCTOR_H
#define CATEGORIES_FUNCTORS_FUNCTOR_H

namespace Categories
{
namespace Functors
{

namespace Details
{

// F : Obj(C) \to Obj(D) 
// F : X \to Y; F(X) \equiv FX
template <typename X, typename FX, FX ObjectMap(X)>
FX object_map(X x)
{
  return ObjectMap(x);
}

} // namespace Details

} // namespace Functors
} // namespace Categories

#endif // CATEGORIES_FUNCTORS_FUNCTOR_H
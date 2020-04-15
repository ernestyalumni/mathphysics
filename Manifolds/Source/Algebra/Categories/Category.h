//------------------------------------------------------------------------------
/// \file Category.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Categories of Category Theory.
/// \ref https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
/// \details group (group operations, group axioms), following CTRP for static
/// polymorphism.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_CATEGORIES_CATEGORY_H
#define ALGEBRA_CATEGORIES_CATEGORY_H


namespace Algebra
{
namespace Categories
{

namespace Category
{

template <typename Object>
Object& identity_morphism(Object& A)
{
  return A;
}

namespace Details
{

// X \in Obj(\mathbf{C}), X is an element of class Obj(\mathbf{C}) of category
// \mathbf{C}

// TODO: Clean up the following.

template <typename X>
struct ObjectElement;

template <typename X>
struct ObjectElement
{
  X identity(X& x) const
  {
    return x;
  }
};

template <typename X, typename Y, typename fMap>
struct Morphism;

template <typename X, typename Y, typename fMap>
struct Morphism
{
  Y operator()(const X& x)
  {
    return fMap(x);
  }

  //template <
  //  typename X,
  //  typename Y,
  //  typename Z,
  //  typename HomXY,
  //  typename HomYZ,
  //  typename HomXZ>
  //friend Z composition()
};

} // namespace Details

} // namespace Category
} // namespace Categories
} // namespace Algebra

#endif // ALGEBRA_CATEGORIES_CATEGORY_H
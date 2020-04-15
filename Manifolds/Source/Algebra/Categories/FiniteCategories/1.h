//------------------------------------------------------------------------------
/// \file 1.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Category with one object and 1 (identity) arrow
/// \ref Saunders Mac Lane. Categories for the Working Mathematician. pp. 10.
/// \details group (group operations, group axioms), following CTRP for static
/// polymorphism.
//------------------------------------------------------------------------------
#ifndef ALGEBRA_CATEGORIES_FINITE_CATEGORIES_1_H
#define ALGEBRA_CATEGORIES_FINITE_CATEGORIES_1_H

#include <string>
#include <type_traits> // std::underlying_type

namespace Algebra
{
namespace Categories
{
namespace FiniteCategories
{

class Category1Object
{
  public:

    Category1Object() = default;

    std::string name() const
    {
      return name_;
    }

  private:

    const std::string name_ {"A"};
};

class Objects1Element
{
  public:

    enum class ObjectName : char
    {
      A = 'A'
    };

    explicit Objects1Element(const ObjectName object_name) :
      name_{object_name}
    {}

    char name() const
    {
      return static_cast<std::underlying_type_t<ObjectName>>(name_);
    }

  private:

    ObjectName name_;
};

} // namespace FiniteCategories
} // namespace Categories
} // namespace Algebra

#endif // ALGEBRA_CATEGORIES_FINITE_CATEGORIES_FINITE_CATEGORIES_H
//------------------------------------------------------------------------------
// \file 1_tests.cpp
//------------------------------------------------------------------------------
#include "Algebra/Categories/Category.h"
#include "Algebra/Categories/FiniteCategories/1.h"

#include <boost/test/unit_test.hpp>
#include <cmath>

using Algebra::Categories::Category::identity_morphism;
using Algebra::Categories::FiniteCategories::Category1Object;
using Algebra::Categories::FiniteCategories::Objects1Element;

BOOST_AUTO_TEST_SUITE(Categories)
BOOST_AUTO_TEST_SUITE(FiniteCategories)
BOOST_AUTO_TEST_SUITE(FiniteCategories_tests)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Category1ObjectAsItsOwnClass)
{
  const Category1Object A;

  BOOST_TEST(A.name() == "A");
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Category1ObjectElementAsAClass)
{
  const Objects1Element A {Objects1Element::ObjectName::A};

  BOOST_TEST(A.name() == 'A');
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(IdentityMorphismWorksOnCategory1Object)
{
  const Category1Object A;

  BOOST_TEST_REQUIRE(A.name() == "A");

  const Category1Object& another_A {identity_morphism(A)};

  BOOST_TEST(another_A.name() == "A");
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(IdentityMorphismWorksOnCategory1ObjectElement)
{
  const Objects1Element A {Objects1Element::ObjectName::A};

  BOOST_TEST_REQUIRE(A.name() == 'A');

  const auto& another_A = identity_morphism(A);

  BOOST_TEST(another_A.name() == 'A');
}

BOOST_AUTO_TEST_SUITE_END() // FiniteCategories_tests
BOOST_AUTO_TEST_SUITE_END() // FiniteCategories
BOOST_AUTO_TEST_SUITE_END() // Categories
//------------------------------------------------------------------------------
/// \file AdditiveIntegers_tests.cpp
/// \ref https://www.cs.utexas.edu/users/fussell/courses/cs429h/lectures/Lecture_2-429h.pdf
//------------------------------------------------------------------------------
#include "Algebra/Groups/AdditiveIntegers.h"

#include <boost/test/unit_test.hpp>
#include <cmath>

using Algebra::Groups::AbelianGroups::AdditiveInteger;

BOOST_AUTO_TEST_SUITE(Algebra)
BOOST_AUTO_TEST_SUITE(Groups)
BOOST_AUTO_TEST_SUITE(AbelianGroups)
BOOST_AUTO_TEST_SUITE(AdditiveIntegers_tests)

BOOST_AUTO_TEST_SUITE(Interface)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(IdentityReturnsZero)
{
  const AdditiveInteger<int> additive_integer_int {2};
  BOOST_TEST(additive_integer_int.identity().value() == 0);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(IdentityReturnsAdditiveInverse)
{
  const AdditiveInteger<int> additive_integer_int {2};
  BOOST_TEST(additive_integer_int.inverse().value() == -2);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GroupLawInheritedAndDoesAbelianGroupAddition)
{
  const AdditiveInteger<int> a {-3};
  const AdditiveInteger<int> b {5};
  BOOST_TEST(a.group_law(b).value() == 2);
  BOOST_TEST(b.group_law(a).value() == 2);
}

BOOST_AUTO_TEST_SUITE_END() // Interface

BOOST_AUTO_TEST_SUITE(Operations)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GroupLawInheritedAndDoesAbelianGroupAddition)
{
  const AdditiveInteger<int> a {1};
  const AdditiveInteger<int> b {2};
  BOOST_TEST((a + b).value() == 3);
}

BOOST_AUTO_TEST_SUITE_END() // Operations

BOOST_AUTO_TEST_SUITE_END() // AdditiveIntegers_tests
BOOST_AUTO_TEST_SUITE_END() // AbelianGroups
BOOST_AUTO_TEST_SUITE_END() // Groups
BOOST_AUTO_TEST_SUITE_END() // Algebra
//------------------------------------------------------------------------------
// \file Rings_tests.cpp
// \ref https://www.cs.utexas.edu/users/fussell/courses/cs429h/lectures/Lecture_2-429h.pdf
//------------------------------------------------------------------------------
#include "Algebra/Rings/Matrices2x2.h"

#include <boost/test/unit_test.hpp>
#include <cmath>

using Algebra::Rings::Matrix2x2;

BOOST_AUTO_TEST_SUITE(Algebra)
BOOST_AUTO_TEST_SUITE(Rings)
BOOST_AUTO_TEST_SUITE(Rings_tests)

BOOST_AUTO_TEST_SUITE(Interface)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(AdditiveIdentityReturnsZero)
{
  BOOST_TEST(false);
}

BOOST_AUTO_TEST_SUITE_END() // Interface

BOOST_AUTO_TEST_SUITE(Operations)


BOOST_AUTO_TEST_SUITE_END() // Operations

BOOST_AUTO_TEST_SUITE_END() // Rings_tests
BOOST_AUTO_TEST_SUITE_END() // Rings
BOOST_AUTO_TEST_SUITE_END() // Algebra
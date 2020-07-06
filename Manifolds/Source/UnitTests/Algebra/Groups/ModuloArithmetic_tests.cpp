//------------------------------------------------------------------------------
/// \file ModuloArithmetic_tests.cpp
/// \ref https://en.cppreference.com/w/cpp/language/operator_arithmetic
/// \brief Demonstrate Modulo arithmetic.
//------------------------------------------------------------------------------
#include "Utilities/SuperBitSet.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Algebra)
BOOST_AUTO_TEST_SUITE(Groups)
BOOST_AUTO_TEST_SUITE(AbelianGroups)
BOOST_AUTO_TEST_SUITE(ModuloArithmetic_tests)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(IdentityReturnsZero)
{
  BOOST_TEST(true);
}

BOOST_AUTO_TEST_SUITE_END() // ModuloArithmetic_tests
BOOST_AUTO_TEST_SUITE_END() // AbelianGroups
BOOST_AUTO_TEST_SUITE_END() // Groups
BOOST_AUTO_TEST_SUITE_END() // Algebra
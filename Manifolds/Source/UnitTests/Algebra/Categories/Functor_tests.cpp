//------------------------------------------------------------------------------
// \file Functor_tests.cpp
//------------------------------------------------------------------------------
#include "Algebra/Categories/Functor.h"

#include <boost/test/unit_test.hpp>
#include <cmath>

using Categories::Functors::Details::object_map;
using std::cos;
using std::sin;
using std::sqrt;

BOOST_AUTO_TEST_SUITE(Categories)
BOOST_AUTO_TEST_SUITE(Functors)
BOOST_AUTO_TEST_SUITE(Functors_tests)

BOOST_AUTO_TEST_SUITE(Details)

BOOST_AUTO_TEST_SUITE(ObjectMap)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(MapAsTemplateParameter)
{
  const auto result = object_map<double, double, &sin>(M_PI_4);
  BOOST_TEST(result == 1.0 / sqrt(2.0));
}

BOOST_AUTO_TEST_SUITE_END() // ObjectMap

BOOST_AUTO_TEST_SUITE_END() // Details

BOOST_AUTO_TEST_SUITE_END() // Functors_tests
BOOST_AUTO_TEST_SUITE_END() // Functors
BOOST_AUTO_TEST_SUITE_END() // Categories
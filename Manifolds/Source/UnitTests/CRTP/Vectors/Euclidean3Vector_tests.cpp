#include "CRTP/Vectors/Euclidean3Vector.h"

#include <boost/test/unit_test.hpp>

using CRTP::Vectors::Euclidean3Vector;

BOOST_AUTO_TEST_SUITE(CRTP)
BOOST_AUTO_TEST_SUITE(Vectors)

BOOST_AUTO_TEST_SUITE(Euclidean3Vectors_tests)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ConstructsFromInitializerList)
{
  Euclidean3Vector<double> x {30, 40, 50};

  BOOST_TEST(true);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(AccessorsRetrieveData)
{
  Euclidean3Vector<double> x {30, 40, 50};

  BOOST_TEST(x.x() == 30);
  BOOST_TEST(x.y() == 40);
  BOOST_TEST(x.z() == 50);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(DimensionRetrievesDimension)
{
  Euclidean3Vector<double> x;

  BOOST_TEST(x.dimension() == 3);
}

BOOST_AUTO_TEST_SUITE_END() // Euclidean3Vectors_tests
BOOST_AUTO_TEST_SUITE_END() // Vectors
BOOST_AUTO_TEST_SUITE_END() // CRTP
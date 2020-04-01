//------------------------------------------------------------------------------
// \file Rings_tests.cpp
// \ref https://www.cs.utexas.edu/users/fussell/courses/cs429h/lectures/Lecture_2-429h.pdf
//------------------------------------------------------------------------------
#include "Algebra/Rings/Matrices2x2.h"

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <sstream>

using Algebra::Rings::Matrix2x2;

BOOST_AUTO_TEST_SUITE(Algebra)
BOOST_AUTO_TEST_SUITE(Rings)
BOOST_AUTO_TEST_SUITE(Rings_tests)

BOOST_AUTO_TEST_SUITE(Matrices2x2)

BOOST_AUTO_TEST_SUITE(Interface)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ConstructsFrom4RealParameters)
{
	{
		const Matrix2x2<int> mat {1, 2, 3, 4};

		const auto result {mat.data()};

		for (int i {0}; i < 4; ++i)
		{
			BOOST_TEST(result.at(i) == (i + 1));
		}
	}
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(OperatorInsertsIntoStringStream)
{
	{
		std::ostringstream oss;

		const Matrix2x2<double> mat {-5.6, 6.2, 7.1, -8.0};

		oss << mat;

		BOOST_TEST(oss.str() == "-5.6 6.2\n7.1 -8");

		oss.flush();
	}
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(AdditiveIdentityReturnsZero)
{
  BOOST_TEST(true);
}

BOOST_AUTO_TEST_SUITE_END() // Interface

BOOST_AUTO_TEST_SUITE(Operations)


BOOST_AUTO_TEST_SUITE_END() // Operations

BOOST_AUTO_TEST_SUITE_END() // Matrices2x2

BOOST_AUTO_TEST_SUITE_END() // Rings_tests
BOOST_AUTO_TEST_SUITE_END() // Rings
BOOST_AUTO_TEST_SUITE_END() // Algebra
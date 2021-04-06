#include "CRTP/Categories/Sets/PosetElement.h"

#include <boost/test/unit_test.hpp>
#include <type_traits>

using CRTP::Categories::Sets::Posets::BoostPartiallyOrderedElement;
using CRTP::Categories::Sets::Posets::PartiallyOrderedSetElement;

BOOST_AUTO_TEST_SUITE(CRTP)
BOOST_AUTO_TEST_SUITE(Categories)
BOOST_AUTO_TEST_SUITE(Sets)

BOOST_AUTO_TEST_SUITE(PosetElement_tests)

BOOST_AUTO_TEST_SUITE(BoostPartiallyOrderedElement_tests)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Reflexivity)
{
  {
    //BoostPartiallyOrderedElement<unsigned int> a {97};
    //BoostPartiallyOrderedElement<unsigned int> b {97};

    //BOOST_TEST((a <= b));
  }
  BOOST_TEST(true);
}

BOOST_AUTO_TEST_SUITE_END() // BoostPartiallyOrderedElement_tests

BOOST_AUTO_TEST_SUITE(PartiallyOrderedSetElement_tests)

template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
class PartiallyOrderedNaturalNumber :
  PartiallyOrderedSetElement<PartiallyOrderedNaturalNumber<T>>
{
  public:

    PartiallyOrderedNaturalNumber():
      data_{static_cast<T>(0)}
    {}

    PartiallyOrderedNaturalNumber(const T n):
      data_{n}
    {}

    bool operator<=(const PartiallyOrderedNaturalNumber& rhs) const
    {
      return this->data_ <= rhs.data_;
    }

    /*
    friend bool operator<=(
      const PartiallyOrderedNaturalNumber& lhs,
      const PartiallyOrderedNaturalNumber& rhs)
    {
      return lhs.data_ <= rhs.data_
    }
    */

  private:

    T data_;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Reflexivity)
{
  PartiallyOrderedNaturalNumber<unsigned int> a {97};
  PartiallyOrderedNaturalNumber<unsigned int> b {97};

  BOOST_TEST((a <= b));
  BOOST_TEST((a <= a));
  BOOST_TEST((b <= b));
}

BOOST_AUTO_TEST_SUITE_END() // PartiallyOrderedSetElement_tests

BOOST_AUTO_TEST_SUITE_END() // PosetElement_tests
BOOST_AUTO_TEST_SUITE_END() // Sets
BOOST_AUTO_TEST_SUITE_END() // Categories
BOOST_AUTO_TEST_SUITE_END() // CRTP
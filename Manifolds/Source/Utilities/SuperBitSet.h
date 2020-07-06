//------------------------------------------------------------------------------
/// \file SuperBitSet.h
/// \author Ernest Yeung
/// \brief std::bitset extension.
/// \ref https://en.cppreference.com/w/cpp/utility/bitset
///-----------------------------------------------------------------------------
#ifndef CPP_UTILITIES_SUPER_BIT_SET_H
#define CPP_UTILITIES_SUPER_BIT_SET_H

#include <bitset>
// cf. https://en.cppreference.com/w/cpp/types/climits
#include <climits> // CHAR_BIT
#include <utility> // std::move

namespace Cpp
{
namespace Utilities
{

namespace
{

//constexpr std::size_t number_of_bits_in_a_byte = 8;
constexpr std::size_t number_of_bits_in_a_byte = CHAR_BIT;

} // anonymous namespace

// N = number of bits
template <std::size_t N = sizeof(unsigned long long) * number_of_bits_in_a_byte>
class SuperBitSet : public std::bitset<N>
{
  public:
    // cf. cppreference.com
    // Includes the following constructors from std::bitset<N>::bitset
    //
    // constexpr bitset() noexcept;
    // Default ctor; constructs a bitset with all bits set to 0
    //
    // constexpr bitset(unsigned long long val) noexcept;
    // Constructs bitset, initializing 1st (rightmost, least significant) M bit
    // positions to corresponding bit values of val, where M is smaller of number
    // of bits in unsigned long long and number of bits N in bitset being
    // constructed.
    // 
    // explicit bitset(const std::basic_string<CharT,Traits,Alloc>& str,
    //   typename std::basic_string<CharT,Traits,Alloc>::size_type pos = 0);
    // Constructs bitset using characters in std::basic_string str. An optional
    // starting position pos and length n can be provided, as well as characters
    // denoting alternative values for set (one) and unset (0) bits.
    using std::bitset<N>::bitset;

    SuperBitSet operator=(const SuperBitSet& bit_set)
    {
      return SuperBitSet{bit_set.to_string()};
    }

    SuperBitSet operator=(const std::bitset<N>& bit_set)
    {
      return std::move(SuperBitSet{bit_set.to_string()});
    }

    // TODO: Consider adding Copy construction.

    // Bit inversion.
    //SuperBitSet operator~() const
    //{
    //  const std::bitset<N>& bits {*this};

    //  return ~bits;
    //}

    // Insert bits to the "front", 0th position
    void insert_at_0(const bool value)
    {
      *this <<= 1;
      (*this)[0] = value;
    }

    // Insert multiple 0 bits at "front"
    void insert_zeroes_from_0(const std::size_t number_of_zeroes)
    {
      *this <<= number_of_zeroes;
    }

    // Remove bit(s) starting from "front"
    void remove_from_0(const std::size_t number_of_bits)
    {
      *this >>= number_of_bits;
    }
};

} // namespace Utilities
} // namespace Cpp

#endif // CPP_UTILITIES_SUPER_BIT_SET_H
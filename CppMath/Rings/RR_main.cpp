//------------------------------------------------------------------------------
/// \file RR_main.cpp
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Integers under addition as a group.
/// \ref Ch. 21 Class Hierarchies, 21.2.Design of Class Hierarchies
///   The C++ Programming Language, 4th Ed., Stroustrup;
/// \details Denote integers as ZZ, group of integers under addition (ZZ, +)
/// denoted throughout.
/// \copyright If you find this code useful, feel free to donate directly
/// (username ernestyalumni or email address above), going directly to:
///
/// paypal.me/ernestyalumni
///
/// which won't go through a 3rd. party like indiegogo, kickstarter, patreon.
/// Otherwise, I receive emails and messages on how all my (free) material on
/// physics, math, and engineering have helped students with their studies, and
/// I know what it's like to not have money as a student, but love physics (or
/// math, sciences, etc.), so I am committed to keeping all my material
/// open-source and free, whether or not sufficiently crowdfunded, under the
/// open-source MIT license: feel free to copy, edit, paste, make your own
/// versions, share, use as you wish.
/// Peace out, never give up! -EY
//------------------------------------------------------------------------------
/// COMPILATION TIPS:
///  g++ --std=c++17 -I ../ RR_main.cpp -o RR_main
//------------------------------------------------------------------------------
#include "RR.h"

#include <iostream>

using Rings::CommutativeRings::Fields::RR;

int main()
{
  // RRDefaultConstructsWithDouble
  {
    std::cout << "\n RRDefaultConstructsWithDouble \n";
    const RR<double> x;
    std::cout << " x default constructs to 0 : " << x.data() << ' ' << 
      (x.data() == 0.0) << '\n';
  }

  // InterfaceTests
  std::cout << "\n InterfaceTests \n";

  // MultiplicativeInverseReturnsMultiplicativeInverseForNonzeroNumbers
  {
    std::cout <<
      "\n MultiplicativeInverseReturnsMultiplicativeInverseForNonzeroNumbers\n";

    const RR<double> x {-5.0};
    std::cout << x.multiplicative_inverse() << '\n';
  }

  // Operator*OverloadedToRingMultiplication
  {
    std::cout << "\n Operator*OverloadedToRingMultiplication \n";
    const RR<double> x {-4.0};
    const RR<double> y {6.0};
    std::cout << " x * y = " << (x * y).data() << '\n';
  } 
}
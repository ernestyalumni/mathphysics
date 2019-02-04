//------------------------------------------------------------------------------
/// \file Integers_main.cpp
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
///  g++ --std=c++17 Integers_main.cpp -o Integers_main
//------------------------------------------------------------------------------
#include "Integers.h"

#include <iostream>

using Groups::AbelianGroups::Integer;

int main()
{
  // IntegerDefaultConstructsWithInt
  {
    std::cout << "\n IntegerDefaultConstructsWithInt \n";
    Integer<int> integer_int;
  }

  // InterfaceTests
  std::cout << "\n InterfaceTests \n";

  // IdentityReturnsZero
  {
    std::cout << "\n IdentityReturnsZero \n";
    const Integer<int> integer_int {2};
    std::cout << integer_int.identity() << '\n';
  }

  // InverseReturnsAdditiveInverse
  {
    std::cout << "\n IdentityReturnsAdditveInverse \n";
    const Integer<int> integer_int {2};
    std::cout << integer_int.inverse() << '\n';
  }

  // GroupLawInheritedAndDoesAbelianGroupAddition
  {
    std::cout << "\n GroupLawInheritedAndDoesAbelianGroupAddition \n";
    const Integer<int> a {-3};
    const Integer<int> b {5};
    std::cout << a.group_law(b) << '\n';
  }

  // OperationsTests

  std::cout << "\n OperationTests \n";

  // Operator+OverloadedToGroupAddition
  {
    std::cout << "\n Operator+OverloadedToGroupAddition \n";
    const Integer<int> a {1};
    const Integer<int> b {2};
    std::cout << " a + b  = " << (a + b).data() << '\n';
  }

}
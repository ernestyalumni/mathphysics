//------------------------------------------------------------------------------
/// \file AbelianGroup.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Abelian group axioms closure, identity, and inverse elements.
/// \ref Ch. 21 Class Hierarchies, 21.2.Design of Class Hierarchies
///   The C++ Programming Language, 4th Ed., Stroustrup;
/// \details Abelian group (group operations, group axioms), following interface
/// implementation. Does operator overload of +.
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
///  g++ -std=c++17 Tuple2_main.cpp -o Tuple2_main
//------------------------------------------------------------------------------
#ifndef _GROUPS_ABELIAN_GROUPS_ABELIAN_GROUP_H_
#define _GROUPS_ABELIAN_GROUPS_ABELIAN_GROUP_H_

#include "Group.h"

namespace Groups
{

namespace AbelianGroups
{

//------------------------------------------------------------------------------
/// \brief Element
/// \details A pure abstract base class for a group Element. Also, use CTRP
/// pattern for the return type.
/// \ref https://en.wikipedia.org/wiki/Group_(mathematics)
/// https://stackoverflow.com/questions/27180342/ \
/// pure-virtual-function-in-abstract-class-with-return-type-of-base-derived-type
/// \tparam AG stands for abelian group AG, the group AG that Element belongs
/// to.
//------------------------------------------------------------------------------
template <typename AG>
class Element : public Groups::Element<AG>
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    // pure virtual function
    AG group_law(const AG& b) const
    {
      return operator+(b);
    }
    
    virtual AG operator+(const AG& b) const = 0;
};

} // namespace AbelianGroups

} // namespace Groups

#endif // _GROUPS_ABELIAN_GROUPS_ABELIAN_GROUP_H_
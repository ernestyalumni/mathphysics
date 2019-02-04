//------------------------------------------------------------------------------
/// \file Group.h
/// \author Ernest Yeung
/// \email  ernestyalumni@gmail.com
/// \brief  Group axioms closure, identity, and inverse elements.
/// \ref Ch. 21 Class Hierarchies, 21.2.Design of Class Hierarchies
///   The C++ Programming Language, 4th Ed., Stroustrup;
/// \details group (group operations, group axioms), following interface
/// implementation.
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
#ifndef _GROUPS_GROUP_H_
#define _GROUPS_GROUP_H_

namespace Groups
{

//------------------------------------------------------------------------------
/// \class Element
/// \details A pure abstract base class for a group Element. Also,
/// use CRTP pattern, for the return type.
/// \ref https://en.wikipedia.org/wiki/Group_(mathematics)
/// https://stackoverflow.com/questions/27180342/ \
/// pure-virtual-function-in-abstract-class-with-return-type-of-base-derived-type
/// \tparam G stands for group G, the group G that Element belongs to.
//------------------------------------------------------------------------------
template <typename G>
class Element
{
  public:

    // Data is gone; ctors gone since there's no data to initialize.    

    // pure virtual function
    virtual G group_law(const G& b) const = 0;

    virtual G identity() const = 0; // pure virtual function

    virtual G inverse() const = 0;    

    //--------------------------------------------------------------------------
    /// \fn (virtual) destructor
    /// \ref Sec. 21.2.2 Interface Inheritance of Ch. 21 Class Hierarchies,
    /// Stroustrup
    /// \details Ensure proper cleanup by defining virtual destructor in base
    /// and overriding it suitably in derived classes.
    //--------------------------------------------------------------------------   
    virtual ~Element()
    {}
};

} // namespace Groups

#endif // _GROUPS_GROUP_H_



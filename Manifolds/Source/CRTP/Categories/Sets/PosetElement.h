#ifndef CRTP_CATEGORIES_SETS_POSET_ELEMENT
#define CRTP_CATEGORIES_SETS_POSET_ELEMENT

#include <boost/operators.hpp>

namespace CRTP
{

namespace Categories
{
namespace Sets
{
namespace Posets
{

//------------------------------------------------------------------------------
/// \class PartiallyOrderedSetElement
/// \details (P, <=)
//------------------------------------------------------------------------------
template <typename Implementation>
class PartiallyOrderedSetElement
{
  public:

    bool operator<=(const Implementation& rhs) const
    {
      return object()->operator<=(rhs);
    }

  private:

    Implementation& object()
    {
      return static_cast<Implementation&>(*this);
    }
};

template <typename SetElement>
class BoostPartiallyOrderedElement :
  private boost::partially_ordered<BoostPartiallyOrderedElement<SetElement>>
{
  public:

    bool operator<=(const SetElement& rhs) const
    {
      return object()->operator<=(rhs);
    }

  private:

    SetElement& object()
    {
      return static_cast<SetElement&>(*this);
    }
};

} // namespace Posets
} // namespace Sets
} // namespace Categories
  
} // namespace CRTP

#endif // CRTP_CATEGORIES_SETS_POSET_ELEMENT
#ifndef CRTP_VECTORS_EUCLIDEAN_VECTOR_H
#define CRTP_VECTORS_EUCLIDEAN_VECTOR_H

#include "Vector.h"

#include <cstddef>

namespace CRTP
{

namespace Vectors
{

template <typename Implementation, std::size_t N, typename Field>
class EuclideanVector :
  public Vector<EuclideanVector<Implementation, N, Field>, Field>
{
  public:

    //--------------------------------------------------------------------------
    /// \brief Length
    //--------------------------------------------------------------------------
    const Field length() const
    {
      return object()->length();
    }

    const std::size_t dimension() const
    {
      return N;
    }

  private:

    Implementation& object()
    {
      return static_cast<Implementation&>(*this);
    }
};

} // namespace Vectors

} // namespace CRTP

#endif // CRTP_VECTORS_EUCLIDEAN_VECTOR_H
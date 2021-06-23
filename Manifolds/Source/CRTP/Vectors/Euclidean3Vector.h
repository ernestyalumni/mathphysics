#ifndef CRTP_VECTORS_EUCLIDEAN_3_VECTOR_H
#define CRTP_VECTORS_EUCLIDEAN_3_VECTOR_H

#include "EuclideanVector.h"

#include <algorithm>
#include <array>
#include <initializer_list>

namespace CRTP
{

namespace Vectors
{

template <typename Field>
class Euclidean3Vector :
  public EuclideanVector<Euclidean3Vector<Field>, 3, Field>
{
  public:

    Euclidean3Vector() = default;

    Euclidean3Vector(const std::initializer_list<Field>& l)
    {
      std::copy(l.begin(), l.end(), data_.begin());
    }

    const Field x() const
    {
      return data_[0];
    }

    const Field y() const
    {
      return data_[1];      
    }

    const Field z() const
    {
      return data_[2];      
    }

  private:

    std::array<Field, 3> data_;
};

} // namespace Vectors

} // namespace CRTP

#endif // CRTP_VECTORS_EUCLIDEAN_3_VECTOR_H
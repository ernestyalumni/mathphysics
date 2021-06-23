#ifndef CRTP_VECTORS_VECTOR_H
#define CRTP_VECTORS_VECTOR_H

namespace CRTP
{

namespace Vectors
{

template <typename Implementation, typename Field>
class Vector
{
  public:

    //--------------------------------------------------------------------------
    /// \brief Vector addition
    //--------------------------------------------------------------------------
    Implementation& operator+(const Implementation& y)
    {
      return object()->operator+(y);
    }

    //--------------------------------------------------------------------------
    /// \brief Scalar Multiplication
    //--------------------------------------------------------------------------
    Implementation& operator*(const Field a)
    {
      return object()->operator*(a);
    }

  private:

    Implementation& object()
    {
      return static_cast<Implementation&>(*this);
    }
};

} // namespace Vectors

} // namespace CRTP

#endif // CRTP_VECTORS_VECTOR_H
#include "ProtoFiles/PhysicalConstants/Fundamental.pb.h"
#include "gtest/gtest.h"

namespace GoogleUnitTests
{
namespace ProtoFiles
{
namespace PhysicalConstants
{

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
TEST(FundamentalTest, DefaultConstructs)
{
  PhysicalConstantsFundamental::FundamentalConstants m {};
  SUCCEED();
}

} // namespace PhysicalConstants
} // namespace Protofiles
} // namespace GoogleUnitTests
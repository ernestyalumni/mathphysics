#include <msgpack.hpp>

#include "gtest/gtest.h"

namespace GoogleUnitTests
{
namespace MessagePack
{
namespace Tutorial
{

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
TEST(TutorialTests, MessagePackVersion)
{
  EXPECT_EQ(MSGPACK_VERSION, "6.0.0");
}

} // namespace Tutorial
} // namespace MessagePack
} // namespace GoogleUnitTests
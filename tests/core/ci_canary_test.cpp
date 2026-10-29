// DO NOT MERGE. Canary proving CI's ASan jobs catch a bug the local dev hook cannot:
// a one-element heap over-read that stays inside malloc padding in a plain build.
#include <gtest/gtest.h>

#include <vector>

namespace {

TEST(CiCanary, HeapOverReadOnlyAsanDetects) {
    const std::vector<int> values(4, 7);
    const volatile int* data = values.data();
    const volatile int past_end = data[values.size()];
    static_cast<void>(past_end);
    SUCCEED();
}

}  // namespace

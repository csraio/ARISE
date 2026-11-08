#include "arise/filter/exchangeable.hpp"

#include "arise/chem/perception.hpp"
#include "arise/core/error.hpp"

#include <string>

namespace arise::filter {

namespace {

// "2", "1 to 3" or "at least 1".
std::string describe(const ConfirmedExchangeable& confirmed) {
    if (!confirmed.max) {
        return "at least " + std::to_string(confirmed.min);
    }
    if (*confirmed.max == confirmed.min) {
        return std::to_string(confirmed.min);
    }
    return std::to_string(confirmed.min) + " to " + std::to_string(*confirmed.max);
}

}  // namespace

core::Result<Verdict> check_confirmed_exchangeable(const chem::Molecule& molecule,
                                                   const ConfirmedExchangeable& confirmed) {
    if (confirmed.min < 0 || (confirmed.max && *confirmed.max < confirmed.min)) {
        return core::Error{
            .code = core::ErrorCode::InvalidArgument,
            .message = "confirmed exchangeable H: the range " + describe(confirmed) +
                       " needs a minimum of 0 or more and a maximum not below it",
        };
    }
    const int count = chem::exchangeable_hydrogens(molecule);
    if (count >= confirmed.min && (!confirmed.max || count <= *confirmed.max)) {
        return Verdict{};
    }
    return Verdict{
        .tier = Tier::Eliminated,
        .reasons =
            {
                "confirmed exchangeable H (on O, N, S, Se): " + describe(confirmed) +
                    "; the candidate has " + std::to_string(count),
            },
    };
}

}  // namespace arise::filter

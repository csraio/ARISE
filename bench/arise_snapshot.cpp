// Development-only: writes tests/golden/equivalence_snapshot.tsv (plan 3d), the results the
// speed work must leave unchanged. Usage: arise_snapshot > tests/golden/equivalence_snapshot.tsv
#include "reference_fixtures.hpp"
#include "snapshot.hpp"

#include <exception>
#include <iostream>

int main() {
    try {
        const auto all = arise::test::join_chunks(arise::test::kPubchemSearchesChunks);
        const auto spectrum = arise::test::join_chunks(arise::test::kBenchmarkSpectrumChunks);
        const auto benchmark_case = arise::bench::parse_benchmark_case(spectrum);
        if (!benchmark_case) {
            std::cerr << benchmark_case.error().message << "\n";
            return 1;
        }
        const auto results =
            arise::bench::results(arise::bench::split_searches(all), benchmark_case.value());
        if (!results) {
            std::cerr << results.error().message << "\n";
            return 1;
        }
        std::cout << arise::bench::snapshot_header() << arise::bench::snapshot(results.value());
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "arise_snapshot: " << error.what() << "\n";
    }
    return 1;
}

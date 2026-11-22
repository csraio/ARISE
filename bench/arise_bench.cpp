// Development-only: times ARISE's cleaning and pruning of the benchmark case (plan 3d;
// ADR 0002 target: 10,000 PubChem rows in under 100 ms, native and wasm). The data are
// compiled in, so the same program runs natively and under Node.
//   usage: arise_bench [--runs N] [--budget-ms X]
// With --budget-ms, exits 1 when the slowest scenario's median total exceeds X; 2 on errors.
#include "prune_pipeline.hpp"
#include "reference_fixtures.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using arise::bench::Scenario;

constexpr std::array kScenarios{Scenario::Worst, Scenario::Typical, Scenario::NoDept};
constexpr int kErrorExit = 2;

double since(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// Milliseconds per run: clean, parse, then prune and rank per scenario.
struct Samples {
    std::vector<double> clean;
    std::vector<double> parse;
    std::array<std::vector<double>, kScenarios.size()> prune;
    std::array<std::vector<double>, kScenarios.size()> rank;
    std::array<std::vector<double>, kScenarios.size()> total;
    std::size_t rows = 0;
    std::size_t candidates = 0;
};

using AllEvidence = std::array<arise::bench::Evidence, kScenarios.size()>;

arise::core::Result<AllEvidence> all_evidence(const arise::bench::BenchmarkCase& c) {
    AllEvidence evidence;
    for (std::size_t s = 0; s < kScenarios.size(); ++s) {
        auto e = arise::bench::evidence_for(kScenarios[s], c);
        if (!e) {
            return e.error();
        }
        evidence[s] = std::move(e).value();
    }
    return evidence;
}

// One run of every stage; samples are kept unless `warm_up`.
std::optional<std::string> run_once(std::string_view csv, const AllEvidence& evidence, bool warm_up,
                                    Samples& samples) {
    auto start = Clock::now();
    const auto imported = arise::bench::clean(csv);
    const double clean_ms = since(start);
    if (!imported) {
        return imported.error().message;
    }
    start = Clock::now();
    const auto parsed = arise::bench::parse(imported.value());
    const double parse_ms = since(start);
    samples.rows = imported.value().rows;
    samples.candidates = parsed.molecules.size();
    if (!warm_up) {
        samples.clean.push_back(clean_ms);
        samples.parse.push_back(parse_ms);
    }
    for (std::size_t s = 0; s < kScenarios.size(); ++s) {
        start = Clock::now();
        const auto fits = arise::bench::prune(parsed.molecules, evidence[s]);
        const double prune_ms = since(start);
        if (!fits) {
            return fits.error().message;
        }
        start = Clock::now();
        const auto order = arise::bench::rank(fits.value());
        const double rank_ms = since(start);
        if (order.size() != samples.candidates) {
            return "the ranking lost candidates";
        }
        if (!warm_up) {
            samples.prune[s].push_back(prune_ms);
            samples.rank[s].push_back(rank_ms);
            samples.total[s].push_back(clean_ms + parse_ms + prune_ms + rank_ms);
        }
    }
    return std::nullopt;
}

void print_row(std::string_view label, const std::vector<double>& samples) {
    const auto summary = arise::bench::summarize(samples);
    if (!summary) {
        std::cerr << label << ": " << summary.error().message << "\n";
        return;
    }
    std::cout << "  " << std::left << std::setw(22) << label << std::right << std::fixed
              << std::setprecision(1) << std::setw(9) << summary.value().min << std::setw(9)
              << summary.value().median << std::setw(9) << summary.value().max << "\n";
}

void print_report(const arise::bench::BenchmarkCase& c, const Samples& samples, std::size_t runs) {
    std::cout << "ARISE benchmark: " << c.formula << ", " << samples.rows << " PubChem rows, "
              << samples.candidates << " candidates; nmrshiftdb2 record " << c.record << "; "
              << runs << " runs after one warm-up\n"
              << "  " << std::left << std::setw(22) << "stage (ms)" << std::right << std::setw(9)
              << "min" << std::setw(9) << "median" << std::setw(9) << "max" << "\n";
    print_row("clean", samples.clean);
    print_row("parse", samples.parse);
    for (std::size_t s = 0; s < kScenarios.size(); ++s) {
        const std::string name{arise::bench::to_string(kScenarios[s])};
        print_row("prune " + name, samples.prune[s]);
        print_row("rank " + name, samples.rank[s]);
        print_row("TOTAL " + name, samples.total[s]);
    }
}

int run(std::span<char*> raw) {
    std::vector<std::string_view> args;
    for (std::size_t i = 1; i < raw.size(); ++i) {
        args.emplace_back(raw[i]);
    }
    const auto options = arise::bench::parse_options(args);
    if (!options) {
        std::cerr << options.error().message << "\nusage: arise_bench [--runs N] [--budget-ms X]\n";
        return kErrorExit;
    }
    const auto spectrum = arise::test::join_chunks(arise::test::kBenchmarkSpectrumChunks);
    const auto csv = arise::test::join_chunks(arise::test::kBenchmarkPubchemChunks);
    const auto benchmark_case = arise::bench::parse_benchmark_case(spectrum);
    const auto evidence = benchmark_case ? all_evidence(benchmark_case.value())
                                         : arise::core::Result<AllEvidence>{benchmark_case.error()};
    if (!evidence) {
        std::cerr << evidence.error().message << "\n";
        return kErrorExit;
    }
    Samples samples;
    for (std::size_t r = 0; r <= options.value().runs; ++r) {  // run 0 warms up
        if (const auto error = run_once(csv, evidence.value(), r == 0, samples)) {
            std::cerr << *error << "\n";
            return kErrorExit;
        }
    }
    print_report(benchmark_case.value(), samples, options.value().runs);
    double slowest = 0;
    for (const auto& totals : samples.total) {
        const auto summary = arise::bench::summarize(totals);
        if (!summary) {
            std::cerr << summary.error().message << "\n";
            return kErrorExit;
        }
        slowest = std::max(slowest, summary.value().median);
    }
    const auto code = arise::bench::exit_code(slowest, options.value().budget_ms);
    if (options.value().budget_ms) {
        std::cout << "budget " << *options.value().budget_ms << " ms for the slowest median total ("
                  << slowest << " ms): " << (code == 0 ? "within" : "OVER") << "\n";
    }
    return code;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(std::span{argv, static_cast<std::size_t>(argc)});
    } catch (const std::exception& error) {
        std::cerr << "arise_bench: " << error.what() << "\n";
    } catch (...) {
        std::cerr << "arise_bench: unknown error\n";
    }
    return kErrorExit;
}

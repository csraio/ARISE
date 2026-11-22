// Development-only: writes data/hybridization_shifts.tsv (shipped with the engine) from the
// nmrshiftdb2 extract. Usage: arise_measure_hybridization_shifts <nmrshiftdb2_13c_all.tsv> >
// hybridization_shifts.tsv
#include "hybridization_shifts.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <span>
#include <sstream>

int main(int argc, char** argv) {
    const std::span args{argv, static_cast<std::size_t>(argc)};
    if (args.size() != 2) {
        std::cerr << "usage: arise_measure_hybridization_shifts <nmrshiftdb2_13c_all.tsv>\n";
        return 2;
    }
    std::ifstream in(args[1]);
    if (!in) {
        std::cerr << "cannot read " << args[1] << "\n";
        return 1;
    }
    std::ostringstream text;
    text << in.rdbuf();
    const auto records = arise::tools::parse_shift_extract(text.str());
    std::cout << arise::tools::hybridization_shifts_header()
              << arise::tools::to_tsv(arise::tools::measure_hybridization_shifts(records));
    return 0;
}

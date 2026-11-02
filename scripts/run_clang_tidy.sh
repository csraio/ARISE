#!/bin/sh
# Runs clang-tidy (config: .clang-tidy) on every project source file in the
# compile database. Headers are checked through HeaderFilterRegex.
#   usage: scripts/run_clang_tidy.sh [build-dir]   (default: build/dev; configure it first)
set -eu

repo_root=$(cd "$(dirname "$0")/.." && pwd)
build_dir=${1:-$repo_root/build/dev}
database="$build_dir/compile_commands.json"

if [ ! -f "$database" ]; then
    echo "run_clang_tidy: $database not found; run 'cmake --preset dev' first." >&2
    exit 1
fi

# The standalone clang-tidy does not know where Apple's SDK lives, and CMake
# does not record a sysroot in the compile commands.
extra_args=""
if [ "$(uname -s)" = "Darwin" ]; then
    extra_args="--extra-arg=-isysroot$(xcrun --show-sdk-path)"
fi

cd "$repo_root"
# Tracked plus new (untracked, not ignored) files, so new code is checked before `git add`.
files=$(git ls-files --cached --others --exclude-standard 'src/*.cpp' 'tests/*.cpp' 'apps/*.cpp' 'tools/*.cpp')
if [ -z "$files" ]; then
    echo "run_clang_tidy: no source files."
    exit 0
fi

tidy_version=$(clang-tidy --version)  # captured so a broken binary fails here
echo "run_clang_tidy: $(echo "$tidy_version" | grep -i 'version' | sed -n 1p)"
failures=0
for file in $files; do
    # shellcheck disable=SC2086
    if ! clang-tidy --quiet -p "$build_dir" $extra_args "$file"; then
        failures=$((failures + 1))
    fi
done

if [ "$failures" -ne 0 ]; then
    echo "run_clang_tidy: $failures file(s) with findings."
    exit 1
fi
echo "run_clang_tidy: clean."

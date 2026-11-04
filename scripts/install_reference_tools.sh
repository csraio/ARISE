#!/bin/sh
# Installs the pinned reference-data tools (scripts/reference-requirements.txt: RDKit,
# numpy, Pillow, pdfminer.six and its requirements) into a private Python 3.13 venv.
# Dev-time only: ARISE never needs them.
#   ARISE_REFERENCE_PYTHON  a Python 3.13 interpreter
#                           (default: the one install_emsdk.sh bootstrapped on macOS,
#                            else python3.13 on PATH)
#   ARISE_REFERENCE_DIR     venv location (default: ~/.local/share/arise-reference)
set -eu

repo_root=$(cd "$(dirname "$0")/.." && pwd)
requirements="$repo_root/scripts/reference-requirements.txt"
venv=${ARISE_REFERENCE_DIR:-$HOME/.local/share/arise-reference}

python=${ARISE_REFERENCE_PYTHON:-}
if [ -z "$python" ]; then
    bootstrap="$HOME/.local/share/arise-tools/python-3.13.3/bin/python3"
    if [ -x "$bootstrap" ]; then
        python=$bootstrap
    elif command -v python3.13 >/dev/null; then
        python=$(command -v python3.13)
    else
        echo "install_reference_tools: need Python 3.13 (set ARISE_REFERENCE_PYTHON)." >&2
        exit 1
    fi
fi
version=$("$python" -c 'import sys; print("%d.%d" % sys.version_info[:2])')
if [ "$version" != "3.13" ]; then
    echo "install_reference_tools: $python is Python $version; the pinned wheels need 3.13." >&2
    exit 1
fi

echo "install_reference_tools: venv at $venv (Python $version)"
"$python" -m venv "$venv"
"$venv/bin/python" -m pip install --quiet --upgrade pip
"$venv/bin/python" -m pip install --quiet --require-hashes --only-binary :all: -r "$requirements"
# Capture first: sh has no pipefail, so a crash must not hide behind a pipe.
rdkit_version=$("$venv/bin/python" -c 'import rdkit; print(rdkit.__version__)')
echo "install_reference_tools: RDKit $rdkit_version"
echo "install_reference_tools: done. Build the reference data with:"
echo "    $venv/bin/python scripts/build_reference_sets.py"
echo "    $venv/bin/python scripts/build_solvent_peaks.py"

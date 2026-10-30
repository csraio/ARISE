#!/bin/sh
# Installs the pinned Emscripten SDK outside the repo, then activates it.
# emsdk needs Python >= 3.10. If python3 is older (macOS 12 ships 3.9), this first
# bootstraps the relocatable Python that emsdk itself bundles, checked by SHA-256.
#   ARISE_EMSDK_DIR  emsdk location           (default: ~/.local/share/emsdk)
#   ARISE_TOOLS_DIR  bootstrap Python parent  (default: ~/.local/share/arise-tools)
# Afterwards, in each shell that builds wasm:  . ~/.local/share/emsdk/emsdk_env.sh
set -eu

emsdk_version=6.0.13
# The SDK bundles Node 24, which needs macOS >= 13.5; Node 22 LTS runs on macOS 11+.
# Pinned on every OS so local and CI test runs use the same JavaScript engine.
node_tool=node-22.16.0-64bit
python_archive=python-3.13.3-0-macos-x86_64.tar.gz
python_sha256=d3f30e1b78354d6f74fae951175dc2bd22b386aafeab13edd6ac7c9613b0d33c
python_url=https://storage.googleapis.com/webassembly/emscripten-releases-builds/deps/$python_archive

emsdk_dir=${ARISE_EMSDK_DIR:-$HOME/.local/share/emsdk}
tools_dir=${ARISE_TOOLS_DIR:-$HOME/.local/share/arise-tools}

tmp_dir=""
cleanup() {
    if [ -n "$tmp_dir" ]; then
        rm -rf "$tmp_dir"
    fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM

python_is_new_enough() {
    "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'
}

# 1. A Python that can run emsdk.
if command -v python3 >/dev/null && python_is_new_enough python3; then
    EMSDK_PYTHON=$(command -v python3)
elif [ "$(uname -s)-$(uname -m)" = "Darwin-x86_64" ]; then
    bootstrap="$tools_dir/python-3.13.3"
    if [ ! -x "$bootstrap/bin/python3" ]; then
        echo "install_emsdk: bootstrapping Python 3.13.3 into $bootstrap"
        tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/arise-emsdk.XXXXXX")
        curl -fsSL -o "$tmp_dir/$python_archive" "$python_url"
        echo "$python_sha256  $tmp_dir/$python_archive" | shasum -a 256 -c -
        mkdir -p "$bootstrap"
        tar -xzf "$tmp_dir/$python_archive" -C "$bootstrap" --strip-components=1
    fi
    EMSDK_PYTHON="$bootstrap/bin/python3"
    python_is_new_enough "$EMSDK_PYTHON"
else
    echo "install_emsdk: need python3 >= 3.10 on PATH." >&2
    exit 1
fi
export EMSDK_PYTHON
echo "install_emsdk: using $("$EMSDK_PYTHON" --version) at $EMSDK_PYTHON"

# 2. emsdk itself, at the pinned tag.
if [ -d "$emsdk_dir/.git" ]; then
    git -C "$emsdk_dir" fetch --depth 1 origin "refs/tags/$emsdk_version:refs/tags/$emsdk_version"
    git -C "$emsdk_dir" checkout --quiet "$emsdk_version"
else
    git clone --quiet --depth 1 --branch "$emsdk_version" \
        https://github.com/emscripten-core/emsdk.git "$emsdk_dir"
fi

# 3. The SDK (clang, binaryen, emscripten) for that version, plus the pinned Node.
"$emsdk_dir/emsdk" install "$emsdk_version" "$node_tool"
"$emsdk_dir/emsdk" activate "$emsdk_version" "$node_tool"

# 4. Prove the activated toolchain runs on this machine.
# shellcheck disable=SC1091
. "$emsdk_dir/emsdk_env.sh" >/dev/null
# Capture before printing: `sh` has no pipefail, so `emcc ... | head` would hide a crash.
emcc_version=$(emcc --version)
# emsdk does not put node on PATH; Emscripten and CMake use $EMSDK_NODE.
node_version=$("$EMSDK_NODE" --version)
echo "install_emsdk: $(echo "$emcc_version" | sed -n 1p)"
echo "install_emsdk: node $node_version at $EMSDK_NODE"

echo "install_emsdk: done. In each shell that builds wasm, run:"
echo "    . $emsdk_dir/emsdk_env.sh"

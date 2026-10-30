#!/bin/sh
# Installs the pinned clang-format and clang-tidy (scripts/dev-tools-requirements.txt)
# into a private venv and links them into a bin directory on PATH.
# Needs only python3 with venv; works on macOS 10.9+ (x86_64/arm64) and Linux.
#   ARISE_TOOLS_DIR  venv location   (default: ~/.local/share/arise-tools)
#   ARISE_BIN_DIR    link directory  (default: ~/.local/bin)
set -eu

repo_root=$(cd "$(dirname "$0")/.." && pwd)
requirements="$repo_root/scripts/dev-tools-requirements.txt"
tools_dir=${ARISE_TOOLS_DIR:-$HOME/.local/share/arise-tools}
bin_dir=${ARISE_BIN_DIR:-$HOME/.local/bin}

echo "install_dev_tools: venv at $tools_dir"
python3 -m venv "$tools_dir"
"$tools_dir/bin/python" -m pip install --quiet --upgrade pip
"$tools_dir/bin/python" -m pip install --quiet \
    --require-hashes --only-binary :all: -r "$requirements"

mkdir -p "$bin_dir"
for tool in clang-format clang-tidy; do
    target="$tools_dir/bin/$tool"
    link="$bin_dir/$tool"
    if [ -e "$link" ] && [ ! -L "$link" ]; then
        echo "install_dev_tools: $link exists and is not a symlink; refusing to replace it." >&2
        exit 1
    fi
    ln -sf "$target" "$link"
    echo "install_dev_tools: $link -> $target"
    "$link" --version | grep -i version
done

case ":$PATH:" in
    *":$bin_dir:"*) ;;
    *) echo "install_dev_tools: NOTE $bin_dir is not on PATH; add it to use the tools." >&2 ;;
esac
echo "install_dev_tools: done."

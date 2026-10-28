#!/bin/sh
# Self-test of the pre-commit gate. Runs entirely inside a throwaway git worktree of
# HEAD, so the caller's working tree is never modified.
#   Check A: the hook accepts clean HEAD.
#   Check B: the hook rejects HEAD plus a deliberately failing canary test.
set -eu

repo_root=$(git rev-parse --show-toplevel)
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/arise-gate.XXXXXX")
worktree="$tmp_dir/worktree"
log="$tmp_dir/hook.log"

cleanup() {
    if [ -d "$worktree" ]; then
        git -C "$repo_root" worktree remove --force "$worktree"
    fi
    rm -rf "$tmp_dir"
    git -C "$repo_root" worktree prune
}
trap cleanup EXIT
trap 'exit 130' INT TERM

git -C "$repo_root" worktree add --detach --quiet "$worktree" HEAD
cd "$worktree"

failures=0

echo "Check A: gate accepts clean HEAD ..."
if .githooks/pre-commit >"$log" 2>&1; then
    echo "  PASS"
else
    echo "  FAIL: hook rejected clean HEAD. Hook output:"
    sed 's/^/    /' "$log"
    failures=$((failures + 1))
fi

echo "Check B: gate rejects a failing test ..."
cat >>tests/core/result_test.cpp <<'EOF'

TEST(GateCanary, DeliberateFailure) { FAIL() << "canary: the commit gate must reject this"; }
EOF
if .githooks/pre-commit >"$log" 2>&1; then
    echo "  FAIL: hook accepted a tree with a failing test."
    failures=$((failures + 1))
elif grep -q "GateCanary.DeliberateFailure" "$log"; then
    echo "  PASS"
else
    echo "  FAIL: hook failed, but not because of the canary test. Hook output:"
    sed 's/^/    /' "$log"
    failures=$((failures + 1))
fi

if [ "$failures" -ne 0 ]; then
    echo "verify_commit_gate: $failures check(s) failed."
    exit 1
fi
echo "verify_commit_gate: all checks passed."

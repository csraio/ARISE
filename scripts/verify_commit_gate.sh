#!/bin/sh
# Self-test of the pre-commit gate. Runs entirely inside a throwaway git worktree of
# HEAD, so the caller's working tree is never modified.
#   Check A: the hook accepts clean HEAD.
#   Check B: the hook rejects a staged, deliberately failing canary test.
#   Check C: an untracked file cannot rescue a commit (the gate tests the index).
#   Check D: an unstaged breaking edit cannot fail a clean commit (same reason).
# ARISE_GATE_HOOK overrides the hook under test (default: the worktree's own
# .githooks/pre-commit), so a new hook can be tested before it is committed.
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
hook=${ARISE_GATE_HOOK:-$worktree/.githooks/pre-commit}
echo "verify_commit_gate: testing $hook"

failures=0
fail() {
    echo "  FAIL: $1"
    if [ "${2:-}" = "show-log" ]; then
        sed 's/^/    /' "$log"
    fi
    failures=$((failures + 1))
}

# Restores the worktree and its index to HEAD between checks.
reset_worktree() {
    git reset --quiet --hard HEAD
    git clean --quiet -d --force
}

# Formatted the way clang-format wants it, so Check B fails on the test, not on style.
append_canary_test() {
    printf '\n%s\n%s\n%s\n' 'TEST(GateCanary, DeliberateFailure) {' \
        '    FAIL() << "canary: the commit gate must reject this";' '}' >>tests/core/result_test.cpp
}

echo "Check A: gate accepts clean HEAD ..."
if "$hook" >"$log" 2>&1; then
    echo "  PASS"
else
    fail "hook rejected clean HEAD. Hook output:" show-log
fi

echo "Check B: gate rejects a staged failing test ..."
append_canary_test
git add tests/core/result_test.cpp
if "$hook" >"$log" 2>&1; then
    fail "hook accepted a staged failing test."
elif grep -q "GateCanary.DeliberateFailure" "$log"; then
    echo "  PASS"
else
    fail "hook failed, but not because of the canary test. Hook output:" show-log
fi
reset_worktree

echo "Check C: an untracked file cannot rescue a commit ..."
# Stage a build reference to a source file that exists only in the working tree.
printf '%s\n' '#include <gtest/gtest.h>' 'TEST(GateCanary, Untracked) { SUCCEED(); }' \
    >tests/core/gate_untracked_test.cpp
sed 's/    result_test.cpp/    result_test.cpp\
    gate_untracked_test.cpp/' tests/core/CMakeLists.txt >"$tmp_dir/CMakeLists.txt"
cp "$tmp_dir/CMakeLists.txt" tests/core/CMakeLists.txt
git add tests/core/CMakeLists.txt
if "$hook" >"$log" 2>&1; then
    fail "hook accepted a commit that only builds thanks to an untracked file."
elif grep -q "gate_untracked_test.cpp" "$log"; then
    echo "  PASS"
else
    fail "hook failed, but not because of the missing file. Hook output:" show-log
fi
reset_worktree

echo "Check D: an unstaged breaking edit cannot fail a clean commit ..."
append_canary_test
if "$hook" >"$log" 2>&1; then
    echo "  PASS"
else
    fail "hook rejected a clean index because of an unstaged edit. Hook output:" show-log
fi
reset_worktree

if [ "$failures" -ne 0 ]; then
    echo "verify_commit_gate: $failures check(s) failed."
    exit 1
fi
echo "verify_commit_gate: all checks passed."

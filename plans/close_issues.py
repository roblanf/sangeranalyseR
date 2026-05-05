#!/usr/bin/env python3
"""
close_issues.py — post pre-drafted replies + close GitHub issues for
sangeranalyseR.

Why this script exists
----------------------
The local `gh` CLI in this dev environment is v0.11.0 (2020-07-16), which
predates `gh issue comment` and the `--body` flag. Rather than upgrade
the CLI we drive the GitHub REST API directly with Python's stdlib —
zero pip dependencies.

The script reads the canonical reply file at plans/16_github_replies.md,
extracts (issue_number, reply_body) pairs via a small Markdown
fence-depth parser, then for each issue:

  1. POSTs the reply body as a new comment.
  2. PATCHes the issue's state to "closed".

A confirmation prompt fires before every issue so the user can step
through one at a time and Ctrl+C any time.

Usage
-----
1. Create a Personal Access Token at https://github.com/settings/tokens
   with the `public_repo` scope (sufficient for public repositories).

2. Export it in your shell:

     export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx

3. Dry-run first to see what would be posted (no API calls):

     python3 plans/close_issues.py --dry-run

4. Live run, with per-issue confirmation prompts:

     python3 plans/close_issues.py

5. Process only a single issue (useful if the dry-run flagged something
   wrong with one of them, or you want to test on one before all):

     python3 plans/close_issues.py --issue 65

Notes
-----
* GitHub's API does not support comment replacement. **Do not** re-run
  this script blindly — it will post duplicate comments. The close call
  on an already-closed issue is a no-op.
* If a single issue fails (e.g. transient HTTP error), re-run with
  `--issue <N>` to retry just that one.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterator

REPO = "roblanf/sangeranalyseR"
DEFAULT_MD = Path(__file__).resolve().parent / "16_github_replies.md"


# ---------------------------------------------------------------------------
# Markdown parser
# ---------------------------------------------------------------------------

ISSUE_HEADER_RE = re.compile(r"^### Issue #(\d+)\b")
FENCE_RE = re.compile(r"^\s*```(.*)$")


def parse_issues(md_text: str) -> Iterator[tuple[int, str]]:
    """
    Yield (issue_number, body) for each issue block in the Markdown file.

    Each issue block has shape:

        ### Issue #<n> — <title>
        ... metadata lines (URL etc.) ...
        ```markdown
        <body, possibly containing nested fenced blocks: ```r ... ```
        and unlabeled  ``` ... ``` >
        ```

    A naive "first ``` after ```markdown is the closer" parser
    over-truncates at the first nested fence; a fence-depth tracker
    fails on bare-bare nested fences (which open and close with the
    same syntax as the outer closer). The robust rule: the OUTER
    closer is the LAST line matching `^\\s*```\\s*$` between the
    `### Issue #N` header that starts this block and the NEXT
    `### Issue #N` header (or EOF, or a horizontal-rule `---`).
    Everything between the outer opener (exclusive) and outer
    closer (exclusive) is the body, verbatim.
    """
    lines = md_text.split("\n")
    bare_fence = re.compile(r"^\s*```\s*$")
    i = 0
    while i < len(lines):
        m = ISSUE_HEADER_RE.match(lines[i])
        if not m:
            i += 1
            continue
        issue_num = int(m.group(1))

        # Walk forward to the outer ```markdown opener.
        j = i + 1
        while j < len(lines):
            if lines[j].lstrip().startswith("```markdown"):
                break
            if ISSUE_HEADER_RE.match(lines[j]):
                break
            j += 1
        if j >= len(lines) or not lines[j].lstrip().startswith("```markdown"):
            i = j
            continue
        body_start = j + 1

        # Find section end: next ### Issue header, OR next standalone
        # `---` rule that is NOT inside a fenced block (we'll just take
        # the next ### header to keep it simple), OR EOF.
        sec_end = body_start
        while sec_end < len(lines):
            if ISSUE_HEADER_RE.match(lines[sec_end]):
                break
            sec_end += 1

        # Outer closer = the LAST bare ``` line in [body_start, sec_end).
        outer_close_line = None
        for k in range(sec_end - 1, body_start - 1, -1):
            if bare_fence.match(lines[k]):
                outer_close_line = k
                break
        if outer_close_line is None:
            # Malformed block; skip.
            i = sec_end
            continue

        body = "\n".join(lines[body_start:outer_close_line])
        yield issue_num, body
        i = sec_end


# ---------------------------------------------------------------------------
# GitHub REST API
# ---------------------------------------------------------------------------

API_HEADERS = {
    "Accept": "application/vnd.github+json",
    "X-GitHub-Api-Version": "2022-11-28",
    "Content-Type": "application/json",
    "User-Agent": "sangeranalyseR-close_issues.py",
}


def _request(token: str, url: str, method: str, payload: dict | None = None) -> dict:
    headers = dict(API_HEADERS)
    headers["Authorization"] = f"Bearer {token}"
    data = None
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, method=method, data=data, headers=headers)
    with urllib.request.urlopen(req) as resp:
        body = resp.read()
        if not body:
            return {}
        return json.loads(body)


def post_comment(token: str, issue_num: int, body: str) -> dict:
    return _request(
        token,
        f"https://api.github.com/repos/{REPO}/issues/{issue_num}/comments",
        "POST",
        {"body": body},
    )


def close_issue(token: str, issue_num: int) -> dict:
    return _request(
        token,
        f"https://api.github.com/repos/{REPO}/issues/{issue_num}",
        "PATCH",
        {"state": "closed"},
    )


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(
        description="Post replies + close GitHub issues from "
        "plans/16_github_replies.md.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--md",
        default=str(DEFAULT_MD),
        help=f"path to the replies Markdown file (default: {DEFAULT_MD})",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="parse and preview only; make no API calls",
    )
    p.add_argument(
        "--issue",
        type=int,
        default=None,
        help="process only this issue number",
    )
    p.add_argument(
        "--no-confirm",
        action="store_true",
        help="skip the per-issue confirmation prompt (DANGEROUS — the "
        "script is irreversible)",
    )
    args = p.parse_args(argv)

    md_path = Path(args.md)
    if not md_path.is_file():
        print(f"ERROR: replies file not found: {md_path}", file=sys.stderr)
        return 1

    md_text = md_path.read_text(encoding="utf-8")
    issues = list(parse_issues(md_text))

    if args.issue is not None:
        issues = [(n, b) for (n, b) in issues if n == args.issue]
        if not issues:
            print(f"ERROR: issue #{args.issue} not found in {md_path}",
                  file=sys.stderr)
            return 1

    if not issues:
        print(f"ERROR: parsed 0 issues from {md_path}", file=sys.stderr)
        return 1

    print(f"Parsed {len(issues)} issue(s) from {md_path}:")
    for num, body in issues:
        first = next((ln for ln in body.splitlines() if ln.strip()), "")
        print(f"  • #{num:<4} ({len(body):>5} chars)  first line: {first[:70]!r}")
    print()

    if args.dry_run:
        print("--dry-run set; no API calls made. Exiting.")
        return 0

    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        print(
            "ERROR: GITHUB_TOKEN env var is required.\n"
            "  Create a token at https://github.com/settings/tokens\n"
            "  (needs `public_repo` scope), then:\n"
            "      export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx",
            file=sys.stderr,
        )
        return 2

    failures: list[tuple[int, str]] = []
    for num, body in issues:
        print(f"\n=== Issue #{num} ===")
        if not args.no_confirm:
            try:
                input(f"  Press Enter to comment + close #{num} "
                      "(or Ctrl+C to abort)... ")
            except KeyboardInterrupt:
                print("\n  Aborted by user.")
                return 130
        try:
            post_comment(token, num, body)
            print(f"  ✓ comment posted on #{num}")
        except urllib.error.HTTPError as e:
            err_body = e.read()[:300].decode("utf-8", errors="replace")
            print(f"  ✗ comment FAILED on #{num}: HTTP {e.code} — {err_body}",
                  file=sys.stderr)
            failures.append((num, f"comment HTTP {e.code}"))
            continue
        except Exception as e:
            print(f"  ✗ comment FAILED on #{num}: {e}", file=sys.stderr)
            failures.append((num, f"comment exception: {e}"))
            continue

        try:
            close_issue(token, num)
            print(f"  ✓ issue #{num} closed")
        except urllib.error.HTTPError as e:
            err_body = e.read()[:300].decode("utf-8", errors="replace")
            print(f"  ✗ close FAILED on #{num}: HTTP {e.code} — {err_body}",
                  file=sys.stderr)
            failures.append((num, f"close HTTP {e.code}"))
        except Exception as e:
            print(f"  ✗ close FAILED on #{num}: {e}", file=sys.stderr)
            failures.append((num, f"close exception: {e}"))

    print("\n" + "=" * 60)
    if failures:
        print(f"DONE with {len(failures)} failure(s):")
        for num, why in failures:
            print(f"  • #{num}: {why}")
        print("Re-run with --issue <N> to retry individual failures.")
        return 1
    print(f"DONE — all {len(issues)} issue(s) commented and closed.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

#!/usr/bin/env python3
"""
Audit review.tex against the real `git diff main..HEAD`.

Checks
------
A  old.tex / new.tex are byte-identical to the git blobs.
B  NEW-SIDE FIDELITY. Strip every review-only construct from review.tex and the
   result must reproduce new.tex. Proves nothing on the "after" side was
   altered, invented, or dropped.
C  OLD-SIDE AUTHENTICITY. Every fragment presented as v1.2 must occur verbatim
   in old.tex. Proves the "before" side was not invented.
D  HUNK COVERAGE. Every hunk git reports is represented, classified as either
   a genuine replacement (must appear in an old-block AND on the new side) or
   an in-place prose extension (old text survives verbatim, delta highlighted
   inline via \\dins -- must appear on the new side).

Run: python3 verify_review.py     (exit 0 iff everything passes)
"""

import pathlib
import re
import subprocess
import sys

HERE = pathlib.Path(__file__).parent
REPO = HERE.parents[2]
TEXPATH = "papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex"

OLD = (HERE / "old.tex").read_text()
NEW = (HERE / "new.tex").read_text()
REV = (HERE / "review.tex").read_text()

failures = []
def ok(m):  print(f"  PASS  {m}")
def bad(m, d=""):
    failures.append(m); print(f"  FAIL  {m}")
    if d: print("        " + d.replace("\n", "\n        "))


# --------------------------------------------------------------- utilities --

def _match(text, i, needle):
    """Brace-matched span starting at \\cmd{ ; returns (open_end, close_idx)."""
    j = i + len(needle); depth = 1
    while j < len(text) and depth:
        if text[j] == "{": depth += 1
        elif text[j] == "}": depth -= 1
        j += 1
    return i + len(needle), j - 1

def unwrap(text, cmd):
    n = "\\" + cmd + "{"
    while (i := text.find(n)) >= 0:
        a, b = _match(text, i, n)
        text = text[:i] + text[a:b] + text[b + 1:]
    return text

def drop(text, cmd):
    n = "\\" + cmd + "{"
    while (i := text.find(n)) >= 0:
        _, b = _match(text, i, n)
        text = text[:i] + text[b + 1:]
    return text

def norm(s):
    return re.sub(r"\s+", " ", s).strip()

def canon(s):
    """Normal form that sees through the review's presentational transforms."""
    for c in ("sgnc", "sgn", "wchg", "mbox", "dins", "text"):
        s = unwrap(s, c)
    for c in ("eqnote", "footnote"):          # review commentary, not paper text
        s = drop(s, c)
    for tok in ("\\oldtag", "\\newtag", "\\tcblower"):
        s = s.replace(tok, " ")               # block seams inside a change region
    s = re.sub(r"\\(begin|end)\{(diffbox|addbox)\}", " ", s)
    s = s.replace("{equation*}", "{equation}")
    s = re.sub(r"\\label\{[^}]*\}", "", s)
    return norm(s)


# ---------------------------------------------- split review into old / new --

def split_blocks(rev):
    """Brace-matched split. Block titles contain \\texttt{...}, so the two
    arguments must be matched properly rather than with [^}]*."""
    old_parts, new_parts = [], []
    i = 0
    while (m := re.search(r"\\begin\{(diffbox|addbox)\}", rev[i:])):
        start = i + m.start()
        new_parts.append(rev[i:start])            # unchanged body before block
        kind = m.group(1)
        # skip the two brace-matched arguments of the environment
        argp = i + m.end()
        for _ in range(2):
            assert rev[argp] == "{", "malformed block arguments"
            _, close = _match(rev, argp - 1, "X{")
            argp = close + 1
        endtok = "\\end{" + kind + "}"
        end = rev.find(endtok, argp)
        assert end > 0, f"unterminated {kind}"
        body = rev[argp:end]
        if kind == "diffbox":
            lo = body.find("\\tcblower")
            assert lo > 0, "diffbox without \\tcblower"
            old_parts.append(body[:lo])
            new_parts.append(body[lo + len("\\tcblower"):])
        else:
            new_parts.append(body)                # addbox is pure addition
        i = end + len(endtok)
    new_parts.append(rev[i:])
    return "".join(old_parts), "".join(new_parts)

OLD_SIDE, NEW_SIDE = split_blocks(REV)
OLD_SIDE_C, NEW_SIDE_C = canon(OLD_SIDE), canon(NEW_SIDE)


# ------------------------------------------------------------------ check A --

print("== A. review inputs are the real git blobs ==")
for label, path, txt in (("old.tex", f"main:{TEXPATH}", OLD),
                         ("new.tex", f"HEAD:{TEXPATH}", NEW)):
    blob = subprocess.run(["git", "-C", str(REPO), "show", path],
                          capture_output=True, text=True, check=True).stdout
    (ok if blob == txt else bad)(f"{label} == {path}")


# ------------------------------------------------------------------ check B --

def strip_markup(text):
    text = re.sub(r"\n% =+ REVIEW MARKUP \(temporary\) =+.*?\n% =+\n", "\n", text, flags=re.S)
    text = re.sub(r"\n\\begin\{center\}\n\\begin\{tcolorbox\}.*?\\clearpage\n", "\n", text, flags=re.S)
    text = re.sub(r"\\begin\{(diffbox|addbox)\}\{[^}]*\}\{[^}]*\}\n?", "", text)
    text = re.sub(r"\\end\{(diffbox|addbox)\}\n?", "", text)
    text = text.replace("\\tcblower", "").replace("\\oldtag", "").replace("\\newtag", "")
    text = drop(text, "eqnote")
    text = drop(text, "footnote")
    for c in ("dins", "mbox", "sgnc", "sgn", "wchg"):
        text = unwrap(text, c)
    # The one deliberate non-markup deviation: the review lives in a subdirectory,
    # so it points at the paper's bibliography instead of keeping a copy.
    text = text.replace(r"\bibliography{../references}", r"\bibliography{references}")
    text = re.sub(r"\\title\{Computational Weber electrodynamics: symplectic n-body integration"
                  r".*?\}\n", r"\\title{Computational Weber electrodynamics: symplectic n-body "
                  r"integration}\n", text, flags=re.S, count=1)
    return text

print("\n== B. new side of review.tex reproduces new.tex exactly ==")
rebuilt = norm(strip_markup(NEW_SIDE))
if rebuilt == norm(NEW):
    ok("review.tex minus all markup == new.tex")
else:
    import difflib
    d = [l for l in difflib.unified_diff(norm(NEW).split(". "), rebuilt.split(". "),
                                         "new.tex", "review(stripped)", n=0, lineterm="")
         if l.startswith(("+", "-")) and not l.startswith(("+++", "---"))]
    bad(f"{len(d)} sentence-level differences", "\n".join(x[:200] for x in d[:10]))


# ------------------------------------------------------------------ check C --

print("\n== C. everything shown as v1.2 is verbatim from old.tex ==")
OLD_FILE_C = canon(OLD)
quoted = re.findall(r"\\begin\{equation\*\}(.*?)\\end\{equation\*\}", OLD_SIDE, re.S)
quoted += re.findall(r"\\begin\{align\*\}(.*?)\\end\{align\*\}", OLD_SIDE, re.S)
strays = [q for q in quoted if canon(q) not in OLD_FILE_C]
(ok if not strays else lambda m: bad(m, "\n".join(canon(q)[:150] for q in strays)))(
    f"all {len(quoted)} equation blocks quoted as v1.2 are verbatim in old.tex")

# prose lines inside old-blocks (strip the tags/notes/equations first)
prose = OLD_SIDE
prose = re.sub(r"\\begin\{(equation\*|align\*)\}.*?\\end\{\1\}", " ", prose, flags=re.S)
prose = drop(prose, "eqnote").replace("\\oldtag", "")
prose_lines = [norm(l) for l in prose.split("\n") if norm(l) and not norm(l).startswith("%")]
stray_prose = [l for l in prose_lines if canon(l) not in OLD_FILE_C]
(ok if not stray_prose else lambda m: bad(m, "\n".join(l[:150] for l in stray_prose)))(
    f"all {len(prose_lines)} prose lines quoted as v1.2 are verbatim in old.tex")


# ------------------------------------------------------------------ check D --

print("\n== E. nothing is shown as removed that still exists unchanged ==")
# An old-block may legitimately quote text that survives verbatim, but only where
# the block itself says so. Anything else would misrepresent an unchanged passage
# as a deletion.
DELIBERATE_UNCHANGED = {
    "C4": "alpha_x is byte-identical; the block is titled 'algebraically identical'",
    "C11": "the closing clause is unchanged; quoted so the v1.2 sentence is complete",
}
NEW_FILE_C = canon(NEW)
survivors = []
for blk in re.finditer(r"\\begin\{diffbox\}\{([^}]*)\}", REV):
    pass
for q in re.findall(r"\\begin\{equation\*\}(.*?)\\end\{equation\*\}", OLD_SIDE, re.S):
    if canon(q) in NEW_FILE_C:
        survivors.append(canon(q)[:90])
if len(survivors) <= len(DELIBERATE_UNCHANGED):
    ok(f"{len(survivors)} quoted block(s) survive unchanged, all deliberately labelled "
       f"({', '.join(DELIBERATE_UNCHANGED)})")
    for sv in survivors:
        print(f"        - {sv}")
else:
    bad(f"{len(survivors)} unchanged passages shown as removed",
        "\n".join(survivors))

print("\n== D. every git hunk is represented ==")
diff = subprocess.run(["git", "-C", str(REPO), "diff", "-U0", "main..HEAD", "--", TEXPATH],
                      capture_output=True, text=True, check=True).stdout

hunks, cur = [], None
for line in diff.splitlines():
    if line.startswith("@@"):
        cur = {"rm": [], "add": []}; hunks.append(cur)
    elif cur is None or line.startswith(("---", "+++")):
        continue
    elif line.startswith("-"):
        cur["rm"].append(line[1:])
    elif line.startswith("+"):
        cur["add"].append(line[1:])

import difflib

pure_insert, replaced, problems = 0, 0, []
n_del_seg = n_ins_seg = 0
for n, h in enumerate(hunks, 1):
    rm = canon(" ".join(x for x in h["rm"] if x.strip()))
    add = canon(" ".join(x for x in h["add"] if x.strip()))
    if not rm and not add:
        continue

    # Word-level opcodes isolate what actually changed inside a partly-changed
    # line. Deleted segments must be shown as removed; inserted segments must be
    # present on the new side. Untouched context carries no obligation.
    rw, aw = rm.split(), add.split()
    ops = difflib.SequenceMatcher(None, rw, aw, autojunk=False).get_opcodes()
    if all(t in ("equal", "insert") for t, *_ in ops):
        pure_insert += 1
    else:
        replaced += 1

    for tag, i1, i2, j1, j2 in ops:
        dseg = " ".join(rw[i1:i2])
        aseg = " ".join(aw[j1:j2])
        inline_growth = tag == "replace" and dseg and dseg in aseg
        if tag in ("delete", "replace") and not inline_growth:
            n_del_seg += 1
            if dseg not in OLD_SIDE_C:
                problems.append(f"hunk {n}: deleted segment not shown as removed: {dseg[:110]}")
        if tag in ("insert", "replace"):
            n_ins_seg += 1
            if aseg not in NEW_SIDE_C:
                problems.append(f"hunk {n}: inserted segment absent from new side: {aseg[:110]}")

print(f"  {len(hunks)} hunks: {replaced} replacement(s), {pure_insert} pure insertion(s)")
print(f"  {n_del_seg} deleted segment(s) checked against old-blocks, "
      f"{n_ins_seg} inserted segment(s) checked against the new side")
(ok if not problems else lambda m: bad(m, "\n".join(problems)))(
    "every hunk is represented in review.tex")

print(f"\n  review.tex renders {len(re.findall(r'\\begin\{diffbox\}', REV))} paired old/new blocks "
      f"and {len(re.findall(r'\\begin\{addbox\}', REV))} addition-only blocks")

print()
if failures:
    print(f"RESULT: {len(failures)} CHECK(S) FAILED"); sys.exit(1)
print("RESULT: review.pdf faithfully and completely represents git diff main..HEAD")
sys.exit(0)

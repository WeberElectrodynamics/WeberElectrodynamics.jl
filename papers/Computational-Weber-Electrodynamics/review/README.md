# Review build — TEMPORARY

Renders the paper with every v1.2 → v1.3 change marked up inline, for review of
PR #92 (canonical Weber Hamiltonian correction).

**Delete this whole directory once the changes are accepted.**

## Reading it

Open [`review.pdf`](review.pdf).

- Dark header bar = one change site, labelled `C1`…`C12`. Page 2 indexes them all.
- **Red field, top half** = removed v1.2 content, deliberately unnumbered so it
  cannot be mistaken for a live equation.
- **Green field, bottom half** = added v1.3 content, with real labels and numbering.
- **Yellow token** = a single symbol that flipped sign or changed name. This is
  the one to scan for — all six sign flips are marked this way.
- Green inline highlight = prose inserted in v1.3.
- Grey sans-serif note under a block = why the change was needed.

Everything outside a marked block is byte-identical to the paper, so the
document still reads as the paper itself.

## Files

| File | Role |
|---|---|
| `review.pdf` | the deliverable |
| `review.tex` | generated source of the PDF |
| `make_review.py` | builds `review.tex` from the two git blobs |
| `verify_review.py` | audits `review.tex` against the real `git diff` |
| `old.tex`, `new.tex` | build inputs, extracted from git — **not tracked** |

## Rebuilding

Needs the **full** TeX Live 2025, not `2025basic` — `tcolorbox` and `tikzfill`
are absent from the basic install:

```bash
export PATH="/usr/local/texlive/2025/bin/universal-darwin:$PATH"
python3 make_review.py && latexmk -pdf review.tex
```

`make_review.py` extracts both inputs straight from git (`main` and `HEAD` by
default, overridable via `WEBER_REVIEW_BASE` / `WEBER_REVIEW_HEAD`), then splices
the old blocks into the new paper at 12 anchored sites. It hard-fails on any
anchor it cannot find, so it cannot silently drift out of sync with the paper.

## Verifying

```bash
python3 verify_review.py     # exit 0 iff the PDF matches the real diff
```

Five checks against `git diff main..HEAD`:

- **A** — the inputs are byte-identical to the `main` and `HEAD` git blobs.
- **B** — stripping every review construct from `review.tex` reproduces `new.tex`
  exactly, proving nothing on the "after" side was altered, invented, or dropped.
- **C** — every fragment presented as v1.2 occurs verbatim in the `main` blob.
- **D** — all 23 hunks are represented, checked at word-segment granularity so
  that partly-changed lines are handled correctly.
- **E** — the converse: nothing is shown as removed that actually survives
  unchanged, apart from one deliberately labelled case (`alpha_x` in C4, which is
  byte-identical but re-contextualised).

The only intentional deviation from the paper source is `\bibliography{../references}`,
since the review lives in a subdirectory; `verify_review.py` whitelists it explicitly.

# Wolfram Language sandbox

A Jupyter notebook-based sandbox for Wolfram Language computations alongside
the Julia package. Runs on the free [Wolfram Engine](https://www.wolfram.com/engine/)
via a Jupyter kernel.

## Quickstart

Open [tour.ipynb](tour.ipynb) in either:

- **Jupyter Lab** — `~/.julia/conda/3/aarch64/bin/jupyter lab wolfram/tour.ipynb`
- **VS Code** — with the official [Wolfram Language extension](https://marketplace.visualstudio.com/items?itemName=WolframResearch.wolfram) installed, just open the file and pick the `Wolfram Language 14.3` kernel.

The notebook is ready-to-run: typeset math renders inline, plots are embedded,
and a final "Exporting to the paper" section writes LaTeX snippets and PNGs
to `wolfram/out/` (gitignored).

## First-time setup

```bash
# 1. Activate the free Wolfram Engine license (once, interactive).
wolframscript -activate

# 2. Register the Wolfram kernel with Jupyter.
git clone https://github.com/WolframResearch/WolframLanguageForJupyter.git \
    ~/.local/share/WolframLanguageForJupyter
~/.local/share/WolframLanguageForJupyter/configure-jupyter.wls add
```

Sign in with a free [Wolfram ID](https://www.wolfram.com/id/) during activation.
After setup, `jupyter kernelspec list` should include `wolframlanguage14.3`
alongside the existing `julia-1.12` and `python3` kernels.

## What the notebook covers

1. **Symbolic computation** — integration, closed-form sums, high-precision numerics.
2. **ODEs** — `DSolve` for a simple harmonic oscillator, plot of solution + derivative.
3. **Linear algebra** — `Eigensystem` of a symmetric tridiagonal matrix.
4. **Weber pair potential** — `U_W(r, ṙ)`, `∂U_W/∂r`, `∂U_W/∂ṙ` as typeset expressions and as LaTeX for the paper.
5. **Numeric exploration** — `ListPlot` of sampled trajectories.
6. **Paper export** — `TeXForm` → `.tex` snippet, `Plot` → PNG, `Rasterize[TraditionalForm[...]]` → PNG.

## Using it as a CAS for the paper

`TeXForm[expr]` converts any symbolic expression to a LaTeX string ready to paste
into [papers/Computational-Weber-Electrodynamics/](../papers/Computational-Weber-Electrodynamics/).
For structured exports the notebook writes a multi-line `align` snippet to
`out/snippets.tex` that can be `\input{}`'d directly.

For quick one-offs from the shell:

```bash
wolframscript -code 'TeXForm[Integrate[Sin[x]^2, x]]'
# → \frac{x}{2}-\frac{1}{4} \sin (2 x)
```

## Files

- [tour.ipynb](tour.ipynb) — the notebook (runs on `wolframlanguage14.3`).
- `out/` — generated `.tex` and `.png` artifacts (gitignored).

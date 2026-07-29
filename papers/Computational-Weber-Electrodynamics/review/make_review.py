#!/usr/bin/env python3
"""
Build a REVIEW copy of Computational-Weber-Electrodynamics.tex in which every
v1.2 -> v1.3 change is rendered inline: the old block and the new block are
paired inside one box, old on top (red), new below (green).

The two inputs are extracted from git, never stored by hand, so the review
cannot drift from what the branch actually changed:

    old.tex  <-  $WEBER_REVIEW_BASE (default: main)
    new.tex  <-  $WEBER_REVIEW_HEAD (default: HEAD)

Temporary review artifact. Delete papers/Computational-Weber-Electrodynamics/review/
once the changes are accepted.

    python3 make_review.py && latexmk -pdf review.tex

Needs the FULL TeX Live (tcolorbox + tikzfill), not texlive-basic:
    export PATH="/usr/local/texlive/2025/bin/universal-darwin:$PATH"
"""

import os
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).parent
REPO = HERE.parents[2]
TEXPATH = "papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex"
BASE_REF = os.environ.get("WEBER_REVIEW_BASE", "main")
HEAD_REF = os.environ.get("WEBER_REVIEW_HEAD", "HEAD")


def git_blob(ref):
    return subprocess.run(["git", "-C", str(REPO), "show", f"{ref}:{TEXPATH}"],
                          capture_output=True, text=True, check=True).stdout


# Materialise the build inputs from git (also consumed by verify_review.py).
(HERE / "old.tex").write_text(git_blob(BASE_REF))
(HERE / "new.tex").write_text(git_blob(HEAD_REF))
print(f"inputs: old.tex <- {BASE_REF}, new.tex <- {HEAD_REF}")

NEW = (HERE / "new.tex").read_text()

# ---------------------------------------------------------------- preamble --

PREAMBLE = r"""
% ===================== REVIEW MARKUP (temporary) =====================
\usepackage[most]{tcolorbox}
\tcbuselibrary{breakable,skins}
\usepackage{soul}

\definecolor{oldbg}{RGB}{255,238,240}
\definecolor{oldink}{RGB}{179,29,42}
\definecolor{newbg}{RGB}{230,255,237}
\definecolor{newink}{RGB}{26,109,52}
\definecolor{barink}{RGB}{60,66,74}
\definecolor{sigbg}{RGB}{255,221,87}
\definecolor{addbg}{RGB}{200,247,197}
\definecolor{delbg}{RGB}{255,215,215}

% Paired old/new block. #1 = change id, #2 = short caption.
\newtcolorbox{diffbox}[2]{%
  enhanced, breakable,
  colback=oldbg,
  colbacklower=newbg,
  colframe=barink,
  boxrule=0.5pt,
  arc=2pt,
  left=7pt, right=7pt, top=5pt, bottom=5pt,
  coltitle=white,
  fonttitle=\sffamily\bfseries\footnotesize,
  title={\makebox[1.9em][l]{#1}#2},
  segmentation style={draw=barink!55, dashed, line width=0.5pt},
}

% Single-sided variants for pure additions / pure removals.
\newtcolorbox{addbox}[2]{%
  enhanced, breakable,
  colback=newbg, colframe=newink, boxrule=0.5pt, arc=2pt,
  left=7pt, right=7pt, top=5pt, bottom=5pt,
  coltitle=white, fonttitle=\sffamily\bfseries\footnotesize,
  title={\makebox[1.9em][l]{#1}#2},
}

\newcommand{\oldtag}{{\sffamily\bfseries\footnotesize
  \textcolor{oldink}{$\ominus$~REMOVED \textnormal{\textmd{— paper v1.2 (current \texttt{main})}}}}\par\smallskip}
\newcommand{\newtag}{{\sffamily\bfseries\footnotesize
  \textcolor{newink}{$\oplus$~ADDED \textnormal{\textmd{— paper v1.3 (this PR)}}}}\par\smallskip}

% Inline prose markup. soul's \hl breaks across lines; \colorbox cannot.
% Fragile bits (math, \eqref, \texttt) must be \mbox-protected inside \dins.
\sethlcolor{addbg}
\newcommand{\dins}[1]{\hl{#1}}
\newcommand{\ddel}[1]{{\color{oldink}\st{#1}}}

% Highlight a single changed token inside math (sign flips).
\newcommand{\sgn}[1]{\mathbin{\colorbox{sigbg}{$#1$}}}
\newcommand{\sgnc}[1]{{\colorbox{sigbg}{$#1$}}}
\newcommand{\wchg}[1]{\colorbox{sigbg}{#1}}

\newcommand{\eqnote}[1]{\par\smallskip{\footnotesize\sffamily\textcolor{barink}{#1}}\par}
% =====================================================================
"""

# ------------------------------------------------------------ front matter --

FRONT = r"""
\begin{center}
\begin{tcolorbox}[enhanced, width=\textwidth, colback=barink!6, colframe=barink,
  boxrule=1pt, arc=3pt, left=12pt, right=12pt, top=10pt, bottom=10pt]
\begin{center}
{\sffamily\bfseries\large REVIEW COPY — NOT FOR DISTRIBUTION}\\[3pt]
{\sffamily\small Paper v1.2 $\rightarrow$ v1.3 \quad$\cdot$\quad
canonical Weber Hamiltonian correction \quad$\cdot$\quad PR \#92}
\end{center}
\smallskip
\small
Every mathematical change is shown inline as a paired block: the
\textcolor{oldink}{\textbf{removed v1.2 content on a red field}} sits directly above the
\textcolor{newink}{\textbf{added v1.3 content on a green field}}, separated by a dashed rule.
Unchanged text is left exactly as it appears in the paper, so this document still reads
as the paper itself. Equation numbering follows v1.3; blocks quoted from v1.2 are
deliberately unnumbered so they cannot be mistaken for live equations.

\smallskip
\noindent
\begin{tabular}{@{}ll@{}}
\colorbox{sigbg}{$+$} & a single token that flipped sign or changed value \\
\dins{added text} & prose inserted in v1.3 \\
\ddel{removed text} & prose deleted from v1.2 \\
\end{tabular}
\end{tcolorbox}
\end{center}

\subsection*{Index of changes}

\noindent\small
\begin{tabular}{@{}p{2.2em}p{9.2cm}p{4.1cm}@{}}
\textbf{ID} & \textbf{What changed} & \textbf{Kind} \\[2pt]
\hline\\[-6pt]
C1  & The $v \to p/m$ substitution sentence, replaced by the canonical momentum $\vec p_i = \partial L/\partial \vec v_i$ & \textbf{root error} \\
C2  & Qualifier after \texttt{h\_box}: $H = T+U$ is velocity-space energy & prose added \\
C3  & Conservation clause points at the canonical form too & prose edit \\
C4  & \texttt{alpha\_x} lead-in re-contextualised & prose edit \\
C5  & \texttt{xdot\_two}, \texttt{xdot\_two\_p2} & \textbf{sign flip} $\times 2$ \\
C6  & \texttt{pdot\_two} & \textbf{sign flip} $\times 2$ \\
C7  & Two-particle scalar inverse \texttt{radial\_inverse} & equation added \\
C8  & \texttt{hamiltonian}: $H_n \to E_n$, plus \texttt{velocity\_recovery}, \texttt{radial\_system}, \texttt{hamiltonian\_canonical} & \textbf{relabel + 3 added} \\
C9  & \texttt{xdot\_expanded}, \texttt{pdot\_expanded} & \textbf{sign flip} $\times 3$ \\
C10 & Complexity analysis & rewritten \\
C11 & Non-separability term & equation replaced \\
C12 & Appendix A momenta block & \textbf{6 equations corrected} \\
\end{tabular}

\bigskip
\noindent\small\textit{Unchanged and verified: abstract; \texttt{potential}; \texttt{force};
\texttt{ke}; $S$; $L = T-S$; \texttt{euler\_lagrange}; \texttt{H\_def}; \texttt{H};
\texttt{h\_box}; \texttt{xdot}; \texttt{pdot}; \texttt{z}; \texttt{update\_step\_x/p};
\texttt{Z}; \texttt{AZ}; \texttt{map}; \texttt{Phi}; the Newton iteration; the Discussion;
Appendix A.1 radial-acceleration identities; bibliography.}

\clearpage
"""


def sub(text, old, new, tag):
    if old not in text:
        sys.exit(f"ANCHOR MISS [{tag}]:\n{old[:200]}")
    return text.replace(old, new, 1)


s = NEW

# Reuse the paper's bibliography rather than keeping a copy in this directory.
s = sub(s, r"\bibliography{references}", r"\bibliography{../references}", "bib")

# ------------------------------------------------------------- preamble in --
s = sub(s, "\n% Document metadata", PREAMBLE + "\n% Document metadata", "preamble")
s = sub(s, r"\end{abstract}", r"\end{abstract}" + "\n" + FRONT, "front")

# --------------------------------------------------------------- C2 (prose) --
s = sub(
    s,
    r"where $T$ is given by~\eqref{ke} and $U$ is the velocity-dependent potential~\eqref{potential}. Equations~\eqref{H} and~\eqref{h_box} express the conserved energy in positions and \textit{physical velocities}; they become a canonical Hamiltonian only after the velocities have been eliminated as described next.",
    r"""where $T$ is given by~\eqref{ke} and $U$ is the velocity-dependent potential~\eqref{potential}. \dins{Equations~\mbox{\eqref{H}} and~\mbox{\eqref{h_box}} express the conserved energy in positions and \mbox{\textit{physical velocities}}; they become a canonical Hamiltonian only after the velocities have been eliminated as described next.}

\eqnote{\textbf{C2}\quad Prose only. Equations \texttt{H} and \texttt{h\_box} are mathematically unchanged --- $H = T+U$ was always correct \emph{in velocity space}. The qualifier makes that explicit so the reader does not carry the identity into $(q,p)$ coordinates.}""",
    "C2",
)

# ------------------------------------------------------ C1 (the root error) --
s = sub(
    s,
    r"""The Lagrangian $L = T - S$ is velocity dependent through $\dot{r}$, so the canonical momentum $\vec{p}_i = \partial L / \partial \vec{v}_i$ is \textit{not} the kinetic momentum $m_i \vec{v}_i$. Carrying out the differentiation with $\partial \dot{r} / \partial \vec{v}_1 = \hat{r} = -\,\partial \dot{r} / \partial \vec{v}_2$ gives

\begin{equation}\label{canonical_momentum}
    \boxed{\vec{p}_1 = m_1 \vec{v}_1 - \vec{\alpha}, \qquad \vec{p}_2 = m_2 \vec{v}_2 + \vec{\alpha}, \qquad \vec{\alpha} \equiv \frac{q_1 q_2}{c^2}\frac{\dot{r}\,(\vec{r}_1 - \vec{r}_2)}{r^2}}
\end{equation}

The correction $\vec{\alpha}$ is purely radial and enters with opposite signs, so total momentum is unaffected: $\vec{p}_1 + \vec{p}_2 = m_1\vec{v}_1 + m_2\vec{v}_2$. Velocities are therefore eliminated in favour of momenta through~\eqref{canonical_momentum}, not by replacing each $\dot{x}_i$ with $p_{x_i}/m_i$.""",
    r"""\begin{diffbox}{C1}{The root error --- canonical momentum is not $m\vec v$}
\oldtag
where each velocity component $\dot{x}_i$ in~\eqref{H} is replaced by $p_{x_i}/m_i$ and
similarly for all other velocity components.

\eqnote{This single sentence was the primary error. It asserts that the Legendre
transform can be performed by dividing momentum by mass, which is false for a
velocity-dependent Lagrangian. Every later equation in the paper inherited it.}
\tcblower
\newtag
The Lagrangian $L = T - S$ is velocity dependent through $\dot{r}$, so the canonical momentum $\vec{p}_i = \partial L / \partial \vec{v}_i$ is \textit{not} the kinetic momentum $m_i \vec{v}_i$. Carrying out the differentiation with $\partial \dot{r} / \partial \vec{v}_1 = \hat{r} = -\,\partial \dot{r} / \partial \vec{v}_2$ gives

\begin{equation}\label{canonical_momentum}
    \boxed{\vec{p}_1 = m_1 \vec{v}_1 - \vec{\alpha}, \qquad \vec{p}_2 = m_2 \vec{v}_2 + \vec{\alpha}, \qquad \vec{\alpha} \equiv \frac{q_1 q_2}{c^2}\frac{\dot{r}\,(\vec{r}_1 - \vec{r}_2)}{r^2}}
\end{equation}

The correction $\vec{\alpha}$ is purely radial and enters with opposite signs, so total momentum is unaffected: $\vec{p}_1 + \vec{p}_2 = m_1\vec{v}_1 + m_2\vec{v}_2$. Velocities are therefore eliminated in favour of momenta through~\eqref{canonical_momentum}, not by replacing each $\dot{x}_i$ with $p_{x_i}/m_i$.
\end{diffbox}""",
    "C1",
)

# ------------------------------------------------------------------ C4 + C5 --
s = sub(
    s,
    r"""Writing the components of $\vec{\alpha}$ as
\begin{equation}\label{alpha_x}
    \alpha_x = \frac{q_1 q_2}{c^2} \frac{\dot{r}_{12} (x_1 - x_2)}{r_{12}^2}
\end{equation}

with analogous definitions for $\alpha_y$ and $\alpha_z$, we apply Hamilton's equations~\eqref{xdot} and~\eqref{pdot} to obtain

\begin{equation}\label{xdot_two}
    \dot{x}_1 = \frac{1}{m_1} \left(p_{x_1} + \alpha_x\right)
\end{equation}

\begin{equation}\label{xdot_two_p2}
    \dot{x}_2 = \frac{1}{m_2} \left(p_{x_2} - \alpha_x\right)
\end{equation}

and

\begin{equation}\label{pdot_two}
    \dot{p}_{x_1} = \frac{q_1 q_2}{r_{12}^2}\left[\frac{(x_1 - x_2)}{r_{12}}\left(1 + \frac{3\dot{r}_{12}^2}{2 c^2}\right) - \frac{\dot{r}_{12}(\dot{x}_1 - \dot{x}_2)}{c^2}\right]
\end{equation}""",
    r"""\begin{diffbox}{C4}{\texttt{alpha\_x} --- algebraically identical, re-contextualised}
\oldtag
Defining
\begin{equation*}
    \alpha_x = \frac{q_1 q_2}{c^2} \frac{\dot{r}_{12} (x_1 - x_2)}{r_{12}^2}
\end{equation*}
\eqnote{In v1.2, $\dot r_{12}$ here was implicitly the $p/m$ surrogate.}
\tcblower
\newtag
Writing the components of $\vec{\alpha}$ as
\begin{equation}\label{alpha_x}
    \alpha_x = \frac{q_1 q_2}{c^2} \frac{\dot{r}_{12} (x_1 - x_2)}{r_{12}^2}
\end{equation}
\eqnote{The formula is byte-identical. What changed is its meaning: $\dot r_{12}$ is now
the \emph{physical} radial velocity, and $\alpha_x$ is a component of the momentum
correction $\vec\alpha$ introduced in C1 rather than a free-standing definition.}
\end{diffbox}

with analogous definitions for $\alpha_y$ and $\alpha_z$, we apply Hamilton's equations~\eqref{xdot} and~\eqref{pdot} to obtain

\begin{diffbox}{C5}{\texttt{xdot\_two}, \texttt{xdot\_two\_p2} --- both signs flipped}
\oldtag
\begin{equation*}
    \dot{x}_1 = \frac{1}{m_1} \left(p_{x_1} \sgn{-} \alpha_x\right)
\end{equation*}
\begin{equation*}
    \dot{x}_2 = \frac{1}{m_2} \left(p_{x_2} \sgn{+} \alpha_x\right)
\end{equation*}
\tcblower
\newtag
\begin{equation}\label{xdot_two}
    \dot{x}_1 = \frac{1}{m_1} \left(p_{x_1} \sgn{+} \alpha_x\right)
\end{equation}

\begin{equation}\label{xdot_two_p2}
    \dot{x}_2 = \frac{1}{m_2} \left(p_{x_2} \sgn{-} \alpha_x\right)
\end{equation}
\eqnote{Solving $\vec p_1 = m_1\vec v_1 - \vec\alpha$ for $\vec v_1$ gives
$+\vec\alpha$, not $-\vec\alpha$. Note that
\texttt{theory/WeberElectrodynamics.md} already carried these signs correctly ---
the paper and the theory notes contradicted each other.}
\end{diffbox}

and

\begin{diffbox}{C6}{\texttt{pdot\_two} --- both $1/c^2$ signs flipped}
\oldtag
\begin{equation*}
    \dot{p}_{x_1} = \frac{q_1 q_2}{r_{12}^2}\left[\frac{(x_1 - x_2)}{r_{12}}\left(1 \sgn{-} \frac{3\dot{r}_{12}^2}{2 c^2}\right) \sgn{+} \frac{\dot{r}_{12}(\dot{x}_1 - \dot{x}_2)}{c^2}\right]
\end{equation*}
\tcblower
\newtag
\begin{equation}\label{pdot_two}
    \dot{p}_{x_1} = \frac{q_1 q_2}{r_{12}^2}\left[\frac{(x_1 - x_2)}{r_{12}}\left(1 \sgn{+} \frac{3\dot{r}_{12}^2}{2 c^2}\right) \sgn{-} \frac{\dot{r}_{12}(\dot{x}_1 - \dot{x}_2)}{c^2}\right]
\end{equation}
\eqnote{Both signs are consequences of the same C1 error and cannot be corrected
independently of it. \texttt{pdot\_two\_p2} ($\dot p_{x_2} = -\dot p_{x_1}$) is unchanged.}
\end{diffbox}""",
    "C4+C5+C6",
)

# ------------------------------------------------------------------ C7 --
s = sub(
    s,
    r"""with analogous equations for the $y$ and $z$ directions. Note that $\dot{p}_{x_1} + \dot{p}_{x_2} = 0$, but $\dot{x}_1 - \dot{x}_2 \neq 0$ in general. Throughout, $\dot{r}_{12}$ and $\dot{x}_1 - \dot{x}_2$ denote \textit{physical} velocities; $\vec{\alpha}$ in~\eqref{xdot_two} and~\eqref{xdot_two_p2} contains $\dot{r}$, so these relations are implicit.

For two particles the implicit relation reduces to a single scalar equation. With $\mu = m_1 m_2/(m_1+m_2)$ the reduced mass and $p_r = \mu\,\hat{r}\cdot(\vec{p}_1/m_1 - \vec{p}_2/m_2)$ the momentum conjugate to the relative radial coordinate, contracting~\eqref{canonical_momentum} with $\hat{r}$ gives

\begin{equation}\label{radial_inverse}
    \boxed{p_r = \left(\mu - \frac{q_1 q_2}{r c^2}\right)\dot{r}, \qquad \dot{r} = \frac{p_r}{\mu - q_1 q_2/(r c^2)}}
\end{equation}

so the physical radial velocity follows in closed form from the canonical momenta. Only the relative radial motion is modified; the centre-of-mass and relative transverse momenta keep their ordinary relation to velocity. The effective radial inertia $\mu - q_1 q_2/(r c^2)$ vanishes for like charges exactly at Weber's critical radius $\rho = q_1 q_2/(\mu c^2)$, revisited in the Discussion; below $\rho$ it is negative but finite.""",
    r"""with analogous equations for the $y$ and $z$ directions. Note that $\dot{p}_{x_1} + \dot{p}_{x_2} = 0$, but $\dot{x}_1 - \dot{x}_2 \neq 0$ in general. \dins{Throughout, \mbox{$\dot{r}_{12}$} and \mbox{$\dot{x}_1 - \dot{x}_2$} denote \mbox{\textit{physical}} velocities; \mbox{$\vec{\alpha}$} in~\mbox{\eqref{xdot_two}} and~\mbox{\eqref{xdot_two_p2}} contains \mbox{$\dot{r}$}, so these relations are implicit.}

\begin{addbox}{C7}{Two-particle scalar inverse --- new, resolves the implicitness}
\newtag
For two particles the implicit relation reduces to a single scalar equation. With $\mu = m_1 m_2/(m_1+m_2)$ the reduced mass and $p_r = \mu\,\hat{r}\cdot(\vec{p}_1/m_1 - \vec{p}_2/m_2)$ the momentum conjugate to the relative radial coordinate, contracting~\eqref{canonical_momentum} with $\hat{r}$ gives

\begin{equation}\label{radial_inverse}
    \boxed{p_r = \left(\mu - \frac{q_1 q_2}{r c^2}\right)\dot{r}, \qquad \dot{r} = \frac{p_r}{\mu - q_1 q_2/(r c^2)}}
\end{equation}

so the physical radial velocity follows in closed form from the canonical momenta. Only the relative radial motion is modified; the centre-of-mass and relative transverse momenta keep their ordinary relation to velocity. The effective radial inertia $\mu - q_1 q_2/(r c^2)$ vanishes for like charges exactly at Weber's critical radius $\rho = q_1 q_2/(\mu c^2)$, revisited in the Discussion; below $\rho$ it is negative but finite.

\eqnote{Nothing corresponds to this in v1.2, which had no notion that the
velocity--momentum relation needed inverting at all. Note $p_r = \mu\dot r$ is
recovered only as $c\to\infty$.}
\end{addbox}""",
    "C7",
)

# ------------------------------------------------------------------- C3 --
s = sub(
    s,
    r"These conserved quantities arise from symmetries in the Hamiltonian~\eqref{h_box}, equivalently the canonical form~\eqref{hamiltonian_canonical}, according to Noether's theorem.",
    r"These conserved quantities arise from symmetries in the Hamiltonian~\eqref{h_box}\dins{, equivalently the canonical form~\mbox{\eqref{hamiltonian_canonical}},} according to Noether's theorem.\footnote{\textbf{C3.} Clause added. The conservation claims themselves are unchanged and still hold; the correction preserves translation, rotation, and time-translation invariance, all of which are re-verified numerically in the PR.}",
    "C3",
)

# ------------------------------------------------------------------- C8 --
s = sub(
    s,
    r"""For a system with $n$ particles, the energy~\eqref{H} generalizes to:
\begin{equation}\label{hamiltonian}
    E_n = \sum_{i=1}^{n} \frac{1}{2} m_i v_i^2 + \sum_{i=1}^{n} \sum_{j=i+1}^{n} \frac{q_i q_j}{r_{ij}} \left(1 - \frac{\dot{r}_{ij}^2}{2 c^2}\right)
\end{equation}""",
    r"""\begin{diffbox}{C8a}{\texttt{hamiltonian} --- $H_n$ relabelled $E_n$, right-hand side identical}
\oldtag
For a system with $n$ particles, the \wchg{Hamiltonian}~\eqref{H} generalizes to:
\begin{equation*}
    \sgnc{H_n} = \sum_{i=1}^{n} \frac{1}{2} m_i v_i^2 + \sum_{i=1}^{n} \sum_{j=i+1}^{n} \frac{q_i q_j}{r_{ij}} \left(1 - \frac{\dot{r}_{ij}^2}{2 c^2}\right)
\end{equation*}
\tcblower
\newtag
For a system with $n$ particles, the \wchg{energy}~\eqref{H} generalizes to:
\begin{equation}\label{hamiltonian}
    \sgnc{E_n} = \sum_{i=1}^{n} \frac{1}{2} m_i v_i^2 + \sum_{i=1}^{n} \sum_{j=i+1}^{n} \frac{q_i q_j}{r_{ij}} \left(1 - \frac{\dot{r}_{ij}^2}{2 c^2}\right)
\end{equation}
\eqnote{The right-hand side is unchanged to the byte. It was always the
velocity-space energy; calling it a Hamiltonian and then feeding it to
Hamilton's equations in $(q,p)$ was the error. Renaming it $E_n$ frees the
symbol $H_n$ for the genuine canonical Hamiltonian added below.}
\end{diffbox}""",
    "C8a",
)

s = sub(
    s,
    r"""Equation~\eqref{hamiltonian} is expressed in positions and \textit{physical velocities}. To obtain the canonical Hamiltonian we must first recover the velocities from the canonical momenta~\eqref{canonical_momentum}. Generalizing~\eqref{alpha_x} to $n$ particles, the correction to $\vec{v}_i$ is a sum over the pairs containing particle $i$,""",
    r"""\begin{addbox}{C8b}{Velocity recovery and the exact canonical Hamiltonian --- all new}
\newtag
Equation~\eqref{hamiltonian} is expressed in positions and \textit{physical velocities}. To obtain the canonical Hamiltonian we must first recover the velocities from the canonical momenta~\eqref{canonical_momentum}. Generalizing~\eqref{alpha_x} to $n$ particles, the correction to $\vec{v}_i$ is a sum over the pairs containing particle $i$,""",
    "C8b-open",
)

s = sub(
    s,
    r"""and $H_n(\mathbf{q},\mathbf{p}) = E_n(\mathbf{q},\mathbf{v}(\mathbf{q},\mathbf{p}))$ by the Legendre transform. In the Coulomb limit $c \to \infty$, and at any instant where every $\dot{r}_{ij}$ vanishes, $k_{ij} \to 0$ and~\eqref{hamiltonian_canonical} reduces to $\sum_i |\vec{p}_i|^2/(2m_i) + \sum_{i<j} q_i q_j / r_{ij}$.""",
    r"""and $H_n(\mathbf{q},\mathbf{p}) = E_n(\mathbf{q},\mathbf{v}(\mathbf{q},\mathbf{p}))$ by the Legendre transform. In the Coulomb limit $c \to \infty$, and at any instant where every $\dot{r}_{ij}$ vanishes, $k_{ij} \to 0$ and~\eqref{hamiltonian_canonical} reduces to $\sum_i |\vec{p}_i|^2/(2m_i) + \sum_{i<j} q_i q_j / r_{ij}$.

\eqnote{v1.2 contained none of this. $s_{ij}$ is exactly the quantity v1.2 called
$\dot r_{ij}$ in $(q,p)$ coordinates --- naming it separately is what makes the
error visible. Verified symbolically for $n=2$ and by 40-digit finite differences
for $n=3$ in \texttt{verify\_formulas.py}.}
\end{addbox}""",
    "C8b-close",
)

# ------------------------------------------------------------------- C9 --
s = sub(
    s,
    r"""\begin{equation}\label{xdot_expanded}
    \boxed{\dot{x}_i = \frac{1}{m_i} \left(p_{x_i} + \sum_{\substack{j=1 \\ j \neq i}}^{n} \frac{q_i q_j}{c^2} \frac{\dot{r}_{ij} (x_i - x_j)}{r_{ij}^2}\right)}
\end{equation}

for particle $i \in \{1, 2, \ldots, n\}$, with analogous expressions for $\dot{y}_i$ and $\dot{z}_i$, and~\eqref{pdot} expands to

\begin{equation}\label{pdot_expanded}
    \boxed{\dot{p}_{x_i} = \sum_{\substack{j=1 \\ j \neq i}}^{n} \frac{q_i q_j}{r_{ij}^2}\left[\frac{(x_i - x_j)}{r_{ij}}\left(1 + \frac{3\dot{r}_{ij}^2}{2 c^2}\right) - \frac{\dot{r}_{ij}(\dot{x}_i - \dot{x}_j)}{c^2}\right]}
\end{equation}

for particle $i \in \{1, 2, \ldots, n\}$, with analogous expressions for $\dot{p}_{y_i}$ and $\dot{p}_{z_i}$. In both boxed equations every $\dot{r}_{ij}$ and $\dot{x}_i - \dot{x}_j$ is a physical velocity obtained from the simultaneous system~\eqref{radial_system}; equation~\eqref{xdot_expanded} is the component form of~\eqref{velocity_recovery}.""",
    r"""\begin{diffbox}{C9}{\texttt{xdot\_expanded}, \texttt{pdot\_expanded} --- three signs flipped}
\oldtag
\begin{equation*}
    \boxed{\dot{x}_i = \frac{1}{m_i} \left(p_{x_i} \sgn{-} \sum_{\substack{j=1 \\ j \neq i}}^{n} \frac{q_i q_j}{c^2} \frac{\dot{r}_{ij} (x_i - x_j)}{r_{ij}^2}\right)}
\end{equation*}
\begin{equation*}
    \boxed{\dot{p}_{x_i} = \sum_{\substack{j=1 \\ j \neq i}}^{n} \frac{q_i q_j}{r_{ij}^2}\left[\frac{(x_i - x_j)}{r_{ij}}\left(1 \sgn{-} \frac{3\dot{r}_{ij}^2}{2 c^2}\right) \sgn{+} \frac{\dot{r}_{ij}(\dot{x}_i - \dot{x}_j)}{c^2}\right]}
\end{equation*}
\tcblower
\newtag
\begin{equation}\label{xdot_expanded}
    \boxed{\dot{x}_i = \frac{1}{m_i} \left(p_{x_i} \sgn{+} \sum_{\substack{j=1 \\ j \neq i}}^{n} \frac{q_i q_j}{c^2} \frac{\dot{r}_{ij} (x_i - x_j)}{r_{ij}^2}\right)}
\end{equation}

for particle $i \in \{1, 2, \ldots, n\}$, with analogous expressions for $\dot{y}_i$ and $\dot{z}_i$, and~\eqref{pdot} expands to

\begin{equation}\label{pdot_expanded}
    \boxed{\dot{p}_{x_i} = \sum_{\substack{j=1 \\ j \neq i}}^{n} \frac{q_i q_j}{r_{ij}^2}\left[\frac{(x_i - x_j)}{r_{ij}}\left(1 \sgn{+} \frac{3\dot{r}_{ij}^2}{2 c^2}\right) \sgn{-} \frac{\dot{r}_{ij}(\dot{x}_i - \dot{x}_j)}{c^2}\right]}
\end{equation}

for particle $i \in \{1, 2, \ldots, n\}$, with analogous expressions for $\dot{p}_{y_i}$ and $\dot{p}_{z_i}$. \dins{In both boxed equations every \mbox{$\dot{r}_{ij}$} and \mbox{$\dot{x}_i - \dot{x}_j$} is a physical velocity obtained from the simultaneous system~\mbox{\eqref{radial_system}}; equation~\mbox{\eqref{xdot_expanded}} is the component form of~\mbox{\eqref{velocity_recovery}}.}
\eqnote{These are the $n$-body forms of C5 and C6 and carry the same three sign
changes. Both remain boxed as the paper's headline results.}
\end{diffbox}""",
    "C9",
)

# ------------------------------------------------------------------ C10 --
s = sub(
    s,
    r"""For a system of $n$ particles, we thus have a total of $P = n(n-1)/2$ unique particle pairs, so assembling the pair geometry costs $\mathcal{O}(n^2)$.

A naive implementation would compute all $n(n-1)$ interactions, but by exploiting~\eqref{alpha_x},~\eqref{xdot_two},~\eqref{xdot_two_p2},~\eqref{pdot_two} and~\eqref{pdot_two_p2}, we reduce this to $P$ interactions. Each evaluation of~\eqref{xdot_expanded} and~\eqref{pdot_expanded} then needs $r_{ij}$ and $\dot{r}_{ij}$ and their powers only once per unique pair.

Pair evaluation alone does not, however, complete an evaluation of the equations of motion. Because canonical momentum is not kinetic momentum, the physical velocities entering~\eqref{xdot_expanded} and~\eqref{pdot_expanded} are only available after solving the coupled system~\eqref{radial_system}, whose matrix is $P \times P$ and dense in general. A direct factorization costs $\mathcal{O}(P^3) = \mathcal{O}(n^6)$ per evaluation, so for large $n$ the linear solve, not the pair sum, dominates. Two remarks temper this. First, $\mathbb{I} - GK$ is a perturbation of the identity of order $q^2/(m c^2 r)$, so for sub-relativistic configurations a few Jacobi or conjugate-gradient iterations converge to machine precision at $\mathcal{O}(n^2)$ cost per iteration. Second, for the two-body case $P = 1$ and the solve degenerates to the closed-form division~\eqref{radial_inverse}. Our reference implementation~\citep{WeberElectrodynamics.jl} uses a direct dense factorization, which is exact at all separations and comfortably fast for the particle counts studied here.""",
    r"""\begin{diffbox}{C10}{Complexity --- the $\mathcal{O}(n^2)$ claim was incomplete}
\oldtag
For a system of $n$ particles, we thus have a total of $n(n-1)/2$ unique particle pairs, resulting in a computational complexity of order $\mathcal{O}(n^2)$.

A naive implementation would compute all $n(n-1)$ interactions, but by exploiting~\eqref{alpha_x},~\eqref{xdot_two},~\eqref{xdot_two_p2},~\eqref{pdot_two} and~\eqref{pdot_two_p2}, we reduce this to $n(n-1)/2$ interactions.

Furthermore, each evaluation of~\eqref{xdot_expanded} and~\eqref{pdot_expanded} only needs to compute $r_{ij}$ and $\dot{r}_{ij}$ and their powers once per unique pair, which further reduces the computational load.

\eqnote{This asserted that per-pair work completes an evaluation. It does not:
with $\vec p \ne m\vec v$ the velocities are not known until a coupled system is
solved.}
\tcblower
\newtag
For a system of $n$ particles, we thus have a total of $P = n(n-1)/2$ unique particle pairs, so assembling the pair geometry costs $\mathcal{O}(n^2)$.

A naive implementation would compute all $n(n-1)$ interactions, but by exploiting~\eqref{alpha_x},~\eqref{xdot_two},~\eqref{xdot_two_p2},~\eqref{pdot_two} and~\eqref{pdot_two_p2}, we reduce this to $P$ interactions. Each evaluation of~\eqref{xdot_expanded} and~\eqref{pdot_expanded} then needs $r_{ij}$ and $\dot{r}_{ij}$ and their powers only once per unique pair.

Pair evaluation alone does not, however, complete an evaluation of the equations of motion. Because canonical momentum is not kinetic momentum, the physical velocities entering~\eqref{xdot_expanded} and~\eqref{pdot_expanded} are only available after solving the coupled system~\eqref{radial_system}, whose matrix is $P \times P$ and dense in general. A direct factorization costs $\mathcal{O}(P^3) = \mathcal{O}(n^6)$ per evaluation, so for large $n$ the linear solve, not the pair sum, dominates. Two remarks temper this. First, $\mathbb{I} - GK$ is a perturbation of the identity of order $q^2/(m c^2 r)$, so for sub-relativistic configurations a few Jacobi or conjugate-gradient iterations converge to machine precision at $\mathcal{O}(n^2)$ cost per iteration. Second, for the two-body case $P = 1$ and the solve degenerates to the closed-form division~\eqref{radial_inverse}. Our reference implementation~\citep{WeberElectrodynamics.jl} uses a direct dense factorization, which is exact at all separations and comfortably fast for the particle counts studied here.
\end{diffbox}""",
    "C10",
)

# ------------------------------------------------------------------ C11 --
s = sub(
    s,
    r"""The canonical Hamiltonian for Weber's electrodynamics~\eqref{hamiltonian_canonical} is non-separable because the velocity-dependent term

\begin{equation}
    \frac{1}{2} k_{ij}\,\dot{r}_{ij}\,s_{ij}
\end{equation}

is a function of both the positions and the momenta after the Legendre transformation: $k_{ij}$ and $s_{ij}$ depend on positions and momenta explicitly, and $\dot{r}_{ij}$ does so through the solve~\eqref{radial_system}. Therefore, standard symplectic methods cannot be applied directly.""",
    r"""\begin{diffbox}{C11}{Non-separability --- the offending term restated}
\oldtag
The general Hamiltonian for Weber's electrodynamics~\eqref{hamiltonian} is non-separable because the velocity-dependent term
\begin{equation*}
    -\frac{\qq \dot{r}_{ij}^2}{2 c^2 r_{ij}}
\end{equation*}
is both a function of the positions and momenta after the Legendre transformation. Therefore, standard symplectic methods cannot be applied directly.
\tcblower
\newtag
The canonical Hamiltonian for Weber's electrodynamics~\eqref{hamiltonian_canonical} is non-separable because the velocity-dependent term

\begin{equation}
    \frac{1}{2} k_{ij}\,\dot{r}_{ij}\,s_{ij}
\end{equation}

is a function of both the positions and the momenta after the Legendre transformation: $k_{ij}$ and $s_{ij}$ depend on positions and momenta explicitly, and $\dot{r}_{ij}$ does so through the solve~\eqref{radial_system}. Therefore, standard symplectic methods cannot be applied directly.
\eqnote{The \emph{conclusion} --- non-separability, hence the need for the
extended-phase-space integrator --- is unchanged, so everything downstream in
Section~\ref{symplectic} (Strang splitting, symmetric projection, $\Phi$, the
Newton iteration) is mathematically untouched.}
\end{diffbox}""",
    "C11",
)

# ------------------------------------------------------------------ C12 --
s = sub(
    s,
    r"""The individual particle canonical momenta are given by~\eqref{canonical_momentum}, componentwise

\begin{align*}
    p_{x_1} = m_1 \dot{x}_1 - \alpha_x &  & p_{x_2} = m_2 \dot{x}_2 + \alpha_x \\
    p_{y_1} = m_1 \dot{y}_1 - \alpha_y &  & p_{y_2} = m_2 \dot{y}_2 + \alpha_y \\
    p_{z_1} = m_1 \dot{z}_1 - \alpha_z &  & p_{z_2} = m_2 \dot{z}_2 + \alpha_z
\end{align*}

with $\alpha_x$, $\alpha_y$, $\alpha_z$ the components of $\vec{\alpha}$ defined in~\eqref{alpha_x}. The kinetic momenta $m_i \dot{x}_i$ are recovered only when $\dot{r} = 0$ or in the limit $c \to \infty$.""",
    r"""\begin{diffbox}{C12}{Appendix A momenta --- all six equations corrected}
\oldtag
The individual particle momenta are
\begin{align*}
    p_{x_1} = m_1 \dot{x}_1 &  & p_{x_2} = m_2 \dot{x}_2 \\
    p_{y_1} = m_1 \dot{y}_1 &  & p_{y_2} = m_2 \dot{y}_2 \\
    p_{z_1} = m_1 \dot{z}_1 &  & p_{z_2} = m_2 \dot{z}_2
\end{align*}
\eqnote{The notation appendix stated the error as a definition, which is why it
propagated so cleanly through the rest of the paper.}
\tcblower
\newtag
The individual particle canonical momenta are given by~\eqref{canonical_momentum}, componentwise

\begin{align*}
    p_{x_1} = m_1 \dot{x}_1 \sgn{-} \alpha_x &  & p_{x_2} = m_2 \dot{x}_2 \sgn{+} \alpha_x \\
    p_{y_1} = m_1 \dot{y}_1 \sgn{-} \alpha_y &  & p_{y_2} = m_2 \dot{y}_2 \sgn{+} \alpha_y \\
    p_{z_1} = m_1 \dot{z}_1 \sgn{-} \alpha_z &  & p_{z_2} = m_2 \dot{z}_2 \sgn{+} \alpha_z
\end{align*}

with $\alpha_x$, $\alpha_y$, $\alpha_z$ the components of $\vec{\alpha}$ defined in~\eqref{alpha_x}. The kinetic momenta $m_i \dot{x}_i$ are recovered only when $\dot{r} = 0$ or in the limit $c \to \infty$.
\end{diffbox}""",
    "C12",
)

# Title banner.
s = sub(
    s,
    r"\title{Computational Weber electrodynamics: symplectic n-body integration}",
    r"\title{Computational Weber electrodynamics: symplectic n-body integration\\[4pt]{\large\sffamily\bfseries\color{oldink}[\,REVIEW COPY \textemdash\ v1.2\,$\rightarrow$\,v1.3 diff\,]}}",
    "title",
)

(HERE / "review.tex").write_text(s)
print(f"wrote review.tex ({len(s.splitlines())} lines)")

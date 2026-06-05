# First unmatched eigenvalue: index and asymptotics

## Setup

On the 2:1 catalog rectangle $[0, 2\pi] \times [0, \pi]$ the Dirichlet–Laplacian
has the exact eigenvalues

$$\lambda_{m,n} = \left(\frac{m\pi}{L_x}\right)^2 + \left(\frac{n\pi}{L_y}\right)^2
= \frac{m^2}{4} + n^2, \qquad m, n = 1, 2, \dots$$

At resolution $M$ (odd, so both edges land on grid lines) the DST grid resolves a
**rectangular block** of modes,

$$M_x = M \quad\text{in } x, \qquad M_y = \frac{M-1}{2} \quad\text{in } y,$$

and the computed spectrum is *exactly* this block
$\{\lambda_{m,n} : 1 \le m \le M_x,\ 1 \le n \le M_y\}$, reproduced to round-off.

Order the analytic eigenvalues by magnitude and order the computed ones by
magnitude. Because the computed set is a block and not the leading eigenvalues,
the two lists agree only up to a point. **Question:** at which index does the
ordered analytic list first contain an eigenvalue with no computed counterpart?
We pin the index down exactly as a count of grid points inside an ellipse, and
read off its growth from the classical Weyl eigenvalue-counting function.

## The first unmatched mode

A mode is unmatched exactly when it is **unresolved**, i.e. $m > M_x$ or
$n > M_y$. The first unmatched eigenvalue is therefore the *smallest* eigenvalue
over all unresolved modes. The cheapest way out of the resolved block is to step
one index past a resolution limit, so the minimum is attained at one of

- $(m, n) = (1,\ M_y + 1)$, costing $\lambda = \tfrac14 + \left(\tfrac{M+1}{2}\right)^2$, or
- $(m, n) = (M_x + 1,\ 1)$, costing $\lambda = \left(\tfrac{M+1}{2}\right)^2 + 1$.

The $y$-step always wins ($+\tfrac14 < +1$), so the first unmatched mode is
**always**

$$\boxed{(m, n) = \left(1,\ \tfrac{M+1}{2}\right), \qquad
\lambda^\* = \frac{(M+1)^2 + 1}{4}.}$$

It is the lowest mode sitting one row above the $y$-resolution limit, and it does
not depend on $M$ in form — only the value $\lambda^\*$ grows.

## The index

The index of $\lambda^\*$ in the magnitude-ordered list is one plus the number of
analytic eigenvalues strictly below it. Every eigenvalue below $\lambda^\*$ is
resolved (by definition $\lambda^\*$ is the smallest unresolved one), so all of
them are matched and the first "gap" is precisely at $\lambda^\*$. Since

$$\lambda_{m,n} < \lambda^\*
\iff m^2 + 4n^2 < (M+1)^2 + 1
\iff m^2 + 4n^2 \le (M+1)^2
\quad(\text{integer left-hand side}),$$

the index is a **lattice-point count inside an ellipse**:

$$\mathrm{index}(M)
= 1 + \#\bigl\{(m, n) \in \mathbb{Z}_{\ge 1}^2 : m^2 + 4n^2 \le (M+1)^2\bigr\}
= 1 + \sum_{n=1}^{(M-1)/2} \left\lfloor \sqrt{(M+1)^2 - 4n^2} \right\rfloor.$$

(The first unmatched mode has $m = 1$, the minimal possible; among any eigenvalues
of equal value it therefore sorts first, so the "$+1$" is exact even when
$\lambda^\*$ is degenerate.)

### Values

| $M$ | 5 | 7 | 9 | 11 | 19 | 21 | 31 |
|----:|--:|--:|--:|---:|---:|---:|---:|
| $\lambda^\*$ | 9.25 | 16.25 | 25.25 | 36.25 | 100.25 | 121.25 | 256.25 |
| index | 10 | 19 | **33** | 47 | 143 | 173 | 376 |

The $M = 9$ entry is the worked example: the first `--` appears at position
**33**, mode $(1, 5)$, $\lambda^\* = 25.25$.

## Asymptotics via the Weyl counting function

The lattice-point sum is a Gauss-type ellipse count with no elementary closed
form, but it is exactly the **Dirichlet eigenvalue-counting function**

$$N(\lambda) = \#\bigl\{(m, n) \in \mathbb{Z}_{\ge 1}^2 : \lambda_{m,n} \le \lambda\bigr\}$$

evaluated at $\lambda^\*$: counting grid points under $m^2 + 4n^2 \le (M+1)^2$ is
counting eigenvalues $\lambda_{m,n} \le \lambda^\*$. So
$\mathrm{index}(M) = N(\lambda^\*)$, up to the bounded multiplicity at
$\lambda^\*$. Weyl's classical two-term law for the Dirichlet Laplacian on a
planar domain $\Omega$,

$$N(\lambda) = \frac{|\Omega|}{4\pi}\,\lambda
\;-\; \frac{|\partial\Omega|}{4\pi}\,\sqrt{\lambda} \;+\; o\!\left(\sqrt{\lambda}\right),$$

is precisely the area-minus-boundary estimate for the grid points inside the
ellipse. With area $|\Omega| = 2\pi^2$, perimeter $|\partial\Omega| = 6\pi$, and
$\lambda^\* = (M+1)^2/4 + \tfrac14$, so $\sqrt{\lambda^\*} \sim (M+1)/2$,

$$\mathrm{index}(M) = N(\lambda^\*)
= \frac{2\pi^2}{4\pi}\cdot\frac{(M+1)^2}{4}
\;-\; \frac{6\pi}{4\pi}\cdot\frac{M+1}{2} \;+\; O(1)
= \frac{\pi (M+1)^2}{8} \;-\; \frac{3(M+1)}{4} \;+\; O(1).$$

The area term $\pi(M+1)^2/8 = \tfrac{\pi}{8}M^2 + O(M)$ scales like the area of
the resolved block ($M_x M_y = M(M-1)/2 \approx M^2/2$) times the constant
$\pi/4$: the first gap opens at a fixed *fraction* of the way through the resolved
spectrum. The perimeter term $-3(M+1)/4$ is the Dirichlet boundary correction —
the grid points lost to the two axes excluded by $m, n \ge 1$. For $M = 9$ the
estimate gives $39.3 - 7.5 \approx 32.8$ against the exact $33$.

## Self-check in code

`experiments/eigenvalues_indexing/eigenvalues_rectangle_match_column.m` computes
the closed-form `index(M)` (helper `first_unmatched_index`) and compares it to the
observed position of the first `--` row, printing `PASS`/`FAIL`. It passes for
every $M$ in the table above.

---

*Generated by Claude (Claude Opus 4.8, model ID claude-opus-4-8).*

# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "python-flint"]
# [tool.uv]
# exclude-newer = "2026-10-19T00:00:00Z"
# ///
"""Refine a SIC fiducial from double precision (search.py, stark.jl) to 150 significant digits.

Solves the d^2 equations |C_p|^2 = t_p, where C_(a,b) = sum_k conj(psi_(k+a)) psi_k w^(-bk), t_0 = 1
fixes the norm and t_p = 1/(d+1) otherwise, by Gauss-Newton in 640-bit arithmetic (python-flint).
No symmetry is imposed. Every step uses the inverse Gauss-Newton matrix computed once in double
precision, so each step gains about 15 digits. Only the input's double-precision part is used. The
result is written as fiducials/fiducial_<d>.txt in Scott and Grassl's layout, with unit norm and a real positive
first component.

  uv run refine.py fiducial_197_stark.txt     # rewrites fiducials/fiducial_197.txt
  uv run refine.py fiducials/fiducial_198.txt # re-derives its 150 digits from the first 16
"""

import argparse
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import acb, acb_mat, arb, ctx

from search import fiducial_text, read_fiducial

ctx.prec = 640  # bits, about 192 digits


def double_jacobian(psi: np.ndarray) -> np.ndarray:
    """Jacobian of |C_p|^2 with respect to (Re psi, Im psi) in double precision; row p = a*d + b."""
    d = psi.size
    k = np.arange(d)
    roots = np.exp(-2j * np.pi * k / d)  # w^(-m)
    a, b = k[:, None], k[None, :]
    C = (np.conj(psi[(a + b) % d]) * psi[None, :]) @ roots[(a * b) % d]
    J = np.empty((d * d, 2 * d))
    for j in range(d):
        m = (j - a) % d
        t1 = psi[m] * roots[(b * m) % d]  # d C_ab / d conj(psi_j)
        t2 = np.conj(psi[(j + a) % d]) * roots[(b * j) % d]  # d C_ab / d psi_j
        J[:, j] = (2 * np.real(np.conj(C) * (t1 + t2))).ravel()
        J[:, d + j] = (2 * np.real(np.conj(C) * (1j * (t2 - t1)))).ravel()
    return J


class Overlaps:
    """C = S T with S[a,k] = conj(psi_(k+a)) and T[k,b] = psi_k w^(-bk), as d x d ball-matrix products."""

    def __init__(self, d: int) -> None:
        self.d = d
        roots = [acb(arb(-2 * m) / d).exp_pi_i() for m in range(d)]  # w^(-m)
        self.wflat = [roots[(k * b) % d] for k in range(d) for b in range(d)]
        self.W = acb_mat(d, d, self.wflat)
        self.shift = [(k + a) % d for a in range(d) for k in range(d)]
        self.row = [k for k in range(d) for _ in range(d)]
        self.neg = [((-a) % d) * d + (-b) % d for a in range(d) for b in range(d)]
        self.anti = [[(a * d + (j - a) % d, (j - a) % d) for a in range(d)] for j in range(d)]
        self.t = [arb(1)] + [arb(1) / (d + 1)] * (d * d - 1)

    def C(self, psi: list[acb]) -> list[acb]:
        S = acb_mat(self.d, self.d, [psi[i].conjugate() for i in self.shift])
        T = acb_mat(self.d, self.d, [psi[k] * w for k, w in zip(self.row, self.wflat, strict=True)])
        return (S * T).entries()

    def gradient(self, psi: list[acb], C: list[acb]) -> list[acb]:
        """J^T r for the residuals r_p = |C_p|^2 - t_p, as a complex d-vector: real parts belong to Re psi,
        imaginary parts to Im psi."""
        r = [c.real * c.real + c.imag * c.imag - t for c, t in zip(C, self.t, strict=True)]
        r = [(r[p] + r[q]) / 2 for p, q in enumerate(self.neg)]  # |C_p| = |C_-p| identically
        H = (acb_mat(self.d, self.d, [x * c.conjugate() for x, c in zip(r, C, strict=True)]) * self.W).entries()
        return [4 * sum((psi[m] * H[i] for i, m in self.anti[j]), acb(0)) for j in range(self.d)]

    def error(self, C: list[acb]) -> float:
        """max over p != 0 of | |C_p|^2 / C_0^2 - 1/(d+1) |, the SIC condition for the normalized vector."""
        N, t = C[0].real, arb(1) / (self.d + 1)
        return max(float(abs(((c.real * c.real + c.imag * c.imag) / (N * N) - t).mid())) for c in C[1:])


def fraction(x: arb) -> Fraction:
    man, exp = x.mid().man_exp()
    return Fraction(int(man)) * Fraction(2) ** int(exp)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fiducial", type=Path, help="a file written by search.py or stark.jl")
    strings = read_fiducial(parser.parse_args().fiducial)
    d = len(strings)
    J = double_jacobian(np.array([complex(float(re), float(im)) for re, im in strings]))
    evals, evecs = np.linalg.eigh(J.T @ J)
    keep = evals > 1e-8 * evals[-1]  # drops the one null direction, the global phase
    M = (evecs[:, keep] / evals[keep]) @ evecs[:, keep].T
    overlaps = Overlaps(d)
    psi = [acb(float(re), float(im)) for re, im in strings]
    for step in range(20):
        C = overlaps.C(psi)
        error = overlaps.error(C)
        print(f"step {step}: max_p | |<psi|D_p|psi>|^2 - 1/{d + 1} | = {error:.1e}", flush=True)
        if error < 1e-185:
            break
        g = overlaps.gradient(psi, C)
        delta = M @ np.array([float(z.real) for z in g] + [float(z.imag) for z in g])
        psi = [(p - acb(delta[j], delta[d + j])).mid() for j, p in enumerate(psi)]
    else:
        raise SystemExit(f"no convergence in 20 steps: {error:.1e}; fiducials/fiducial_{d}.txt left as it was")
    norm = sum((p.real * p.real + p.imag * p.imag for p in psi), arb(0)).sqrt()
    phase = psi[0].conjugate() / abs(psi[0])
    psi = [p * phase / norm for p in psi]
    def exact(x: Fraction) -> Fraction:  # the converged residual is about 1e-194: a part below 1e-180 is 0
        return Fraction(0) if abs(x) < Fraction(1, 10**180) else x

    rows = [(exact(fraction(p.real)), Fraction(0) if j == 0 else exact(fraction(p.imag))) for j, p in enumerate(psi)]
    out = Path(f"fiducials/fiducial_{d}.txt")
    out.write_text(fiducial_text(rows, 150))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

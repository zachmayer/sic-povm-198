# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy"]
# ///
"""Search for SIC-POVM fiducials in dimension d with Scott's numerical method (arXiv:1703.03993).

Each trial draws a random vector in an eigenspace of Zauner's order-3 unitary and minimizes the
frame potential F(psi) = sum_p |<psi|D_p|psi>|^4 / <psi|psi>^4, summed over the d^2 displacements
D_p, with L-BFGS. F >= 2d/(d+1), with equality exactly at SIC fiducials. A trial that ends within
1e-5 of the bound is polished by Gauss-Newton on the SIC equations and kept if every overlap
|<psi|D_p|psi>|^2 matches 1/(d+1) to within 1e-10.

The search stops at its first SIC and writes it to fiducial_<d>_seed<seed>_trial<trial>.txt, to
double precision, in the file layout and conventions of Scott and Grassl (arXiv:0910.5784);
refine.py takes it to 150 digits.

  uv run search.py 19                 # a SIC in d = 19 in under a second
  uv run search.py 198 --minutes 60   # what one core of the cloud search did
"""

import argparse
import re
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy.fft as sfft
from scipy.optimize import minimize


def tau_power(d: int, e: np.ndarray) -> np.ndarray:
    """tau^e for tau = -exp(i pi/d) = exp(i pi (d+1)/d), with exponents taken mod 2d."""
    return np.exp(1j * np.pi * ((((d + 1) * e) % (2 * d)) / d))


def zauner_eigenspace(d: int) -> np.ndarray:
    """Orthonormal basis of the largest eigenspace of U_Z = exp(i pi (d-1)/12)/sqrt(d) sum_rs tau^(2rs + (d+1)r^2) |r><s|.

    For odd d, U_Z is Scott and Grassl's Zauner unitary. For even d it is D_(0,d/2) times theirs,
    so its eigenvectors become theirs under D_(d/2,d/2); see scott_convention.
    """
    r = np.arange(d)
    U = np.exp(1j * np.pi * (d - 1) / 12) / np.sqrt(d) * tau_power(d, 2 * np.outer(r, r) + (d + 1) * (r**2)[:, None])
    root = (U @ U @ U)[0, 0] ** (1 / 3)  # U_Z^3 is a scalar
    bases = []
    for m in range(3):
        lam = root * np.exp(2j * np.pi * m / 3)
        u, s, _ = np.linalg.svd((np.eye(d) + U / lam + (U @ U) / lam**2) / 3)  # projector onto eigenvalue lam
        bases.append(u[:, : int(np.sum(s > 0.5))])
    return max(bases, key=lambda B: B.shape[1])


class Problem:
    """The frame potential on one Zauner eigenspace, in real coordinates x = (Re c, Im c), psi = B c."""

    def __init__(self, d: int) -> None:
        self.d = d
        self.B = zauner_eigenspace(d)
        self.n = self.B.shape[1]
        k = np.arange(d)
        self.shift = (k[None, :] + k[:, None]) % d  # [a, k] -> k + a
        self.unshift = (k[None, :] - k[:, None]) % d  # [a, j] -> j - a
        self.rows = k[:, None]

    def psi(self, x: np.ndarray) -> np.ndarray:
        return self.B @ (x[: self.n] + 1j * x[self.n :])

    def overlaps(self, psi: np.ndarray) -> np.ndarray:
        """|<psi|D_(a,b)|psi>|^2 for all (a, b), unnormalized: one length-d FFT per shift a."""
        C = sfft.fft(np.conj(psi[self.shift]) * psi[None, :], axis=1, workers=1)
        return C.real**2 + C.imag**2

    def potential(self, x: np.ndarray) -> tuple[float, np.ndarray]:
        """F and its gradient in x."""
        c = x[: self.n] + 1j * x[self.n :]
        psi = self.B @ c
        N = float(np.vdot(c, c).real)
        C = sfft.fft(np.conj(psi[self.shift]) * psi[None, :], axis=1, workers=1)
        A2 = C.real**2 + C.imag**2
        f = float(np.sum(A2 * A2))
        H = sfft.fft(A2 * np.conj(C), axis=1, workers=1)
        G = 4.0 * (psi[None, :] * H)[self.rows, self.unshift].sum(axis=0)  # df/d conj(psi)
        F = f / N**4
        Gc = self.B.conj().T @ (G / N**4 - 4.0 * F / N * psi)
        return F, np.concatenate([2.0 * Gc.real, 2.0 * Gc.imag])

    def random_start(self, rng: np.random.Generator) -> np.ndarray:
        x = rng.standard_normal(2 * self.n)
        return x / np.linalg.norm(x)

    def residuals(self, x: np.ndarray) -> np.ndarray:
        """|<psi|D_p|psi>|^2 / <psi|psi>^2 - 1/(d+1) for every p != 0."""
        psi = self.psi(x)
        return (self.overlaps(psi) / np.vdot(psi, psi).real ** 2).ravel()[1:] - 1.0 / (self.d + 1)

    def jacobian(self, x: np.ndarray) -> np.ndarray:
        """d residuals / dx, one batch of FFTs per coordinate of x."""
        psi = self.psi(x)
        N = np.vdot(psi, psi).real
        C = sfft.fft(np.conj(psi[self.shift]) * psi[None, :], axis=1, workers=1)
        A2 = C.real**2 + C.imag**2
        columns = []
        for e in np.eye(2 * self.n):
            dpsi = self.psi(e)
            dV = np.conj(dpsi[self.shift]) * psi[None, :] + np.conj(psi[self.shift]) * dpsi[None, :]
            dC = sfft.fft(dV, axis=1, workers=1)
            dA2 = 2.0 * (C.real * dC.real + C.imag * dC.imag)
            dN = 2.0 * np.vdot(psi, dpsi).real
            columns.append((dA2 / N**2 - 2.0 * A2 * dN / N**3).ravel()[1:])
        return np.stack(columns, axis=1)

    def trial(self, x0: np.ndarray) -> np.ndarray | None:
        """Minimize F from x0; return the polished unit SIC fiducial, or None if this trial missed."""
        x = minimize(
            self.potential,
            x0,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 20000, "maxcor": 10, "ftol": 1e-16, "gtol": 1e-11, "maxfun": 80000},
        ).x
        psi = self.psi(x)
        A = self.overlaps(psi) / np.vdot(psi, psi).real ** 2
        if np.sum(A * A) - 2.0 * self.d / (self.d + 1) > 1e-5:
            return None
        x = x / np.linalg.norm(self.psi(x))
        for _ in range(8):  # Gauss-Newton: quadratic convergence near a SIC
            r = self.residuals(x)
            if np.max(np.abs(r)) < 1e-15:
                break
            x = x + np.linalg.lstsq(self.jacobian(x), -r, rcond=1e-12)[0]
            x = x / np.linalg.norm(self.psi(x))
        if np.max(np.abs(self.residuals(x))) > 1e-10:
            return None
        return self.psi(x)


def scott_convention(psi: np.ndarray) -> np.ndarray:
    """Map an eigenvector of our U_Z to Scott and Grassl's (D_(d/2,d/2) for even d), then fix the
    gauge: unit norm, first component real and positive."""
    d = psi.size
    if d % 2 == 0:
        psi = np.roll(psi * (-1.0) ** np.arange(d), d // 2)
    psi = psi * np.conj(psi[0]) / abs(psi[0]) / np.linalg.norm(psi)
    psi[0] = psi[0].real
    return psi


def search(d: int, seed: int, minutes: float = 1.0) -> tuple[int, tuple[int, str] | None]:
    """Run trials until one finds a SIC or time runs out.

    Returns the number of trials run and, for a SIC, its trial number and fiducial-file text.
    """
    problem = Problem(d)
    rng = np.random.default_rng(seed)
    deadline, k = time.time() + 60 * minutes, 0
    while time.time() < deadline:
        k += 1
        psi = problem.trial(problem.random_start(rng))
        if psi is not None:
            return k, (k, double_precision_file(psi))
    return k, None


def double_precision_file(psi: np.ndarray) -> str:
    return fiducial_text([(Fraction(z.real), Fraction(z.imag)) for z in scott_convention(psi)], 17)


def scott_number(x: Fraction, digits: int) -> str:
    """x as '+.<digits>e<exponent>', correctly rounded (half-even), as in Scott and Grassl's files."""
    if x == 0:
        return "+." + "0" * digits + "e0"
    a = abs(x)
    e = len(str(a.numerator)) - len(str(a.denominator))
    while a >= Fraction(10) ** e:
        e += 1
    while a < Fraction(10) ** (e - 1):
        e -= 1
    m = round(a * Fraction(10) ** (digits - e))
    if m == 10**digits:
        m, e = m // 10, e + 1
    return f"{'+' if x > 0 else '-'}.{m:0{digits}d}e{e}"


def fiducial_text(components: list[tuple[Fraction, Fraction]], digits: int) -> str:
    """One line per component: real and imaginary parts run together, ending in 'i'."""
    return "".join(scott_number(re, digits) + scott_number(im, digits) + "i\n" for re, im in components)


def read_fiducial(path: Path) -> list[tuple[str, str]]:
    """The (real, imaginary) parts of each component as plain decimal strings."""
    number = re.compile(r"([+-])\.(\d+e-?\d+)")
    return [tuple(f"{s}0.{x}" for s, x in number.findall(line)) for line in path.read_text().split()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("d", type=int, help="dimension")
    parser.add_argument("--seed", type=int, default=0, help="random seed (default 0)")
    parser.add_argument("--minutes", type=float, default=1.0, help="time budget (default 1)")
    args = parser.parse_args()
    count, found = search(args.d, args.seed, args.minutes)
    if found is None:
        print(f"d = {args.d}, seed {args.seed}: no SIC in {count} trials")
        return
    k, text = found
    path = Path(f"fiducial_{args.d}_seed{args.seed}_trial{k}.txt")
    path.write_text(text)
    print(f"d = {args.d}, seed {args.seed}: SIC in trial {k}, written to {path}")


if __name__ == "__main__":
    main()

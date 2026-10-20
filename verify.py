# /// script
# requires-python = ">=3.10"
# dependencies = ["mpmath"]
# [tool.uv]
# exclude-newer = "2026-10-19T00:00:00Z"
# ///
"""Check that files hold SIC-POVM fiducials in Scott and Grassl's conventions; by default every fiducial_<d>.txt here.

A unit vector psi in C^d is a SIC fiducial iff |<psi|X^a Z^b|psi>|^2 = 1/(d+1) for every (a, b) != (0, 0),
where X|k> = |k+1 mod d> and Z|k> = exp(2 pi i k/d)|k>. Scott and Grassl's conventions (arXiv:0910.5784):
psi_0 real and positive, and psi an eigenvector of their order-3 Zauner unitary
U_Z[r,s] = exp(i pi (d-1)/12)/sqrt(d) exp(i pi (2rs + (d+1) r^2)/d), whose eigenvectors Scott's Zauner-type
fiducials in arXiv:1703.03993 are (checked: d = 52, 55, 60, 67, 70, 88). This check shares no code with the other scripts.

python verify.py [FILE ...]             # double precision, standard library only, about 1 s per file
uv run verify.py [FILE ...] --precise   # 160 digits, about 3 min per file
"""

import cmath
import re
import sys
from pathlib import Path

precise = "--precise" in sys.argv
if precise:
    import mpmath

    mpmath.mp.dps = 160
    one, cplx, root = mpmath.mpf(1), mpmath.mpc, lambda t: mpmath.expjpi(2 * mpmath.mpf(t))
else:
    one, cplx, root = 1.0, lambda re_, im: complex(float(re_), float(im)), lambda t: cmath.exp(2j * cmath.pi * t)

files = [Path(arg) for arg in sys.argv[1:] if arg != "--precise"]
here = sorted(Path(__file__).parent.glob("fiducial_*.txt"))
failed = False
for path in files or [p for p in here if re.fullmatch(r"fiducial_\d+\.txt", p.name)]:
    rows = [re.findall(r"[+-]\.\d+e-?\d+", line) for line in path.read_text().splitlines() if line.strip()]
    psi = [cplx(re_.replace(".", "0.", 1), im.replace(".", "0.", 1)) for re_, im in rows]
    d = len(psi)
    w = [root(one * e / (2 * d)) for e in range(2 * d)]  # exp(i pi e/d)
    worst = abs(sum(abs(z) ** 2 for z in psi) - 1)
    for a in range(d):
        v = [psi[(k + a) % d].conjugate() * psi[k] for k in range(d)]
        for b in range(d):
            if (a, b) != (0, 0):
                overlap = sum(v[k] * w[2 * b * k % (2 * d)] for k in range(d))
                worst = max(worst, abs(abs(overlap) ** 2 - one / (d + 1)))
    u = root(one * (d - 1) / 24) / (one * d) ** (one / 2)
    U_psi = [u * sum(w[(2 * r * s + (d + 1) * r * r) % (2 * d)] * psi[s] for s in range(d)) for r in range(d)]
    lam = sum(p.conjugate() * q for p, q in zip(psi, U_psi, strict=True))
    zauner = max(abs(q - lam * p) for p, q in zip(psi, U_psi, strict=True))
    named = re.fullmatch(r"fiducial_(\d+)\.txt", path.name)
    ok = d > 1 and (not named or int(named[1]) == d) and psi[0].imag == 0 and psi[0].real > 0
    ok = ok and max(worst, zauner) < (1e-147 if precise else 1e-12)
    failed |= not ok
    k = round(cmath.phase(complex(lam)) / (2 * cmath.pi / 3)) % 3  # U_Z's eigenvalues are the cube roots of unity
    print(f"{path.name}: d = {d}, max deviation {float(worst):.1e}, Zauner residual {float(zauner):.1e}, eigenvalue ω^{k} {'PASS' if ok else 'FAIL'}")
sys.exit(failed)

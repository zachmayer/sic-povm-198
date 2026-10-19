"""Check that files hold SIC-POVM fiducials; by default every fiducial_<d>.txt beside this script.

A vector psi in C^d is a SIC fiducial iff |<psi|X^a Z^b|psi>|^2 = |psi|^4/(d+1) for every
(a, b) != (0, 0), where X|k> = |k+1 mod d> and Z|k> = exp(2 pi i k/d)|k>.
File layout and conventions follow Scott and Grassl (arXiv:0910.5784, arXiv:1703.03993).
This check shares no code with the other scripts.

python verify.py [FILE ...]             # double precision, standard library only, about 1 s per file
python verify.py [FILE ...] --precise   # 160 digits, needs mpmath, 1-3 min per file
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
    norm2 = sum(abs(z) ** 2 for z in psi)
    w = [root(one * e / d) for e in range(d)]
    worst = 0 * one
    for a in range(d):
        v = [psi[(k + a) % d].conjugate() * psi[k] for k in range(d)]
        for b in range(d):
            if (a, b) != (0, 0):
                overlap = sum(v[k] * w[b * k % d] for k in range(d))
                worst = max(worst, abs(abs(overlap) ** 2 / norm2**2 - one / (d + 1)))
    ok = worst < (1e-140 if precise else 1e-12)
    failed |= not ok
    print(f"{path.name}: max over {d * d - 1} displacements of | |<psi|D_p|psi>|^2 - 1/{d + 1} | = {float(worst):.2e} {'PASS' if ok else 'FAIL'}")
sys.exit(failed)

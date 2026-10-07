# SIC-POVM in dimension 198

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for every pair i ≠ j. Zauner conjectured in 1999 that one exists in every dimension, and that's still open. SICs are known in every dimension up to 196 ([Bengtsson and Grassl](https://arxiv.org/abs/2512.13201)) and in some higher ones. As far as I can tell, this is the first in dimension 198.

The usual SICs are the orbit of one fiducial vector under the d² shifts X^a Z^b, so that one vector is the whole certificate. `fiducial_198.txt` has its 198 components to 150 digits, in the file format and conventions of Scott and Grassl ([arXiv:0910.5784](https://arxiv.org/abs/0910.5784), [arXiv:1703.03993](https://arxiv.org/abs/1703.03993)). It is an eigenvector of their order-3 Zauner unitary.

```sh
python verify.py             # max deviation of |<psi|D_p|psi>|^2 from 1/199: about 1e-17 (standard library, ~1 s)
python verify.py --precise   # about 1e-151 (needs mpmath, ~1 min)
```

Found with Scott's numerical method on rented cloud CPUs.

The problem comes from Zauner's thesis ([1999, in German](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English translation, 2011](https://doi.org/10.1142/S0219749911006776)) and from Renes, Blume-Kohout, Scott and Caves, [Symmetric informationally complete quantum measurements](https://arxiv.org/abs/quant-ph/0310075) (2004).

*Mostly AI-generated: the search code, the checker and this README are by Claude Opus 5.5.*

Public domain (CC0).

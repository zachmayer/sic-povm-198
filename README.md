# SIC-POVMs in dimensions 197 and 198

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for every pair i ≠ j. Zauner conjectured in 1999 that one exists in every dimension. SICs are known in every dimension up to 196 ([Bengtsson and Grassl](https://arxiv.org/abs/2512.13201)) and in some higher ones; as far as I know, these are the first written down in dimensions 197 and 198.

The usual SICs are the orbit of one fiducial vector under the d² shifts X^a Z^b, so that one vector is the whole certificate. `fiducial_197.txt` and `fiducial_198.txt` have their components to 150 digits, in the file format and conventions of Scott and Grassl ([arXiv:0910.5784](https://arxiv.org/abs/0910.5784), [arXiv:1703.03993](https://arxiv.org/abs/1703.03993)). `fiducial_198.txt` is an eigenvector of their order-3 Zauner unitary. `fiducial_197.txt` is not; it is an eigenvector of a conjugate one, the Clifford unitary of the order-3 symplectic matrix (0 −3; 66 −1) mod 197 in place of Zauner's (0 −1; 1 −1).

```sh
python verify.py                     # max deviation of |<psi|D_p|psi>|^2 from 1/199: about 1e-17 (standard library, ~1 s)
python verify.py --precise           # about 1e-151 (needs mpmath, ~1 min)
python verify.py fiducial_197.txt    # the same checks for d = 197
```

## How they were found

**d = 198** by Scott's numerical method: random starting vectors in an eigenspace of Zauner's unitary, each minimized with L-BFGS. The search ran on 800 cloud cores with [Modal](https://modal.com) and found this SIC in its first hour, after 10–17 million trials (500–800 core-hours, about 0.2 core-seconds per trial). Gauss–Newton in 640-bit arithmetic then took it to 150 digits.

**d = 197** by the algebraic method of Appleby, Flammia and Kopp ([arXiv:2501.03970](https://arxiv.org/abs/2501.03970)), which builds a fiducial from Stark units instead of searching for one. It ran on a laptop in hours, using Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl) with modifications; the same random search at d = 197 ran 46.5 million trials without finding a SIC. The steps: take the admissible form ⟨22, −66, 1⟩ of conductor 1 (`AdmissibleTuple(197)` defaults to conductor 3, a field four times larger); compute its ghost overlaps and raise their precision; identify them exactly in the ring class field of ℚ(√1067); and apply the Galois automorphism that sends √1067 to −√1067, which turns the ghost into a SIC. At d = 197 the stock package's integer-relation step fails and several others take hours, so those needed faster replacements.

Each script lists its dependencies, so [uv](https://docs.astral.sh/uv/) runs it directly:

```sh
uv run search.py 19                            # a SIC in d = 19, written to fiducial_19_seed0_trial<k>.txt
uv run refine.py fiducial_19_seed0_trial*.txt  # takes it to 150 digits: fiducial_19.txt
uv run refine.py fiducial_198.txt              # re-derives this file's 150 digits; it comes back unchanged
uvx modal run modal_search.py --d 198 --containers 25 --minutes 60  # the search that found it; needs Modal
```

Floating-point differences between machines send each trial down its own path, so a repeat of the search finds SICs at its own pace rather than retracing this one.

| File | |
| --- | --- |
| `search.py` | Scott's method: random starts, L-BFGS, a Gauss–Newton polish |
| `modal_search.py` | runs `search.py` on many cloud cores |
| `refine.py` | takes a fiducial from double precision to 150 digits |
| `verify.py` | checks a fiducial; shares no code with the other scripts |
| `fiducial_197.txt`, `fiducial_198.txt` | the results |

The problem comes from Zauner's thesis ([1999, in German](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English translation, 2011](https://doi.org/10.1142/S0219749911006776)) and from Renes, Blume-Kohout, Scott and Caves, [Symmetric informationally complete quantum measurements](https://arxiv.org/abs/quant-ph/0310075) (2004).

*Mostly AI-generated: Gemini 4 Argon found the d = 197 SIC; the code and this README are by Claude Opus 5.5.*

Public domain (CC0).

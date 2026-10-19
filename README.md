# SIC-POVMs in 15 new dimensions

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for every pair i ≠ j. Zauner conjectured in 1999 that one exists in every dimension. SICs are known in every dimension up to 196 ([Bengtsson and Grassl](https://arxiv.org/abs/2512.13201)) and in some higher ones; as far as we know, these are the first written down in dimensions 197, 198, 203, 207, 211, 215, 223, 227, 239, 241, 243, 245, 249, 251 and 253.

Each SIC is the orbit of one fiducial vector under the d² shifts X^a Z^b, so that one vector is the whole certificate. `fiducial_<d>.txt` has its components to 150 digits, in the file format and conventions of Scott and Grassl ([arXiv:0910.5784](https://arxiv.org/abs/0910.5784), [arXiv:1703.03993](https://arxiv.org/abs/1703.03993)): unit norm, first component real and positive, and an eigenvector of their order-3 Zauner unitary.

```sh
python verify.py                    # each file's max deviation of |<psi|D_p|psi>|^2 from 1/(d+1): about 1e-17 (standard library, 11 s)
python verify.py --precise          # about 1e-151 (needs mpmath, 1-3 min per file)
python verify.py fiducial_197.txt   # one file
```

## How they were found

**d = 198** by Scott's numerical method: random starting vectors in an eigenspace of Zauner's unitary, each minimized with L-BFGS. The search ran on 800 cloud cores with [Modal](https://modal.com) and found this SIC in its first hour, after 10–17 million trials (500–800 core-hours, about 0.2 core-seconds per trial). Gauss–Newton in 640-bit arithmetic then took it to 150 digits.

**The other 14** by the algebraic method of Appleby, Flammia and Kopp ([arXiv:2501.03970](https://arxiv.org/abs/2501.03970)), which builds a fiducial from Stark units instead of searching for one: it identifies the overlaps of a "ghost" exactly in a ring class field of ℚ(√((d+1)(d−3))), and the Galois automorphism that negates that square root turns them into the overlaps of a SIC. `stark.jl` runs it with Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl), pinned in `Manifest.toml`, replacing the steps that fail or take hours at these sizes. Each dimension takes 1–4 minutes and about 6 GB on a laptop; the random search at d = 197 had run 46.5 million trials without finding a SIC.

The Python scripts list their dependencies, so [uv](https://docs.astral.sh/uv/) runs them directly:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once, a few minutes
julia --project=. stark.jl 197                        # fiducial_197_stark.txt, to double precision
uv run refine.py fiducial_197_stark.txt               # takes it to 150 digits: rewrites fiducial_197.txt, unchanged
uv run search.py 19                                   # a SIC in d = 19 by search: fiducial_19_seed0_trial<k>.txt
uv run refine.py fiducial_19_seed0_trial*.txt         # fiducial_19.txt
uvx modal run modal_search.py --d 198 --containers 25 --minutes 60  # the search that found d = 198; needs Modal
```

`stark.jl` is deterministic, so it and `refine.py` reproduce each algebraic file digit for digit. The search is not: floating-point differences between machines send each trial down its own path, so a repeat finds SICs at its own pace rather than retracing this one, and `uv run refine.py fiducial_198.txt` re-derives that file instead.

| File | |
| --- | --- |
| `stark.jl`, `Project.toml`, `Manifest.toml` | the algebraic method |
| `search.py` | Scott's method: random starts, L-BFGS, a Gauss–Newton polish |
| `modal_search.py` | runs `search.py` on many cloud cores |
| `refine.py` | takes a fiducial from double precision to 150 digits |
| `verify.py` | checks fiducials; shares no code with the other scripts |

The problem comes from Zauner's thesis ([1999, in German](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English translation, 2011](https://doi.org/10.1142/S0219749911006776)) and from Renes, Blume-Kohout, Scott and Caves, [Symmetric informationally complete quantum measurements](https://arxiv.org/abs/quant-ph/0310075) (2004).

*Mostly AI-generated: Gemini 4 Argon first found the algebraic SICs; Claude Opus 5.5 found d = 198 and wrote the code (`stark.jl` from Argon's) and this README.*

Public domain (CC0).

# SIC-POVMs in 15 new dimensions

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for every pair i ≠ j. Zauner conjectured in 1999 ([thesis](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English](https://doi.org/10.1142/S0219749911006776)) that one exists in every dimension. SICs are known in every dimension up to 196 ([Bengtsson and Grassl](https://arxiv.org/abs/2512.13201)) and in some higher ones; as far as we know, these are the first written down in dimensions 197, 198, 203, 207, 211, 215, 223, 227, 239, 241, 243, 245, 249, 251 and 253.

Each SIC here is the orbit of a fiducial vector ψ under the d² displacements D_p = X^a Z^b (up to phase), so ψ determines it. `fiducial_<d>.txt` gives ψ numerically, to 150 digits, in the file format and conventions of Scott and Grassl ([arXiv:0910.5784](https://arxiv.org/abs/0910.5784)) and Scott ([arXiv:1703.03993](https://arxiv.org/abs/1703.03993)): unit norm, ψ_0 real and positive, and an eigenvector of their order-3 Zauner unitary. These are numerical solutions, not exact ones. `verify.py` checks all of it and shares no code with the other scripts:

```sh
python3 verify.py            # every file: max deviation about 1e-17 (Python 3.10+, standard library, 15 s)
uv run verify.py --precise   # 160 digits: about 1e-151 (about 3 min per file)
```

## How they were found

**d = 198** by Scott's numerical method: random starts in the eigenvalue-1 eigenspace of Zauner's unitary, each minimized with L-BFGS. It ran on 800 cloud cores with [Modal](https://modal.com), about 0.2 core-seconds per trial, and found this SIC within its first hour.

**The other 14** by the algebraic method of Appleby, Flammia and Kopp ([arXiv:2501.03970](https://arxiv.org/abs/2501.03970)), which builds a fiducial from Stark units instead of searching for one. `stark.jl` computes a ghost's overlaps on one Galois orbit, recognizes them by LLL as integers of a ring class field H of ℚ(√((d+1)(d−3))), applies an automorphism of H that negates the square root, and completes the SIC numerically by alternating projections. It runs Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl), pinned by commit, and replaces the steps that fail or take hours at these sizes; each dimension takes a few minutes and about 5 GB. Its choices are canonical (see its header): of a dimension's Galois-conjugate SICs it returns a fixed one, and of that SIC's images under the Clifford unitaries that commute with U_Z, a fixed fiducial.

`refine.py` takes each fiducial from double precision to 150 digits by Gauss–Newton in 640-bit arithmetic.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once; Julia 1.13.2
julia --project=. stark.jl 197                        # fiducial_197_stark.txt, to double precision
uv run refine.py fiducial_197_stark.txt               # rewrites fiducial_197.txt, unchanged
uv run refine.py fiducial_198.txt                     # re-derives its 150 digits from the first 16, unchanged
uv run search.py 19                                   # Scott's search: a SIC in d = 19 in a second
uvx modal run modal_search.py --d 198 --containers 25 --minutes 60  # a search the size of the one that found d = 198
```

`stark.jl` then `refine.py` rewrite each algebraic file byte for byte, whatever the random seeds or the order of dimensions (tested with Julia 1.13.2 on macOS, Apple silicon). Searches are not repeatable: each trial's path depends on floating-point details, so a repeat finds SICs at its own pace.

*Mostly AI-generated: Gemini 4 Argon first found the algebraic SICs; Claude Opus 5.5 found d = 198 and wrote the code (`stark.jl` from Argon's) and this README.*

Public domain (CC0).

# SIC-POVMs in every dimension from 197 to 300

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for every pair i ≠ j. Zauner conjectured in 1999 ([thesis](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English](https://doi.org/10.1142/S0219749911006776)) that one exists in every dimension. SICs are known in every dimension up to 196 ([Bengtsson and Grassl](https://arxiv.org/abs/2512.13201)) and in some higher ones.

This repo has a numerical Weyl–Heisenberg SIC fiducial for every d from 197 to 300 (104 files), so with Bengtsson and Grassl's, every d ≤ 300 has one. Which of these dimensions are new rests on a literature search (Grassl's tables, Scott 2017, Bengtsson–Grassl, Appleby–Flammia–Kopp, recent arXiv), and our view of unpublished and talk-only results is spotty:

- **No SIC found published before (95 dimensions):** every d in 197..300 except the nine below.
- **SIC known before (9):** in papers, 199, 228 and 259 ([Scott](https://arxiv.org/abs/1703.03993)) and 292 ([Bengtsson, Grassl and McConnell](https://arxiv.org/abs/2403.02872)); only in talks or surveys, 201 ([Grassl, 2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), p. 21) and 204, 224, 255 and 288 ([Horodecki, Rudnicki and Życzkowski](https://arxiv.org/abs/2002.03233), crediting Grassl). Ours may lie in different orbits from those; we have not checked orbit equivalence.

This is about explicit fiducials. A recent manuscript, [*Stark units for real quadratic fields and reciprocity laws*](https://www.danrad.net/papers/stark-reciprocity.pdf), claims the results needed for Zauner's conjecture in general; we have not checked it.

Each SIC is the orbit of a fiducial ψ under the d² displacements D_p = X^a Z^b (up to phase), so ψ determines it. `fiducial_<d>.txt` gives ψ to 150 digits in the file format and conventions of Scott and Grassl ([arXiv:0910.5784](https://arxiv.org/abs/0910.5784)) and Scott (arXiv:1703.03993): unit norm, ψ_0 real and positive, and an eigenvector of their order-3 Zauner unitary U_Z. These are numerical solutions, not exact ones. `verify.py` checks all of it and shares no code with the other scripts:

```sh
python3 verify.py            # every file, double precision; max deviation at most about 5e-17 (Python 3.10+, standard library)
uv run verify.py --precise   # 160 digits; PASS means deviation and Zauner residual below 1e-147 (here at most 1e-150); minutes per file
```

## How they were found

**d = 198** by Scott's numerical method: random starts in the eigenvalue-1 eigenspace of U_Z, each minimized with L-BFGS (`search.py`). It ran on 800 cloud cores with [Modal](https://modal.com) (`modal_search.py`) and found this SIC within its first hour.

**The other 103** by the algebraic method of Appleby, Flammia and Kopp ([arXiv:2501.03970](https://arxiv.org/abs/2501.03970)), which builds a fiducial from Stark units instead of searching for one. `stark.jl` computes a ghost's overlaps on one Galois orbit to thousands of bits (Gauss–Newton), recognizes integral trace forms of them by LLL in a ring class field H of ℚ(√((d+1)(d−3))), applies an automorphism of H that negates the square root, reconstructs the SIC's overlaps on that orbit, and completes the fiducial by alternating projections. It runs Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl), pinned by commit, handles odd and even d and forms beyond Zauner.jl's table, and replaces the steps that fail or take hours at these sizes (each replacement says why). One thread per dimension: from 4 minutes to about a day at class number h ≈ 40, under 25 GB. Its choices are canonical (see its header): of a dimension's Galois-conjugate SICs it returns a fixed one, and of that SIC's images under the Clifford unitaries that commute with U_Z, a fixed fiducial. The rule goes by Galois data, not eigenvalue, so some files lie in U_Z's ω or ω² eigenspace (`verify.py` prints each file's).

`refine.py` takes each fiducial from double precision to 150 digits by Gauss–Newton in 640-bit arithmetic, and writes a real or imaginary part below 1e-180 as 0 (the converged residual is about 1e-194).

[RESULTS.md](RESULTS.md) gives, for the 89 dimensions in 199..300 added last, the form, conductor, class number, eigenvalue, machine, run time and memory.

## Reproducing

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once; Julia 1.13.2
julia --project=. stark.jl 197                        # fiducial_197_stark.txt, to double precision
uv run refine.py fiducial_197_stark.txt               # rewrites fiducial_197.txt, unchanged
uv run refine.py fiducial_198.txt                     # re-derives its 150 digits from the first 16, unchanged
uv run search.py 19                                   # Scott's search: a SIC in d = 19 in a second
uvx modal run modal_search.py --d 198 --containers 25 --minutes 60  # a search the size of the one that found d = 198
```

With `STARK_CHECKPOINTS=<dir>` (an existing directory, fresh for each version of `stark.jl`), `stark.jl` saves each Gauss–Newton step and its other expensive intermediates there, and a rerun resumes from them.

`stark.jl` then `refine.py` rewrite the 14 algebraic files from 197 to 253 byte for byte, whatever the random seeds or the order of dimensions (Julia 1.13.2, macOS on Apple silicon). The early ones among the other 89 ran on earlier versions of `stark.jl` (each of which regenerated those 14 byte for byte); of those, only 206 and 214 were rerun with this one, and their refined files came out byte-identical. Of d = 199..300, 16 ran on Modal (Linux x86, same Project and Manifest), 294 on Modal then the Mac, and the rest on the Mac. The double-precision `stark.jl` output differs between the two in the 17th digit, but the refined files agreed byte for byte in both dimensions run on both (197 and 224). Searches are not repeatable: each trial's path depends on floating-point details, so a repeat finds SICs at its own pace.

*Mostly AI-generated: Gemini 4 Argon first found the algebraic SICs; Claude Opus 5.5 found d = 198, ran the sweep to 300, and wrote the code (`stark.jl` from Argon's) and this README.*

Public domain (CC0).
